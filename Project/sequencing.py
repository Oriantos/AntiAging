"""
Compare sequencing reads to a reference using sliding overlap Levenshtein
(see ``docs/levenshtein_distance.md``): slide the read along the reference, score each
overlap with classic Levenshtein (Wagner–Fischer / optional ``python-Levenshtein``),
normalize by overlap length, and take the minimum score.

Run as a script after setting ``.env`` values (paths, read limits, logging, tolerances).
``MAX_READS`` caps how many FASTQ reads (records) to process; ``None`` means the whole file.
``LEVENSHTEIN_PROCESSES`` sets the process-pool size for read scoring; ``1`` keeps it serial.
``RUN_ALL_RECORDS`` forces unlimited FASTQ processing and disables sampling limits.

Logging: ``main`` calls :func:`configure_logging`, which writes only to :data:`LOG_FILE`
(no console). If you import this module, call it yourself. Set ``LOG_LEVEL`` to
``logging.DEBUG`` for per-read and sliding-window detail. ``READ_PROGRESS_INTERVAL`` sets
how often INFO progress lines appear while scoring (``None`` disables).
"""

from __future__ import annotations

import csv
from concurrent.futures import ProcessPoolExecutor
import logging
import math
import os
import random
from pathlib import Path
from typing import Iterator, List, Optional, Tuple

from Project.post_process_csv import post_process_results_csv

_ROOT = Path(__file__).resolve().parent.parent

logger = logging.getLogger("Project.sequencing")
_WORKER_REFERENCE_SEQUENCE: Optional[str] = None


def _load_dotenv(dotenv_path: Path) -> None:
    """Load simple KEY=VALUE pairs from a .env file into os.environ."""
    if not dotenv_path.exists():
        return

    for raw_line in dotenv_path.read_text(encoding="utf-8").splitlines():
        line = raw_line.strip()
        if not line or line.startswith("#") or "=" not in line:
            continue
        key, value = line.split("=", 1)
        key = key.strip()
        if not key:
            continue
        value = value.strip()
        if value and value[0] == value[-1] and value[0] in {'"', "'"}:
            value = value[1:-1]
        os.environ.setdefault(key, value)


def _env_str(name: str, default: str) -> str:
    return os.getenv(name, default)


def _env_path(name: str, default: Path) -> Path:
    raw = os.getenv(name)
    candidate = Path(raw) if raw else default
    return candidate if candidate.is_absolute() else (_ROOT / candidate)


def _env_int(name: str, default: int) -> int:
    raw = os.getenv(name)
    if raw is None:
        return default
    return int(raw)


def _env_float(name: str, default: float) -> float:
    raw = os.getenv(name)
    if raw is None:
        return default
    return float(raw)


def _env_optional_int(name: str, default: Optional[int]) -> Optional[int]:
    raw = os.getenv(name)
    if raw is None:
        return default
    if raw.strip().lower() in {"", "none", "null"}:
        return None
    return int(raw)


def _env_log_level(name: str, default: int) -> int:
    raw = os.getenv(name)
    if raw is None:
        return default
    level_text = raw.strip().upper()
    if level_text.isdigit():
        return int(level_text)
    return int(getattr(logging, level_text, default))


def _env_bool(name: str, default: bool) -> bool:
    raw = os.getenv(name)
    if raw is None:
        return default
    return raw.strip().lower() in {"1", "true", "yes", "on"}


def _initialize_scoring_worker(reference_sequence: str) -> None:
    global _WORKER_REFERENCE_SEQUENCE
    _WORKER_REFERENCE_SEQUENCE = reference_sequence


_DOTENV_PATH = _ROOT / ".env"
_load_dotenv(_DOTENV_PATH)

# Edit these before running as a script (``python sequencing.py``).
SEQUENCES_FILE: str | Path = _env_path(
    "SEQUENCES_FILE",
    Path("data") / "reads" / "SRR31234567_1.fastq",
)
REFERENCE_FASTA: str | Path = _env_path(
    "REFERENCE_FASTA",
    Path("data") / "reference" / "reference.fasta",
)
RESULT_FILE: str | Path = _env_path(
    "RESULT_FILE",
    Path("output") / "levenshtein_best_matches.csv",
)
SEQUENCE_FORMAT: str = _env_str("SEQUENCE_FORMAT", "fastq")
MAX_READS: Optional[int] = _env_optional_int("MAX_READS", 10000)
SAMPLE_READ_POOL: Optional[int] = _env_optional_int("SAMPLE_READ_POOL", 10000)
SAMPLE_READ_COUNT: Optional[int] = _env_optional_int("SAMPLE_READ_COUNT", 100)
RANDOM_SEED: Optional[int] = _env_optional_int("RANDOM_SEED", 42)
LEVENSHTEIN_PROCESSES: Optional[int] = _env_optional_int("LEVENSHTEIN_PROCESSES", 1)
RUN_ALL_RECORDS: bool = _env_bool("RUN_ALL_RECORDS", False)

if LEVENSHTEIN_PROCESSES is not None and LEVENSHTEIN_PROCESSES < 1:
    raise ValueError("LEVENSHTEIN_PROCESSES must be a positive integer when set")

# Logging (call :func:`configure_logging` from ``main``; output goes to the file only).
LOG_FILE: str | Path = _env_path(
    "LOG_FILE",
    Path("output") / "sequencing.log",
)
LOG_LEVEL: int = _env_log_level("LOG_LEVEL", logging.INFO)
READ_PROGRESS_INTERVAL: Optional[int] = _env_optional_int("READ_PROGRESS_INTERVAL", 25)
LEVENSHTEIN_RESULT_LOG_EVERY: int = _env_int("LEVENSHTEIN_RESULT_LOG_EVERY", 10)
BEST_SCORE_EPSILON: float = _env_float("BEST_SCORE_EPSILON", 1e-15)
TIE_SCORE_ABS_TOL: float = _env_float("TIE_SCORE_ABS_TOL", 1e-12)

try:
    import Levenshtein as _Lev  # type: ignore

    def _levenshtein_fast(first: str, second: str) -> int:
        return int(_Lev.distance(first, second))

    _HAS_C_LEV = True
except ImportError:
    _HAS_C_LEV = False


def configure_logging(log_file: str | Path, log_level: int) -> None:
    """
    Send all ``Project.sequencing`` log records to *log_file* only
    (``propagate=False``, so nothing is printed to stderr/stdout).
    """
    log_path = Path(log_file).resolve()
    log_path.parent.mkdir(parents=True, exist_ok=True)

    for existing in logger.handlers[:]:
        logger.removeHandler(existing)
        existing.close()

    logger.setLevel(log_level)
    file_handler = logging.FileHandler(log_path, mode="w", encoding="utf-8")
    file_handler.setFormatter(
        logging.Formatter(
            "%(asctime)s | %(levelname)-8s | %(name)s | %(message)s",
            datefmt="%H:%M:%S",
        ),
    )
    logger.addHandler(file_handler)
    logger.propagate = False

    logger.debug(
        "logging configured: path=%s level=%s",
        log_path,
        logging.getLevelName(log_level),
    )


def levenshtein_distance(first_string: str, second_string: str) -> int:
    """
    Minimum number of single-character insertions, deletions, and substitutions
    needed to transform *first_string* into *second_string* (symmetric). Matches
    the recurrence in levenshtein_distance_markdown.md via the standard iterative DP.
    """
    if _HAS_C_LEV:
        return _levenshtein_fast(first_string, second_string)

    shorter, longer = (
        (first_string, second_string)
        if len(first_string) <= len(second_string)
        else (second_string, first_string)
    )
    short_len, long_len = len(shorter), len(longer)
    previous_row = list(range(long_len + 1))
    current_row = [0] * (long_len + 1)
    for row_index in range(1, short_len + 1):
        current_row[0] = row_index
        char_from_shorter = shorter[row_index - 1]
        for col_index in range(1, long_len + 1):
            substitution_cost = 0 if char_from_shorter == longer[col_index - 1] else 1
            current_row[col_index] = min(
                previous_row[col_index] + 1,
                current_row[col_index - 1] + 1,
                previous_row[col_index - 1] + substitution_cost,
            )
        previous_row, current_row = current_row, previous_row
    return int(previous_row[long_len])


def sliding_levenshtein(reference: str, read: str) -> tuple[float, float]:
    """
    Implements the overlap-based ``sliding_levenshtein`` pseudocode in
    ``levenshtein_distance.md``: the *read* slides along the *reference*
    (``reference_length = n``, ``read_length = m``, ``shift in range(0, n-m+1)``),
    overlapping fragments use ``read_index = reference_index - alignment_shift``.
    Return ``(best_score, worst_score)`` across all valid overlaps where the read is fully contained.
    Lists plus ``join`` replace ``+=`` for speed.
    """
    reference_length = len(reference)
    read_length = len(read)
    if reference_length == 0 or read_length == 0:
        return float("inf"), float("inf")

    best_score = float("inf")
    worst_score = float("-inf")
    for alignment_shift in range(0, reference_length - read_length + 1):
        overlap_reference_chars: List[str] = []
        overlap_read_chars: List[str] = []
        for reference_index in range(reference_length):
            read_index = reference_index - alignment_shift
            if 0 <= read_index < read_length:
                overlap_reference_chars.append(reference[reference_index])
                overlap_read_chars.append(read[read_index])
        if not overlap_reference_chars:
            continue
        overlap_from_reference = "".join(overlap_reference_chars)
        overlap_from_read = "".join(overlap_read_chars)
        edit_distance = levenshtein_distance(overlap_from_reference, overlap_from_read)
        best_score = min(best_score, edit_distance)
        worst_score = max(worst_score, edit_distance)

    if worst_score == float("-inf"):
        worst_score = float("inf")

    if logger.isEnabledFor(logging.DEBUG):
        num_shifts = reference_length - read_length + 1
        logger.debug(
            "sliding_levenshtein done: ref_len=%d read_len=%d alignment_shifts=%d -> "
            "best_score=%.8g worst_score=%.8g",
            reference_length,
            read_length,
            num_shifts,
            best_score,
            worst_score,
        )
    return best_score, worst_score


def _score_read_against_reference(
    reference_sequence: str,
    read_id: str,
    read_sequence: str,
) -> Tuple[str, float, float, int, float, str]:
    sliding_score, worst_score = sliding_levenshtein(reference_sequence, read_sequence)
    record_length = len(read_sequence)
    matching_precentage = sliding_score / record_length if record_length > 0 else float("inf")
    return (
        read_id,
        sliding_score,
        worst_score,
        record_length,
        matching_precentage,
        read_sequence,
    )


def _score_read_in_worker(read_item: Tuple[str, str]) -> Tuple[str, float, float, int, float, str]:
    if _WORKER_REFERENCE_SEQUENCE is None:
        raise RuntimeError("worker reference sequence was not initialized")
    read_id, read_sequence = read_item
    return _score_read_against_reference(_WORKER_REFERENCE_SEQUENCE, read_id, read_sequence)


def load_reference_fasta(path: str | Path) -> str:
    """Concatenate all non-header sequence lines from a FASTA file."""
    path = Path(path)
    logger.debug("loading reference FASTA: %s", path)
    chunks: List[str] = []
    sequence_lines = 0
    with open(path, encoding="utf-8", errors="replace") as fasta_file:
        for line in fasta_file:
            line = line.strip()
            if not line or line.startswith(">"):
                continue
            chunks.append(line)
            sequence_lines += 1
    reference = "".join(chunks)
    logger.info(
        "loaded reference: %s (%d sequence lines, %d bp)",
        path.name,
        sequence_lines,
        len(reference),
    )
    return reference


def iter_fastq_sequences(
    path: str | Path,
    *,
    max_reads: Optional[int] = None,
) -> Iterator[Tuple[str, str]]:
    """
    Yield (read_id, sequence) for each FASTQ record (4 lines per read).

    The file human_sequences/SRR31234567_1.fastq is standard FASTQ; records
    follow one another without blank lines between them.

    If *max_reads* is given, yield at most that many reads from the start of the
    file (each read is one 4-line FASTQ record). ``None`` means read until end of file.

    Note: each read occupies **four lines** in the file (header, sequence, ``+``,
    quality). A limit on *lines* would not match one row per read in the output;
    *max_reads* is therefore a **read count**, not a line count.
    """
    if max_reads is not None and max_reads < 0:
        raise ValueError("max_reads must be non-negative when given")

    path = Path(path)
    logger.debug(
        "opening FASTQ for streaming: %s (max_reads=%s)",
        path,
        max_reads if max_reads is not None else "unlimited",
    )
    reads_yielded = 0
    with open(path, encoding="utf-8", errors="replace") as fastq_file:
        while True:
            if max_reads is not None and reads_yielded >= max_reads:
                break
            header_line = fastq_file.readline()
            if not header_line:
                break
            header_line = header_line.rstrip("\r\n")
            sequence_line = fastq_file.readline()
            plus_line = fastq_file.readline()
            quality_line = fastq_file.readline()
            if not sequence_line or not plus_line or not quality_line:
                break
            sequence = sequence_line.rstrip("\r\n")
            read_id = (
                header_line[1:].split()[0]
                if header_line.startswith("@") and len(header_line) > 1
                else header_line
            )
            reads_yielded += 1
            yield read_id, sequence


def iter_blank_block_sequences(path: str | Path) -> Iterator[Tuple[str, str]]:
    """
    Yield (synthetic_id, sequence) for blocks separated by one or more empty lines.

    Within each block, sequence lines are stripped and concatenated. Lines starting
    with ``@`` (FASTQ header) are skipped; if the block looks like FASTQ, the line
    after ``@`` is taken as the sequence for that block.
    """
    path = Path(path)
    logger.debug("loading blank-block sequences: %s", path)
    file_text = path.read_text(encoding="utf-8", errors="replace")
    paragraph_blocks = file_text.split("\n\n")
    block_index = 0
    for raw_block in paragraph_blocks:
        non_empty_lines = [ln.rstrip("\r\n") for ln in raw_block.splitlines()]
        non_empty_lines = [ln for ln in non_empty_lines if ln.strip()]
        if not non_empty_lines:
            continue
        block_index += 1
        if non_empty_lines[0].startswith("@") and len(non_empty_lines) >= 2:
            yield f"block_{block_index}", non_empty_lines[1].strip()
        else:
            yield f"block_{block_index}", "".join(ln.strip() for ln in non_empty_lines)

    logger.info("parsed %d non-empty blank-block sequences from %s", block_index, path.name)


def find_best_matching_sequences(
    sequences_file: str | Path,
    reference_fasta: str | Path,
    result_file: str | Path,
    *,
    sequence_format: str = "fastq",
    max_reads: Optional[int] = None,
    sample_read_pool: Optional[int] = 10000,
    sample_read_count: Optional[int] = 100,
    run_all_records: bool = False,
    random_seed: Optional[int] = 42,
    levenshtein_processes: Optional[int] = 1,
    read_progress_interval: Optional[int] = 25,
    levenshtein_result_log_every: int = 10,
    best_score_epsilon: float = 1e-15,
    tie_score_abs_tol: float = 1e-12,
) -> Tuple[float, int]:
    """
    For each sequence in *sequences_file*, compute the sliding overlap Levenshtein
    score (normalized, see :func:`sliding_levenshtein`) with the read sliding along
    the reference. Write *result_file* as CSV: summary key/value rows, a blank line,
    then columns ``read_id``, ``sliding_score``, ``sequence`` for every read.

    sequence_format: ``"fastq"`` (default) or ``"blank_blocks"`` (paragraphs
    separated by empty lines).

    If *max_reads* is set, only that many FASTQ **reads** from the start are used
    when ``sequence_format`` is ``"fastq"`` (see :func:`iter_fastq_sequences`).

    Returns ``(minimum_sliding_score, number_of_reads_at_that_minimum)``.
    """
    logger.info(
        "starting run: sequences_file=%s reference_fasta=%s sequence_format=%s "
        "max_reads=%s result_file=%s levenshtein_backend=%s",
        sequences_file,
        reference_fasta,
        sequence_format,
        max_reads if max_reads is not None else "unlimited",
        result_file,
        "python-Levenshtein" if _HAS_C_LEV else "pure Python DP",
    )

    reference_sequence = load_reference_fasta(reference_fasta)
    if not reference_sequence:
        raise ValueError(f"No sequence found in reference FASTA: {reference_fasta}")

    if max_reads is not None and sequence_format != "fastq":
        raise ValueError("max_reads applies only when sequence_format is 'fastq'")

    if run_all_records and sequence_format == "fastq":
        max_reads = None
        sample_read_pool = None
        sample_read_count = None
        logger.info("run_all_records enabled: reading all FASTQ records with no sampling limits")

    if levenshtein_processes is not None and levenshtein_processes < 1:
        raise ValueError("levenshtein_processes must be a positive integer when set")

    if levenshtein_result_log_every <= 0:
        raise ValueError("levenshtein_result_log_every must be greater than zero")

    rng = random.Random(random_seed) if random_seed is not None else random

    if sequence_format == "fastq":
        pool_size = max_reads if max_reads is not None else sample_read_pool
        fastq_iterator = iter_fastq_sequences(
            sequences_file,
            max_reads=pool_size,
        )
        if sample_read_pool is not None and sample_read_count is not None:
            reads = list(fastq_iterator)
            if len(reads) > sample_read_count:
                sequence_iterator = iter(rng.sample(reads, sample_read_count))
                logger.info(
                    "sampled %d reads randomly from first %d FASTQ reads (seed=%s)",
                    sample_read_count,
                    len(reads),
                    random_seed if random_seed is not None else "system-random",
                )
            else:
                sequence_iterator = iter(reads)
        else:
            sequence_iterator = fastq_iterator
    elif sequence_format == "blank_blocks":
        sequence_iterator = iter_blank_block_sequences(sequences_file)
    else:
        raise ValueError("sequence_format must be 'fastq' or 'blank_blocks'")

    logger.info(
        "scoring reads (sliding overlap Levenshtein); progress every %s reads at INFO",
        read_progress_interval if read_progress_interval else "no intermediate",
    )

    reads_to_score: List[Tuple[str, str]] = []
    for read_id, read_sequence in sequence_iterator:
        if not read_sequence:
            logger.warning("skipping empty sequence for read_id=%s", read_id)
            continue
        reads_to_score.append((read_id, read_sequence))

    if not reads_to_score:
        raise ValueError(f"No sequences read from: {sequences_file}")

    use_processes = levenshtein_processes is not None and levenshtein_processes > 1
    if use_processes:
        logger.info("scoring with %d worker processes", levenshtein_processes)
    else:
        logger.info("scoring sequentially")

    best_sliding_score: Optional[float] = None
    reads_tied_for_best = 0
    per_read_results: List[Tuple[str, float, float, int, float, str]] = []
    reads_scored = 0

    if use_processes:
        with ProcessPoolExecutor(
            max_workers=levenshtein_processes,
            initializer=_initialize_scoring_worker,
            initargs=(reference_sequence,),
        ) as executor:
            scored_reads = executor.map(_score_read_in_worker, reads_to_score, chunksize=1)
            for (
                read_id,
                sliding_score,
                worst_score,
                record_length,
                matching_precentage,
                read_sequence,
            ) in scored_reads:
                per_read_results.append(
                    (
                        read_id,
                        sliding_score,
                        worst_score,
                        record_length,
                        matching_precentage,
                        read_sequence,
                    )
                )
                reads_scored += 1

                if reads_scored % levenshtein_result_log_every == 0:
                    logger.info(
                        "every-%d-reads levenshtein result: read #%d id=%s sliding_score=%.8g worst_score=%.8g length=%d matching_precentage=%.8g",
                        levenshtein_result_log_every,
                        reads_scored,
                        read_id,
                        sliding_score,
                        worst_score,
                        record_length,
                        matching_precentage,
                    )

                if best_sliding_score is None or sliding_score < best_sliding_score - best_score_epsilon:
                    best_sliding_score = sliding_score
                    reads_tied_for_best = 1
                    logger.debug(
                        "new best sliding score: %.8g (read_id=%s len=%d)",
                        best_sliding_score,
                        read_id,
                        len(read_sequence),
                    )
                elif math.isclose(
                    sliding_score,
                    best_sliding_score,
                    rel_tol=0.0,
                    abs_tol=tie_score_abs_tol,
                ):
                    reads_tied_for_best += 1

                logger.debug(
                    "read #%d id=%s bp=%d sliding_score=%.8g worst_score=%.8g matching_precentage=%.8g best_so_far=%.8g",
                    reads_scored,
                    read_id,
                    record_length,
                    sliding_score,
                    worst_score,
                    matching_precentage,
                    best_sliding_score if best_sliding_score is not None else float("nan"),
                )

                if read_progress_interval and reads_scored % read_progress_interval == 0:
                    logger.info(
                        "progress: %d reads scored (last read_id=%s score=%.6g best_so_far=%.6g)",
                        reads_scored,
                        read_id,
                        sliding_score,
                        best_sliding_score,
                    )
    else:
        for read_id, read_sequence in reads_to_score:
            sliding_score, worst_score = sliding_levenshtein(reference_sequence, read_sequence)
            record_length = len(read_sequence)
            matching_precentage = (
                sliding_score / record_length if record_length > 0 else float("inf")
            )
            per_read_results.append(
                (
                    read_id,
                    sliding_score,
                    worst_score,
                    record_length,
                    matching_precentage,
                    read_sequence,
                )
            )
            reads_scored += 1

            if reads_scored % levenshtein_result_log_every == 0:
                logger.info(
                    "every-%d-reads levenshtein result: read #%d id=%s sliding_score=%.8g worst_score=%.8g length=%d matching_precentage=%.8g",
                    levenshtein_result_log_every,
                    reads_scored,
                    read_id,
                    sliding_score,
                    worst_score,
                    record_length,
                    matching_precentage,
                )

            if best_sliding_score is None or sliding_score < best_sliding_score - best_score_epsilon:
                best_sliding_score = sliding_score
                reads_tied_for_best = 1
                logger.debug(
                    "new best sliding score: %.8g (read_id=%s len=%d)",
                    best_sliding_score,
                    read_id,
                    len(read_sequence),
                )
            elif math.isclose(
                sliding_score,
                best_sliding_score,
                rel_tol=0.0,
                abs_tol=tie_score_abs_tol,
            ):
                reads_tied_for_best += 1

            logger.debug(
                "read #%d id=%s bp=%d sliding_score=%.8g worst_score=%.8g matching_precentage=%.8g best_so_far=%.8g",
                reads_scored,
                read_id,
                record_length,
                sliding_score,
                worst_score,
                matching_precentage,
                best_sliding_score if best_sliding_score is not None else float("nan"),
            )

            if read_progress_interval and reads_scored % read_progress_interval == 0:
                logger.info(
                    "progress: %d reads scored (last read_id=%s score=%.6g best_so_far=%.6g)",
                    reads_scored,
                    read_id,
                    sliding_score,
                    best_sliding_score,
                )

    logger.info(
        "scoring finished: %d reads; best_sliding_score=%.8g (%d reads tied)",
        reads_scored,
        best_sliding_score,
        reads_tied_for_best,
    )

    output_path = Path(result_file)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    logger.info("writing results CSV: %s", output_path.resolve())
    with open(output_path, "w", encoding="utf-8", newline="") as csv_file:
        csv_writer = csv.writer(csv_file)
        csv_writer.writerow(["minimum_sliding_score", best_sliding_score])
        csv_writer.writerow(["reads_at_minimum", reads_tied_for_best])
        csv_writer.writerow(["total_reads", len(per_read_results)])
        csv_writer.writerow([])
        csv_writer.writerow(
            [
                "read_id",
                "sliding_score",
                "worst_score",
                "record_length",
                "matching_precentage",
                "sequence",
            ]
        )
        for (
            read_id,
            sliding_score,
            worst_score,
            record_length,
            matching_precentage,
            read_sequence,
        ) in per_read_results:
            csv_writer.writerow(
                [
                    read_id,
                    sliding_score,
                    worst_score,
                    record_length,
                    matching_precentage,
                    read_sequence,
                ]
            )

    post_process_results_csv(output_path)

    logger.info("wrote %d data rows to %s", len(per_read_results), output_path.resolve())
    return best_sliding_score, reads_tied_for_best


def run_sequencing(
    *,
    sequences_file: str | Path,
    reference_fasta: str | Path,
    result_file: str | Path,
    sequence_format: str,
    max_reads: Optional[int],
    sample_read_pool: Optional[int],
    sample_read_count: Optional[int],
    run_all_records: bool,
    random_seed: Optional[int],
    levenshtein_processes: Optional[int],
    log_file: str | Path,
    log_level: int,
    read_progress_interval: Optional[int],
    levenshtein_result_log_every: int,
    best_score_epsilon: float,
    tie_score_abs_tol: float,
) -> Tuple[float, int]:
    """Run sequencing pipeline using explicit parameters supplied by caller."""
    configure_logging(log_file, log_level)
    minimum_score, reads_at_minimum = find_best_matching_sequences(
        sequences_file,
        reference_fasta,
        result_file,
        sequence_format=sequence_format,
        max_reads=max_reads,
        sample_read_pool=sample_read_pool,
        sample_read_count=sample_read_count,
        run_all_records=run_all_records,
        random_seed=random_seed,
        levenshtein_processes=levenshtein_processes,
        read_progress_interval=read_progress_interval,
        levenshtein_result_log_every=levenshtein_result_log_every,
        best_score_epsilon=best_score_epsilon,
        tie_score_abs_tol=tie_score_abs_tol,
    )
    logger.info(
        "done: minimum_sliding_score=%s reads_at_minimum=%s output=%s",
        minimum_score,
        reads_at_minimum,
        Path(result_file).resolve(),
    )
    if not _HAS_C_LEV:
        logger.warning(
            "install python-Levenshtein for faster runs on large FASTQ files "
            "(pip install python-Levenshtein)",
        )
    return minimum_score, reads_at_minimum


def run_sequencing_from_env() -> Tuple[float, int]:
    """
    Run the sequencing pipeline using values loaded from the root .env file.

    Returns ``(minimum_sliding_score, reads_at_minimum)``.
    """
    return run_sequencing(
        sequences_file=SEQUENCES_FILE,
        reference_fasta=REFERENCE_FASTA,
        result_file=RESULT_FILE,
        sequence_format=SEQUENCE_FORMAT,
        max_reads=MAX_READS,
        sample_read_pool=SAMPLE_READ_POOL,
        sample_read_count=SAMPLE_READ_COUNT,
        run_all_records=RUN_ALL_RECORDS,
        random_seed=RANDOM_SEED,
        levenshtein_processes=LEVENSHTEIN_PROCESSES,
        log_file=LOG_FILE,
        log_level=LOG_LEVEL,
        read_progress_interval=READ_PROGRESS_INTERVAL,
        levenshtein_result_log_every=LEVENSHTEIN_RESULT_LOG_EVERY,
        best_score_epsilon=BEST_SCORE_EPSILON,
        tie_score_abs_tol=TIE_SCORE_ABS_TOL,
    )


def run_post_processing_from_env(result_file: Optional[str | Path] = None) -> Path:
    """Run only CSV post-processing using .env RESULT_FILE by default."""
    configure_logging(LOG_FILE, LOG_LEVEL)
    target_file = result_file if result_file is not None else RESULT_FILE
    return post_process_results_csv(target_file)


def main() -> None:
    """Backward-compatible CLI entrypoint for this module."""
    run_sequencing_from_env()


if __name__ == "__main__":
    main()
