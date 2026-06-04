"""
Sliding Group Levenshtein Distance implementation.
See docs/k_mer_levinshtein_distance.md for the full algorithm description.

Partitions the query into fixed-size k-mer groups, finds exact-match anchor positions
for each group in the reference, then enumerates anchor combinations (cartesian product)
to find the globally optimal alignment score.

New .env variables:
  K_MER_GROUP_SIZE              Fixed group size (default 15)
  K_MER_RESULT_FILE             Output CSV path
  K_MER_MAX_ANCHOR_COMBINATIONS Cartesian-product cap (default 10000)

Reuses SEQUENCES_FILE, REFERENCE_FASTA, MAX_READS, LOG_FILE, and other shared keys.
"""

from __future__ import annotations

import csv
import itertools
import logging
import math
import os
import random
from pathlib import Path
from typing import List, Optional, Tuple

from post_process_csv import post_process_results_csv
from sequencing import (
    _HAS_C_LEV,
    configure_logging,
    iter_blank_block_sequences,
    iter_fastq_sequences,
    levenshtein_distance,
    load_reference_fasta,
)

_ROOT = Path(__file__).resolve().parent.parent
logger = logging.getLogger("Project.k_mer_levenshtein")


# ---------------------------------------------------------------------------
# .env helpers (identical pattern to sequencing.py)
# ---------------------------------------------------------------------------

def _load_dotenv(dotenv_path: Path) -> None:
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


def _env_path(name: str, default: Path) -> Path:
    raw = os.getenv(name)
    candidate = Path(raw) if raw else default
    return candidate if candidate.is_absolute() else (_ROOT / candidate)


def _env_int(name: str, default: int) -> int:
    raw = os.getenv(name)
    return int(raw) if raw is not None else default


def _env_str(name: str, default: str) -> str:
    return os.getenv(name, default)


def _env_float(name: str, default: float) -> float:
    raw = os.getenv(name)
    return float(raw) if raw is not None else default


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


_DOTENV_PATH = _ROOT / ".env"
_load_dotenv(_DOTENV_PATH)

# ---------------------------------------------------------------------------
# Configuration (shared keys re-read here so the module is self-contained)
# ---------------------------------------------------------------------------

SEQUENCES_FILE: Path = _env_path(
    "SEQUENCES_FILE", Path("data") / "reads" / "RNA" / "young_small.fastq"
)
REFERENCE_FASTA: Path = _env_path(
    "REFERENCE_FASTA", Path("data") / "reference" / "DNA" / "reference.fasta"
)
SEQUENCE_FORMAT: str = _env_str("SEQUENCE_FORMAT", "fastq")
MAX_READS: Optional[int] = _env_optional_int("MAX_READS", 10000)
SAMPLE_READ_POOL: Optional[int] = _env_optional_int("SAMPLE_READ_POOL", 10000)
SAMPLE_READ_COUNT: Optional[int] = _env_optional_int("SAMPLE_READ_COUNT", 100)
RANDOM_SEED: Optional[int] = _env_optional_int("RANDOM_SEED", 42)
RUN_ALL_RECORDS: bool = _env_bool("RUN_ALL_RECORDS", False)
LOG_FILE: Path = _env_path("LOG_FILE", Path("output") / "sequencing.log")
LOG_LEVEL: int = _env_log_level("LOG_LEVEL", logging.INFO)
READ_PROGRESS_INTERVAL: Optional[int] = _env_optional_int("READ_PROGRESS_INTERVAL", 25)
LEVENSHTEIN_RESULT_LOG_EVERY: int = _env_int("LEVENSHTEIN_RESULT_LOG_EVERY", 10)
BEST_SCORE_EPSILON: float = _env_float("BEST_SCORE_EPSILON", 1e-15)
TIE_SCORE_ABS_TOL: float = _env_float("TIE_SCORE_ABS_TOL", 1e-12)

# K-mer-specific
K_MER_GROUP_SIZE: int = _env_int("K_MER_GROUP_SIZE", 15)
K_MER_RESULT_FILE: Path = _env_path(
    "K_MER_RESULT_FILE",
    Path("output") / "k_mer_levenshtein_best_matches.csv",
)
K_MER_MAX_ANCHOR_COMBINATIONS: int = _env_int("K_MER_MAX_ANCHOR_COMBINATIONS", 10_000)

if K_MER_GROUP_SIZE < 1:
    raise ValueError("K_MER_GROUP_SIZE must be at least 1")


# ---------------------------------------------------------------------------
# Core algorithm
# ---------------------------------------------------------------------------

def partition_query(query: str, group_size: int) -> List[str]:
    """Step 1: split query into consecutive non-overlapping groups of size group_size."""
    return [query[i : i + group_size] for i in range(0, len(query), group_size)]


def find_anchor_set(group: str, reference: str) -> List[int]:
    """Step 2: all positions in reference where group appears as an exact match."""
    positions: List[int] = []
    g = len(group)
    if g == 0:
        return positions
    start = 0
    while True:
        pos = reference.find(group, start)
        if pos == -1:
            break
        positions.append(pos)
        start = pos + 1  # allow overlapping occurrences
    return positions


def _score_combination(
    query: str,
    reference: str,
    groups: List[str],
    anchor_positions: List[Optional[int]],
) -> int:
    """
    Step 4: global alignment score for one anchor combination.

    Anchored groups contribute 0 (exact match). Each segment of the query between
    (and around) anchored groups is scored against the corresponding reference
    window using levenshtein_distance. The reference window is determined by the
    first anchored group's implied query-start position in the reference.
    """
    m = len(reference)
    anchored_indices = [i for i in range(len(groups)) if anchor_positions[i] is not None]

    if not anchored_indices:
        return m

    # Cumulative start positions of each group within the query string.
    group_starts: List[int] = []
    acc = 0
    for g in groups:
        group_starts.append(acc)
        acc += len(g)

    first_gi = anchored_indices[0]
    first_pos = anchor_positions[first_gi]  # type: ignore[index]
    start_in_r = first_pos - group_starts[first_gi]  # implied start of Q in R (may be < 0)

    total_score = 0
    prev_q_end = 0
    prev_r_end = start_in_r  # tracks R cursor; may be negative

    for gi in anchored_indices:
        pos: int = anchor_positions[gi]  # type: ignore[assignment]
        q_start = group_starts[gi]

        # Score the query segment before this anchor against the corresponding R window.
        q_seg = query[prev_q_end:q_start]
        r_seg_start = max(0, prev_r_end)
        r_seg_end = max(0, min(m, pos))
        r_seg = reference[r_seg_start:r_seg_end]
        total_score += levenshtein_distance(q_seg, r_seg)

        # Anchored group: exact match in R → 0 cost.
        prev_q_end = q_start + len(groups[gi])
        prev_r_end = pos + len(groups[gi])

    # Score the trailing query segment.
    q_seg = query[prev_q_end:]
    if q_seg:
        r_seg_start = max(0, prev_r_end)
        r_seg_end = min(m, r_seg_start + len(q_seg))
        r_seg = reference[r_seg_start:r_seg_end]
        total_score += levenshtein_distance(q_seg, r_seg)

    return total_score


def sliding_group_levenshtein(
    query: str,
    reference: str,
    group_size: int = K_MER_GROUP_SIZE,
    max_combinations: int = K_MER_MAX_ANCHOR_COMBINATIONS,
) -> int:
    """
    Sliding Group Levenshtein distance between query and reference.

    Returns the minimum global alignment score across all valid anchor combinations.
    Returns len(reference) when no group has any perfect match in the reference,
    representing a maximally penalised unaligned comparison (as specified in the doc).

    max_combinations caps the cartesian product enumeration to avoid exponential blow-up
    on highly repetitive references; the first max_combinations combos are evaluated.
    """
    n = len(query)
    m = len(reference)
    if n == 0:
        return m
    if m == 0:
        return n

    # Step 1 — partition
    groups = partition_query(query, group_size)
    k = len(groups)

    # Step 2 — anchor sets via fast exact-match scan
    anchor_sets = [find_anchor_set(g, reference) for g in groups]

    logger.debug(
        "anchor sets: %d groups, anchor counts=%s",
        k,
        [len(a) for a in anchor_sets],
    )

    # Step 3 — check whether any group has anchors
    anchored_group_indices = [i for i, a in enumerate(anchor_sets) if a]
    if not anchored_group_indices:
        logger.debug("no anchors found; returning max penalty m=%d", m)
        return m

    anchored_sets = [anchor_sets[i] for i in anchored_group_indices]
    num_combinations = math.prod(len(s) for s in anchored_sets)
    logger.debug(
        "%d anchor combinations (anchored group indices: %s)",
        num_combinations,
        anchored_group_indices,
    )

    combo_iter = (
        itertools.islice(itertools.product(*anchored_sets), max_combinations)
        if num_combinations > max_combinations
        else itertools.product(*anchored_sets)
    )
    if num_combinations > max_combinations:
        logger.warning(
            "anchor combinations (%d) exceed K_MER_MAX_ANCHOR_COMBINATIONS (%d); "
            "evaluating first %d only",
            num_combinations,
            max_combinations,
            max_combinations,
        )

    best_score = m  # pessimistic upper bound

    for combo in combo_iter:
        # Map each anchored group to its candidate position in R.
        anchor_positions: List[Optional[int]] = [None] * k
        for j, gi in enumerate(anchored_group_indices):
            anchor_positions[gi] = combo[j]

        # Validate: anchor positions must be non-overlapping and ordered in R.
        valid = True
        prev_r_end = -1
        for gi in anchored_group_indices:
            pos = anchor_positions[gi]
            if pos < prev_r_end:  # type: ignore[operator]
                valid = False
                break
            prev_r_end = pos + len(groups[gi])  # type: ignore[operator]
        if not valid:
            continue

        # Step 4 — score this combination
        score = _score_combination(query, reference, groups, anchor_positions)
        if score < best_score:
            best_score = score
            if best_score == 0:
                break  # perfect alignment found; no combination can do better

    return best_score


# ---------------------------------------------------------------------------
# Pipeline (mirrors sequencing.py structure)
# ---------------------------------------------------------------------------

def find_best_matching_sequences_k_mer(
    sequences_file: str | Path,
    reference_fasta: str | Path,
    result_file: str | Path,
    *,
    sequence_format: str = SEQUENCE_FORMAT,
    max_reads: Optional[int] = MAX_READS,
    sample_read_pool: Optional[int] = SAMPLE_READ_POOL,
    sample_read_count: Optional[int] = SAMPLE_READ_COUNT,
    run_all_records: bool = RUN_ALL_RECORDS,
    random_seed: Optional[int] = RANDOM_SEED,
    read_progress_interval: Optional[int] = READ_PROGRESS_INTERVAL,
    levenshtein_result_log_every: int = LEVENSHTEIN_RESULT_LOG_EVERY,
    best_score_epsilon: float = BEST_SCORE_EPSILON,
    tie_score_abs_tol: float = TIE_SCORE_ABS_TOL,
    group_size: int = K_MER_GROUP_SIZE,
    max_anchor_combinations: int = K_MER_MAX_ANCHOR_COMBINATIONS,
) -> Tuple[int, int]:
    """
    Score all reads against the reference using Sliding Group Levenshtein.
    Writes a CSV to result_file (same layout as sequencing.py) and runs
    post_process_results_csv.  Returns (minimum_k_mer_score, reads_at_minimum).
    """
    logger.info(
        "k-mer run: sequences_file=%s reference_fasta=%s format=%s "
        "max_reads=%s group_size=%d max_combinations=%d result_file=%s backend=%s",
        sequences_file,
        reference_fasta,
        sequence_format,
        max_reads if max_reads is not None else "unlimited",
        group_size,
        max_anchor_combinations,
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
        logger.info("run_all_records=true: reading all FASTQ records without sampling")

    rng = random.Random(random_seed) if random_seed is not None else random

    if sequence_format == "fastq":
        pool_size = max_reads if max_reads is not None else sample_read_pool
        fastq_iter = iter_fastq_sequences(sequences_file, max_reads=pool_size)
        if sample_read_pool is not None and sample_read_count is not None:
            reads_pool = list(fastq_iter)
            if len(reads_pool) > sample_read_count:
                sequence_iterator = iter(rng.sample(reads_pool, sample_read_count))
                logger.info(
                    "sampled %d reads from first %d FASTQ reads (seed=%s)",
                    sample_read_count,
                    len(reads_pool),
                    random_seed if random_seed is not None else "system-random",
                )
            else:
                sequence_iterator = iter(reads_pool)
        else:
            sequence_iterator = fastq_iter
    elif sequence_format == "blank_blocks":
        sequence_iterator = iter_blank_block_sequences(sequences_file)
    else:
        raise ValueError("sequence_format must be 'fastq' or 'blank_blocks'")

    reads_to_score = [
        (rid, seq)
        for rid, seq in sequence_iterator
        if seq or logger.warning("skipping empty sequence for read_id=%s", rid) or False  # type: ignore[func-returns-value]
    ]
    if not reads_to_score:
        raise ValueError(f"No sequences read from: {sequences_file}")

    logger.info(
        "scoring %d reads with Sliding Group Levenshtein (group_size=%d)",
        len(reads_to_score),
        group_size,
    )

    best_score: Optional[int] = None
    reads_tied_for_best = 0
    per_read_results: List[Tuple[str, int, int, float, str]] = []
    reads_scored = 0

    for read_id, read_seq in reads_to_score:
        k_mer_score = sliding_group_levenshtein(
            read_seq, reference_sequence, group_size, max_anchor_combinations
        )
        record_length = len(read_seq)
        matching_percentage = (
            k_mer_score / record_length if record_length > 0 else float("inf")
        )

        per_read_results.append(
            (read_id, k_mer_score, record_length, matching_percentage, read_seq)
        )
        reads_scored += 1

        if reads_scored % levenshtein_result_log_every == 0:
            logger.info(
                "every-%d-reads k-mer result: read #%d id=%s k_mer_score=%d "
                "length=%d matching_pct=%.8g",
                levenshtein_result_log_every,
                reads_scored,
                read_id,
                k_mer_score,
                record_length,
                matching_percentage,
            )

        if best_score is None or k_mer_score < best_score - best_score_epsilon:
            best_score = k_mer_score
            reads_tied_for_best = 1
            logger.debug(
                "new best k-mer score: %d (read_id=%s len=%d)",
                best_score,
                read_id,
                record_length,
            )
        elif math.isclose(k_mer_score, best_score, rel_tol=0.0, abs_tol=tie_score_abs_tol):
            reads_tied_for_best += 1

        logger.debug(
            "read #%d id=%s bp=%d k_mer_score=%d matching_pct=%.8g best_so_far=%s",
            reads_scored,
            read_id,
            record_length,
            k_mer_score,
            matching_percentage,
            best_score,
        )

        if read_progress_interval and reads_scored % read_progress_interval == 0:
            logger.info(
                "progress: %d reads scored (last id=%s score=%d best_so_far=%d)",
                reads_scored,
                read_id,
                k_mer_score,
                best_score,
            )

    logger.info(
        "scoring finished: %d reads; best_k_mer_score=%d (%d reads tied)",
        reads_scored,
        best_score,
        reads_tied_for_best,
    )

    output_path = Path(result_file)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    logger.info("writing k-mer results CSV: %s", output_path.resolve())

    with open(output_path, "w", encoding="utf-8", newline="") as csv_file:
        writer = csv.writer(csv_file)
        writer.writerow(["minimum_k_mer_score", best_score])
        writer.writerow(["reads_at_minimum", reads_tied_for_best])
        writer.writerow(["total_reads", len(per_read_results)])
        writer.writerow(["group_size", group_size])
        writer.writerow([])
        writer.writerow(
            ["read_id", "k_mer_score", "record_length", "matching_percentage", "sequence"]
        )
        for rid, score, rlen, mpct, seq in per_read_results:
            writer.writerow([rid, score, rlen, mpct, seq])

    post_process_results_csv(output_path)

    logger.info("wrote %d data rows to %s", len(per_read_results), output_path.resolve())
    return best_score, reads_tied_for_best  # type: ignore[return-value]


def run_k_mer_levenshtein(
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
    log_file: str | Path,
    log_level: int,
    read_progress_interval: Optional[int],
    levenshtein_result_log_every: int,
    best_score_epsilon: float,
    tie_score_abs_tol: float,
    group_size: int,
    max_anchor_combinations: int,
) -> Tuple[int, int]:
    """Run the k-mer pipeline using explicit parameters supplied by caller."""
    configure_logging(log_file, log_level)
    minimum_score, reads_at_minimum = find_best_matching_sequences_k_mer(
        sequences_file,
        reference_fasta,
        result_file,
        sequence_format=sequence_format,
        max_reads=max_reads,
        sample_read_pool=sample_read_pool,
        sample_read_count=sample_read_count,
        run_all_records=run_all_records,
        random_seed=random_seed,
        read_progress_interval=read_progress_interval,
        levenshtein_result_log_every=levenshtein_result_log_every,
        best_score_epsilon=best_score_epsilon,
        tie_score_abs_tol=tie_score_abs_tol,
        group_size=group_size,
        max_anchor_combinations=max_anchor_combinations,
    )
    logger.info(
        "done: minimum_k_mer_score=%s reads_at_minimum=%s output=%s",
        minimum_score,
        reads_at_minimum,
        Path(result_file).resolve(),
    )
    return minimum_score, reads_at_minimum


def run_k_mer_from_env() -> Tuple[int, int]:
    """Run the k-mer pipeline using values loaded from the root .env file."""
    return run_k_mer_levenshtein(
        sequences_file=SEQUENCES_FILE,
        reference_fasta=REFERENCE_FASTA,
        result_file=K_MER_RESULT_FILE,
        sequence_format=SEQUENCE_FORMAT,
        max_reads=MAX_READS,
        sample_read_pool=SAMPLE_READ_POOL,
        sample_read_count=SAMPLE_READ_COUNT,
        run_all_records=RUN_ALL_RECORDS,
        random_seed=RANDOM_SEED,
        log_file=LOG_FILE,
        log_level=LOG_LEVEL,
        read_progress_interval=READ_PROGRESS_INTERVAL,
        levenshtein_result_log_every=LEVENSHTEIN_RESULT_LOG_EVERY,
        best_score_epsilon=BEST_SCORE_EPSILON,
        tie_score_abs_tol=TIE_SCORE_ABS_TOL,
        group_size=K_MER_GROUP_SIZE,
        max_anchor_combinations=K_MER_MAX_ANCHOR_COMBINATIONS,
    )


def main() -> None:
    """CLI entrypoint — runs k-mer pipeline from .env configuration."""
    run_k_mer_from_env()


if __name__ == "__main__":
    main()
