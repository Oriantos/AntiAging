"""Root entrypoint for DNA sequencing operations."""

from __future__ import annotations

import logging
import os
from multiprocessing import freeze_support
from pathlib import Path
from typing import Optional

from Project.sequencing import run_sequencing

_ROOT = Path(__file__).resolve().parent


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


def start_dna_operations() -> None:
    """Start DNA operations using root .env runtime configuration."""
    _load_dotenv(_ROOT / ".env")

    run_all_records = _env_bool("RUN_ALL_RECORDS", False)
    max_reads = _env_optional_int("MAX_READS", 10000)
    sample_read_pool = _env_optional_int("SAMPLE_READ_POOL", 10000)
    sample_read_count = _env_optional_int("SAMPLE_READ_COUNT", 100)

    if run_all_records:
        max_reads = None
        sample_read_pool = None
        sample_read_count = None

    run_sequencing(
        sequences_file=_env_path("SEQUENCES_FILE", Path("data") / "reads" / "SRR31234567_1.fastq"),
        reference_fasta=_env_path("REFERENCE_FASTA", Path("data") / "reference" / "reference.fasta"),
        result_file=_env_path("RESULT_FILE", Path("output") / "levenshtein_best_matches.csv"),
        sequence_format=_env_str("SEQUENCE_FORMAT", "fastq"),
        max_reads=max_reads,
        sample_read_pool=sample_read_pool,
        sample_read_count=sample_read_count,
        run_all_records=run_all_records,
        random_seed=_env_optional_int("RANDOM_SEED", 42),
        levenshtein_processes=_env_optional_int("LEVENSHTEIN_PROCESSES", 1),
        log_file=_env_path("LOG_FILE", Path("output") / "sequencing.log"),
        log_level=_env_log_level("LOG_LEVEL", logging.INFO),
        read_progress_interval=_env_optional_int("READ_PROGRESS_INTERVAL", 25),
        levenshtein_result_log_every=_env_int("LEVENSHTEIN_RESULT_LOG_EVERY", 10),
        best_score_epsilon=_env_float("BEST_SCORE_EPSILON", 1e-15),
        tie_score_abs_tol=_env_float("TIE_SCORE_ABS_TOL", 1e-12),
    )


def main() -> None:
    freeze_support()
    start_dna_operations()


if __name__ == "__main__":
    main()
