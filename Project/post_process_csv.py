"""Post-process sequencing CSV output with derived columns."""

from __future__ import annotations

import argparse
import csv
import os
from pathlib import Path
from typing import Optional

_ROOT = Path(__file__).resolve().parent.parent


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


def post_process_results_csv(result_file: str | Path) -> Path:
    """Add read_length and score_percentage columns to the sequencing CSV."""
    output_path = Path(result_file)

    with open(output_path, encoding="utf-8", newline="") as csv_file:
        rows = list(csv.reader(csv_file))

    if not rows:
        raise ValueError(f"CSV is empty: {output_path}")

    blank_row_index: Optional[int] = None
    for index, row in enumerate(rows):
        if not row or all(cell.strip() == "" for cell in row):
            blank_row_index = index
            break

    if blank_row_index is None or blank_row_index + 1 >= len(rows):
        raise ValueError(
            f"CSV does not contain expected summary/data split: {output_path}",
        )

    header_index = blank_row_index + 1
    header = rows[header_index]
    data_rows = rows[header_index + 1 :]

    has_read_length = "read_length" in header
    has_score_percentage = "score_percentage" in header
    if has_read_length and has_score_percentage:
        return output_path

    if "sequence" not in header:
        raise ValueError(f"CSV is missing required 'sequence' column: {output_path}")

    sequence_index = header.index("sequence")
    sliding_score_index = header.index("sliding_score") if "sliding_score" in header else None
    record_length_index = header.index("record_length") if "record_length" in header else None
    matching_precentage_index = (
        header.index("matching_precentage") if "matching_precentage" in header else None
    )

    if sliding_score_index is None and not has_score_percentage:
        raise ValueError(
            f"CSV is missing required 'sliding_score' column for score_percentage: {output_path}",
        )

    new_header = list(header)
    if not has_read_length:
        new_header.append("read_length")
    if not has_score_percentage:
        new_header.append("score_percentage")

    new_rows = rows[:header_index]
    new_rows.append(new_header)

    for original_row in data_rows:
        if not original_row or all(cell.strip() == "" for cell in original_row):
            continue

        row = list(original_row)
        if len(row) < len(header):
            row.extend([""] * (len(header) - len(row)))

        sequence_value = row[sequence_index]
        computed_read_length = len(sequence_value)

        if record_length_index is not None:
            source_record_length = row[record_length_index].strip()
            read_length_value = source_record_length if source_record_length else str(computed_read_length)
        else:
            read_length_value = str(computed_read_length)

        if matching_precentage_index is not None:
            source_matching = row[matching_precentage_index].strip()
            if source_matching:
                score_percentage_value = source_matching
            else:
                sliding_score = float(row[sliding_score_index]) if sliding_score_index is not None else float("inf")
                score_percentage_value = (
                    str(sliding_score / computed_read_length)
                    if computed_read_length > 0
                    else str(float("inf"))
                )
        else:
            sliding_score = float(row[sliding_score_index]) if sliding_score_index is not None else float("inf")
            score_percentage_value = (
                str(sliding_score / computed_read_length)
                if computed_read_length > 0
                else str(float("inf"))
            )

        if not has_read_length:
            row.append(read_length_value)
        if not has_score_percentage:
            row.append(score_percentage_value)

        new_rows.append(row)

    with open(output_path, "w", encoding="utf-8", newline="") as csv_file:
        csv.writer(csv_file).writerows(new_rows)

    return output_path


def main() -> None:
    _load_dotenv(_ROOT / ".env")
    default_csv = _env_path("RESULT_FILE", Path("output") / "levenshtein_best_matches.csv")

    parser = argparse.ArgumentParser(
        description="Post-process sequencing CSV and append derived columns",
    )
    parser.add_argument(
        "--csv-file",
        default=str(default_csv),
        help="CSV file to post-process (defaults to RESULT_FILE from .env)",
    )
    args = parser.parse_args()

    result_path = post_process_results_csv(args.csv_file)
    print(f"Post-processing complete: {result_path}")


if __name__ == "__main__":
    main()
