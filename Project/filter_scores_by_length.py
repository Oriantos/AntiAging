"""Remove CSV rows where score is larger than the available data length."""

from __future__ import annotations

import argparse
import csv
import tempfile
from pathlib import Path
from typing import Iterable, Optional

ROOT = Path(__file__).resolve().parent.parent

# Preferred column names in priority order.
SCORE_COLUMNS = ("sliding_score", "k_mer_score", "score")
LENGTH_COLUMNS = ("record_length", "read_length", "data_length")
SEQUENCE_COLUMNS = ("sequence", "data")


def _is_blank_row(row: list[str]) -> bool:
    return not row or all(cell.strip() == "" for cell in row)


def _find_first_index(header: list[str], names: Iterable[str]) -> Optional[int]:
    for name in names:
        if name in header:
            return header.index(name)
    return None


def _parse_length(value: str) -> Optional[float]:
    text = value.strip()
    if not text:
        return None
    return float(text)


def _process_csv_file(input_path: Path, output_path: Path, dry_run: bool = False) -> tuple[int, int]:
    """Filter one CSV and return (kept_rows, removed_rows)."""
    kept_rows = 0
    removed_rows = 0

    # Keep summary/preamble rows before the data header unchanged.
    preamble: list[list[str]] = []
    in_data_section = False
    header: Optional[list[str]] = None
    score_index: Optional[int] = None
    length_index: Optional[int] = None
    sequence_index: Optional[int] = None

    if dry_run:
        writer = None
        tmp_path = None
    else:
        output_path.parent.mkdir(parents=True, exist_ok=True)
        tmp_file = tempfile.NamedTemporaryFile(
            mode="w",
            encoding="utf-8",
            newline="",
            delete=False,
            dir=str(output_path.parent),
            prefix=f".{output_path.stem}.",
            suffix=".tmp",
        )
        tmp_path = Path(tmp_file.name)
        tmp_file.close()

        out_handle = tmp_path.open("w", encoding="utf-8", newline="")
        writer = csv.writer(out_handle)

    try:
        with input_path.open("r", encoding="utf-8", newline="") as in_handle:
            reader = csv.reader(in_handle)

            for row in reader:
                if not in_data_section:
                    preamble.append(row)
                    if _is_blank_row(row):
                        in_data_section = True
                    continue

                if header is None:
                    header = row
                    score_index = _find_first_index(header, SCORE_COLUMNS)
                    length_index = _find_first_index(header, LENGTH_COLUMNS)
                    sequence_index = _find_first_index(header, SEQUENCE_COLUMNS)

                    if score_index is None:
                        raise ValueError(
                            f"Could not find score column in {input_path}. "
                            f"Expected one of: {', '.join(SCORE_COLUMNS)}"
                        )
                    if length_index is None and sequence_index is None:
                        raise ValueError(
                            f"Could not find a length source column in {input_path}. "
                            f"Expected one of: {', '.join(LENGTH_COLUMNS)} or sequence/data"
                        )

                    if writer is not None:
                        writer.writerows(preamble)
                        writer.writerow(header)
                    continue

                if _is_blank_row(row):
                    continue

                if len(row) < len(header):
                    row = row + [""] * (len(header) - len(row))

                score_value = float(row[score_index])

                if length_index is not None:
                    length_value = _parse_length(row[length_index])
                else:
                    length_value = None

                if length_value is None:
                    length_value = float(len(row[sequence_index])) if sequence_index is not None else None

                if length_value is None:
                    raise ValueError(f"Unable to determine length value for a row in {input_path}")

                if score_value > length_value:
                    removed_rows += 1
                    continue

                kept_rows += 1
                if writer is not None:
                    writer.writerow(row)

        if writer is not None:
            out_handle.close()
            tmp_path.replace(output_path)

    except Exception:
        if writer is not None:
            out_handle.close()
            tmp_path.unlink(missing_ok=True)
        raise

    return kept_rows, removed_rows


def _default_csv_targets() -> list[Path]:
    candidates = [
        ROOT / "output" / "levenshtein_best_matches.csv",
        ROOT / "output" / "k_mer_levenshtein_best_matches.csv",
    ]
    return [path for path in candidates if path.exists()]


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Delete rows from output CSV files where score is larger than data length"
        )
    )
    parser.add_argument(
        "--csv-files",
        nargs="+",
        help="One or more CSV files to process. Defaults to known output CSV files.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Analyze and print how many rows would be removed without modifying files.",
    )
    args = parser.parse_args()

    if args.csv_files:
        targets = [Path(p) for p in args.csv_files]
    else:
        targets = _default_csv_targets()

    if not targets:
        raise SystemExit("No target CSV files found.")

    for target in targets:
        input_path = target if target.is_absolute() else (ROOT / target)
        if not input_path.exists():
            print(f"Skipping missing file: {input_path}")
            continue

        kept, removed = _process_csv_file(
            input_path=input_path,
            output_path=input_path,
            dry_run=args.dry_run,
        )
        action = "Would remove" if args.dry_run else "Removed"
        print(f"{action} {removed} rows from {input_path} (kept {kept})")


if __name__ == "__main__":
    main()
