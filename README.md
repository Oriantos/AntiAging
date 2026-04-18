# AntiAging Sequence Matching

This project compares sequencing reads to a reference using sliding-window Levenshtein distance.

## Organized structure

- `data/reads/` input FASTQ reads
- `data/reference/` reference FASTA
- `docs/` project notes and algorithm description
- `output/` generated CSV and logs
- `Project/sequencing.py` main script
- `.env` runtime configuration

## Configure

Edit `.env` to control:

- input/output file paths
- sampling and read limits
- run-all mode for processing every FASTQ record
- multiprocessing worker count for Levenshtein scoring
- logging intervals and level
- score tolerance settings

## Run

```powershell
python DNA_main.py
```
