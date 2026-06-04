# Sliding Group Levenshtein Distance

## Overview

The **Sliding Group Levenshtein Distance** is a variant of the classical Levenshtein (edit) distance algorithm, designed to find the optimal alignment between a **query string** and a **reference string** by sliding fixed-size groups of the query over the reference and computing a globally coherent edit distance score.

---

## Motivation

Standard Levenshtein distance compares two strings globally and symmetrically. In many real-world scenarios, however, one string serves as a structured **reference** (e.g., a long sequence, a captured stream, or a corpus), while the other is a **query** that is partitioned and searched within the reference at different offsets. The Sliding Group Levenshtein algorithm exploits this asymmetry to find the best positional alignment of the query relative to the reference.

---

## Definitions

| Term | Description |
|---|---|
| **Query string** `Q` | The structured string of length `n`, divided into fixed-size groups |
| **Reference string** `R` | The string of length `m` that the query groups are slid across |
| **Group size** `g` | Fixed size of each group; default is **15 characters** |
| **Group** `G_i` | The `i`-th consecutive substring of `Q` of length `g`: `Q[i·g : (i+1)·g]` |
| **Anchor position** | A position in `R` where a group `G_i` has a **perfect match** (zero edit distance). A group may have multiple anchor positions. |
| **Anchor set** `A_i` | The set of all positions in `R` where group `G_i` matches perfectly: `A_i = { p : levenshtein(G_i, R[p : p + |G_i|]) = 0 }` |

---

## Algorithm

### Step 1 — Partition the Query

Divide the query string `Q` into consecutive, non-overlapping groups of size `g = 15`:

```
G_0 = Q[0:15]
G_1 = Q[15:30]
G_2 = Q[30:45]
...
G_k = Q[k·g : min((k+1)·g, n)]
```

The last group may be shorter than `g` if `n` is not a multiple of 15.

### Step 2 — Slide Each Group Over the Reference and Collect Perfect Matches

For each group `G_i`, slide it across the entire reference string `R` and record every position where it matches **perfectly** (zero edit distance):

```
A_i = {}
for each position p in R (0 ≤ p ≤ m - |G_i|):
    if G_i == R[p : p + |G_i|]:   // exact match, i.e. levenshtein = 0
        A_i.add(p)
```

`A_i` is the **anchor set** for group `G_i` — the collection of all positions in `R` where `G_i` appears verbatim. If `A_i` is empty for **all** groups, the algorithm terminates immediately with a score of `m` (the length of `R`). If `A_i` is empty for only some groups, those groups are excluded from anchoring and aligned by the global step around the anchored groups.

### Step 3 — Enumerate All Anchor Combinations

For each group that has at least one anchor, every position in its anchor set is a candidate. If **no group** has any anchor at all, the algorithm stops and returns `m` as the final score. Otherwise, the algorithm enumerates **all combinations** of one anchor per anchored group:

```
anchor_combinations = cartesian_product(A_0, A_1, ..., A_k)

for each combination C = (p_0, p_1, ..., p_k) in anchor_combinations:
    compute global alignment score using C
    keep track of the combination with the lowest score
```

This exhaustive search ensures no valid alignment is overlooked.

### Step 4 — Compute the Global Alignment Score per Combination

For each anchor combination `C = (p_0, p_1, ..., p_k)`, reconstruct the alignment of the full query `Q` against `R` such that each group `G_i` is positioned at `p_i` in `R`, with all other groups maintaining their relative order around the anchored positions.

The alignment score for combination `C` is:

```
score(C) = levenshtein(Q_aligned(C), R)
```

### Step 5 — Select the Best Alignment

The final distance is the minimum score across all anchor combinations:

```
best_score = min( score(C) for C in anchor_combinations )
```

Where `Q_aligned(C)` is the version of `Q` re-positioned according to combination `C`, and all non-anchored groups fill in around the anchored ones in their original relative order.

---

## Scoring

The final distance score reflects:

- How well the full query aligns with the reference when each group is pinned to one of its perfect-match positions
- The cost of fitting non-anchored groups around the anchored ones using standard Levenshtein operations (insertion, deletion, substitution), each with a cost of 1
- The best possible alignment across **all valid anchor combinations**

A **score of 0** means the query, under some anchor combination, is identical to the reference. Higher scores indicate greater dissimilarity. If no group has any perfect match in `R`, the score is set to `m` — the length of the reference string — representing a maximally penalized, unaligned comparison.

---

## Complexity

| Phase | Complexity |
|---|---|
| Partitioning `Q` | `O(n)` |
| Sliding each group over `R` (exact match scan) | `O(k · m · g)` |
| Enumerating anchor combinations | `O(∏ |A_i|)` — product of anchor set sizes |
| Computing global alignment per combination | `O(n · m)` |
| **Total** | `O(k · m · g + n · m · ∏ |A_i|)` |

Where `k = ⌈n/g⌉` is the number of groups and `|A_i|` is the number of perfect-match positions found for group `G_i`. In the worst case (e.g., highly repetitive references), the number of combinations can be large; in practice it is bounded by the frequency of each group pattern in `R`.

---

## Example

Given:
- Query: `"ABCDEFGHIJKLMNO_PQRSTUVWXYZ"` (n = 27)
- Reference: `"XYZPQRSTUVWXYZABCDEFGHIJKLMNOABCDEFGHIJKLMNO"` (m = 44)
- Group size: `g = 15`

**Groups:**
- `G_0 = "ABCDEFGHIJKLMNO"`
- `G_1 = "_PQRSTUVWXYZ"` (remainder)

**Perfect-match anchor sets:**
- `A_0 = { 14, 29 }` → `G_0` appears verbatim at positions 14 and 29 in `R`
- `A_1 = {}` → `G_1` has no perfect match in `R`; excluded from anchoring

**Anchor combinations (only `G_0` anchors):**
- Combination 1: `G_0` anchored at position 14
- Combination 2: `G_0` anchored at position 29

**Global scores:**
- `score(C1) = levenshtein(Q_aligned(14), R)`
- `score(C2) = levenshtein(Q_aligned(29), R)`

**Final score:** `min(score(C1), score(C2))`

---

## Use Cases

- **Sequence matching** in bioinformatics where a known motif is searched within a longer genome
- **Log pattern matching** where a template log is compared against a noisy log stream
- **Network flow fingerprinting** where a reference traffic pattern is matched against captured flows
- **Approximate substring search** with structured reference templates

---

## Relation to Classical Levenshtein

| Property | Classical Levenshtein | Sliding Group Levenshtein |
|---|---|---|
| Symmetry | Symmetric | Asymmetric (query vs. reference) |
| Alignment | Global | Local per group, global score |
| Query structure | None | Partitioned into groups of size `g` |
| Best for | Full string comparison | Pattern matching with known structure |