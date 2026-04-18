### Levenshtein Distance (Recursive Definition)

$$
\text{lev}(a, b) =
\begin{cases}
|a| & \text{if } |b| = 0, \\
|b| & \text{if } |a| = 0, \\
\text{lev}(\text{tail}(a), \text{tail}(b)) & \text{if } \text{head}(a) = \text{head}(b), \\
1 + \min \begin{cases}
\text{lev}(\text{tail}(a), b), \\
\text{lev}(a, \text{tail}(b)), \\
\text{lev}(\text{tail}(a), \text{tail}(b))
\end{cases} & \text{otherwise}.
\end{cases}
$$

---

### Pseudocode Version

```pseudo
function lev(a, b):
    if length(b) == 0:
        return length(a)

    if length(a) == 0:
        return length(b)

    if head(a) == head(b):
        return lev(tail(a), tail(b))

    return 1 + min(
        lev(tail(a), b),        // deletion
        lev(a, tail(b)),        // insertion
        lev(tail(a), tail(b))   // substitution
    )
```

---

### Notes (Quick Intuition)

- **Deletion**: remove a character from `a`
- **Insertion**: add a character to `a`
- **Substitution**: replace a character in `a`
- The algorithm finds the *minimum number of edits* needed to transform one string into another.

