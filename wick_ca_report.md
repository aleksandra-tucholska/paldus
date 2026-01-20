# Report on `wick_ca` Subroutine

## 1. Overview of `wick_ca`

The `wick_ca` subroutine is a method of the `cas` (creator/annihilator string) class in `paldus_cas.py`. Its primary purpose is to apply **Wick's Theorem** to a string of creation and annihilation operators.

When `wick_ca` is called on a `cas` object, it performs the following steps:

1.  **Generate Contractions**: It calls `all_contractions()` to generate all possible pairings of operators.
2.  **Filter Contractions**: It calls `remove_null_contractions()` to keep only non-zero contractions (typically annihilator-creator pairs).
3.  **Process Terms**: It creates a list of resulting terms.
    *   **Term 0**: The **normal ordered** form of the original operator string (without any contractions). This involves moving all creation operators to the left of all annihilation operators, applying the necessary sign changes based on the number of swaps (fermionic statistics).
    *   **Subsequent Terms**: For each valid contraction, it creates a new term where the contracted operators are replaced by a Kronecker delta (and removed from the operator list), and the remaining operators are normal ordered.
4.  **Return**: It returns a list of `cas` objects representing the sum of the normal ordered term and all contracted terms.

In essence, `x.wick_ca()` computes the Wick expansion of the operator string `x`:
$$ x = \{x\} + \sum \{x \text{ with contractions}\} $$
where $\{ \dots \}$ denotes normal ordering.

## 2. Walkthrough: `KL = evaluate(h1cas, AAqp, AArs)`

This section analyzes the execution of the example provided:
```python
KL = evaluate(h1cas, AAqp, AArs)
res = []
for x in KL:
    result_list = x.wick_ca()
    res = res + result_list
```

### Definitions
From `templates.py`:
*   `h1cas`: Represents the one-electron Hamiltonian $\sum_{tu} h_{tu} t^\dagger u$.
*   `AAqp`: Represents the operator $q^\dagger p$. (Indices $q, p$ are fixed/external).
*   `AArs`: Represents the operator $r^\dagger s$. (Indices $r, s$ are fixed/external).

### Step 2.1: `M = evaluate(h1cas, AAqp)`

The `evaluate` function calculates the commutator $[h1cas, AAqp]$.
Internally, `commute_cas` computes this by forming two products:
1.  **L1** $= h1cas \times AAqp = \sum_{tu} h_{tu} t^\dagger u q^\dagger p$.
2.  **L2** $= AAqp \times h1cas = \sum_{tu} h_{tu} q^\dagger p t^\dagger u$.

It calls `wick_ca()` on both L1 and L2.

**For L1 ($t^\dagger u q^\dagger p$):**
*   **Normal Order**: $-t^\dagger q^\dagger u p$ (1 swap needed to move $q^\dagger$ past $u$).
*   **Contraction**: $u$ (ANI) with $q^\dagger$ (CRE) $\rightarrow \delta_{uq}$.
    *   Result: $\delta_{uq} t^\dagger p$.
*   **L1 Wick Result**: $\{-t^\dagger q^\dagger u p\} + \{\delta_{uq} t^\dagger p\}$.

**For L2 ($q^\dagger p t^\dagger u$):**
*   **Normal Order**: $-q^\dagger t^\dagger p u$ (1 swap needed to move $p$ past $t^\dagger$).
    *   Note: $-q^\dagger t^\dagger p u = -(-t^\dagger q^\dagger)(-u p) = -t^\dagger q^\dagger u p$. This matches the normal ordered term of L1.
*   **Contraction**: $p$ (ANI) with $t^\dagger$ (CRE) $\rightarrow \delta_{pt}$.
    *   Result: $\delta_{pt} q^\dagger u$.
*   **L2 Wick Result**: $\{-t^\dagger q^\dagger u p\} + \{\delta_{pt} q^\dagger u\}$.

**Commutator ($L1 - L2$):**
The normal ordered terms $\{-t^\dagger q^\dagger u p\}$ cancel out perfectly.
The remaining terms are:
$$ M = \sum_{tu} h_{tu} (\delta_{uq} t^\dagger p - \delta_{pt} q^\dagger u) $$
Simplifying deltas:
$$ M = \sum_t h_{tq} t^\dagger p - \sum_u h_{pu} q^\dagger u $$

### Step 2.2: `KL = evaluate(M, AArs)`

Now `evaluate` calculates $[M, AArs]$. This is linear, so we compute $[M_1, AArs]$ and $[M_2, AArs]$.

**Term 1**: $M_1 = \sum_t h_{tq} t^\dagger p$.
*   **L1** $= M_1 \times r^\dagger s = \sum_t h_{tq} t^\dagger p r^\dagger s$.
    *   Contraction: $p$ with $r^\dagger \rightarrow \delta_{pr}$.
    *   Result: $\{t^\dagger p r^\dagger s\} + \delta_{pr} \sum_t h_{tq} t^\dagger s$.
*   **L2** $= r^\dagger s \times M_1 = \sum_t h_{tq} r^\dagger s t^\dagger p$.
    *   Contraction: $s$ with $t^\dagger \rightarrow \delta_{st}$.
    *   Result: $\{r^\dagger s t^\dagger p\} + \delta_{st} \sum_t h_{tq} r^\dagger p$.
*   **Commutator**: The normal ordered terms cancel.
    *   Result: $\delta_{pr} \sum_t h_{tq} t^\dagger s - h_{sq} r^\dagger p$ (using $\delta_{st}$ to sum over t).

**Term 2**: $M_2 = - \sum_u h_{pu} q^\dagger u$.
*   **L1** $= M_2 \times r^\dagger s = - \sum_u h_{pu} q^\dagger u r^\dagger s$.
    *   Contraction: $u$ with $r^\dagger \rightarrow \delta_{ur}$.
    *   Result: $\{ \dots \} - \sum_u h_{pu} \delta_{ur} q^\dagger s = - h_{pr} q^\dagger s$.
*   **L2** $= r^\dagger s \times M_2 = - \sum_u h_{pu} r^\dagger s q^\dagger u$.
    *   Contraction: $s$ with $q^\dagger \rightarrow \delta_{sq}$.
    *   Result: $\{ \dots \} - \sum_u h_{pu} \delta_{sq} r^\dagger u$.
*   **Commutator ($L1 - L2$):**
    *   Result: $- h_{pr} q^\dagger s - (- \sum_u h_{pu} \delta_{sq} r^\dagger u) = - h_{pr} q^\dagger s + \delta_{sq} \sum_u h_{pu} r^\dagger u$.

### Final KL Content
The variable `KL` is a list containing the following 4 terms (representing the double commutator):
1.  $+\delta_{pr} \sum_t h_{tq} t^\dagger s$
2.  $-h_{sq} r^\dagger p$
3.  $-h_{pr} q^\dagger s$
4.  $+\delta_{sq} \sum_u h_{pu} r^\dagger u$

### Step 2.3: Loop and Final `wick_ca` Call

```python
res = [] 
for x in KL: 
    result_list = x.wick_ca() 
    res = res + result_list 
```

The code iterates over each term in `KL`. Each `x` is a `cas` object representing a one-body operator (e.g., $t^\dagger s$).
*   **x**: $t^\dagger s$ (Creation operator followed by Annihilation operator).
*   **x.wick_ca()**:
    *   Checks for contractions. Since the string is already normal ordered ($CRE \cdot ANI$), there are **no valid contractions** between an annihilator on the left and a creator on the right.
    *   Returns a list containing only the normal ordered term, which is `x` itself.
    *   `result_list` = `[x]`.

**Conclusion**: The loop effectively just copies the terms from `KL` into `res`, as they are already fully simplified and normal-ordered one-body operators.

The final `res` list contains exactly the same 4 terms derived above.
