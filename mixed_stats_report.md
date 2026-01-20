# Report on Mixed Statistics Implementation in `wick_ca`

## Overview of Changes

The `wick_ca` subroutine and associated helper methods in `paldus_cas.py` have been modified to support mixed fermion-boson systems. The code now respects the statistical properties of operators defined in the `operator_stat` attribute of the `cas` class.

### Key Modifications

1.  **`cas` Class Initialization**:
    - The `cas` class now includes `operator_stat`, a list parallel to `operator_idx` and `operator_type`, storing the particle type (`ferm` or `boso`).

2.  **`contraction` Method**:
    - Updated to account for particle statistics when determining the sign of the contraction term.
    - Previously, the method calculated `num` (the number of operators crossed by the contracted pair) and applied a sign change if `num` was odd. This implicitly assumed all operators were fermions.
    - **New Logic**: The method now iterates through the operators between the contracted pair. It decrements the count `num` if the crossed operator is a boson (`boso`). This ensures that only crossings of fermions contribute to the sign change (fermionic parity).

3.  **`normal_order` Method**:
    - This method reorders operators (creators to the left, annihilators to the right) and sorts them by index. It calculates the sign change based on the number of swaps required.
    - **New Logic**:
        - It extracts the statistics (`operator_stat`) for the creators and annihilators.
        - It uses a new helper function `swap_count_stats` instead of the generic `swap_count`.
        - It uses an updated `count_anlator_pass` which is now aware of statistics.

4.  **`count_anlator_pass` Function**:
    - Previously, this counted how many annihilators each creator passed to reach the left side, changing sign for every pass (assuming fermions).
    - **New Logic**: It accepts `crean_stat` (list of statistics). It iterates through the operators and increments `pass_count` **only if both** the moving creator and the passed annihilator are fermions. If either is a boson, they commute without a sign change.

5.  **`swap_count_stats` Function (New)**:
    - A new function replacing `swap_count` for mixed statistics.
    - It performs a bubble sort simulation on the list of operators (zipped with their statistics).
    - During the sort, whenever a swap occurs, it checks the statistics of the two swapped operators.
    - A sign change (incrementing `swaps`) is recorded **only if both** swapped operators are fermions. Swaps involving bosons do not contribute to the sign count.

6.  **`remove_null_contractions` Method**:
    - Updated to enforce that contractions can only occur between operators of the same statistic.
    - It checks if `operator_stat[i] == operator_stat[j]`. If they differ (e.g., trying to contract a fermion with a boson), the contraction is treated as zero and removed.

## Integration Points for Fermionic Rules

The following specific locations in the integration/Wick subroutines apply fermionic rules and have been updated:

1.  **`paldus_cas.py: contraction`**:
    - The loop checking `range(list_of_c[i][0]+1, list_of_c[i][1])` now checks `self.operator_stat[j] == boso` to skip sign changes for bosons.

2.  **`paldus_cas.py: normal_order`**:
    - Calls `swap_count_stats` (new) instead of `swap_count` to handle reordering of creators and annihilators amongst themselves.
    - Calls `count_anlator_pass` (updated) to handle moving creators past annihilators.

3.  **`paldus_cas.py: count_anlator_pass`**:
    - The condition `if not crean_stat or (crean_stat[i] == ferm and crean_stat[j] == ferm):` ensures sign changes only for fermion-fermion crossings.

## Verification

A test script `test_mixed_stats.py` was created to verify the changes:
- **Pure Fermions**: Behaves as before (anticommutation relations applied).
- **Pure Bosons**: Validated that swapping operators does not introduce a sign change ($[b_i, b_j^\dagger] = \delta_{ij}$).
- **Mixed Systems**: Validated that fermions and bosons commute with each other without sign changes.

The implementation correctly defaults to fermionic behavior if `operator_stat` is empty, preserving backward compatibility.
