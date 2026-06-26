#!/usr/bin/env python3
"""Targeted C3 x C3 residual counterexample search for WS-0001.

This checks the first genuinely noncyclic quotient with no C2 x C2 plane:

    G = (Z/9Z)^2,    H = 3G,    G/H ~= C3 x C3.

The pair is residual because H contains the 3-torsion of G.  The script uses a
backtracking search over normalized transversals and proves that no transversal
has |T-T| <= 23.  The standard box transversal has |T-T| = 25, so the exact
value is 25.
"""

from __future__ import annotations

import argparse
import json
from itertools import product
from pathlib import Path
from typing import Any

from ws0001_quotient_fibre_sharpness_check import (
    Elt,
    all_subgroups,
    difference_set,
    elems,
    m_value,
)

Q_ELTS: list[Elt] = list(product(range(3), repeat=2))
SEARCH_ORDER: list[Elt] = [
    (0, 0),
    (1, 0),
    (0, 1),
    (1, 1),
    (1, 2),
    (2, 0),
    (0, 2),
    (2, 1),
    (2, 2),
]
H_COORDS: list[Elt] = list(product(range(3), repeat=2))


def q_sub(x: Elt, y: Elt) -> Elt:
    return ((x[0] - y[0]) % 3, (x[1] - y[1]) % 3)


def q_neg(x: Elt) -> Elt:
    return ((-x[0]) % 3, (-x[1]) % 3)


def lift(q: Elt, h_coord: Elt) -> Elt:
    return ((q[0] + 3 * h_coord[0]) % 9, (q[1] + 3 * h_coord[1]) % 9)


def g_sub(x: Elt, y: Elt) -> Elt:
    return ((x[0] - y[0]) % 9, (x[1] - y[1]) % 9)


def direction_pairs() -> list[tuple[Elt, Elt]]:
    pairs = []
    seen = {(0, 0)}
    for a in Q_ELTS:
        if a in seen:
            continue
        b = q_neg(a)
        pairs.append((a, b))
        seen.add(a)
        seen.add(b)
    return pairs


def lower_bound_from_partial_fibres(
    fibres: dict[Elt, set[Elt]], pairs: list[tuple[Elt, Elt]]
) -> int:
    """Lower bound for any residual completion extending current fibres.

    In this odd-order example, F_{-a}=-F_a, so opposite fibres have equal final
    size.  Since m(G,H)=1, every nonzero final fibre has size at least two.
    """

    total = 1
    for a, b in pairs:
        total += 2 * max(len(fibres.get(a, set())), len(fibres.get(b, set())), 2)
    return total


def exists_transversal_at_most(limit: int) -> dict[str, Any]:
    """Search for a normalized transversal with |T-T| <= limit."""

    pairs = direction_pairs()
    assigned: dict[Elt, Elt] = {(0, 0): (0, 0)}
    points: dict[Elt, Elt] = {(0, 0): lift((0, 0), (0, 0))}
    fibres: dict[Elt, set[Elt]] = {(0, 0): {(0, 0)}}
    nodes = 0
    solution: dict[Elt, Elt] | None = None
    solution_fibre_sizes: dict[str, int] | None = None
    solution_size: int | None = None

    def recurse(i: int) -> bool:
        nonlocal nodes, solution, solution_fibre_sizes, solution_size
        nodes += 1
        if lower_bound_from_partial_fibres(fibres, pairs) > limit:
            return False

        if i == len(SEARCH_ORDER):
            sizes = {a: len(fibres.get(a, set())) for a in Q_ELTS}
            total = sum(sizes.values())
            if total <= limit and all(sizes[a] >= 2 for a in Q_ELTS if a != (0, 0)):
                solution = dict(assigned)
                solution_fibre_sizes = {str(k): v for k, v in sorted(sizes.items())}
                solution_size = total
                return True
            return False

        q = SEARCH_ORDER[i]
        for h_coord in H_COORDS:
            p = lift(q, h_coord)
            additions: list[tuple[Elt, Elt]] = []
            ok = True

            for q2, p2 in points.items():
                for a, d in ((q_sub(q, q2), g_sub(p, p2)), (q_sub(q2, q), g_sub(p2, p))):
                    fibre = fibres.setdefault(a, set())
                    if d in fibre:
                        continue
                    fibre.add(d)
                    additions.append((a, d))
                    if a == (0, 0) and len(fibre) > 1:
                        ok = False
                        break
                if not ok:
                    break

            if ok and lower_bound_from_partial_fibres(fibres, pairs) <= limit:
                assigned[q] = h_coord
                points[q] = p
                if recurse(i + 1):
                    return True
                del assigned[q]
                del points[q]

            for a, d in reversed(additions):
                fibres[a].remove(d)
                if not fibres[a]:
                    del fibres[a]

        return False

    found = recurse(1)
    return {
        "limit": limit,
        "exists": found,
        "nodes_visited": nodes,
        "solution_h_coordinates": {str(k): v for k, v in sorted(solution.items())}
        if solution
        else None,
        "solution_diffset_size": solution_size,
        "solution_fibre_sizes": solution_fibre_sizes,
    }


def box_transversal_summary() -> dict[str, Any]:
    T = tuple((x, y) for x in range(3) for y in range(3))
    D = difference_set((9, 9), T)
    return {
        "transversal": T,
        "diffset_size": len(D),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--json-out", type=Path, required=True)
    args = parser.parse_args()

    mods = (9, 9)
    H = frozenset((3 * x, 3 * y) for x, y in product(range(3), repeat=2))
    subgroups = all_subgroups(mods)
    m = m_value(mods, H, subgroups)
    quotient_order = len(elems(mods)) // len(H)
    lower = 2 * quotient_order - m
    threshold_results = [exists_transversal_at_most(limit) for limit in (17, 19, 21, 23, 25)]
    box = box_transversal_summary()

    result = {
        "script": "experiments/python/ws0001_c3_plane_counterexample.py",
        "case": {
            "group_moduli": mods,
            "group_order": len(elems(mods)),
            "subgroup": sorted(H),
            "subgroup_order": len(H),
            "quotient_invariants": [3, 3],
            "quotient_order": quotient_order,
            "quotient_has_c2x_c2": False,
            "m_value": m,
            "quotient_fibre_lower_bound": lower,
            "normalized_transversal_count": len(H) ** (quotient_order - 1),
        },
        "threshold_searches": threshold_results,
        "upper_bound": box,
        "conclusion": {
            "exact_delta": box["diffset_size"],
            "gap_over_quotient_fibre_lower_bound": box["diffset_size"] - lower,
            "reason": (
                "No normalized transversal has |T-T| <= 23, while the box "
                "transversal has |T-T| = 25. Since T-T is negation-stable in "
                "the odd-order group (Z/9Z)^2, its cardinality is odd."
            ),
        },
    }

    args.json_out.parent.mkdir(parents=True, exist_ok=True)
    args.json_out.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps(result["conclusion"], indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
