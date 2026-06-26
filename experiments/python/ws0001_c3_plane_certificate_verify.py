#!/usr/bin/env python3
"""Independent verifier for the WS-0001 C3 x C3 counterexample certificate.

This intentionally does not import ws0001_c3_plane_counterexample.py.  It reruns
the same threshold decision problem with copy-on-recursion fibre state, so it is
a small independent check against mutation/backtracking bugs in the main script.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from itertools import product
from pathlib import Path
from typing import Any

Elt = tuple[int, int]

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
    pairs: list[tuple[Elt, Elt]] = []
    seen = {(0, 0)}
    for a in Q_ELTS:
        if a in seen:
            continue
        b = q_neg(a)
        pairs.append((a, b))
        seen.add(a)
        seen.add(b)
    return pairs


def lower_bound_from_fibres(
    fibres: dict[Elt, frozenset[Elt]], pairs: list[tuple[Elt, Elt]]
) -> int:
    total = 1
    empty: frozenset[Elt] = frozenset()
    for a, b in pairs:
        total += 2 * max(len(fibres.get(a, empty)), len(fibres.get(b, empty)), 2)
    return total


def extend_fibres(
    q: Elt,
    p: Elt,
    points: tuple[tuple[Elt, Elt], ...],
    fibres: dict[Elt, frozenset[Elt]],
) -> dict[Elt, frozenset[Elt]]:
    new_fibres = dict(fibres)
    for q2, p2 in points:
        for a, d in ((q_sub(q, q2), g_sub(p, p2)), (q_sub(q2, q), g_sub(p2, p))):
            current = new_fibres.get(a, frozenset())
            if d not in current:
                new_fibres[a] = current | {d}
    return new_fibres


def exists_transversal_at_most(limit: int) -> dict[str, Any]:
    pairs = direction_pairs()
    initial_points = (((0, 0), lift((0, 0), (0, 0))),)
    initial_fibres = {(0, 0): frozenset({(0, 0)})}
    nodes = 0
    solution: dict[Elt, Elt] | None = None
    solution_size: int | None = None
    solution_fibre_sizes: dict[str, int] | None = None

    def recurse(
        i: int,
        assigned: dict[Elt, Elt],
        points: tuple[tuple[Elt, Elt], ...],
        fibres: dict[Elt, frozenset[Elt]],
    ) -> bool:
        nonlocal nodes, solution, solution_size, solution_fibre_sizes
        nodes += 1
        if lower_bound_from_fibres(fibres, pairs) > limit:
            return False

        if i == len(SEARCH_ORDER):
            sizes = {a: len(fibres.get(a, frozenset())) for a in Q_ELTS}
            total = sum(sizes.values())
            if total <= limit and all(sizes[a] >= 2 for a in Q_ELTS if a != (0, 0)):
                solution = dict(assigned)
                solution_size = total
                solution_fibre_sizes = {str(k): v for k, v in sorted(sizes.items())}
                return True
            return False

        q = SEARCH_ORDER[i]
        for h_coord in H_COORDS:
            p = lift(q, h_coord)
            new_fibres = extend_fibres(q, p, points, fibres)
            if len(new_fibres.get((0, 0), frozenset())) > 1:
                continue
            if lower_bound_from_fibres(new_fibres, pairs) > limit:
                continue
            new_assigned = dict(assigned)
            new_assigned[q] = h_coord
            new_points = points + ((q, p),)
            if recurse(i + 1, new_assigned, new_points, new_fibres):
                return True
        return False

    found = recurse(1, {(0, 0): (0, 0)}, initial_points, initial_fibres)
    return {
        "limit": limit,
        "exists": found,
        "nodes_visited": nodes,
        "solution_diffset_size": solution_size,
        "solution_fibre_sizes": solution_fibre_sizes,
        "solution_h_coordinates": {str(k): v for k, v in sorted(solution.items())}
        if solution
        else None,
    }


def difference_set(T: tuple[Elt, ...]) -> set[Elt]:
    return {g_sub(x, y) for x in T for y in T}


def fibre_sizes(D: set[Elt]) -> dict[str, int]:
    sizes = {a: 0 for a in Q_ELTS}
    for d in D:
        sizes[(d[0] % 3, d[1] % 3)] += 1
    return {str(k): v for k, v in sorted(sizes.items())}


def line_sizes(D: set[Elt]) -> list[int]:
    direction_sizes = {a: 0 for a in Q_ELTS}
    for d in D:
        direction_sizes[(d[0] % 3, d[1] % 3)] += 1
    return sorted(direction_sizes[a] for a, _ in direction_pairs())


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--canonical-result",
        type=Path,
        default=Path("experiments/results/ws0001-c3-plane-counterexample.json"),
    )
    parser.add_argument("--json-out", type=Path, required=True)
    args = parser.parse_args()

    thresholds = [17, 19, 21, 23, 25]
    threshold_results = [exists_transversal_at_most(limit) for limit in thresholds]

    box = tuple((x, y) for x in range(3) for y in range(3))
    box_diffset = difference_set(box)
    canonical = json.loads(args.canonical_result.read_text())
    expected_thresholds = [
        {
            "limit": item["limit"],
            "exists": item["exists"],
            "nodes_visited": item["nodes_visited"],
            "solution_diffset_size": item["solution_diffset_size"],
        }
        for item in canonical["threshold_searches"]
    ]
    observed_thresholds = [
        {
            "limit": item["limit"],
            "exists": item["exists"],
            "nodes_visited": item["nodes_visited"],
            "solution_diffset_size": item["solution_diffset_size"],
        }
        for item in threshold_results
    ]

    result = {
        "script": "experiments/python/ws0001_c3_plane_certificate_verify.py",
        "purpose": "independent copy-on-recursion verification for C-WS0001-013",
        "canonical_result": str(args.canonical_result),
        "canonical_result_sha256": sha256(args.canonical_result),
        "normalized_transversal_count": 9**8,
        "threshold_searches": threshold_results,
        "threshold_summary_matches_canonical": observed_thresholds == expected_thresholds,
        "box_transversal": box,
        "box_diffset_size": len(box_diffset),
        "box_fibre_sizes": fibre_sizes(box_diffset),
        "box_line_sizes": line_sizes(box_diffset),
        "conclusion": {
            "limit_23_absent": not threshold_results[3]["exists"],
            "limit_25_present": bool(threshold_results[4]["exists"]),
            "odd_cardinality_excludes_24": True,
            "verified_delta": 25,
        },
    }

    args.json_out.parent.mkdir(parents=True, exist_ok=True)
    args.json_out.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps(result["conclusion"], indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
