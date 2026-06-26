#!/usr/bin/env python3
"""Exact p=5 low-direction classification for the square-plane problem.

For G=(Z/25Z)^2 and H=5G, write a normalized section as

    s(x)=x+5 f(x),        x in F_5^2, f(0)=0.

For a nonzero quotient direction u, define

    F_u={s(x+u)-s(x): x in F_5^2}.

The command-line --pairs parameter chooses how many antipodal direction-pairs
are forced to be low.  The original run with --pairs 5 refuted the auxiliary
"at most 8 oriented low directions" target by finding a 10-oriented-direction
example.  The later --pairs 7 run is the fixed-p=5 certificate used by
C-WS0001-021: it enumerates GL_2(F_5)-orbits of 7-subsets of antipodal
direction-pairs and decides that no normalized section has all corresponding
14 oriented directions with |F_u|<=3.

The search is not random.  It is a finite proof-reduction/classification
routine with affine-lift normalization f(0)=f(e1)=f(e2)=0.
"""

from __future__ import annotations

import argparse
import json
from dataclasses import dataclass, field
from itertools import combinations, product
from functools import lru_cache
from pathlib import Path
from typing import Any

P = 5
Elt = tuple[int, int]
Matrix = tuple[int, int, int, int]

V: list[Elt] = list(product(range(P), repeat=2))
NONZERO: list[Elt] = [x for x in V if x != (0, 0)]
VALUE_ORDER: list[Elt] = V[:]
FIXED_VALUES: dict[Elt, Elt] = {
    (0, 0): (0, 0),
    (1, 0): (0, 0),
    (0, 1): (0, 0),
}


def add(x: Elt, y: Elt) -> Elt:
    return ((x[0] + y[0]) % P, (x[1] + y[1]) % P)


def sub(x: Elt, y: Elt) -> Elt:
    return ((x[0] - y[0]) % P, (x[1] - y[1]) % P)


def neg(x: Elt) -> Elt:
    return ((-x[0]) % P, (-x[1]) % P)


def sum_elts(values: list[Elt]) -> Elt:
    total = (0, 0)
    for value in values:
        total = add(total, value)
    return total


def mat_apply(m: Matrix, x: Elt) -> Elt:
    a, b, c, d = m
    return ((a * x[0] + b * x[1]) % P, (c * x[0] + d * x[1]) % P)


def det(m: Matrix) -> int:
    a, b, c, d = m
    return (a * d - b * c) % P


def gl2() -> list[Matrix]:
    return [m for m in product(range(P), repeat=4) if det(m) != 0]  # type: ignore[list-item]


def pair_rep(u: Elt) -> Elt:
    v = neg(u)
    return min(u, v)


PAIR_REPS: list[Elt] = sorted({pair_rep(u) for u in NONZERO})
PAIR_INDEX: dict[Elt, int] = {u: i for i, u in enumerate(PAIR_REPS)}


def pair_index(u: Elt) -> int:
    return PAIR_INDEX[pair_rep(u)]


def projective_rep(u: Elt) -> Elt:
    for c in range(1, P):
        v = ((c * u[0]) % P, (c * u[1]) % P)
        if v[0] == 1 or (v[0] == 0 and v[1] == 1):
            return v
    raise AssertionError(f"bad nonzero vector {u}")


PROJECTIVE_REPS: list[Elt] = sorted({projective_rep(u) for u in NONZERO})


def orbit_representatives(k: int) -> list[tuple[int, ...]]:
    matrices = gl2()
    actions: list[list[int]] = []
    for m in matrices:
        actions.append([pair_index(mat_apply(m, u)) for u in PAIR_REPS])

    seen: set[tuple[int, ...]] = set()
    reps: list[tuple[int, ...]] = []
    for subset in combinations(range(len(PAIR_REPS)), k):
        subset = tuple(subset)
        if subset in seen:
            continue
        images = {
            tuple(sorted(action[i] for i in subset))
            for action in actions
        }
        canonical = min(images)
        seen.update(images)
        reps.append(canonical)
    return sorted(set(reps))


def subset_profile(subset: tuple[int, ...]) -> dict[str, Any]:
    per_projective: dict[str, int] = {str(rep): 0 for rep in PROJECTIVE_REPS}
    directions = [PAIR_REPS[i] for i in subset]
    for u in directions:
        per_projective[str(projective_rep(u))] += 1
    counts = sorted(per_projective.values(), reverse=True)
    return {
        "pair_reps": directions,
        "projective_line_counts": per_projective,
        "projective_count_partition": counts,
    }


def carry(u: Elt, x: Elt) -> Elt:
    """Return ([x+u]-[x]-[u])/5 in F_5^2, coordinatewise."""

    return (
        (-1 if x[0] + u[0] >= P else 0) % P,
        (-1 if x[1] + u[1] >= P else 0) % P,
    )


@lru_cache(maxsize=None)
def repeated_sumset(values: tuple[Elt, ...], count: int) -> frozenset[Elt]:
    if count == 0:
        return frozenset({(0, 0)})
    previous = repeated_sumset(values, count - 1)
    return frozenset(add(x, y) for x in previous for y in values)


def affine_lines(direction: Elt) -> list[list[Elt]]:
    seen: set[Elt] = set()
    lines: list[list[Elt]] = []
    for x in V:
        if x in seen:
            continue
        line = []
        cur = x
        for _ in range(P):
            line.append(cur)
            seen.add(cur)
            cur = add(cur, direction)
        lines.append(line)
    return lines


def default_search_order(selected: set[Elt]) -> list[Elt]:
    fixed = set(FIXED_VALUES)
    order = list(FIXED_VALUES)
    assigned = set(order)
    remaining = [x for x in V if x not in fixed]
    while remaining:
        best = max(
            remaining,
            key=lambda x: (
                sum(1 for y in assigned if sub(x, y) in selected or sub(y, x) in selected),
                -x[0] - x[1],
            ),
        )
        order.append(best)
        assigned.add(best)
        remaining.remove(best)
    return order


@dataclass
class Solver:
    selected: set[Elt]
    limit: int = 3
    node_limit: int | None = None
    nodes: int = 0
    timed_out_by_node_limit: bool = False
    fibres: dict[Elt, set[Elt]] = field(default_factory=dict)
    lines_by_direction: dict[Elt, list[list[Elt]]] = field(default_factory=dict)
    assigned: dict[Elt, Elt] = field(default_factory=dict)
    unassigned: set[Elt] = field(default_factory=set)
    solution: dict[Elt, Elt] | None = None

    def __post_init__(self) -> None:
        self.fibres = {u: set() for u in self.selected}
        self.lines_by_direction = {u: affine_lines(u) for u in self.selected}
        self.unassigned = set(V)

    def add_assignment(self, q: Elt, value: Elt) -> list[tuple[Elt, Elt]]:
        additions: list[tuple[Elt, Elt]] = []
        for y, y_value in self.assigned.items():
            u = sub(q, y)
            if u in self.selected:
                h = add(sub(value, y_value), carry(u, y))
                fibre = self.fibres[u]
                if h not in fibre:
                    fibre.add(h)
                    additions.append((u, h))
            v = sub(y, q)
            if v in self.selected:
                h = add(sub(y_value, value), carry(v, q))
                fibre = self.fibres[v]
                if h not in fibre:
                    fibre.add(h)
                    additions.append((v, h))
        self.assigned[q] = value
        return additions

    def rollback(self, q: Elt, additions: list[tuple[Elt, Elt]]) -> None:
        del self.assigned[q]
        for u, h in reversed(additions):
            self.fibres[u].remove(h)

    def line_requirements(self, u: Elt) -> list[tuple[int, Elt]]:
        requirements: list[tuple[int, Elt]] = []
        for line in self.lines_by_direction[u]:
            known: list[Elt] = []
            unknown = 0
            for x in line:
                y = add(x, u)
                if x in self.assigned and y in self.assigned:
                    known.append(add(sub(self.assigned[y], self.assigned[x]), carry(u, x)))
                else:
                    unknown += 1
            needed = sub(neg(u), sum_elts(known))
            requirements.append((unknown, needed))
        return requirements

    def requirements_feasible_with_values(
        self, values: tuple[Elt, ...], requirements: list[tuple[int, Elt]]
    ) -> bool:
        return all(needed in repeated_sumset(values, unknown) for unknown, needed in requirements)

    def direction_cycle_feasible(self, u: Elt) -> bool:
        current = tuple(sorted(self.fibres[u]))
        if len(current) < 2:
            return True
        requirements = self.line_requirements(u)
        if len(current) == self.limit:
            return self.requirements_feasible_with_values(current, requirements)
        if len(current) == self.limit - 1:
            if self.requirements_feasible_with_values(current, requirements):
                return True
            current_set = set(current)
            for extra in V:
                if extra in current_set:
                    continue
                values = tuple(sorted((*current, extra)))
                if self.requirements_feasible_with_values(values, requirements):
                    return True
            return False
        # This branch is only relevant if limit is raised.  It keeps the
        # current p=5 proof search conservative rather than exponential.
        return True

    def feasible_current_fibres(self) -> bool:
        return all(len(values) <= self.limit for values in self.fibres.values()) and all(
            self.direction_cycle_feasible(u) for u in self.selected
        )

    def candidate_values(self, q: Elt) -> list[Elt]:
        candidates: list[Elt] = []
        for value in VALUE_ORDER:
            additions = self.add_assignment(q, value)
            if self.feasible_current_fibres():
                candidates.append(value)
            self.rollback(q, additions)
        return candidates

    def choose_next_variable(self) -> tuple[Elt, list[Elt]] | None:
        best_q: Elt | None = None
        best_candidates: list[Elt] | None = None
        best_assigned_neighbours = -1

        for q in sorted(self.unassigned):
            candidates = self.candidate_values(q)
            if not candidates:
                return q, []
            assigned_neighbours = sum(
                1
                for y in self.assigned
                if sub(q, y) in self.selected or sub(y, q) in self.selected
            )
            if (
                best_candidates is None
                or len(candidates) < len(best_candidates)
                or (
                    len(candidates) == len(best_candidates)
                    and assigned_neighbours > best_assigned_neighbours
                )
            ):
                best_q = q
                best_candidates = candidates
                best_assigned_neighbours = assigned_neighbours

        if best_q is None or best_candidates is None:
            return None
        return best_q, best_candidates

    def recurse(self) -> bool:
        if self.node_limit is not None and self.nodes >= self.node_limit:
            self.timed_out_by_node_limit = True
            return False
        self.nodes += 1
        if not self.unassigned:
            self.solution = dict(self.assigned)
            return True

        choice = self.choose_next_variable()
        if choice is None:
            self.solution = dict(self.assigned)
            return True
        q, candidates = choice
        if not candidates:
            return False

        self.unassigned.remove(q)
        for value in candidates:
            additions = self.add_assignment(q, value)
            if self.feasible_current_fibres() and self.recurse():
                return True
            self.rollback(q, additions)
            if self.timed_out_by_node_limit:
                self.unassigned.add(q)
                return False
        self.unassigned.add(q)
        return False

    def solve(self) -> dict[str, Any]:
        for q in sorted(FIXED_VALUES):
            self.unassigned.remove(q)
            additions = self.add_assignment(q, FIXED_VALUES[q])
            if not self.feasible_current_fibres():
                raise AssertionError(f"fixed normalization infeasible at {q}: {additions}")
        exists = self.recurse()
        return {
            "exists": exists,
            "nodes": self.nodes,
            "timed_out_by_node_limit": self.timed_out_by_node_limit,
            "solution": {str(k): v for k, v in sorted(self.solution.items())}
            if self.solution
            else None,
            "fibre_sizes_at_solution": {
                str(u): len(values) for u, values in sorted(self.fibres.items())
            }
            if self.solution
            else None,
        }


def selected_from_pairs(pair_subset: tuple[int, ...]) -> set[Elt]:
    result: set[Elt] = set()
    for i in pair_subset:
        u = PAIR_REPS[i]
        result.add(u)
        result.add(neg(u))
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--json-out", type=Path, required=True)
    parser.add_argument("--pairs", type=int, default=5)
    parser.add_argument("--node-limit", type=int, default=None)
    parser.add_argument("--only-orbit", type=int, default=None)
    args = parser.parse_args()

    reps = orbit_representatives(args.pairs)
    selected_reps = (
        [(args.only_orbit, reps[args.only_orbit])]
        if args.only_orbit is not None
        else list(enumerate(reps))
    )

    orbit_results: list[dict[str, Any]] = []
    for orbit_index, rep in selected_reps:
        solver = Solver(selected_from_pairs(rep), node_limit=args.node_limit)
        solved = solver.solve()
        orbit_results.append(
            {
                "orbit_index": orbit_index,
                "representative_pair_indices": rep,
                "profile": subset_profile(rep),
                "search": solved,
            }
        )

    any_exists = any(item["search"]["exists"] for item in orbit_results)
    any_timeout = any(item["search"]["timed_out_by_node_limit"] for item in orbit_results)
    result = {
        "script": "experiments/python/ws0001_p5_low_direction_classification.py",
        "case": {
            "prime": P,
            "group": "(Z/25Z)^2",
            "subgroup": "5G",
            "quotient": "C5 x C5",
            "low_direction_limit": 3,
            "affine_lift_normalization": "f(0,0)=f(1,0)=f(0,1)=0",
            "antipodal_pair_representatives": PAIR_REPS,
            "gl2_order": len(gl2()),
        },
        "classification": {
            "antipodal_pairs_requested": args.pairs,
            "orbit_count": len(reps),
            "orbits_checked": len(orbit_results),
            "any_representative_exists": any_exists,
            "any_node_limit_timeout": any_timeout,
            "interpretation": (
                f"If all {args.pairs}-pair orbit representatives are infeasible "
                "without node-limit timeout, then no normalized section over "
                f"F_5^2 has at least {2 * args.pairs} oriented nonzero "
                "directions u with |F_u|<=3."
            ),
        },
        "orbit_results": orbit_results,
    }

    args.json_out.parent.mkdir(parents=True, exist_ok=True)
    args.json_out.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps(result["classification"], indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
