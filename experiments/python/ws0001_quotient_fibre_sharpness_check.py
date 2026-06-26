#!/usr/bin/env python3
"""Exact checks for quotient-fibre lower-bound sharpness.

This script is intentionally dependency-free.  It verifies the smallest
non-sharp example suggested by the expert answer and records fibre histograms
for small primitive residual examples.
"""

from __future__ import annotations

import argparse
import json
from collections import Counter, deque
from itertools import product
from math import prod
from pathlib import Path
from typing import Iterable

Elt = tuple[int, ...]


def invariant_tuples_upto(nmax: int) -> list[tuple[int, ...]]:
    out: list[tuple[int, ...]] = []

    def rec(start: int, cur_prod: int, arr: list[int]) -> None:
        if arr:
            out.append(tuple(arr))
        for n in range(start, nmax // cur_prod + 1):
            if arr and n % arr[-1] != 0:
                continue
            if cur_prod * n <= nmax:
                rec(n, cur_prod * n, arr + [n])

    rec(2, 1, [])
    return sorted(out, key=lambda t: (prod(t), len(t), t))


def elems(mods: tuple[int, ...]) -> list[Elt]:
    return list(product(*[range(n) for n in mods]))


def zero(mods: tuple[int, ...]) -> Elt:
    return tuple(0 for _ in mods)


def add(x: Elt, y: Elt, mods: tuple[int, ...]) -> Elt:
    return tuple((x[i] + y[i]) % mods[i] for i in range(len(mods)))


def sub(x: Elt, y: Elt, mods: tuple[int, ...]) -> Elt:
    return tuple((x[i] - y[i]) % mods[i] for i in range(len(mods)))


def neg(x: Elt, mods: tuple[int, ...]) -> Elt:
    return tuple((-x[i]) % mods[i] for i in range(len(mods)))


def generated_subgroup(mods: tuple[int, ...], gens: Iterable[Elt]) -> frozenset[Elt]:
    z = zero(mods)
    H = {z}
    q: deque[Elt] = deque([z])
    genlist = []
    for g in gens:
        if g != z:
            genlist.append(g)
            genlist.append(neg(g, mods))

    while q:
        h = q.popleft()
        for g in genlist:
            u = add(h, g, mods)
            if u not in H:
                H.add(u)
                q.append(u)

    return frozenset(H)


def all_subgroups(mods: tuple[int, ...]) -> list[frozenset[Elt]]:
    G = elems(mods)
    z = zero(mods)
    seen = {frozenset([z])}
    q: deque[frozenset[Elt]] = deque([frozenset([z])])

    while q:
        H = q.popleft()
        for g in G:
            if g not in H:
                K = generated_subgroup(mods, list(H) + [g])
                if K not in seen:
                    seen.add(K)
                    q.append(K)

    return sorted(seen, key=lambda H: (len(H), sorted(H)))


def cosets_of_subgroup(
    mods: tuple[int, ...], H: frozenset[Elt]
) -> tuple[list[Elt], list[frozenset[Elt]]]:
    unused = set(elems(mods))
    reps: list[Elt] = []
    cosets: list[frozenset[Elt]] = []
    z = zero(mods)

    while unused:
        r = min(unused)
        C = frozenset(add(r, h, mods) for h in H)
        reps.append(r)
        cosets.append(C)
        unused -= set(C)

    idx = next(i for i, C in enumerate(cosets) if z in C)
    reps[0], reps[idx] = reps[idx], reps[0]
    cosets[0], cosets[idx] = cosets[idx], cosets[0]
    return reps, cosets


def difference_set(mods: tuple[int, ...], T: Iterable[Elt]) -> set[Elt]:
    T_tuple = tuple(T)
    return {sub(x, y, mods) for x in T_tuple for y in T_tuple}


def m_value(
    mods: tuple[int, ...], H: frozenset[Elt], subgroups: list[frozenset[Elt]]
) -> int:
    z = zero(mods)
    Hset = set(H)
    return max(len(K) for K in subgroups if set(K) & Hset == {z})


def exact_delta(
    mods: tuple[int, ...], H: frozenset[Elt], max_transversals: int
) -> dict[str, object]:
    reps, cosets = cosets_of_subgroup(mods, H)
    total = len(H) ** (len(cosets) - 1)
    if total > max_transversals:
        return {
            "complete": False,
            "total_transversals_after_translation": total,
            "reason": f"exceeds cap {max_transversals}",
        }

    Hlist = list(H)
    z = zero(mods)
    best: int | None = None
    best_T: tuple[Elt, ...] | None = None
    best_count = 0
    diff_hist: Counter[int] = Counter()
    fibre_hist: Counter[tuple[int, ...]] = Counter()

    for choices in product(Hlist, repeat=len(cosets) - 1):
        T = tuple([z] + [add(r, h, mods) for r, h in zip(reps[1:], choices)])
        D = difference_set(mods, T)
        d = len(D)
        diff_hist[d] += 1
        fibre_sizes = tuple(sorted(len(D & set(C)) for C in cosets))
        fibre_hist[fibre_sizes] += 1
        if best is None or d < best:
            best = d
            best_T = T
            best_count = 1
        elif d == best:
            best_count += 1

    return {
        "complete": True,
        "total_transversals_after_translation": total,
        "delta": best,
        "best_T": best_T,
        "best_count": best_count,
        "diffset_histogram": {str(k): diff_hist[k] for k in sorted(diff_hist)},
        "fibre_histogram": {
            ",".join(map(str, k)): fibre_hist[k] for k in sorted(fibre_hist)
        },
    }


def summarize_pair(
    mods: tuple[int, ...],
    H: frozenset[Elt],
    max_transversals: int,
    case_id: str,
) -> dict[str, object]:
    subgroups = all_subgroups(mods)
    G = elems(mods)
    quotient_size = len(G) // len(H)
    m = m_value(mods, H, subgroups)
    lower = 2 * quotient_size - m
    exact = exact_delta(mods, H, max_transversals)
    delta = exact.get("delta")
    return {
        "case_id": case_id,
        "group_moduli": mods,
        "group_order": len(G),
        "subgroup_order": len(H),
        "quotient_size": quotient_size,
        "m_value": m,
        "quotient_fibre_lower_bound": lower,
        "delta": delta,
        "gap_over_quotient_fibre_lower_bound": (
            delta - lower if isinstance(delta, int) else None
        ),
        "subgroups_examined": len(subgroups),
        "exact_search": exact,
    }


def check_groups_below(
    max_order_exclusive: int, max_transversals: int
) -> dict[str, object]:
    checked = 0
    skipped = []
    gaps = []
    for mods in invariant_tuples_upto(max_order_exclusive - 1):
        G = elems(mods)
        subgroups = all_subgroups(mods)
        for H in subgroups:
            if len(H) == 1 or len(H) == len(G):
                continue
            total = len(H) ** (len(G) // len(H) - 1)
            if total > max_transversals:
                skipped.append(
                    {
                        "group_moduli": mods,
                        "subgroup_order": len(H),
                        "quotient_size": len(G) // len(H),
                        "total_transversals_after_translation": total,
                    }
                )
                continue
            checked += 1
            m = m_value(mods, H, subgroups)
            lower = 2 * (len(G) // len(H)) - m
            exact = exact_delta(mods, H, max_transversals)
            delta = exact["delta"]
            if delta > lower:
                gaps.append(
                    {
                        "group_moduli": mods,
                        "subgroup": sorted(H),
                        "subgroup_order": len(H),
                        "quotient_size": len(G) // len(H),
                        "m_value": m,
                        "lower_bound": lower,
                        "delta": delta,
                    }
                )
    return {
        "max_order_exclusive": max_order_exclusive,
        "checked_nontrivial_pairs": checked,
        "skipped_by_transversal_cap": skipped,
        "gap_count": len(gaps),
        "gaps": gaps,
    }


def c4_power_pair(r: int) -> tuple[tuple[int, ...], frozenset[Elt]]:
    mods = tuple(4 for _ in range(r))
    H = frozenset(
        tuple(2 * bit for bit in bits) for bits in product([0, 1], repeat=r)
    )
    return mods, H


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--max-order-exclusive", type=int, default=16)
    parser.add_argument("--max-transversals", type=int, default=3_000_000)
    parser.add_argument("--c4-family-rmax", type=int, default=3)
    parser.add_argument("--json-out", type=Path, required=True)
    args = parser.parse_args()

    below = check_groups_below(args.max_order_exclusive, args.max_transversals)
    c4_family = []
    for r in range(2, args.c4_family_rmax + 1):
        mods, H = c4_power_pair(r)
        c4_family.append(
            summarize_pair(
                mods,
                H,
                args.max_transversals,
                case_id=f"C4^{r}-2G",
            )
        )

    result = {
        "script": "experiments/python/ws0001_quotient_fibre_sharpness_check.py",
        "parameters": {
            "max_order_exclusive": args.max_order_exclusive,
            "max_transversals": args.max_transversals,
            "c4_family_rmax": args.c4_family_rmax,
        },
        "groups_below_threshold": below,
        "c4_power_family": c4_family,
    }
    args.json_out.parent.mkdir(parents=True, exist_ok=True)
    args.json_out.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps(result["groups_below_threshold"], indent=2, sort_keys=True))
    for case in c4_family:
        print(
            case["case_id"],
            "delta=",
            case["delta"],
            "lower=",
            case["quotient_fibre_lower_bound"],
            "gap=",
            case["gap_over_quotient_fibre_lower_bound"],
        )


if __name__ == "__main__":
    main()
