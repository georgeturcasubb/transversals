#!/usr/bin/env python3
"""Profile-level p=5 weighted triangle-deficit classification.

This is not a solver for the full square-plane problem.  It records the finite
profile reduction left after the reviewed parallel-type and triangle-count
criteria:

* In a counterexample, all collinear positive-deficit fibres must lie on one
  projective line.
* Every remaining positive-deficit fibre is triangle-type of size 3.
* For p=5, a counterexample to the box bound needs signed deficit at least 18
  from the baseline size 4 over the 24 oriented nonzero directions.

The script enumerates the minimal labelled antipodal-pair profiles that can
reach that deficit, modulo GL_2(F_5), and audits the known stress-test
sections against the same profile vocabulary.
"""

from __future__ import annotations

import argparse
import ast
import json
from dataclasses import dataclass
from itertools import combinations, product
from pathlib import Path
from typing import Any

P = 5
Elt = tuple[int, int]
Matrix = tuple[int, int, int, int]

V: list[Elt] = list(product(range(P), repeat=2))
NONZERO: list[Elt] = [x for x in V if x != (0, 0)]

NEUTRAL = 0
PARALLEL_SIZE2 = 1
PARALLEL_SIZE3 = 2
TRIANGLE_SIZE3 = 3
CATEGORY_NAMES = {
    NEUTRAL: "neutral_or_large",
    PARALLEL_SIZE2: "parallel_size2",
    PARALLEL_SIZE3: "parallel_size3",
    TRIANGLE_SIZE3: "triangle_size3",
}


def add(x: Elt, y: Elt) -> Elt:
    return ((x[0] + y[0]) % P, (x[1] + y[1]) % P)


def sub(x: Elt, y: Elt) -> Elt:
    return ((x[0] - y[0]) % P, (x[1] - y[1]) % P)


def neg(x: Elt) -> Elt:
    return ((-x[0]) % P, (-x[1]) % P)


def scalar_mul(c: int, x: Elt) -> Elt:
    return ((c * x[0]) % P, (c * x[1]) % P)


def pair_rep(u: Elt) -> Elt:
    return min(u, neg(u))


PAIR_REPS: list[Elt] = sorted({pair_rep(u) for u in NONZERO})
PAIR_INDEX: dict[Elt, int] = {u: i for i, u in enumerate(PAIR_REPS)}


def pair_index(u: Elt) -> int:
    return PAIR_INDEX[pair_rep(u)]


def projective_rep(u: Elt) -> Elt:
    if u == (0, 0):
        raise ValueError("zero vector has no projective representative")
    for c in range(1, P):
        v = scalar_mul(c, u)
        if v[0] == 1 or (v[0] == 0 and v[1] == 1):
            return v
    raise AssertionError(f"bad nonzero vector {u}")


PROJECTIVE_REPS: list[Elt] = sorted({projective_rep(u) for u in NONZERO})
PAIR_PROJECTIVE_REP: dict[int, Elt] = {
    i: projective_rep(u) for i, u in enumerate(PAIR_REPS)
}
PAIRS_BY_PROJECTIVE: dict[Elt, list[int]] = {
    line: [i for i, u in enumerate(PAIR_REPS) if projective_rep(u) == line]
    for line in PROJECTIVE_REPS
}


def mat_apply(m: Matrix, x: Elt) -> Elt:
    a, b, c, d = m
    return ((a * x[0] + b * x[1]) % P, (c * x[0] + d * x[1]) % P)


def det(m: Matrix) -> int:
    a, b, c, d = m
    return (a * d - b * c) % P


def gl2() -> list[Matrix]:
    return [m for m in product(range(P), repeat=4) if det(m) != 0]  # type: ignore[list-item]


def gl2_pair_actions() -> list[list[int]]:
    actions: list[list[int]] = []
    for m in gl2():
        actions.append([pair_index(mat_apply(m, u)) for u in PAIR_REPS])
    return actions


GL2_PAIR_ACTIONS = gl2_pair_actions()


def canonical_assignment(categories: tuple[int, ...]) -> tuple[int, ...]:
    images = []
    for action in GL2_PAIR_ACTIONS:
        image = [NEUTRAL] * len(PAIR_REPS)
        for old_i, new_i in enumerate(action):
            image[new_i] = categories[old_i]
        images.append(tuple(image))
    return min(images)


def category_deficit(category: int) -> int:
    if category == PARALLEL_SIZE2:
        return 4
    if category in (PARALLEL_SIZE3, TRIANGLE_SIZE3):
        return 2
    return 0


def line_profile(categories: tuple[int, ...]) -> list[dict[str, Any]]:
    result = []
    for line in PROJECTIVE_REPS:
        counts = {name: 0 for name in CATEGORY_NAMES.values()}
        for pair_i in PAIRS_BY_PROJECTIVE[line]:
            counts[CATEGORY_NAMES[categories[pair_i]]] += 1
        result.append(
            {
                "projective_line": line,
                "counts": counts,
                "positive_pair_deficit": sum(
                    category_deficit(categories[pair_i])
                    for pair_i in PAIRS_BY_PROJECTIVE[line]
                ),
            }
        )
    return result


def assignment_record(categories: tuple[int, ...]) -> dict[str, Any]:
    nonneutral = [
        {
            "pair_index": i,
            "pair_rep": PAIR_REPS[i],
            "projective_line": PAIR_PROJECTIVE_REP[i],
            "category": CATEGORY_NAMES[category],
            "oriented_deficit_contribution": category_deficit(category),
        }
        for i, category in enumerate(categories)
        if category != NEUTRAL
    ]
    return {
        "nonneutral_pairs": nonneutral,
        "line_profile": line_profile(categories),
        "total_oriented_positive_deficit": sum(category_deficit(c) for c in categories),
        "triangle_pair_count": sum(1 for c in categories if c == TRIANGLE_SIZE3),
        "parallel_positive_pair_count": sum(
            1 for c in categories if c in (PARALLEL_SIZE2, PARALLEL_SIZE3)
        ),
    }


@dataclass(frozen=True)
class Template:
    parallel_size2_pairs: int
    parallel_size3_pairs: int
    triangle_size3_pairs: int

    @property
    def oriented_positive_deficit(self) -> int:
        return 4 * self.parallel_size2_pairs + 2 * self.parallel_size3_pairs + 2 * self.triangle_size3_pairs

    @property
    def positive_low_pair_count(self) -> int:
        return self.parallel_size2_pairs + self.parallel_size3_pairs + self.triangle_size3_pairs


def minimal_templates() -> list[Template]:
    templates = []
    for p2 in range(3):
        for p3 in range(3 - p2):
            tri = 9 - 2 * p2 - p3
            if tri >= 5 and p2 + p3 <= 2 and tri + p2 + p3 <= len(PAIR_REPS):
                template = Template(p2, p3, tri)
                if template.oriented_positive_deficit == 18:
                    templates.append(template)
    return templates


def enumerate_template_orbits(template: Template) -> dict[str, Any]:
    canonical: dict[tuple[int, ...], int] = {}
    raw_assignments = 0
    all_pairs = set(range(len(PAIR_REPS)))

    for exceptional_line in PROJECTIVE_REPS:
        exceptional_pairs = PAIRS_BY_PROJECTIVE[exceptional_line]
        for p2_pairs in combinations(exceptional_pairs, template.parallel_size2_pairs):
            remaining_exceptional = [i for i in exceptional_pairs if i not in p2_pairs]
            for p3_pairs in combinations(remaining_exceptional, template.parallel_size3_pairs):
                blocked = set(p2_pairs) | set(p3_pairs)
                triangle_available = sorted(all_pairs - blocked)
                for tri_pairs in combinations(triangle_available, template.triangle_size3_pairs):
                    categories = [NEUTRAL] * len(PAIR_REPS)
                    for i in p2_pairs:
                        categories[i] = PARALLEL_SIZE2
                    for i in p3_pairs:
                        categories[i] = PARALLEL_SIZE3
                    for i in tri_pairs:
                        categories[i] = TRIANGLE_SIZE3
                    raw_assignments += 1
                    key = canonical_assignment(tuple(categories))
                    canonical[key] = canonical.get(key, 0) + 1

    orbit_records = []
    for orbit_index, (categories, multiplicity) in enumerate(sorted(canonical.items()), start=1):
        record = assignment_record(categories)
        record["orbit_index"] = orbit_index
        record["raw_assignment_multiplicity"] = multiplicity
        orbit_records.append(record)

    return {
        "template": {
            "parallel_size2_pairs": template.parallel_size2_pairs,
            "parallel_size3_pairs": template.parallel_size3_pairs,
            "triangle_size3_pairs": template.triangle_size3_pairs,
            "positive_low_pair_count": template.positive_low_pair_count,
            "oriented_positive_deficit": template.oriented_positive_deficit,
        },
        "raw_assignments": raw_assignments,
        "gl2_orbit_count": len(canonical),
        "orbits": orbit_records,
    }


def carry(u: Elt, x: Elt) -> Elt:
    return (
        (-1 if x[0] + u[0] >= P else 0) % P,
        (-1 if x[1] + u[1] >= P else 0) % P,
    )


def lift(q: Elt, h: Elt) -> Elt:
    return ((q[0] + P * h[0]) % (P * P), (q[1] + P * h[1]) % (P * P))


def g_sub(x: Elt, y: Elt) -> Elt:
    return ((x[0] - y[0]) % (P * P), (x[1] - y[1]) % (P * P))


def derivative_image(section: dict[Elt, Elt], u: Elt) -> set[Elt]:
    return {
        add(sub(section[add(x, u)], section[x]), carry(u, x))
        for x in V
    }


def is_collinear(points: set[Elt]) -> bool:
    if len(points) <= 2:
        return True
    pts = list(points)
    base = pts[0]
    directions = [sub(pt, base) for pt in pts[1:] if pt != base]
    nonzero = [d for d in directions if d != (0, 0)]
    if not nonzero:
        return True
    line = projective_rep(nonzero[0])
    return all(projective_rep(d) == line for d in nonzero)


def affine_line_direction(points: set[Elt]) -> Elt | None:
    if len(points) <= 1:
        return None
    pts = list(points)
    for i, x in enumerate(pts):
        for y in pts[i + 1 :]:
            diff = sub(y, x)
            if diff != (0, 0):
                return projective_rep(diff)
    return None


def section_metrics(section: dict[Elt, Elt], label: str) -> dict[str, Any]:
    images = {u: derivative_image(section, u) for u in NONZERO}
    sizes = {u: len(image) for u, image in images.items()}
    points = [lift(q, h) for q, h in section.items()]
    diffset_size = len({g_sub(x, y) for x in points for y in points})

    direction_records = []
    for u in sorted(NONZERO):
        image = images[u]
        size = sizes[u]
        collinear = is_collinear(image)
        line_dir = affine_line_direction(image) if collinear else None
        parallel_type = collinear and line_dir == projective_rep(u)
        triangle_type = size == 3 and not collinear
        direction_records.append(
            {
                "direction": u,
                "projective_line": projective_rep(u),
                "size": size,
                "signed_deficit_from_4": 4 - size,
                "collinear": collinear,
                "parallel_type": parallel_type,
                "triangle_type": triangle_type,
                "image": sorted(image),
            }
        )

    profile: dict[str, int] = {}
    for size in sizes.values():
        profile[str(size)] = profile.get(str(size), 0) + 1

    line_sizes: dict[str, list[int]] = {str(line): [] for line in PROJECTIVE_REPS}
    line_triangle_counts: dict[str, int] = {str(line): 0 for line in PROJECTIVE_REPS}
    line_parallel_positive_counts: dict[str, int] = {str(line): 0 for line in PROJECTIVE_REPS}
    for record in direction_records:
        line_key = str(record["projective_line"])
        line_sizes[line_key].append(record["size"])
        if record["triangle_type"]:
            line_triangle_counts[line_key] += 1
        if record["parallel_type"] and record["size"] < 4:
            line_parallel_positive_counts[line_key] += 1

    signed_deficit = sum(4 - size for size in sizes.values())
    positive_deficit = sum(max(0, 4 - size) for size in sizes.values())
    large_penalty = sum(max(0, size - 4) for size in sizes.values())

    return {
        "label": label,
        "diffset_size": diffset_size,
        "fibre_sum": sum(sizes.values()),
        "fibre_sum_identity_holds": diffset_size == 1 + sum(sizes.values()),
        "box_diffset_size": (2 * P - 1) ** 2,
        "signed_deficit_from_4": signed_deficit,
        "positive_deficit_from_4": positive_deficit,
        "large_fibre_penalty_over_4": large_penalty,
        "profile": dict(sorted(profile.items(), key=lambda kv: int(kv[0]))),
        "triangle_type_direction_count": sum(1 for r in direction_records if r["triangle_type"]),
        "parallel_positive_direction_count": sum(
            1 for r in direction_records if r["parallel_type"] and r["size"] < 4
        ),
        "line_size_profiles": {k: sorted(v) for k, v in sorted(line_sizes.items())},
        "line_triangle_counts": line_triangle_counts,
        "line_parallel_positive_counts": line_parallel_positive_counts,
        "direction_records": direction_records,
    }


def parse_tuple_key_dict(raw: dict[str, Any]) -> dict[Elt, Elt]:
    return {
        ast.literal_eval(key): (value[0], value[1])
        for key, value in raw.items()
    }


def load_low_direction_section(path: Path) -> dict[Elt, Elt] | None:
    if not path.exists():
        return None
    data = json.loads(path.read_text())
    return parse_tuple_key_dict(data["section_h_coordinates"])


def load_sage_final_section(path: Path) -> dict[Elt, Elt] | None:
    if not path.exists():
        return None
    data = json.loads(path.read_text())
    table = data["final"]["table"]
    section: dict[Elt, Elt] = {}
    for x, row in enumerate(table):
        for y, value in enumerate(row):
            section[(x, y)] = (value[0], value[1])
    return section


def known_section_metrics(args: argparse.Namespace) -> list[dict[str, Any]]:
    sections: list[tuple[str, dict[Elt, Elt] | None]] = [
        (
            "p5_low_direction_counterexample_support_87",
            load_low_direction_section(args.low_direction_verify),
        ),
        (
            "p5_sage_heuristic_stress_support_85",
            load_sage_final_section(args.sage_heuristic_summary),
        ),
    ]
    return [section_metrics(section, label) for label, section in sections if section is not None]


def seven_pair_certificate(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {
            "path": str(path),
            "available": False,
            "usable": False,
            "interpretation": "No seven-low-pair classification result was found.",
        }

    data = json.loads(path.read_text())
    classification = data.get("classification", {})
    requested = classification.get("antipodal_pairs_requested")
    usable = (
        requested == 7
        and classification.get("any_representative_exists") is False
        and classification.get("any_node_limit_timeout") is False
        and classification.get("orbits_checked") == classification.get("orbit_count")
    )
    return {
        "path": str(path),
        "available": True,
        "usable": usable,
        "antipodal_pairs_requested": requested,
        "orbit_count": classification.get("orbit_count"),
        "orbits_checked": classification.get("orbits_checked"),
        "any_representative_exists": classification.get("any_representative_exists"),
        "any_node_limit_timeout": classification.get("any_node_limit_timeout"),
        "total_search_nodes": sum(
            item.get("search", {}).get("nodes", 0)
            for item in data.get("orbit_results", [])
        ),
        "interpretation": classification.get("interpretation"),
    }


def p5_box_consequence(certificate: dict[str, Any]) -> dict[str, Any]:
    consequence = {
        "counterexample_signed_deficit_threshold": 18,
        "counterexample_requires_at_least_low_antipodal_pairs": 7,
        "reason": (
            "After C-WS0001-018, all parallel-type positive-deficit fibres in a "
            "counterexample must lie on one projective line, hence in at most "
            "two antipodal pairs.  With at most six low antipodal pairs total, "
            "the maximum signed deficit is 2 size-2 pairs contributing 8 plus "
            "4 remaining size-3 pairs contributing 8, namely 16.  Therefore a "
            "p=5 counterexample needs at least seven antipodal pairs with "
            "|A_u|<=3."
        ),
        "certificate_rules_out_required_low_pair_count": bool(certificate.get("usable")),
        "conditional_conclusion": "pending-certificate-review",
    }
    if certificate.get("usable"):
        consequence["conditional_conclusion"] = (
            "No p=5 counterexample to the box bound exists, conditional on the "
            "seven-low-pair classification code and its GL2 reduction being "
            "accepted on review."
        )
    return consequence


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--json-out", type=Path, required=True)
    parser.add_argument(
        "--low-direction-verify",
        type=Path,
        default=Path("experiments/results/ws0001-p5-low-direction-counterexample-verify.json"),
    )
    parser.add_argument(
        "--sage-heuristic-summary",
        type=Path,
        default=Path("experiments/results/ws0001-p5-sage-heuristic-summary.json"),
    )
    parser.add_argument(
        "--seven-pair-classification",
        type=Path,
        default=Path("experiments/results/ws0001-p5-seven-low-pair-classification.json"),
    )
    args = parser.parse_args()

    template_records = [enumerate_template_orbits(template) for template in minimal_templates()]
    seven_pair = seven_pair_certificate(args.seven_pair_classification)
    result = {
        "script": "experiments/python/ws0001_p5_weighted_triangle_profile.py",
        "case": {
            "prime": P,
            "group": "(Z/25Z)^2",
            "subgroup": "5G",
            "quotient": "C5 x C5",
            "oriented_nonzero_direction_count": len(NONZERO),
            "antipodal_pair_count": len(PAIR_REPS),
            "projective_line_count": len(PROJECTIVE_REPS),
            "box_diffset_size": (2 * P - 1) ** 2,
            "box_fibre_sum": (2 * P - 1) ** 2 - 1,
            "baseline_size_4_fibre_sum": 4 * len(NONZERO),
            "counterexample_signed_deficit_threshold": 18,
        },
        "reviewed_reduction_used": {
            "parallel_type_subcase": "C-WS0001-018",
            "triangle_count_criterion": "C-WS0001-019",
            "assumption_for_remaining_counterexamples": (
                "all collinear positive-deficit fibres are confined to one "
                "projective line; all other positive-deficit fibres are "
                "triangle-type size-3 fibres"
            ),
            "minimal_deficit_equation": (
                "4*p2_pairs + 2*p3_parallel_pairs + 2*triangle_pairs = 18"
            ),
        },
        "minimal_dangerous_templates": template_records,
        "seven_low_pair_certificate": seven_pair,
        "p5_box_consequence": p5_box_consequence(seven_pair),
        "known_stress_sections": known_section_metrics(args),
        "interpretation": (
            "The profile enumeration is necessary-only by itself.  When the "
            "seven-low-pair certificate is usable, it discharges all minimal "
            "dangerous profiles and gives a computation-dependent p=5 box "
            "optimality conclusion pending proof/code review."
        ),
    }

    args.json_out.parent.mkdir(parents=True, exist_ok=True)
    args.json_out.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(
        json.dumps(
            {
                "templates": [
                    {
                        "template": item["template"],
                        "gl2_orbit_count": item["gl2_orbit_count"],
                    }
                    for item in template_records
                ],
                "known_sections": [
                    {
                        "label": item["label"],
                        "diffset_size": item["diffset_size"],
                        "signed_deficit_from_4": item["signed_deficit_from_4"],
                        "profile": item["profile"],
                        "triangle_type_direction_count": item["triangle_type_direction_count"],
                        "parallel_positive_direction_count": item["parallel_positive_direction_count"],
                    }
                    for item in result["known_stress_sections"]
                ],
                "p5_box_consequence": result["p5_box_consequence"],
                "seven_low_pair_certificate_usable": seven_pair["usable"],
            },
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
