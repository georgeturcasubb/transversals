# Section 6 computation checks

This repository contains the Python scripts and recorded JSON artifacts for the
small-prime computational claims in Section 6 of the paper
`Transversal Difference Numbers in Finite Abelian Quotients`.

The computations support the two claims stated in the paper:

- for `p = 3`, the square-plane value is `delta((Z/9Z)^2, 3(Z/9Z)^2) = 25`;
- for `p = 5`, the square-plane value is `delta((Z/25Z)^2, 5(Z/25Z)^2) = 81`.


## What is included

The repository is intentionally small.  It should contain only:

```text
experiments/
  python/
    ws0001_quotient_fibre_sharpness_check.py
    ws0001_c3_plane_counterexample.py
    ws0001_c3_plane_certificate_verify.py
    ws0001_p5_low_direction_classification.py
    ws0001_p5_weighted_triangle_profile.py
  results/
    ws0001-c3-plane-counterexample.json
    ws0001-c3-plane-certificate-verify.json
    ws0001-p5-seven-low-pair-classification.json
    ws0001-p5-weighted-triangle-profile.json
    ws0001-p5-low-direction-counterexample-verify.json
    ws0001-p5-sage-heuristic-summary.json
```

The two `p=5` files with `counterexample-verify` and `sage-heuristic` in their
names are diagnostic inputs used by the weighted-profile script.  They are not
the main proof certificates, but they are included so the recorded
weighted-profile artifact can be reproduced exactly.

## Requirements

Use Python 3.  The scripts use only the Python standard library; no Sage,
Magma, PARI/GP, or Python package installation is needed for these reviewer
checks.

The commands below should be run from the repository root.

## Quick integrity check

First check that Python can parse the scripts:

```bash
python3 -m py_compile \
  experiments/python/ws0001_quotient_fibre_sharpness_check.py \
  experiments/python/ws0001_c3_plane_counterexample.py \
  experiments/python/ws0001_c3_plane_certificate_verify.py \
  experiments/python/ws0001_p5_low_direction_classification.py \
  experiments/python/ws0001_p5_weighted_triangle_profile.py
```

## Recheck the `p = 3` certificate

This regenerates the `p=3` search output and then verifies it with a separate
verifier script.

```bash
python3 experiments/python/ws0001_c3_plane_counterexample.py \
  --json-out experiments/results/ws0001-c3-plane-counterexample.json

python3 experiments/python/ws0001_c3_plane_certificate_verify.py \
  --canonical-result experiments/results/ws0001-c3-plane-counterexample.json \
  --json-out experiments/results/ws0001-c3-plane-certificate-verify.json
```

The reviewer should look in
`experiments/results/ws0001-c3-plane-certificate-verify.json` for:

```json
{
  "verified_delta": 25,
  "limit_23_absent": true,
  "limit_25_present": true,
  "odd_cardinality_excludes_24": true
}
```

This is the computation behind the `p=3` value in Section 6.

## Recheck the `p = 5` certificate

First rerun the seven-low-pair classification:

```bash
python3 experiments/python/ws0001_p5_low_direction_classification.py \
  --pairs 7 \
  --json-out experiments/results/ws0001-p5-seven-low-pair-classification.json
```

The reviewer should look in
`experiments/results/ws0001-p5-seven-low-pair-classification.json` for:

```json
{
  "antipodal_pairs_requested": 7,
  "orbit_count": 7,
  "orbits_checked": 7,
  "any_representative_exists": false,
  "any_node_limit_timeout": false
}
```

Then rerun the weighted-profile summary, using the seven-low-pair file just
created:

```bash
python3 experiments/python/ws0001_p5_weighted_triangle_profile.py \
  --seven-pair-classification experiments/results/ws0001-p5-seven-low-pair-classification.json \
  --json-out experiments/results/ws0001-p5-weighted-triangle-profile.json
```

The reviewer should look in
`experiments/results/ws0001-p5-weighted-triangle-profile.json` for:

```json
{
  "seven_low_pair_certificate": {
    "usable": true,
    "orbits_checked": 7,
    "any_representative_exists": false,
    "any_node_limit_timeout": false
  }
}
```

and for the `p5_box_consequence` block, which explains how the seven-low-pair
certificate rules out the remaining dangerous `p=5` profiles.

## A note on timing

The `p=3` commands are small.  The `p=5` classification is more substantial but
still intended to be a finite foreground Python check on an ordinary machine.
The recorded run checked seven projective low-pair orbits and used 25017 search
nodes in total.

## Provenance

This repository was prepared with the aid of Codex and model GPT-5.5.  The
mathematical responsibility for the claims, scripts, and interpretation remains
with the paper authors.
