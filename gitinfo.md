# Git Info

## Repository Scope

- Repository name: `cifplane`
- Goal: enumerate atomic planes and representative surface terminations for a given `cif` and `hkl`
- Tracking policy: source code and minimal example CIF inputs are versioned; generated outputs are not versioned and should be reproduced from the scripts

## Versioning Policy

- `v1`, `v2`, ... are milestone snapshots of the workflow
- Each version entry should record:
  - problem scope
  - algorithm-level behavior
  - input/output contract
  - validation status

## Versions

### v1

- Date: `2026-04-15`
- Status: initial standalone Git snapshot
- Scope:
  - standardize a CIF to a conventional bulk cell
  - build an oriented unit cell for a target `(hkl)`
  - cluster atoms along the slab normal into atomic planes
  - canonicalize each plane motif in in-plane fractional coordinates
  - enumerate representative terminations from the cyclic plane stack
  - export summary tables and optional surface-only visualization files
- Equivalent-termination rule:
  - duplicate removal is based on the cyclic sequence of plane motifs in one oriented bulk repeat
  - the code keeps the first encountered `top_plane_id` whose full cyclic stack signature is unique
  - this is not a full symmetry-equivalence reduction over surface rotations or mirrors
- Main files:
  - `analyze_planes.py`: main workflow for plane construction, termination enumeration, reporting, and optional surface export
  - `export_vesta_pngs.py`: helper for VESTA-oriented PNG export
  - `run_surface_pipeline.sh`: shell entry point for batch-style execution
  - `environment.yml`: reproducible package environment
  - `cif/`: example CIF inputs
- Inputs:
  - CIF file path
  - Miller index `hkl`
  - tolerances and slab-generation options from CLI flags
- Outputs:
  - `plane_stack.csv`
  - `terminations.json`
  - `report.md`
  - optional surface-only CIF/PNG artifacts under `results_plane/`
- Validation:
  - smoke-test on a small example such as `SrTiO3 (100)` should complete without error
  - generated `terminations.json` should contain one representative per unique cyclic stack signature

## Update Rule

- When a new version is cut, append a new section instead of rewriting the old one
- If the algorithm changes, explicitly state whether the termination deduplication rule changed
