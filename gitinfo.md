# Git Info

## Purpose

This file is the version ledger for `cifplane`.

It should answer four questions for each version:

1. What problem does this version solve?
2. What algorithmic rule does it use?
3. What are the input and output contracts?
4. How was the version validated?

## Repository Contract

- Repository name: `cifplane`
- Core goal: enumerate atomic planes and representative terminations for a given `cif` and `hkl`
- Versioned content: source code, environment definition, example CIF inputs, and milestone documentation
- Non-versioned content: generated outputs such as `results_plane/`, cache files, and local visualization artifacts

## Versioning Rule

- `v1`, `v2`, ... are milestone snapshots
- New versions must be appended, not rewritten over older entries
- If the termination deduplication logic changes, the change must be stated explicitly in the version entry

## Entry Template

Use the following fields when adding a new version:

- Version
- Date
- Status
- Objective
- Algorithm Summary
- Equivalent-Termination Rule
- Inputs
- Outputs
- Main Files
- Validation
- Known Limits
- Notes for Next Version

## Version Ledger

### Version: `v1`

**Date**

`2026-04-15`

**Status**

Initial standalone Git snapshot and first public repository baseline.

**Objective**

Create a reproducible workflow that accepts one CIF and one Miller index, then reports the atomic planes and representative terminations associated with the corresponding oriented slab repeat.

**Algorithm Summary**

1. Parse the CIF and standardize the structure to a conventional bulk cell.
2. Build an oriented unit cell using `pymatgen.core.surface.SlabGenerator`.
3. Cluster atomic sites into unique atomic planes along the slab normal.
4. Canonicalize each plane motif in in-plane fractional coordinates.
5. Build cyclic plane-stack signatures from every candidate `top_plane_id`.
6. Keep one representative for each unique cyclic stack signature.
7. Export CSV, JSON, markdown, and optional surface-only visualization artifacts.

**Equivalent-Termination Rule**

- Duplicate removal is based on the cyclic sequence of plane motifs in one oriented bulk repeat.
- The workflow keeps the first encountered `top_plane_id` whose full cyclic stack signature has not appeared before.
- This rule removes periodic duplicates inside the repeat stack.
- This rule does not perform a full symmetry-equivalence reduction over all surface rotations, mirrors, or other in-plane operations.

**Inputs**

- CIF file path
- Miller index `hkl`
- CLI parameters such as `symprec`, `plane_tol_ang`, `motif_round_digits`, slab thickness, and vacuum thickness

**Outputs**

- `plane_stack.csv`
- `terminations.json`
- `report.md`
- optional `surface_only/` CIF and PNG artifacts when `--surface-only-export` is enabled

**Main Files**

- `analyze_planes.py`
- `export_vesta_pngs.py`
- `run_surface_pipeline.sh`
- `environment.yml`
- `cif/`

**Validation**

- Smoke test command:

```bash
python analyze_planes.py \
  --cif cif/SrTiO3.cif \
  --hkl 1 0 0 \
  --output-dir /tmp/cifplane_smoke \
  --termination-only
```

- Expected pass criteria:
  - execution completes without error
  - `/tmp/cifplane_smoke/SrTiO3/hkl_1_0_0/plane_stack.csv` exists
  - `/tmp/cifplane_smoke/SrTiO3/hkl_1_0_0/terminations.json` exists
  - `/tmp/cifplane_smoke/SrTiO3/hkl_1_0_0/report.md` exists
  - `terminations.json` contains one representative per unique cyclic stack signature

**Known Limits**

- Equivalent termination removal is sequence-based rather than full surface-symmetry-based.
- The batch shell script expects the user to edit `CIF_PATH` manually before running.
- VESTA PNG export assumes a local VESTA installation compatible with the helper script.
- Generated outputs are intentionally excluded from version control.

**Notes for Next Version**

- Add richer README examples with representative output snippets and screenshots.
- Consider making the batch pipeline configurable by CLI arguments instead of editing shell variables.
- Consider exposing duplicate-to-representative mappings directly in the JSON output.
