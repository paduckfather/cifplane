# cifplane

`cifplane` is a small surface-analysis workflow for enumerating atomic planes and representative terminations for a given crystal structure `cif` and Miller index `hkl`.

The current implementation is oriented toward inspection and reproducible reporting rather than high-throughput screening. It generates plane-level summaries, termination summaries, markdown reports, and optional surface-only visualization files.

## What This Repository Does

For one input CIF and one target `(hkl)`, the workflow:

1. standardizes the bulk structure to a conventional cell
2. constructs an oriented unit cell normal to the target plane
3. clusters atoms into atomic planes along the slab normal
4. builds a canonical in-plane motif signature for each plane
5. enumerates representative terminations from the cyclic plane stack
6. exports tables, JSON summaries, markdown reports, and optional visualization artifacts

## Equivalent-Termination Rule

The repository currently removes duplicate terminations using a periodic stack rule:

- each candidate termination is represented by the cyclic sequence of plane motifs starting from one `top_plane_id`
- if two candidates produce exactly the same cyclic sequence, only the first encountered one is kept
- this is a periodic stack deduplication rule inside one oriented bulk repeat
- this is not a full symmetry-equivalence reduction over all in-plane rotations or mirrors

## Repository Layout

- `analyze_planes.py`: main CLI for plane construction, termination enumeration, and report export
- `export_vesta_pngs.py`: helper for exporting VESTA-based PNG files from generated surface-only CIFs
- `run_surface_pipeline.sh`: batch-style shell entry point for a fixed CIF and a list of hkls
- `environment.yml`: reproducible conda environment
- `cif/`: example CIF inputs
- `gitinfo.md`: version ledger and milestone notes

## Environment

The repository is designed around a conda environment.

```bash
conda env create -f environment.yml
conda activate cifplane
```

## CLI Usage

Single-direction analysis:

```bash
python analyze_planes.py \
  --cif cif/SrTiO3.cif \
  --hkl 1 0 0 \
  --output-dir results_plane \
  --termination-only
```

Surface-only exports for visualization:

```bash
python analyze_planes.py \
  --cif cif/LaMnO3.cif \
  --hkl 1 1 1 \
  --output-dir results_plane \
  --termination-only \
  --surface-only-export
```

Batch-style shell entry point:

```bash
bash run_surface_pipeline.sh
```

Before running `run_surface_pipeline.sh`, edit the `CIF_PATH` near the top of the script to point to the desired input CIF.

## Main Outputs

Results are written under:

```text
<output-dir>/<cif-stem>/hkl_<h>_<k>_<l>/
```

Main files:

- `plane_stack.csv`: plane-by-plane summary inside one oriented repeat
- `terminations.json`: representative termination records with motif and cut metadata
- `report.md`: compact markdown report

Optional files when `--surface-only-export` is enabled:

- `surface_only/termination_*_top.cif`
- `surface_only/termination_*_top.png`
- `surface_only/termination_*_top_vesta_view.cif`
- `surface_only/open_all_in_vesta.command`

## Representative Smoke Test

```bash
python analyze_planes.py \
  --cif cif/SrTiO3.cif \
  --hkl 1 0 0 \
  --output-dir /tmp/cifplane_smoke \
  --termination-only
```

Expected outputs:

- `/tmp/cifplane_smoke/SrTiO3/hkl_1_0_0/plane_stack.csv`
- `/tmp/cifplane_smoke/SrTiO3/hkl_1_0_0/terminations.json`
- `/tmp/cifplane_smoke/SrTiO3/hkl_1_0_0/report.md`

## Current Scope and Limits

- the workflow is designed for inspection of a single structure and selected hkls
- equivalent termination removal is sequence-based, not full surface-symmetry-based
- generated outputs under `results_plane/` are not versioned by default
- VESTA PNG export depends on the local VESTA installation path expected by the helper script

## Version History

See [gitinfo.md](./gitinfo.md) for milestone-level version notes.
