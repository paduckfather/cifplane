#!/usr/bin/env python3
"""Enumerate unique atomic planes and terminations for a given (hkl)."""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import re
import tempfile
from collections import Counter
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
from typing import Any

import numpy as np
import spglib
from ase.data import atomic_numbers, covalent_radii
from ase.data.colors import jmol_colors
from pymatgen.core import Lattice, Structure
from pymatgen.core.surface import SlabGenerator
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer


DEFAULT_VESTA_ELEMENTS_INI = Path("/Applications/VESTA.app/Contents/Resources/elements.ini")


@dataclass(frozen=True)
class PlaneRecord:
    """One unique atomic plane in the oriented bulk repeat cell."""

    plane_id: int
    z_frac: float
    z_ang: float
    species_counts: tuple[tuple[str, int], ...]
    motif_signature: tuple[tuple[str, float, float], ...]
    site_indices: tuple[int, ...]
    motif_sites: tuple[tuple[str, float, float], ...]


def parse_args() -> argparse.Namespace:
    """Parse CLI arguments for a single-hkl analysis."""
    parser = argparse.ArgumentParser(
        description="Enumerate unique planes and possible terminations for one (hkl).",
    )
    parser.add_argument("--cif", required=True, help="Input CIF path.")
    parser.add_argument(
        "--hkl",
        required=True,
        nargs=3,
        type=int,
        metavar=("H", "K", "L"),
        help="Target Miller index.",
    )
    parser.add_argument(
        "--output-dir",
        default="results_plane",
        help="Base output directory. Results are stored under <output-dir>/<cif>/<hkl>/.",
    )
    parser.add_argument(
        "--symprec",
        type=float,
        default=1e-3,
        help="Symmetry tolerance passed to pymatgen/spglib.",
    )
    parser.add_argument(
        "--plane-tol-ang",
        type=float,
        default=0.05,
        help="Tolerance in angstrom for clustering atomic planes along the normal.",
    )
    parser.add_argument(
        "--motif-round-digits",
        type=int,
        default=8,
        help="Rounding digits for canonical motif signatures.",
    )
    parser.add_argument(
        "--min-slab-size",
        type=float,
        default=12.0,
        help="Representative slab thickness in angstrom.",
    )
    parser.add_argument(
        "--min-vacuum-size",
        type=float,
        default=15.0,
        help="Representative slab vacuum thickness in angstrom.",
    )
    parser.add_argument(
        "--termination-only",
        action="store_true",
        help="Enumerate/report terminations without exporting representative slab CIF files.",
    )
    parser.add_argument(
        "--surface-only-export",
        action="store_true",
        help="Export top-plane-only CIF files for VESTA visualization.",
    )
    parser.add_argument(
        "--surface-only-vacuum-size",
        type=float,
        default=20.0,
        help="Vacuum thickness in angstrom for top-plane-only CIF export.",
    )
    parser.add_argument(
        "--png-boundary-xmax",
        type=int,
        default=3,
        help="Number of in-plane repeats along the x lattice direction for PNG export.",
    )
    parser.add_argument(
        "--png-boundary-ymax",
        type=int,
        default=3,
        help="Number of in-plane repeats along the y lattice direction for PNG export.",
    )
    return parser.parse_args()


def wrap01(value: float, digits: int = 12) -> float:
    """Wrap a scalar into [0, 1) with stable rounding."""
    wrapped = float(value % 1.0)
    if math.isclose(wrapped, 1.0, abs_tol=10 ** (-digits)):
        return 0.0
    if math.isclose(wrapped, 0.0, abs_tol=10 ** (-digits)):
        return 0.0
    return wrapped


def wrap_vector(values: np.ndarray, digits: int = 12) -> np.ndarray:
    """Wrap a vector into [0, 1) component-wise."""
    wrapped = np.mod(values, 1.0)
    wrapped[np.isclose(wrapped, 1.0, atol=10 ** (-digits))] = 0.0
    wrapped[np.isclose(wrapped, 0.0, atol=10 ** (-digits))] = 0.0
    return wrapped


def cluster_periodic_positions(values: np.ndarray, tol_frac: float) -> list[dict[str, Any]]:
    """Cluster positions on a unit circle, merging the wrap-around cluster when needed."""
    wrapped = np.mod(np.asarray(values, dtype=float), 1.0)
    order = np.argsort(wrapped)
    sorted_values = wrapped[order]

    if len(sorted_values) == 0:
        return []

    clusters: list[dict[str, Any]] = [
        {"members": [int(order[0])], "values": [float(sorted_values[0])]},
    ]
    for idx, value in zip(order[1:], sorted_values[1:], strict=True):
        if value - clusters[-1]["values"][-1] <= tol_frac:
            clusters[-1]["members"].append(int(idx))
            clusters[-1]["values"].append(float(value))
        else:
            clusters.append({"members": [int(idx)], "values": [float(value)]})

    if len(clusters) > 1 and (clusters[0]["values"][0] + 1.0 - clusters[-1]["values"][-1] <= tol_frac):
        merged_values = clusters[-1]["values"] + [value + 1.0 for value in clusters[0]["values"]]
        merged_members = clusters[-1]["members"] + clusters[0]["members"]
        clusters = [{"members": merged_members, "values": merged_values}] + clusters[1:-1]

    result: list[dict[str, Any]] = []
    for cluster in clusters:
        angles = np.array(cluster["values"], dtype=float) * 2.0 * np.pi
        center = math.atan2(np.sin(angles).mean(), np.cos(angles).mean()) / (2.0 * np.pi)
        center = wrap01(center)
        result.append(
            {
                "members": sorted(cluster["members"]),
                "center_frac": center,
            },
        )
    result.sort(key=lambda item: item["center_frac"])
    return result


def cluster_linear_positions(values: np.ndarray, tol: float) -> list[dict[str, Any]]:
    """Cluster scalar positions on a line without periodic wrapping."""
    arr = np.asarray(values, dtype=float)
    order = np.argsort(arr)
    sorted_values = arr[order]
    if len(sorted_values) == 0:
        return []

    clusters: list[dict[str, Any]] = [
        {"members": [int(order[0])], "values": [float(sorted_values[0])]},
    ]
    for idx, value in zip(order[1:], sorted_values[1:], strict=True):
        if value - clusters[-1]["values"][-1] <= tol:
            clusters[-1]["members"].append(int(idx))
            clusters[-1]["values"].append(float(value))
        else:
            clusters.append({"members": [int(idx)], "values": [float(value)]})

    result: list[dict[str, Any]] = []
    for cluster in clusters:
        result.append(
            {
                "members": sorted(cluster["members"]),
                "center": float(np.mean(cluster["values"])),
            },
        )
    result.sort(key=lambda item: item["center"])
    return result


def canonicalize_planar_motif(
    species: list[str],
    xy_frac: np.ndarray,
    digits: int,
) -> tuple[tuple[str, float, float], ...]:
    """Make a translation-invariant canonical motif signature in fractional in-plane coordinates."""
    xy = wrap_vector(np.asarray(xy_frac, dtype=float)[:, :2])
    best: tuple[tuple[str, float, float], ...] | None = None

    for anchor in xy:
        shifted = wrap_vector(xy - anchor)
        candidate = tuple(
            sorted(
                (
                    species[idx],
                    round(float(shifted[idx, 0]), digits),
                    round(float(shifted[idx, 1]), digits),
                )
                for idx in range(len(species))
            ),
        )
        if best is None or candidate < best:
            best = candidate

    return best or tuple()


def counter_to_tuple(counter: Counter[str]) -> tuple[tuple[str, int], ...]:
    """Serialize a Counter deterministically."""
    return tuple(sorted((str(key), int(value)) for key, value in counter.items()))


def tuple_to_formula_string(items: tuple[tuple[str, int], ...]) -> str:
    """Compact species counter as a human-readable formula-like string."""
    return " ".join(f"{symbol}{count}" for symbol, count in items)


def canonicalize_species_symbol(species_like: Any) -> str:
    """Normalize pymatgen species labels to plain element symbols for lookup/reporting."""
    element = getattr(species_like, "element", None)
    if element is not None and getattr(element, "symbol", None) is not None:
        return str(element.symbol)

    symbol = getattr(species_like, "symbol", None)
    if symbol is not None:
        return str(symbol)

    match = re.match(r"[A-Z][a-z]?", str(species_like))
    if match is not None:
        return match.group(0)

    raise ValueError(f"Failed to canonicalize species label: {species_like!r}")


def get_standard_conventional_structure(cif_path: Path, symprec: float) -> tuple[Structure, str, int]:
    """Recover symmetry and standardize the CIF to a conventional cell."""
    raw_structure = Structure.from_file(cif_path)
    analyzer = SpacegroupAnalyzer(raw_structure, symprec=symprec)
    conventional = analyzer.get_conventional_standard_structure()
    return conventional, analyzer.get_space_group_symbol(), analyzer.get_space_group_number()


def build_plane_records(
    oriented_unit_cell: Structure,
    plane_tol_ang: float,
    motif_round_digits: int,
) -> list[PlaneRecord]:
    """Enumerate unique atomic planes in one oriented bulk repeat cell."""
    c_repeat = oriented_unit_cell.lattice.c
    tol_frac = plane_tol_ang / c_repeat
    z_frac = wrap_vector(oriented_unit_cell.frac_coords[:, 2:3]).ravel()
    plane_clusters = cluster_periodic_positions(z_frac, tol_frac)

    records: list[PlaneRecord] = []
    for plane_id, cluster in enumerate(plane_clusters):
        indices = cluster["members"]
        species = [canonicalize_species_symbol(oriented_unit_cell[idx].specie) for idx in indices]
        xy = oriented_unit_cell.frac_coords[indices, :2]
        motif_signature = canonicalize_planar_motif(species, xy, digits=motif_round_digits)
        motif_sites = tuple(
            sorted(
                (
                    species[pos],
                    round(wrap01(float(xy[pos, 0])), motif_round_digits),
                    round(wrap01(float(xy[pos, 1])), motif_round_digits),
                )
                for pos in range(len(indices))
            ),
        )
        records.append(
            PlaneRecord(
                plane_id=plane_id,
                z_frac=float(cluster["center_frac"]),
                z_ang=float(cluster["center_frac"] * c_repeat),
                species_counts=counter_to_tuple(Counter(species)),
                motif_signature=motif_signature,
                site_indices=tuple(indices),
                motif_sites=motif_sites,
            ),
        )
    records.sort(key=lambda item: item.z_frac)
    return records


def midpoint_on_unit_interval(left: float, right: float) -> float:
    """Return the midpoint between two periodic positions on [0, 1)."""
    delta = (right - left) % 1.0
    return wrap01(left + 0.5 * delta)


def cyclic_rotated_signatures(
    plane_records: list[PlaneRecord],
    start_index: int,
) -> tuple[tuple[tuple[str, float, float], ...], ...]:
    """Return the cyclic plane-signature sequence starting at a given top plane."""
    signatures = [plane.motif_signature for plane in plane_records]
    return tuple(signatures[start_index:] + signatures[:start_index])


def reverse_cyclic_rotated_signatures(
    plane_records: list[PlaneRecord],
    bottom_index: int,
) -> tuple[tuple[tuple[str, float, float], ...], ...]:
    """Return the reversed cyclic signature sequence for the opposite normal."""
    n_planes = len(plane_records)
    ordered = []
    for offset in range(n_planes):
        ordered.append(plane_records[(bottom_index - offset) % n_planes].motif_signature)
    return tuple(ordered)


def lattice_rows_to_spglib_cell(structure: Structure) -> tuple[np.ndarray, np.ndarray, list[int]]:
    """Convert a pymatgen Structure to the spglib cell tuple."""
    return (
        np.array(structure.lattice.matrix, dtype=float),
        np.array(structure.frac_coords, dtype=float),
        list(map(int, structure.atomic_numbers)),
    )


def project_unique_2d_operations(dataset: Any, digits: int = 8) -> list[tuple[np.ndarray, np.ndarray]]:
    """Project layer-group operations to unique 2D operations."""
    seen: set[tuple[tuple[int, ...], tuple[float, ...]]] = set()
    operations: list[tuple[np.ndarray, np.ndarray]] = []
    for rotation, translation in zip(dataset.rotations, dataset.translations, strict=True):
        rot2 = np.array(rotation[:2, :2], dtype=int)
        trans2 = wrap_vector(np.array(translation[:2], dtype=float), digits=digits)
        key = (
            tuple(int(value) for value in rot2.ravel()),
            tuple(round(float(value), digits) for value in trans2),
        )
        if key not in seen:
            seen.add(key)
            operations.append((rot2, trans2))
    return operations


def matrix_order_2d(rotation: np.ndarray) -> int:
    """Return the smallest positive n with R^n = I for a 2D integer rotation matrix."""
    current = np.eye(2, dtype=int)
    identity = np.eye(2, dtype=int)
    for order in range(1, 7):
        current = current @ rotation
        if np.array_equal(current, identity):
            return order
    return 0


def fractional_linear_solution_exists(matrix: np.ndarray, rhs: np.ndarray, tol: float = 1e-8) -> bool:
    """Check solvability of Mx = rhs + integer_shift over reals for small integer shifts."""
    for n1 in (-1, 0, 1):
        for n2 in (-1, 0, 1):
            target = rhs + np.array([n1, n2], dtype=float)
            solution, residuals, _, _ = np.linalg.lstsq(matrix, target, rcond=None)
            error = np.linalg.norm(matrix @ solution - target)
            if error <= tol:
                return True
    return False


def is_mirror_operation(rotation: np.ndarray, translation: np.ndarray, tol: float = 1e-8) -> bool:
    """Classify det=-1 operations into mirror or glide by fixed-line solvability."""
    system = np.eye(2, dtype=float) - rotation.astype(float)
    return fractional_linear_solution_exists(system, translation, tol=tol)


def lattice_family_2d(lattice_matrix: np.ndarray, atol: float = 1e-6) -> str:
    """Classify the in-plane Bravais family from the 2D metric."""
    a_vec = np.array(lattice_matrix[0], dtype=float)
    b_vec = np.array(lattice_matrix[1], dtype=float)
    a = float(np.linalg.norm(a_vec))
    b = float(np.linalg.norm(b_vec))
    cos_gamma = float(np.dot(a_vec, b_vec) / (a * b))
    cos_gamma = max(-1.0, min(1.0, cos_gamma))
    gamma = math.degrees(math.acos(cos_gamma))

    if math.isclose(a, b, rel_tol=1e-6, abs_tol=atol) and math.isclose(gamma, 90.0, abs_tol=1e-4):
        return "square"
    if math.isclose(a, b, rel_tol=1e-6, abs_tol=atol) and (
        math.isclose(gamma, 60.0, abs_tol=1e-4) or math.isclose(gamma, 120.0, abs_tol=1e-4)
    ):
        return "hexagonal"
    if math.isclose(gamma, 90.0, abs_tol=1e-4):
        return "rectangular"
    return "oblique"


def cartesian_rotation(rotation_frac: np.ndarray, inplane_lattice_rows: np.ndarray) -> np.ndarray:
    """Convert a fractional-basis 2D operation to Cartesian coordinates."""
    basis = np.array(inplane_lattice_rows, dtype=float).T
    return basis @ rotation_frac @ np.linalg.inv(basis)


def mirror_axis_angles(
    operations: list[tuple[np.ndarray, np.ndarray]],
    inplane_lattice_rows: np.ndarray,
) -> list[float]:
    """Return mirror-line angles in degrees modulo 180."""
    angles = []
    for rotation, translation in operations:
        if round(np.linalg.det(rotation)) != -1 or not is_mirror_operation(rotation, translation):
            continue
        rotation_cart = cartesian_rotation(rotation, inplane_lattice_rows)
        eigenvalues, eigenvectors = np.linalg.eig(rotation_cart)
        mirror_idx = int(np.argmin(np.abs(eigenvalues - 1.0)))
        direction = np.real(eigenvectors[:, mirror_idx])
        direction /= np.linalg.norm(direction)
        angle = math.degrees(math.atan2(direction[1], direction[0])) % 180.0
        angles.append(angle)
    return angles


def angle_distance_mod_180(a: float, b: float) -> float:
    """Distance between two line angles modulo 180 degrees."""
    delta = abs(a - b) % 180.0
    return min(delta, 180.0 - delta)


def classify_hexagonal_mirror_group(
    operations: list[tuple[np.ndarray, np.ndarray]],
    inplane_lattice_rows: np.ndarray,
) -> str:
    """Disambiguate p3m1 vs p31m from mirror-line orientation."""
    a_vec = np.array(inplane_lattice_rows[0], dtype=float)
    b_vec = np.array(inplane_lattice_rows[1], dtype=float)
    basis_directions = [
        math.degrees(math.atan2(vec[1], vec[0])) % 180.0
        for vec in (a_vec, b_vec, a_vec - b_vec)
    ]
    bisector_directions = [
        math.degrees(math.atan2(vec[1], vec[0])) % 180.0
        for vec in (a_vec + b_vec, 2.0 * a_vec - b_vec, a_vec - 2.0 * b_vec)
    ]
    mirror_angles = mirror_axis_angles(operations, inplane_lattice_rows)

    if not mirror_angles:
        return "p3"

    basis_score = min(
        angle_distance_mod_180(angle, candidate)
        for angle in mirror_angles
        for candidate in basis_directions
    )
    bisector_score = min(
        angle_distance_mod_180(angle, candidate)
        for angle in mirror_angles
        for candidate in bisector_directions
    )
    return "p3m1" if basis_score <= bisector_score else "p31m"


def classify_plane_group(dataset: Any, inplane_lattice_rows: np.ndarray) -> str | None:
    """Infer the 17 wallpaper-group symbol from layer-group operations."""
    if dataset is None:
        return None

    operations = project_unique_2d_operations(dataset)
    proper_orders = [
        matrix_order_2d(rotation)
        for rotation, _ in operations
        if round(np.linalg.det(rotation)) == 1
    ]
    max_order = max(proper_orders) if proper_orders else 1

    mirror_count = 0
    glide_count = 0
    for rotation, translation in operations:
        if round(np.linalg.det(rotation)) != -1:
            continue
        if is_mirror_operation(rotation, translation):
            mirror_count += 1
        else:
            glide_count += 1

    has_mirror = mirror_count > 0
    has_glide = glide_count > 0
    family = lattice_family_2d(inplane_lattice_rows)
    centered_rectangular = family == "oblique" and (has_mirror or has_glide)

    if max_order == 6:
        return "p6mm" if has_mirror or has_glide else "p6"
    if max_order == 4:
        if not (has_mirror or has_glide):
            return "p4"
        return "p4gm" if has_glide else "p4mm"
    if max_order == 3:
        if not (has_mirror or has_glide):
            return "p3"
        return classify_hexagonal_mirror_group(operations, inplane_lattice_rows)
    if max_order == 2:
        if not (has_mirror or has_glide):
            return "p2"
        if has_mirror and not has_glide:
            return "c2mm" if centered_rectangular else "p2mm"
        if has_glide and not has_mirror:
            return "p2gg"
        return "p2mg"
    if max_order == 1:
        if not (has_mirror or has_glide):
            return "p1"
        if has_glide and not has_mirror:
            return "pg"
        if has_mirror and not has_glide:
            return "cm" if centered_rectangular else "pm"
        return "cm" if centered_rectangular else "pg"
    return None


def build_top_plane_only_structure(
    oriented_unit_cell: Structure,
    top_plane: PlaneRecord,
    vacuum_size: float,
) -> Structure:
    """Build a flat visualization cell containing only the algorithm-defined top termination plane."""
    lattice_matrix = np.array(oriented_unit_cell.lattice.matrix, dtype=float)
    a_vec = lattice_matrix[0]
    b_vec = lattice_matrix[1]
    normal = np.cross(a_vec, b_vec)
    normal_norm = np.linalg.norm(normal)
    if normal_norm == 0:
        raise ValueError("Failed to build a non-singular surface-only cell.")
    lattice_matrix[2] = normal / normal_norm * vacuum_size

    species = []
    frac_coords = []
    for species_label, x_frac, y_frac in top_plane.motif_sites:
        frac_coords.append(
            [
                wrap01(float(x_frac)),
                wrap01(float(y_frac)),
                0.5,
            ],
        )
        species.append(species_label)

    return Structure(
        lattice=Lattice(lattice_matrix),
        species=species,
        coords=frac_coords,
        coords_are_cartesian=False,
        to_unit_cell=True,
    )


def build_vesta_view_structure(
    structure: Structure,
    repeat_x: int,
    repeat_y: int,
) -> Structure:
    """Expand the top-plane-only structure in-plane for immediate VESTA viewing."""
    view_structure = structure.copy()
    view_structure.make_supercell([[repeat_x, 0, 0], [0, repeat_y, 0], [0, 0, 1]])
    return view_structure


@lru_cache(maxsize=1)
def load_vesta_element_colors() -> dict[str, tuple[float, float, float]]:
    """Load VESTA default element RGB colors from elements.ini when available."""
    elements_ini = Path(os.environ.get("VESTA_ELEMENTS_INI", DEFAULT_VESTA_ELEMENTS_INI))
    if not elements_ini.exists():
        return {}

    colors: dict[str, tuple[float, float, float]] = {}
    for line in elements_ini.read_text(encoding="utf-8").splitlines():
        fields = line.split()
        if len(fields) < 8:
            continue
        try:
            symbol = fields[1]
            rgb = tuple(float(component) for component in fields[-3:])
        except ValueError:
            continue
        colors[symbol] = rgb
    return colors


def select_annotation_primitive_cell(
    prim_a_xy: np.ndarray,
    prim_b_xy: np.ndarray,
    min_xy: np.ndarray,
    max_xy: np.ndarray,
    u_min: int,
    u_max: int,
    v_min: int,
    v_max: int,
) -> tuple[int, int]:
    """Choose the fully visible primitive cell closest to the lower-left corner."""
    candidate_cells: list[tuple[float, float, int, int]] = []
    for u_index in range(u_min, u_max):
        for v_index in range(v_min, v_max):
            primitive_outline_candidate = np.array(
                [
                    u_index * prim_a_xy + v_index * prim_b_xy,
                    (u_index + 1) * prim_a_xy + v_index * prim_b_xy,
                    (u_index + 1) * prim_a_xy + (v_index + 1) * prim_b_xy,
                    u_index * prim_a_xy + (v_index + 1) * prim_b_xy,
                ],
                dtype=float,
            )
            if (
                primitive_outline_candidate[:, 0].min() < min_xy[0]
                or primitive_outline_candidate[:, 0].max() > max_xy[0]
                or primitive_outline_candidate[:, 1].min() < min_xy[1]
                or primitive_outline_candidate[:, 1].max() > max_xy[1]
            ):
                continue
            candidate_cells.append(
                (
                    float(primitive_outline_candidate[:, 0].min()),
                    float(primitive_outline_candidate[:, 1].min()),
                    u_index,
                    v_index,
                ),
            )

    if candidate_cells:
        _, _, axis_u, axis_v = min(candidate_cells)
        return axis_u, axis_v
    return 0, 0


def render_structure_top_view_png(
    structure: Structure,
    png_path: Path,
    boundary_xmax: int,
    boundary_ymax: int,
    plane_group: str | None = None,
) -> None:
    """Render a top-view PNG for one surface-only structure with primitive-cell annotations."""
    temp_cache_dir = str(Path(tempfile.gettempdir()) / "codex_plot_cache")
    os.environ.setdefault("MPLCONFIGDIR", temp_cache_dir)
    os.environ.setdefault("XDG_CACHE_HOME", temp_cache_dir)
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    from matplotlib.patches import Circle

    a_vec = np.array(structure.lattice.matrix[0], dtype=float)
    b_vec = np.array(structure.lattice.matrix[1], dtype=float)
    normal_vec = np.cross(a_vec, b_vec)
    normal_unit = normal_vec / np.linalg.norm(normal_vec)
    e1 = a_vec / np.linalg.norm(a_vec)
    b_proj = b_vec - np.dot(b_vec, e1) * e1
    if np.linalg.norm(b_proj) == 0:
        raise ValueError("In-plane lattice vectors are linearly dependent; cannot render PNG.")
    e2 = b_proj / np.linalg.norm(b_proj)

    fig, ax = plt.subplots(figsize=(8.6, 7.2), dpi=220)
    fig.subplots_adjust(right=0.78)

    primitive_structure = structure.get_primitive_structure()
    primitive_rows = np.array(primitive_structure.lattice.matrix, dtype=float)
    inplane_order = np.argsort(np.abs(primitive_rows @ normal_unit))[:2]
    prim_a_vec = primitive_rows[inplane_order[0]]
    prim_b_vec = primitive_rows[inplane_order[1]]
    prim_c_vec = primitive_rows[[idx for idx in range(3) if idx not in inplane_order][0]]
    prim_a_xy = np.array([np.dot(prim_a_vec, e1), np.dot(prim_a_vec, e2)], dtype=float)
    prim_b_xy = np.array([np.dot(prim_b_vec, e1), np.dot(prim_b_vec, e2)], dtype=float)
    primitive_basis = np.column_stack([prim_a_xy, prim_b_xy])
    if abs(np.linalg.det(primitive_basis)) < 1e-10:
        raise ValueError("Primitive in-plane lattice vectors are linearly dependent; cannot render PNG.")
    primitive_lattice = Lattice(np.array([prim_a_vec, prim_b_vec, prim_c_vec], dtype=float))

    atom_points: list[np.ndarray] = []
    legend_entries: dict[str, tuple[float, float, float]] = {}
    vesta_colors = load_vesta_element_colors()
    for ix in range(boundary_xmax):
        for iy in range(boundary_ymax):
            shift = ix * a_vec + iy * b_vec
            for site in structure:
                coord = np.array(site.coords, dtype=float) + shift
                xy = np.array([np.dot(coord, e1), np.dot(coord, e2)], dtype=float)
                atom_points.append(xy)
                symbol = canonicalize_species_symbol(site.specie)
                atomic_number = atomic_numbers[symbol]
                color = vesta_colors.get(symbol, tuple(float(value) for value in jmol_colors[atomic_number]))
                radius = 0.35 * covalent_radii[atomic_number]
                legend_entries.setdefault(symbol, color)
                ax.add_patch(
                    Circle(
                        (float(xy[0]), float(xy[1])),
                        radius=radius,
                        facecolor=color,
                        edgecolor="black",
                        linewidth=0.5,
                    ),
                )

    cell_corners = [
        np.zeros(3, dtype=float),
        boundary_xmax * a_vec,
        boundary_xmax * a_vec + boundary_ymax * b_vec,
        boundary_ymax * b_vec,
        np.zeros(3, dtype=float),
    ]
    outline = np.array([[np.dot(corner, e1), np.dot(corner, e2)] for corner in cell_corners], dtype=float)

    points = np.array(atom_points + [outline_point for outline_point in outline], dtype=float)
    min_xy = points.min(axis=0)
    max_xy = points.max(axis=0)
    span = np.maximum(max_xy - min_xy, 1.0)
    margin = 0.08 * max(span)

    coeff_points = np.linalg.solve(primitive_basis, points.T).T
    u_min = math.floor(float(coeff_points[:, 0].min())) - 1
    u_max = math.ceil(float(coeff_points[:, 0].max())) + 1
    v_min = math.floor(float(coeff_points[:, 1].min())) - 1
    v_max = math.ceil(float(coeff_points[:, 1].max())) + 1

    for u_index in range(u_min, u_max + 1):
        start = u_index * prim_a_xy + v_min * prim_b_xy
        end = u_index * prim_a_xy + v_max * prim_b_xy
        ax.plot([start[0], end[0]], [start[1], end[1]], color="0.78", linewidth=0.6, zorder=0)
    for v_index in range(v_min, v_max + 1):
        start = u_min * prim_a_xy + v_index * prim_b_xy
        end = u_max * prim_a_xy + v_index * prim_b_xy
        ax.plot([start[0], end[0]], [start[1], end[1]], color="0.78", linewidth=0.6, zorder=0)

    ax.plot(outline[:, 0], outline[:, 1], color="black", linewidth=1.2, zorder=1)

    parent_outline = np.array(
        [
            np.zeros(2, dtype=float),
            np.array([np.dot(a_vec, e1), np.dot(a_vec, e2)], dtype=float),
            np.array([np.dot(a_vec + b_vec, e1), np.dot(a_vec + b_vec, e2)], dtype=float),
            np.array([np.dot(b_vec, e1), np.dot(b_vec, e2)], dtype=float),
            np.zeros(2, dtype=float),
        ],
        dtype=float,
    )
    ax.plot(
        parent_outline[:, 0],
        parent_outline[:, 1],
        color="0.25",
        linewidth=1.0,
        linestyle=(0, (2.0, 2.0)),
        zorder=2,
    )

    ax.set_xlim(float(min_xy[0] - margin), float(max_xy[0] + margin))
    ax.set_ylim(float(min_xy[1] - margin), float(max_xy[1] + margin))
    ax.set_aspect("equal")
    ax.set_xlabel("Projected x (Angstrom)")
    ax.set_ylabel("Projected y (Angstrom)")
    ax.tick_params(axis="both", which="major", labelsize=8)
    ax.grid(False)

    atom_coeffs = np.linalg.solve(primitive_basis, np.array(atom_points, dtype=float).T).T
    axis_u, axis_v = select_annotation_primitive_cell(
        prim_a_xy=prim_a_xy,
        prim_b_xy=prim_b_xy,
        min_xy=min_xy,
        max_xy=max_xy,
        u_min=u_min,
        u_max=u_max,
        v_min=v_min,
        v_max=v_max,
    )

    primitive_outline = np.array(
        [
            axis_u * prim_a_xy + axis_v * prim_b_xy,
            (axis_u + 1) * prim_a_xy + axis_v * prim_b_xy,
            (axis_u + 1) * prim_a_xy + (axis_v + 1) * prim_b_xy,
            axis_u * prim_a_xy + (axis_v + 1) * prim_b_xy,
            axis_u * prim_a_xy + axis_v * prim_b_xy,
        ],
        dtype=float,
    )
    ax.plot(primitive_outline[:, 0], primitive_outline[:, 1], color="black", linewidth=1.4, zorder=3)

    lattice_point = axis_u * prim_a_xy + axis_v * prim_b_xy

    rel_coeffs = atom_coeffs - np.array([axis_u, axis_v], dtype=float)
    motif_mask = (
        (rel_coeffs[:, 0] >= -1e-8)
        & (rel_coeffs[:, 0] < 1.0 - 1e-8)
        & (rel_coeffs[:, 1] >= -1e-8)
        & (rel_coeffs[:, 1] < 1.0 - 1e-8)
    )
    motif_points = np.array(atom_points, dtype=float)[motif_mask]
    if len(motif_points) > 0:
        ax.scatter(
            motif_points[:, 0],
            motif_points[:, 1],
            s=140,
            facecolors="none",
            edgecolors="#d62728",
            linewidths=1.2,
            zorder=4,
        )

    ax.scatter(
        [float(lattice_point[0])],
        [float(lattice_point[1])],
        s=46,
        facecolors="white",
        edgecolors="black",
        linewidths=1.1,
        zorder=5,
    )
    ax.scatter(
        [float(lattice_point[0])],
        [float(lattice_point[1])],
        s=58,
        marker="+",
        color="black",
        linewidths=1.2,
        zorder=6,
    )

    arrow_origin = lattice_point + 0.14 * prim_a_xy + 0.14 * prim_b_xy
    a_tip = arrow_origin + 0.42 * prim_a_xy
    b_tip = arrow_origin + 0.42 * prim_b_xy
    ax.annotate("", xy=a_tip, xytext=arrow_origin, arrowprops={"arrowstyle": "->", "lw": 1.2, "color": "black"}, zorder=7)
    ax.annotate("", xy=b_tip, xytext=arrow_origin, arrowprops={"arrowstyle": "->", "lw": 1.2, "color": "black"}, zorder=7)
    ax.text(float(a_tip[0]), float(a_tip[1]), "a", fontsize=9, ha="center", va="center", zorder=8)
    ax.text(float(b_tip[0]), float(b_tip[1]), "b", fontsize=9, ha="center", va="center", zorder=8)

    legend_handles = [
        Line2D(
            [0],
            [0],
            marker="o",
            color="none",
            label=symbol,
            markerfacecolor=color,
            markeredgecolor="black",
            markeredgewidth=0.5,
            markersize=8,
        )
        for symbol, color in sorted(legend_entries.items())
    ]
    legend = ax.legend(
        handles=legend_handles,
        title="Elements",
        loc="upper left",
        bbox_to_anchor=(1.02, 1.0),
        borderaxespad=0.0,
        frameon=True,
        fontsize=8,
        title_fontsize=9,
    )
    fig.canvas.draw()
    legend_bbox = legend.get_window_extent(renderer=fig.canvas.get_renderer()).transformed(fig.transFigure.inverted())
    lattice_text = "\n".join(
        [
            "Primitive cell",
            "Parent cell = dotted overlay",
            f"plane group = {plane_group if plane_group is not None else 'n/a'}",
            f"a = {primitive_lattice.a:.3f} A",
            f"b = {primitive_lattice.b:.3f} A",
            f"c = {primitive_lattice.c:.3f} A",
            f"alpha = {primitive_lattice.alpha:.2f} deg",
            f"beta = {primitive_lattice.beta:.2f} deg",
            f"gamma = {primitive_lattice.gamma:.2f} deg",
        ],
    )
    fig.text(
        legend_bbox.x0,
        legend_bbox.y0 - 0.015,
        lattice_text,
        ha="left",
        va="top",
        fontsize=8,
        bbox={"boxstyle": "round,pad=0.28", "facecolor": "white", "edgecolor": "0.75", "alpha": 0.95},
    )
    fig.savefig(png_path, bbox_inches="tight", pad_inches=0.05)
    plt.close(fig)


def motif_summary_from_plane(plane: PlaneRecord) -> dict[str, Any]:
    """Serialize the plane-level motif summary."""
    return {
        "formula": tuple_to_formula_string(plane.species_counts),
        "species_counts": dict(plane.species_counts),
        "motif_sites_frac_xy": [
            {"species": species, "x_frac": x, "y_frac": y}
            for species, x, y in plane.motif_sites
        ],
        "motif_signature": [
            {"species": species, "x_frac": x, "y_frac": y}
            for species, x, y in plane.motif_signature
        ],
    }


def infer_bulk_oxidation_state_map(structure: Structure) -> tuple[dict[str, float], str]:
    """Infer one oxidation-state map for the bulk structure."""
    explicit_map: dict[str, float] = {}
    explicit_available = False
    for site in structure:
        symbol = site.specie.symbol
        oxi_state = getattr(site.specie, "oxi_state", None)
        if oxi_state is None:
            continue
        explicit_available = True
        explicit_map[symbol] = float(oxi_state)

    structure_symbols = {site.specie.symbol for site in structure}
    if explicit_available and structure_symbols.issubset(explicit_map):
        return explicit_map, "explicit oxidation states from CIF/site species"

    guesses = structure.composition.oxi_state_guesses(max_sites=-1)
    if guesses:
        guessed = {str(symbol): float(value) for symbol, value in guesses[0].items()}
        return guessed, "bulk composition guess via pymatgen Composition.oxi_state_guesses"

    return {}, "unavailable"


def ionic_valence_sum_from_species_counts(
    species_counts: tuple[tuple[str, int], ...],
    oxidation_state_map: dict[str, float],
) -> float | None:
    """Return the total ionic valence sum for one species-count tuple."""
    if not oxidation_state_map:
        return None
    total = 0.0
    for symbol, count in species_counts:
        if symbol not in oxidation_state_map:
            return None
        total += oxidation_state_map[symbol] * count
    return total


def format_signed_charge(value: float | None) -> str:
    """Format one signed charge-like value for tables and markdown."""
    if value is None:
        return "n/a"
    rounded = round(value)
    if math.isclose(value, rounded, abs_tol=1e-8):
        return f"{int(rounded):+d}"
    return f"{value:+.3f}"


def format_oxidation_state_map(oxidation_state_map: dict[str, float]) -> str:
    """Format the oxidation-state map for human-readable report text."""
    if not oxidation_state_map:
        return "n/a"
    items = sorted(oxidation_state_map.items())
    return ", ".join(f"{symbol}:{format_signed_charge(value)}" for symbol, value in items)


def build_termination_records(
    plane_records: list[PlaneRecord],
    slab_generator: SlabGenerator,
    symprec: float,
    oxidation_state_map: dict[str, float],
) -> list[dict[str, Any]]:
    """Enumerate unique top terminations by cyclic plane-stack signatures."""
    if not plane_records:
        return []

    unique_records: list[dict[str, Any]] = []
    seen_top_sequences: dict[tuple[Any, ...], int] = {}
    n_planes = len(plane_records)

    for top_index in range(n_planes):
        top_sequence = cyclic_rotated_signatures(plane_records, top_index)
        if top_sequence in seen_top_sequences:
            continue

        prev_index = (top_index - 1) % n_planes
        cut_shift_frac = midpoint_on_unit_interval(plane_records[prev_index].z_frac, plane_records[top_index].z_frac)
        opposite_sequence = reverse_cyclic_rotated_signatures(plane_records, prev_index)

        # Reuse the already-oriented generator so that cut shifts stay on the same convention.
        representative_slab = slab_generator.get_slab(shift=cut_shift_frac)
        dataset = spglib.get_layergroup(
            lattice_rows_to_spglib_cell(representative_slab),
            aperiodic_dir=2,
            symprec=symprec,
        )
        plane_group = classify_plane_group(dataset, np.array(representative_slab.lattice.matrix[:2], dtype=float))

        unique_id = len(unique_records)
        seen_top_sequences[top_sequence] = unique_id
        unique_records.append(
            {
                "termination_id": unique_id,
                "top_plane_id": plane_records[top_index].plane_id,
                "bottom_plane_id": plane_records[prev_index].plane_id,
                "cut_shift_frac": cut_shift_frac,
                "top_formula": tuple_to_formula_string(plane_records[top_index].species_counts),
                "bottom_formula": tuple_to_formula_string(plane_records[prev_index].species_counts),
                "top_ionic_valence_sum": ionic_valence_sum_from_species_counts(
                    plane_records[top_index].species_counts,
                    oxidation_state_map,
                ),
                "bottom_ionic_valence_sum": ionic_valence_sum_from_species_counts(
                    plane_records[prev_index].species_counts,
                    oxidation_state_map,
                ),
                "top_termination_sequence": [plane.plane_id for plane in plane_records[top_index:] + plane_records[:top_index]],
                "bottom_opposite_normal_sequence": [
                    plane_records[(prev_index - offset) % n_planes].plane_id for offset in range(n_planes)
                ],
                "top_motif": motif_summary_from_plane(plane_records[top_index]),
                "bottom_motif": motif_summary_from_plane(plane_records[prev_index]),
                "surface_layer_group_symbol": None if dataset is None else dataset.international,
                "surface_layer_group_number": None if dataset is None else int(dataset.number),
                "surface_point_group": None if dataset is None else dataset.pointgroup,
                "plane_group": plane_group,
                "paired_opposite_normal_signature_is_identical": top_sequence == opposite_sequence,
                "representative_slab": representative_slab,
            },
        )

    return unique_records


def annotate_termination_intervals(
    termination_records: list[dict[str, Any]],
    c_repeat: float,
) -> None:
    """Annotate cyclic spacing from the previous listed unique termination."""
    if not termination_records:
        return
    for index, record in enumerate(termination_records):
        previous_cut = termination_records[index - 1]["cut_shift_frac"]
        current_cut = record["cut_shift_frac"]
        delta_frac = (current_cut - previous_cut) % 1.0
        record["spacing_from_prev_frac"] = delta_frac
        record["spacing_from_prev_ang"] = delta_frac * c_repeat


def write_plane_stack_csv(
    path: Path,
    plane_records: list[PlaneRecord],
    oxidation_state_map: dict[str, float],
) -> None:
    """Write the plane stack table."""
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "plane_id",
                "z_frac",
                "z_ang",
                "n_sites",
                "species_formula",
                "ionic_valence_sum",
                "site_indices",
                "motif_signature",
            ],
        )
        writer.writeheader()
        for plane in plane_records:
            writer.writerow(
                {
                    "plane_id": plane.plane_id,
                    "z_frac": f"{plane.z_frac:.10f}",
                    "z_ang": f"{plane.z_ang:.10f}",
                    "n_sites": len(plane.site_indices),
                    "species_formula": tuple_to_formula_string(plane.species_counts),
                    "ionic_valence_sum": format_signed_charge(
                        ionic_valence_sum_from_species_counts(plane.species_counts, oxidation_state_map),
                    ),
                    "site_indices": ",".join(map(str, plane.site_indices)),
                    "motif_signature": json.dumps(
                        [
                            {"species": species, "x_frac": x, "y_frac": y}
                            for species, x, y in plane.motif_signature
                        ],
                        ensure_ascii=True,
                    ),
                },
            )


def write_terminations_json(path: Path, termination_records: list[dict[str, Any]]) -> None:
    """Write the termination summary without embedding full slab objects."""
    payload = []
    for record in termination_records:
        serializable = {key: value for key, value in record.items() if key != "representative_slab"}
        payload.append(serializable)
    path.write_text(json.dumps(payload, indent=2, ensure_ascii=True), encoding="utf-8")


def write_report_md(
    path: Path,
    cif_path: Path,
    hkl: tuple[int, int, int],
    termination_only: bool,
    space_group_symbol: str,
    space_group_number: int,
    plane_records: list[PlaneRecord],
    termination_records: list[dict[str, Any]],
    oxidation_state_map: dict[str, float],
    oxidation_state_source: str,
) -> None:
    """Write a compact markdown report for quick inspection."""
    lines = [
        f"# Plane Analysis for `{cif_path.name}`",
        "",
        f"- hkl: `{hkl}`",
        f"- mode: `{'termination-only' if termination_only else 'full'}`",
        f"- recovered bulk space group: `{space_group_symbol}` ({space_group_number})",
        f"- unique atomic planes in one oriented repeat: `{len(plane_records)}`",
        f"- unique top terminations along +normal: `{len(termination_records)}`",
        f"- oxidation-state source for ionic valence sums: `{oxidation_state_source}`",
        f"- oxidation-state map used: `{format_oxidation_state_map(oxidation_state_map)}`",
        "",
        "## Definitions",
        "",
        "- `cut_shift_frac` is the fractional cut position passed to `SlabGenerator.get_slab(shift=...)` along the reoriented slab normal.",
        "- It is chosen as the midpoint between the previous plane and the top plane in wrapped fractional `z`, so the cleavage plane sits between those two atomic planes.",
        "- `spacing_from_prev_frac` and `spacing_from_prev_ang` are the cyclic distances from the previous listed unique termination to the current one along the same slab-normal repeat.",
        "- `top_ionic_valence_sum` is the signed sum of `count x oxidation_state` over atoms in the top termination plane.",
        "",
        "## Libraries Used",
        "",
        "| library | used feature in this workflow |",
        "| --- | --- |",
        "| `pymatgen.core.Structure` | CIF parsing, structure I/O, primitive-cell reduction |",
        "| `pymatgen.symmetry.analyzer.SpacegroupAnalyzer` | conventional-standard bulk cell recovery |",
        "| `pymatgen.core.surface.SlabGenerator` | oriented unit cell construction and slab cutting with `shift=cut_shift_frac` |",
        "| `pymatgen.core.Composition` | oxidation-state guessing for ionic valence sums |",
        "| `spglib` | layer-group symmetry detection for surface/layer symmetry |",
        "| `numpy` | vector projection, in-plane basis transforms, clustering support |",
        "| `matplotlib` | top-view `top.png` rendering and annotation overlays |",
        "| `ase.data` | atomic radii and default element colors for 2D atom drawing |",
        "",
        "## Terminations",
        "",
        "| termination_id | cut_shift_frac | spacing_from_prev_frac | spacing_from_prev_ang | top_formula | top_ionic_valence_sum | plane_group | layer_group | point_group |",
        "| --- | ---: | ---: | ---: | --- | ---: | --- | --- | --- |",
    ]
    for record in termination_records:
        lines.append(
            "| "
            f"{record['termination_id']} | "
            f"{record['cut_shift_frac']:.8f} | "
            f"{record['spacing_from_prev_frac']:.8f} | "
            f"{record['spacing_from_prev_ang']:.6f} | "
            f"{record['top_formula']} | "
            f"{format_signed_charge(record['top_ionic_valence_sum'])} | "
            f"{record['plane_group']} | "
            f"{record['surface_layer_group_symbol']} | "
            f"{record['surface_point_group']} |"
        )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_open_all_in_vesta_script(surface_only_dir: Path, vesta_view_paths: list[Path]) -> Path:
    """Write a launcher that opens all generated VESTA-view CIF files."""
    script_path = surface_only_dir / "open_all_in_vesta.command"
    lines = [
        "#!/bin/zsh",
        "set -euo pipefail",
        "",
    ]
    for cif_path in vesta_view_paths:
        lines.append(f"open -na /Applications/VESTA.app --args -open '{cif_path}'")
    lines.append("")
    script_path.write_text("\n".join(lines), encoding="utf-8")
    script_path.chmod(0o755)
    return script_path


def analyze_direction(
    cif_path: Path,
    hkl: tuple[int, int, int],
    output_base_dir: Path,
    termination_only: bool = False,
    surface_only_export: bool = False,
    symprec: float = 1e-3,
    plane_tol_ang: float = 0.05,
    motif_round_digits: int = 8,
    min_slab_size: float = 12.0,
    min_vacuum_size: float = 15.0,
    surface_only_vacuum_size: float = 20.0,
    png_boundary_xmax: int = 3,
    png_boundary_ymax: int = 3,
) -> dict[str, Any]:
    """Run the full plane-stack and termination analysis for one CIF and one hkl."""
    if tuple(hkl) == (0, 0, 0):
        raise ValueError("hkl=(0, 0, 0) is not a valid plane normal.")

    structure, space_group_symbol, space_group_number = get_standard_conventional_structure(cif_path, symprec=symprec)
    oxidation_state_map, oxidation_state_source = infer_bulk_oxidation_state_map(structure)
    slab_generator = SlabGenerator(
        initial_structure=structure,
        miller_index=hkl,
        min_slab_size=min_slab_size,
        min_vacuum_size=min_vacuum_size,
        lll_reduce=False,
        center_slab=False,
        in_unit_planes=False,
        primitive=False,
        max_normal_search=None,
        reorient_lattice=True,
    )
    oriented_unit_cell = slab_generator.oriented_unit_cell
    plane_records = build_plane_records(
        oriented_unit_cell=oriented_unit_cell,
        plane_tol_ang=plane_tol_ang,
        motif_round_digits=motif_round_digits,
    )
    termination_records = build_termination_records(
        plane_records=plane_records,
        slab_generator=slab_generator,
        symprec=symprec,
        oxidation_state_map=oxidation_state_map,
    )
    annotate_termination_intervals(termination_records, oriented_unit_cell.lattice.c)

    cif_stem = cif_path.stem
    hkl_label = f"hkl_{hkl[0]}_{hkl[1]}_{hkl[2]}"
    result_dir = output_base_dir / cif_stem / hkl_label
    result_dir.mkdir(parents=True, exist_ok=True)
    slab_dir = None if termination_only else result_dir / "slabs"
    surface_only_dir = result_dir / "surface_only" if surface_only_export else None
    if slab_dir is not None:
        slab_dir.mkdir(parents=True, exist_ok=True)
    if surface_only_dir is not None:
        surface_only_dir.mkdir(parents=True, exist_ok=True)

    write_plane_stack_csv(result_dir / "plane_stack.csv", plane_records, oxidation_state_map)

    if slab_dir is not None:
        for record in termination_records:
            slab_path = slab_dir / f"termination_{record['termination_id']:03d}.cif"
            record["representative_slab"].to(filename=str(slab_path))
            record["representative_slab_cif"] = str(slab_path)

    launcher_path = None
    if surface_only_dir is not None:
        vesta_view_paths: list[Path] = []
        for record in termination_records:
            top_plane = plane_records[record["top_plane_id"]]
            top_only_structure = build_top_plane_only_structure(
                oriented_unit_cell=oriented_unit_cell,
                top_plane=top_plane,
                vacuum_size=surface_only_vacuum_size,
            )
            surface_only_path = surface_only_dir / f"termination_{record['termination_id']:03d}_top.cif"
            top_only_structure.to(filename=str(surface_only_path))
            record["top_plane_only_cif"] = str(surface_only_path)
            record["top_plane_only_site_count"] = len(top_plane.motif_sites)

            vesta_view_structure = build_vesta_view_structure(
                structure=top_only_structure,
                repeat_x=png_boundary_xmax,
                repeat_y=png_boundary_ymax,
            )
            vesta_view_path = surface_only_path.with_name(surface_only_path.stem + "_vesta_view.cif")
            vesta_view_structure.to(filename=str(vesta_view_path))
            record["top_plane_only_vesta_view_cif"] = str(vesta_view_path)
            record["top_plane_only_vesta_view_site_count"] = len(vesta_view_structure)
            vesta_view_paths.append(vesta_view_path)

            png_path = surface_only_path.with_suffix(".png")
            render_structure_top_view_png(
                structure=top_only_structure,
                png_path=png_path,
                boundary_xmax=png_boundary_xmax,
                boundary_ymax=png_boundary_ymax,
                plane_group=record["plane_group"],
            )
            record["top_plane_only_png"] = str(png_path)
        launcher_path = write_open_all_in_vesta_script(surface_only_dir=surface_only_dir, vesta_view_paths=vesta_view_paths)

    write_terminations_json(result_dir / "terminations.json", termination_records)
    write_report_md(
        result_dir / "report.md",
        cif_path=cif_path,
        hkl=hkl,
        termination_only=termination_only,
        space_group_symbol=space_group_symbol,
        space_group_number=space_group_number,
        plane_records=plane_records,
        termination_records=termination_records,
        oxidation_state_map=oxidation_state_map,
        oxidation_state_source=oxidation_state_source,
    )

    return {
        "cif_path": str(cif_path),
        "hkl": list(hkl),
        "mode": "termination-only" if termination_only else "full",
        "space_group_symbol": space_group_symbol,
        "space_group_number": space_group_number,
        "n_planes": len(plane_records),
        "n_unique_terminations": len(termination_records),
        "result_dir": str(result_dir),
        "plane_stack_csv": str(result_dir / "plane_stack.csv"),
        "terminations_json": str(result_dir / "terminations.json"),
        "report_md": str(result_dir / "report.md"),
        "slab_dir": None if slab_dir is None else str(slab_dir),
        "surface_only_dir": None if surface_only_dir is None else str(surface_only_dir),
        "surface_only_vesta_launcher": None if launcher_path is None else str(launcher_path),
    }


def main() -> None:
    """CLI entry point."""
    args = parse_args()
    summary = analyze_direction(
        cif_path=Path(args.cif).resolve(),
        hkl=tuple(args.hkl),
        output_base_dir=Path(args.output_dir).resolve(),
        termination_only=args.termination_only,
        surface_only_export=args.surface_only_export,
        symprec=args.symprec,
        plane_tol_ang=args.plane_tol_ang,
        motif_round_digits=args.motif_round_digits,
        min_slab_size=args.min_slab_size,
        min_vacuum_size=args.min_vacuum_size,
        surface_only_vacuum_size=args.surface_only_vacuum_size,
        png_boundary_xmax=args.png_boundary_xmax,
        png_boundary_ymax=args.png_boundary_ymax,
    )
    print(json.dumps(summary, indent=2, ensure_ascii=True))


if __name__ == "__main__":
    main()
