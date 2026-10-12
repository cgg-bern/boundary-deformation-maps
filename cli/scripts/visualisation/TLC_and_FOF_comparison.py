#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Analyze mesh processing results across subdirectories and generate scatter plots.

--- Default mode (star analysis) ---

Directory structure expected:
    <input_directory>/
    +-- remeshing-maps-data/
    |   +-- block_A_star_<idx>_boundary_deformation_map_output.json   (DAF/VE)
    +-- remeshing-maps-data-preremeshing-only/
    |   +-- block_A_star_<idx>_boundary_deformation_map_output.json   (DAF pre-remeshing only)
    +-- remeshing-maps-data-no-remeshing/
    |   +-- block_A_star_<idx>_boundary_deformation_map_output.json   (DAF no remeshing)
    +-- block_A/
    |   +-- star_<idx>_FOF.time
    |   +-- star_<idx>_TLC.time
    |   +-- star_<idx>_FOF_stats.txt
    |   +-- star_<idx>_TLC_stats.txt
    |   +-- expansion_data/
    |   |   +-- star_codomain_<idx>.json                   (SAE)
    |   +-- expansion_data_foliation/
    |       +-- star_domain_<idx>_foliation.json           (FOL)
    +-- block_B/
        +-- ...

--- Spine mode (--analyse-spine) ---

Directory structure expected:
    <input_directory>/
    +-- AIS-1/
    |   +-- FSU-XX_FOF.time
    |   +-- FSU-XX_FOF_stats.txt
    |   +-- FSU-XX_target_vol_TLC.time
    |   +-- FSU-XX_target_vol_TLC_stats.txt
    |   +-- FSU-XX_output.json                             (DAF/VE)
    +-- AIS-2/
        +-- ...

Cases are identified as AIS-N_FSU_XX. Only DAF, FOF and TLC are available;
SAE and FOL are not used.

--- Usage ---

    python analyze_results.py --input-directory /path/to/dir --plot-output-path plot.png
    python analyze_results.py --input-directory /path/to/dir --plot-output-path plot.png --show-plots
    python analyze_results.py --input-directory /path/to/dir --plot-output-path plot.png --analyse-spine
"""


import argparse
import json
import re
import sys
from pathlib import Path

import matplotlib.pyplot as plt


TOTAL_CASES_SPINE = 117
TOTAL_CASES_VOLMAP = 1949
TOTAL_CASES_STARS = 1784

TOTAL_CASES = TOTAL_CASES_SPINE


# ---------------------------------------------------------------------------
# File parsers
# ---------------------------------------------------------------------------

def parse_stats_txt(filepath):
    """Parse a _stats.txt file, returning a dict of the 8 fields."""
    fields = [
        "mesh_name",
        "n_init_vertices",
        "n_init_edges",
        "n_init_faces",
        "n_init_cells",
        "n_degenerate_cells",
        "n_flipped_cells",
        "idx",
    ]
    values = {}
    try:
        with open(filepath) as f:
            lines = [l.strip() for l in f.readlines() if l.strip()]
        for i, field in enumerate(fields):
            if i < len(lines):
                raw = lines[i]
                if field != "mesh_name":
                    try:
                        raw = int(raw)
                    except ValueError:
                        try:
                            raw = float(raw)
                        except ValueError:
                            pass
                values[field] = raw
    except Exception as e:
        print(f"  Warning: could not parse {filepath}: {e}")
    return values


def parse_time_file(filepath):
    """Parse a .time file, returning runtime in seconds as float.

    Expects a value in nanoseconds with a trailing 'ns', e.g. '123412341234ns'.
    """
    try:
        with open(filepath) as f:
            content = f.read().strip()

        if not content:
            print(f"  Warning: empty time file {filepath}")
            return None

        if not content.endswith("ns"):
            print(f"  Warning: unexpected format in {filepath.name}: {content!r}")
            return None

        return int(content[:-2]) / 1e9

    except Exception as e:
        print(f"  Warning: could not parse {filepath}: {e}")
        return None


def parse_json_file(filepath):
    """Parse a JSON file, returning its contents as a dict.

    Tolerates trailing commas before closing braces/brackets, which is
    technically invalid JSON but appears in some generated files.
    """
    try:
        with open(filepath) as f:
            content = f.read()
        # Strip trailing commas before } or ] (e.g. from code generators)
        content = re.sub(r",(\s*[}\]])", r"\1", content)
        return json.loads(content)
    except Exception as e:
        print(f"  Warning: could not parse {filepath}: {e}")
        return None


def extract_idx(filename, pattern):
    """Extract the integer index from a filename using a regex pattern.
    The pattern must contain one capturing group for the index."""
    m = re.search(pattern, filename)
    if m:
        return int(m.group(1))
    return None


# ---------------------------------------------------------------------------
# Data collection
# ---------------------------------------------------------------------------

def collect_ve_from_dir(ve_dir, subdirs):
    """Load VE/DAF data from a flat directory of JSON files.
    Returns a dict keyed by (subdir_name, idx) -> parsed JSON dict."""
    result = {}
    if not ve_dir.is_dir():
        print(f"  Warning: {ve_dir} not found, DAF variant data will be missing.")
        return result
    for subdir in subdirs:
        for f in ve_dir.iterdir():
            if not f.is_file():
                continue
            pattern = rf"^{re.escape(subdir.name)}_star_(\d+)_boundary_deformation_map_output\.json$"
            m = re.match(pattern, f.name)
            if m:
                idx = int(m.group(1))
                raw = parse_json_file(f)
                if raw is not None:
                    result[(subdir.name, idx)] = raw
    return result


def collect_data(input_dir):
    """
    Walk input_dir -> block subdirs -> files.

    Cases are uniquely identified by (subdir, idx).

    Returns five lists of dicts, one per method: FOF, TLC, VE, SAE, FOL.
    Each dict contains at minimum:
        subdir, idx, n_init_vertices, n_init_cells, runtime_s, success
    FOF and TLC additionally contain 'stats'.
    VE, SAE, FOL additionally contain 'raw' (the full JSON).
    """
    fof_data = []
    tlc_data = []
    ve_data  = []
    sae_data = []
    fol_data = []

    input_path = Path(input_dir)
    if not input_path.is_dir():
        print(f"Error: {input_dir} is not a valid directory.")
        sys.exit(1)

    all_subdirs = [p for p in input_path.iterdir() if p.is_dir()]
    subdirs = [p for p in all_subdirs if "block" in p.name]
    if not subdirs:
        print(f"No subdirectories containing 'block' in their name found in {input_dir}.")
        print(f"  (Found {len(all_subdirs)} subdirectory(ies) total, none matched.)")
        sys.exit(1)

    print(f"Found {len(subdirs)} subdirectory(ies) with 'block' in name to process "
          f"(skipped {len(all_subdirs) - len(subdirs)} other(s)).")

    # VE files now live in a flat remeshing-maps-data/ directory at the input root,
    # named <subdir>_star_<idx>_boundary_deformation_map_output.json.
    ve_dir = input_path / "remeshing-maps-data"
    if not ve_dir.is_dir():
        print(f"  Warning: remeshing-maps-data/ directory not found in {input_dir}. "
              f"VE data will be missing.")

    for subdir in sorted(subdirs):
        print(f"\nProcessing subdirectory: {subdir.name}")

        # --- FOF / TLC: files directly in subdir ---
        # Group by idx: star_<idx>_FOF.time etc.
        cases = {}  # idx -> dict of file paths

        for f in subdir.iterdir():
            if not f.is_file():
                continue
            name = f.name
            if name.endswith("_FOF.time"):
                idx = extract_idx(name, r"star_(\d+)_FOF\.time")
                if idx is not None:
                    cases.setdefault(idx, {})["FOF_time"] = f
            elif name.endswith("_TLC.time"):
                idx = extract_idx(name, r"star_(\d+)_TLC\.time")
                if idx is not None:
                    cases.setdefault(idx, {})["TLC_time"] = f
            elif name.endswith("_FOF_stats.txt"):
                idx = extract_idx(name, r"star_(\d+)_FOF_stats\.txt")
                if idx is not None:
                    cases.setdefault(idx, {})["FOF_stats"] = f
            elif name.endswith("_TLC_stats.txt"):
                idx = extract_idx(name, r"star_(\d+)_TLC_stats\.txt")
                if idx is not None:
                    cases.setdefault(idx, {})["TLC_stats"] = f

        # --- SAE: expansion_data/star_codomain_<idx>.json ---
        sae_dir = subdir / "expansion_data"
        if sae_dir.is_dir():
            for f in sae_dir.iterdir():
                if not f.is_file():
                    continue
                idx = extract_idx(f.name, r"star_codomain_(\d+)\.json")
                if idx is not None:
                    cases.setdefault(idx, {})["SAE_json"] = f
        else:
            print(f"  Warning: no expansion_data/ subdirectory found.")

        # --- FOL: expansion_data_foliation/star_domain_<idx>_foliation.json ---
        fol_dir = subdir / "expansion_data_foliation"
        if fol_dir.is_dir():
            for f in fol_dir.iterdir():
                if not f.is_file():
                    continue
                idx = extract_idx(f.name, r"star_domain_(\d+)_foliation\.json")
                if idx is not None:
                    cases.setdefault(idx, {})["FOL_json"] = f
        else:
            print(f"  Warning: no expansion_data_foliation/ subdirectory found.")

        # --- VE: remeshing-maps-data/<subdir>_star_<idx>_boundary_deformation_map_output.json ---
        if ve_dir.is_dir():
            for idx in list(cases.keys()):
                ve_file = ve_dir / f"{subdir.name}_star_{idx}_boundary_deformation_map_output.json"
                if ve_file.is_file():
                    cases[idx]["VE_json"] = ve_file

        print(f"  Found {len(cases)} case(s).")

        for idx, file_map in cases.items():
            # Parse stats for FOF and TLC eagerly; used to extract n_init_vertices
            # which is shared across all five methods and sanity-checked for consistency.
            fof_stats = None
            tlc_stats = None
            ve_raw    = None
            sae_raw   = None
            fol_raw   = None

            if "FOF_stats" in file_map:
                fof_stats = parse_stats_txt(file_map["FOF_stats"])
            if "TLC_stats" in file_map:
                tlc_stats = parse_stats_txt(file_map["TLC_stats"])
            if "VE_json" in file_map:
                ve_raw = parse_json_file(file_map["VE_json"])
            if "SAE_json" in file_map:
                sae_raw = parse_json_file(file_map["SAE_json"])
            if "FOL_json" in file_map:
                fol_raw = parse_json_file(file_map["FOL_json"])

            # Collect n_init_vertices (for storage); not sanity-checked since VE can differ.
            # FOF/TLC/VE use "n_init_vertices"; SAE/FOL use "n_vertices".
            n_init_vertices = None
            if fof_stats and "n_init_vertices" in fof_stats:
                n_init_vertices = fof_stats["n_init_vertices"]
            elif tlc_stats and "n_init_vertices" in tlc_stats:
                n_init_vertices = tlc_stats["n_init_vertices"]
            elif sae_raw and "n_vertices" in sae_raw:
                n_init_vertices = sae_raw["n_vertices"]
            elif fol_raw and "n_vertices" in fol_raw:
                n_init_vertices = fol_raw["n_vertices"]

            # Sanity check n_init_cells across FOF, TLC, SAE, FOL (VE excluded).
            # This shared value is used as the x-axis for all methods including VE.
            cell_sources = {}
            if fof_stats and "n_init_cells" in fof_stats:
                cell_sources["FOF"] = fof_stats["n_init_cells"]
            if tlc_stats and "n_init_cells" in tlc_stats:
                cell_sources["TLC"] = tlc_stats["n_init_cells"]
            if sae_raw and "n_init_cells" in sae_raw:
                cell_sources["SAE"] = sae_raw["n_init_cells"]
            if fol_raw and "n_init_cells" in fol_raw:
                cell_sources["FOL"] = fol_raw["n_init_cells"]

            unique_cell_values = set(cell_sources.values())
            if len(unique_cell_values) > 1:
                print(f"  Warning [idx={idx}] n_init_cells mismatch across methods: "
                      f"{cell_sources}")
            n_init_cells = next(iter(unique_cell_values)) if unique_cell_values else None
            if n_init_cells is None:
                print(f"  Warning [idx={idx}] could not determine n_init_cells, skipping.")
                continue

            # --- FOF ---
            if "FOF_time" in file_map and fof_stats is not None:
                runtime = parse_time_file(file_map["FOF_time"])
                if runtime is not None:
                    fof_data.append({
                        "subdir": subdir.name,
                        "idx": idx,
                        "n_init_vertices": n_init_vertices,
                        "n_init_cells": n_init_cells,
                        "runtime_s": runtime,
                        "success": is_success_stats(fof_stats, n_init_cells),
                        "stats": fof_stats,
                    })
            else:
                missing = [k for k in ("FOF_time", "FOF_stats") if k not in file_map]
                if missing:
                    print(f"  [idx={idx}] Missing FOF files: {missing}")

            # --- TLC ---
            if "TLC_time" in file_map and tlc_stats is not None:
                runtime = parse_time_file(file_map["TLC_time"])
                if runtime is not None:
                    tlc_data.append({
                        "subdir": subdir.name,
                        "idx": idx,
                        "n_init_vertices": n_init_vertices,
                        "n_init_cells": n_init_cells,
                        "runtime_s": runtime,
                        "success": is_success_stats(tlc_stats, n_init_cells),
                        "stats": tlc_stats,
                    })
            else:
                missing = [k for k in ("TLC_time", "TLC_stats") if k not in file_map]
                if missing:
                    print(f"  [idx={idx}] Missing TLC files: {missing}")

            # --- VE ---
            if ve_raw is not None:
                runtime = ve_raw.get("runtime_s")
                if runtime is not None:
                    ve_data.append({
                        "subdir": subdir.name,
                        "idx": idx,
                        "n_init_vertices": n_init_vertices,
                        "n_init_cells": n_init_cells,
                        "n_final_vertices": ve_raw.get("n_final_vertices"),
                        "runtime_s": runtime,
                        "success": ve_raw.get("opti_result") == 0,
                        "raw": ve_raw,
                    })
                else:
                    print(f"  [idx={idx}] VE JSON missing key: runtime_s")
            else:
                print(f"  [idx={idx}] Missing VE JSON file.")

            # --- SAE ---
            if sae_raw is not None:
                runtime = sae_raw.get("expansion_time_s")
                if runtime is not None:
                    sae_data.append({
                        "subdir": subdir.name,
                        "idx": idx,
                        "n_init_vertices": n_init_vertices,
                        "n_init_cells": n_init_cells,
                        "n_final_vertices": sae_raw.get("n_vertices_after_expansion"),
                        "runtime_s": runtime,
                        "success": sae_raw.get("expansion_result") == 0,
                        "raw": sae_raw,
                    })
                else:
                    print(f"  [idx={idx}] SAE JSON missing key: expansion_time_s")
            else:
                print(f"  [idx={idx}] Missing SAE JSON file.")

            # --- FOL ---
            if fol_raw is not None:
                runtime = fol_raw.get("foliation_time_s")
                if runtime is not None:
                    fol_data.append({
                        "subdir": subdir.name,
                        "idx": idx,
                        "n_init_vertices": n_init_vertices,
                        "n_init_cells": n_init_cells,
                        "n_final_vertices": fol_raw.get("n_final_vertices"),
                        "runtime_s": runtime,
                        "success": fol_raw.get("valid_in_double") == 1,
                        "raw": fol_raw,
                    })
                else:
                    print(f"  [idx={idx}] FOL JSON missing key: foliation_time_s")
            else:
                print(f"  [idx={idx}] Missing FOL JSON file.")

    return fof_data, tlc_data, ve_data, sae_data, fol_data


def collect_daf_variant(input_dir, subdir_name):
    """Load one DAF variant from <input_dir>/<subdir_name>/.
    Returns a dict keyed by (subdir, idx) -> entry dict (same shape as ve_data entries),
    where subdir matches the block subdirectory name."""
    result = {}
    input_path = Path(input_dir)
    variant_dir = input_path / subdir_name

    if not variant_dir.is_dir():
        print(f"  Warning: {subdir_name}/ not found, variant data will be missing.")
        return result

    for f in variant_dir.iterdir():
        if not f.is_file():
            continue
        # Pattern: <block_name>_star_<idx>_boundary_deformation_map_output.json
        m = re.match(r"^(.+)_star_(\d+)_boundary_deformation_map_output.*\.json$", f.name)
        if not m:
            continue
        block_name = m.group(1)
        idx = int(m.group(2))
        raw = parse_json_file(f)
        if raw is None:
            continue
        runtime = raw.get("runtime_s")
        if runtime is None:
            print(f"  Warning [{subdir_name}] {f.name} missing runtime_s")
            continue
        result[(block_name, idx)] = {
            "subdir": block_name,
            "idx": idx,
            "n_init_vertices": raw.get("n_init_vertices"),
            "n_init_cells": raw.get("n_init_cells"),
            "n_final_vertices": raw.get("n_final_vertices"),
            "runtime_s": runtime,
            "success": raw.get("opti_result") == 0,
            "raw": raw,
        }
    n_success = sum(1 for e in result.values() if e["success"])
    print(f"  {subdir_name}: {len(result)} cases found, {n_success} successful")
    if result:
        sample_key = next(iter(result))
        print(f"    Sample key: {sample_key}")
    return result


# ---------------------------------------------------------------------------
# Analysis helpers
# ---------------------------------------------------------------------------

def is_success_stats(stats, n_init_cells):
    """A stats result is successful if there are no degenerate cells and
    flipped cells are either zero or equal to n_init_cells (fully flipped,
    which is just a convention difference)."""
    n_degen   = stats.get("n_degenerate_cells", -1)
    n_flipped = stats.get("n_flipped_cells", -1)
    return n_degen == 0 and (n_flipped == 0 or n_flipped == n_init_cells)


def is_success(entry):
    return entry.get("success", False)


def compute_common_successes(fof_data, tlc_data, ve_data, sae_data, fol_data,
                             ignore_tlc=False):
    """Return the subset of cases where all relevant methods succeeded,
    as five parallel filtered lists keyed by (subdir, idx).
    When ignore_tlc=True, TLC is excluded from the intersection."""
    def success_map(data):
        return {(d["subdir"], d["idx"]): d for d in data if is_success(d)}

    fof_s = success_map(fof_data)
    tlc_s = success_map(tlc_data)
    ve_s  = success_map(ve_data)
    sae_s = success_map(sae_data)
    fol_s = success_map(fol_data)

    common_keys = set(fof_s) & set(ve_s) & set(sae_s) & set(fol_s)
    if not ignore_tlc:
        common_keys &= set(tlc_s)

    return (
        [fof_s[k] for k in common_keys],
        [tlc_s.get(k) for k in common_keys],
        [ve_s[k]  for k in common_keys],
        [sae_s[k] for k in common_keys],
        [fol_s[k] for k in common_keys],
    )


def print_summary(fof_data, tlc_data, ve_data, sae_data, fol_data):
    all_methods = [("FOF", fof_data), ("TLC", tlc_data), ("DAF", ve_data),
                   ("SAE", sae_data), ("FOL", fol_data)]

    print("\n--- Summary ---")
    for label, data in all_methods:
        if data:
            runtimes = [d["runtime_s"] for d in data]
            cells = [d["n_init_cells"] for d in data]
            print(
                f"  {label}: {len(data)} cases | "
                f"runtime [{min(runtimes):.3f}s – {max(runtimes):.3f}s] | "
                f"n_init_cells [{min(cells)} – {max(cells)}]"
            )
        else:
            print(f"  {label}: 0 cases")

    # Per-method success rates
    print(f"\n--- Success rates (out of {TOTAL_CASES} total cases) ---")
    for label, data in all_methods:
        n_success = sum(1 for d in data if is_success(d))
        print(f"  {label}: {n_success} / {TOTAL_CASES} ({100 * n_success / TOTAL_CASES:.1f}%)")

    def print_runtime_per_cell(common, label):
        for method, data in zip(["FOF", "TLC", "DAF", "SAE", "FOL"], common):
            if not data or any(d is None for d in data):
                continue
            rates = [d["runtime_s"] / d["n_init_cells"] for d in data]
            print(f"  {method}: avg = {sum(rates)/len(rates):.6f} s/cell | "
                  f"max = {max(rates):.6f} s/cell || avg = {sum(d["runtime_s"])/(len(d["runtime_s"]*60)):.6f} min | "
                  f"max = {max(d["runtime_s"])/60:.6f} min" )

    # Common success rate and runtime/cell — all five methods
    common_5 = compute_common_successes(fof_data, tlc_data, ve_data, sae_data, fol_data,
                                        ignore_tlc=False)
    n_common_5 = len(common_5[0])
    print(f"\n--- Common success rate (all five methods) ---")
    print(f"  {n_common_5} / {TOTAL_CASES} ({100 * n_common_5 / TOTAL_CASES:.1f}%)")
    print(f"\n--- Runtime per input cell (on {n_common_5} common successes, all five methods) ---")
    print_runtime_per_cell(common_5, "all five methods")

    # Common success rate and runtime/cell — excluding TLC
    common_4 = compute_common_successes(fof_data, tlc_data, ve_data, sae_data, fol_data,
                                        ignore_tlc=True)
    n_common_4 = len(common_4[0])
    print(f"\n--- Common success rate (FOF, VE, SAE, FOL only) ---")
    print(f"  {n_common_4} / {TOTAL_CASES} ({100 * n_common_4 / TOTAL_CASES:.1f}%)")
    print(f"\n--- Runtime per input cell (on {n_common_4} common successes, excl. TLC) ---")
    print_runtime_per_cell(common_4, "excl. TLC")

    # 10 shortest-runtime TLC failures
    tlc_failures = sorted(
        (d for d in tlc_data if not is_success(d)),
        key=lambda d: d["runtime_s"]
    )[:10]
    print(f"\n--- 10 shortest-runtime TLC failures ---")
    if not tlc_failures:
        print("  (none)")
    else:
        print(f"  {'subdir':<30}  {'idx':>5}  {'n_init_cells':>12}  {'runtime_s':>10}  "
              f"{'n_degen':>7}  {'n_flipped':>9}")
        for d in tlc_failures:
            stats = d.get("stats", {})
            print(f"  {d['subdir']:<30}  {d['idx']:>5}  {d['n_init_cells']:>12}  "
                  f"{d['runtime_s']:>10.4f}  "
                  f"{stats.get('n_degenerate_cells', '?'):>7}  "
                  f"{stats.get('n_flipped_cells', '?'):>9}")

    # Growth ratio stats (VE, SAE, FOL only) on common successes (excl. TLC)
    print(f"\n--- Growth ratio (final / initial vertices) (on {n_common_4} common successes, excl. TLC) ---")
    for method, data in [("DAF", common_4[2]), ("SAE", common_4[3]), ("FOL", common_4[4])]:
        valid = [d for d in data if d is not None and d.get("n_final_vertices") is not None
                 and d.get("n_init_vertices") not in (None, 0)]
        if not valid:
            print(f"  {method}: no data")
            continue
        ratios = [d["n_final_vertices"] / d["n_init_vertices"] for d in valid]
        print(f"  {method}: avg = {sum(ratios)/len(ratios):.4f} | "
              f"max = {max(ratios):.4f} | "
              f"min = {min(ratios):.4f}")

    # Per-block FOL analysis (only when FOL data exists)
    if fol_data:
        from collections import defaultdict
        fol_by_block = defaultdict(list)
        for d in fol_data:
            fol_by_block[d["subdir"]].append(d)

        RUNTIME_THRESHOLDS_H = [6, 12, 24]
        blocks_with_failure    = []
        blocks_over_runtime    = {t: [] for t in RUNTIME_THRESHOLDS_H}
        blocks_with_large_case = []

        # A block "failed" if any case is missing a FOL JSON (timeout) or has a failure result
        # We detect missing files by comparing FOL case count against FOF/TLC case count per block
        fof_by_block = {}
        for d in fof_data:
            fof_by_block.setdefault(d["subdir"], set()).add(d["idx"])

        for block, entries in sorted(fol_by_block.items()):
            fol_idxs = {e["idx"] for e in entries}
            fof_idxs = fof_by_block.get(block, set())
            has_missing = bool(fof_idxs - fol_idxs)  # cases in FOF but not in FOL
            has_failure = any(not is_success(e) for e in entries)
            if has_missing or has_failure:
                blocks_with_failure.append(block)
            total_runtime_h = sum(e["runtime_s"] for e in entries) / 3600
            for t in RUNTIME_THRESHOLDS_H:
                if total_runtime_h > t:
                    blocks_over_runtime[t].append(block)
            if any(e.get("n_init_cells", 0) > 1000 for e in entries):
                blocks_with_large_case.append(block)

        n_blocks = len(fol_by_block)
        print(f"\n--- FOL per-block analysis ({n_blocks} blocks) ---")

        # List individual FOL failures with their source block
        fol_failures = [(d["subdir"], d["idx"]) for d in fol_data if not is_success(d)]
        print(f"  FOL failures ({len(fol_failures)} total):")
        if fol_failures:
            for subdir, idx in sorted(fol_failures):
                print(f"    block={subdir}, idx={idx}")
        else:
            print("    (none)")

        print(f"  Blocks with at least one FOL failure: "
              f"{len(blocks_with_failure)} / {n_blocks}")
        for t in RUNTIME_THRESHOLDS_H:
            print(f"  Blocks with total FOL runtime > {t}h: "
                  f"{len(blocks_over_runtime[t])} / {n_blocks}")
        print(f"  Blocks with at least one case > 1000 cells: "
              f"{len(blocks_with_large_case)} / {n_blocks}")

        # Block with highest cumulative runtime
        block_runtimes = {block: sum(e["runtime_s"] for e in entries)
                          for block, entries in fol_by_block.items()}
        worst_block = max(block_runtimes, key=block_runtimes.get)
        print(f"  Block with highest cumulative FOL runtime: {worst_block} "
              f"({block_runtimes[worst_block]/3600:.2f}h, "
              f"{len(fol_by_block[worst_block])} cases)")


def print_daf_variants_summary(input_dir, fof_data, tlc_data, ve_data):
    """Load and compare DAF, DAF-preremeshing-only, and DAF-no-remeshing.

    Restricts analysis to cases where TLC or FOF succeeded (no remeshing needed).
    Reports success rate, growth ratio, and runtime/cell for each variant,
    where growth/runtime are computed only on cases where all three DAF variants
    AND (TLC or FOF) succeeded.
    """
    print("\n--- DAF variants analysis ---")

    # Cases where TLC or FOF succeeded
    no_remesh_keys = (
        {(d["subdir"], d["idx"]) for d in fof_data if is_success(d)} |
        {(d["subdir"], d["idx"]) for d in tlc_data if is_success(d)}
    )
    n_no_remesh = len(no_remesh_keys)
    print(f"  Base set: {n_no_remesh} cases where TLC or FOF succeeded")

    # Load the two extra variants
    pre_remesh  = collect_daf_variant(input_dir, "remeshing-maps-data-preremeshing-only")
    no_remesh_v = collect_daf_variant(input_dir, "remeshing-maps-data-no-remeshing")

    # Build lookup for the original DAF
    daf_orig = {(d["subdir"], d["idx"]): d for d in ve_data}

    variants = [
        ("DAF (full)",             daf_orig),
        ("DAF (pre-remesh only)",  pre_remesh),
        ("DAF (no remeshing)",     no_remesh_v),
    ]

    # Success rates within the no-remesh base set
    print(f"\n  Success rates (within {n_no_remesh} cases where TLC or FOF succeeded):")
    for name, lookup in variants:
        n_success = sum(1 for k in no_remesh_keys if lookup.get(k) and is_success(lookup[k]))
        pct = 100 * n_success / n_no_remesh if n_no_remesh else float("nan")
        print(f"    {name}: {n_success} / {n_no_remesh} ({pct:.1f}%)")
    fof_lookup = {(d["subdir"], d["idx"]): d for d in fof_data}
    tlc_lookup = {(d["subdir"], d["idx"]): d for d in tlc_data}
    for name, lookup in [("FOF", fof_lookup), ("TLC", tlc_lookup)]:
        n_success = sum(1 for k in no_remesh_keys if lookup.get(k) and is_success(lookup[k]))
        pct = 100 * n_success / n_no_remesh if n_no_remesh else float("nan")
        print(f"    {name}: {n_success} / {n_no_remesh} ({pct:.1f}%)")

    # Common successes: all three DAF variants succeed AND TLC or FOF succeeded
    common_keys = {
        k for k in no_remesh_keys
        if (daf_orig.get(k) and is_success(daf_orig[k]) and
            pre_remesh.get(k) and is_success(pre_remesh[k]) and
            no_remesh_v.get(k) and is_success(no_remesh_v[k]))
    }
    n_common = len(common_keys)
    print(f"\n  Common successes (all three DAF variants + TLC or FOF): {n_common} cases")

    if n_common == 0:
        print("  No common successes — skipping growth ratio and runtime/cell stats.")
        return

    # FOF and TLC runtime/cell on common successes
    print(f"\n  FOF and TLC runtime/cell (on {n_common} common successes):")
    for name, lookup in [("FOF", fof_lookup), ("TLC", tlc_lookup)]:
        entries = [lookup[k] for k in common_keys if lookup.get(k)]
        valid = [e for e in entries if e.get("n_init_cells") not in (None, 0)]
        if not valid:
            print(f"    {name}: no data")
            continue
        rates = [e["runtime_s"] / e["n_init_cells"] for e in valid]
        print(f"    {name}: avg = {sum(rates)/len(rates):.6f} s/cell | "
              f"max = {max(rates):.6f} s/cell")

    # Growth ratio and runtime/cell on common successes
    print(f"\n  DAF variants growth ratio and runtime/cell (on {n_common} common successes):")
    for name, lookup in variants:
        entries = [lookup[k] for k in common_keys]
        valid_growth = [e for e in entries
                        if e.get("n_final_vertices") is not None
                        and e.get("n_init_vertices") not in (None, 0)
                        and e.get("n_init_cells") not in (None, 0)]
        if not valid_growth:
            print(f"    {name}: no data")
            continue
        ratios = [e["n_final_vertices"] / e["n_init_vertices"] for e in valid_growth]
        rates  = [e["runtime_s"] / e["n_init_cells"] for e in valid_growth]
        print(f"    {name}:")
        print(f"      growth ratio: avg = {sum(ratios)/len(ratios):.4f} | "
              f"max = {max(ratios):.4f} | min = {min(ratios):.4f}")
        print(f"      runtime/cell: avg = {sum(rates)/len(rates):.6f} s/cell | "
              f"max = {max(rates):.6f} s/cell")


# ---------------------------------------------------------------------------
# Spine data collection (AIS-1 / AIS-2, FSU cases, DAF + FOF + TLC only)
# ---------------------------------------------------------------------------


def collect_data_spine(input_dir):
    """
    Walk input_dir -> AIS-N subdirs -> FSU files.

    Directory structure:
        <input_dir>/
        +-- AIS-1/
        |   +-- FSU-XX_FOF.time
        |   +-- FSU-XX_FOF_stats.txt
        |   +-- FSU-XX_target_vol_TLC.time
        |   +-- FSU-XX_target_vol_TLC_stats.txt
        |   +-- FSU-XX_output.json               (DAF/VE)
        +-- AIS-2/
            +-- ...

    Cases are identified by (ais_dir, XX) and labelled "AIS-N_FSU_XX".
    Returns fof_data, tlc_data, ve_data (SAE and FOL are empty lists).
    """
    fof_data = []
    tlc_data = []
    ve_data  = []

    input_path = Path(input_dir)
    if not input_path.is_dir():
        print(f"Error: {input_dir} is not a valid directory.")
        sys.exit(1)

    ais_dirs = sorted([p for p in input_path.iterdir()
                       if p.is_dir() and re.match(r"AIS-\d+", p.name)])
    if not ais_dirs:
        print(f"Error: no AIS-N subdirectories found in {input_dir}.")
        sys.exit(1)

    print(f"Found {len(ais_dirs)} AIS directory(ies): {[p.name for p in ais_dirs]}")

    # DAF files live in remeshing-maps-data/ alongside the AIS dirs
    daf_dir = input_path / "remeshing-maps-data"
    if not daf_dir.is_dir():
        print(f"  Warning: remeshing-maps-data/ not found in {input_dir}. "
              f"DAF data will be missing.")
    # Build a lookup: ID -> Path for all DAF json files
    daf_files = {}
    if daf_dir.is_dir():
        for f in daf_dir.iterdir():
            if f.is_file() and f.name.endswith("_output.json"):
                xx = f.name[:-len("_output.json")]
                daf_files[xx] = f
        print(f"  remeshing-maps-data: {len(daf_files)} DAF file(s) found")

    # Seed all cases from DAF files; FOF/TLC will be added from AIS dirs when available
    cases = {}
    for xx, f in daf_files.items():
        cases.setdefault(xx, {})["VE_json"] = f

    for ais_dir in ais_dirs:
        print(f"\nScanning {ais_dir.name} for FOF/TLC files...")
        for f in ais_dir.iterdir():
            if not f.is_file():
                continue
            name = f.name
            if name.endswith("_FOF.time"):
                xx = name[:-len("_FOF.time")]
                cases.setdefault(xx, {})["FOF_time"] = f
            elif name.endswith("_FOF_stats.txt"):
                xx = name[:-len("_FOF_stats.txt")]
                cases.setdefault(xx, {})["FOF_stats"] = f
            elif name.endswith("_TLC.time"):
                xx = name[:-len("_TLC.time")]
                cases.setdefault(xx, {})["TLC_time"] = f
            elif name.endswith("_TLC_stats.txt"):
                xx = name[:-len("_TLC_stats.txt")]
                cases.setdefault(xx, {})["TLC_stats"] = f

    print(f"  Found {len(cases)} case(s) total.")

    # Report cases that are missing data for one or more methods
    mismatches = []
    for xx, file_map in cases.items():
        missing = []
        if "VE_json" not in file_map:
            missing.append("DAF")
        has_fof = "FOF_time" in file_map and "FOF_stats" in file_map
        has_tlc = "TLC_time" in file_map and "TLC_stats" in file_map
        if not has_fof and ("FOF_time" in file_map or "FOF_stats" in file_map):
            missing.append("FOF (partial)")
        if not has_tlc and ("TLC_time" in file_map or "TLC_stats" in file_map):
            missing.append("TLC (partial)")
        if missing:
            mismatches.append((xx, missing))
    if mismatches:
        print(f"  Cases with missing method data ({len(mismatches)}):")
        for xx, missing in sorted(mismatches):
            print(f"    {xx}: missing {missing}")

    for xx, file_map in cases.items():
        case_id = xx
        fof_stats = None
        tlc_stats = None
        ve_raw    = None

        if "FOF_stats" in file_map:
            fof_stats = parse_stats_txt(file_map["FOF_stats"])
        if "TLC_stats" in file_map:
            tlc_stats = parse_stats_txt(file_map["TLC_stats"])
        if "VE_json" in file_map:
            ve_raw = parse_json_file(file_map["VE_json"])

        # n_init_vertices: from FOF/TLC stats
        n_init_vertices = None
        if fof_stats and "n_init_vertices" in fof_stats:
            n_init_vertices = fof_stats["n_init_vertices"]
        elif tlc_stats and "n_init_vertices" in tlc_stats:
            n_init_vertices = tlc_stats["n_init_vertices"]

        # n_init_cells: sanity-check FOF vs TLC; fall back to VE JSON
        cell_sources = {}
        if fof_stats and "n_init_cells" in fof_stats:
            cell_sources["FOF"] = fof_stats["n_init_cells"]
        if tlc_stats and "n_init_cells" in tlc_stats:
            cell_sources["TLC"] = tlc_stats["n_init_cells"]
        unique_cell_values = set(cell_sources.values())
        if len(unique_cell_values) > 1:
            print(f"  Warning [{case_id}] n_init_cells mismatch: {cell_sources}")
        n_init_cells = next(iter(unique_cell_values)) if unique_cell_values else None
        # Fall back to VE JSON if stats not available
        if n_init_cells is None and ve_raw is not None:
            n_init_cells = ve_raw.get("n_init_cells")
            if n_init_cells is None:
                print(f"  Warning [{case_id}] VE JSON keys: {list(ve_raw.keys())[:10]}")
        if n_init_cells is None:
            print(f"  Warning [{case_id}] could not determine n_init_cells, skipping.")
            continue

        # FOF
        if "FOF_time" in file_map and fof_stats is not None:
            runtime = parse_time_file(file_map["FOF_time"])
            if runtime is not None:
                fof_data.append({
                    "subdir": case_id,
                    "idx": xx,
                    "n_init_vertices": n_init_vertices,
                    "n_init_cells": n_init_cells,
                    "runtime_s": runtime,
                    "success": is_success_stats(fof_stats, n_init_cells),
                    "stats": fof_stats,
                })

        # TLC
        if "TLC_time" in file_map and tlc_stats is not None:
            runtime = parse_time_file(file_map["TLC_time"])
            if runtime is not None:
                tlc_data.append({
                    "subdir": case_id,
                    "idx": xx,
                    "n_init_vertices": n_init_vertices,
                    "n_init_cells": n_init_cells,
                    "runtime_s": runtime,
                    "success": is_success_stats(tlc_stats, n_init_cells),
                    "stats": tlc_stats,
                })

        # DAF/VE
        if ve_raw is not None:
            runtime = ve_raw.get("runtime_s")
            if runtime is not None:
                ve_data.append({
                    "subdir": case_id,
                    "idx": xx,
                    "n_init_vertices": n_init_vertices,
                    "n_init_cells": n_init_cells,
                    "n_final_vertices": ve_raw.get("n_final_vertices"),
                    "runtime_s": runtime,
                    "success": ve_raw.get("opti_result") == 0,
                    "raw": ve_raw,
                })
            else:
                print(f"  [{case_id}] VE JSON missing key: runtime_s")

    return fof_data, tlc_data, ve_data


def compute_common_successes_spine(fof_data, tlc_data, ve_data, ignore_tlc=False):
    """Common successes for spine mode (DAF, FOF, TLC only).
    Returns the same 5-tuple shape as compute_common_successes, with empty SAE/FOL lists."""
    def success_map(data):
        return {(d["subdir"], d["idx"]): d for d in data if is_success(d)}

    fof_s = success_map(fof_data)
    tlc_s = success_map(tlc_data)
    ve_s  = success_map(ve_data)

    common_keys = set(fof_s) & set(ve_s)
    if not ignore_tlc:
        common_keys &= set(tlc_s)

    return (
        [fof_s[k] for k in common_keys],
        [tlc_s.get(k) for k in common_keys],
        [ve_s[k] for k in common_keys],
        [],  # SAE
        [],  # FOL
    )


def print_summary_spine(fof_data, tlc_data, ve_data, total_cases=None, title="spine"):
    n_actual = len({(d["subdir"], d["idx"])
                    for data in (fof_data, tlc_data, ve_data) for d in data})
    if total_cases is None:
        total_cases = n_actual
    all_methods = [("DAF", ve_data), ("FOF", fof_data), ("TLC", tlc_data)]

    print(f"\n--- Summary ({title}) ---")
    for label, data in all_methods:
        if data:
            runtimes = [d["runtime_s"] for d in data]
            cells = [d["n_init_cells"] for d in data]
            avg_runtime_s = (sum(runtimes)/len(runtimes))
            print(
                f"  {label}: {len(data)} cases | "
                f"runtime [{min(runtimes):.3f}s – {max(runtimes):.3f}s] (avg: {avg_runtime_s:.3}s / {avg_runtime_s/60:.3}min) | "
                f"n_init_cells [{min(cells)} – {max(cells)}]"
            )
        else:
            print(f"  {label}: 0 cases")

    print(f"\n--- Success rates (out of {TOTAL_CASES} total cases) ---")
    for label, data in all_methods:
        n_success = sum(1 for d in data if is_success(d))
        print(f"  {label}: {n_success} / {TOTAL_CASES} ({100 * n_success / TOTAL_CASES:.1f}%)")

    def print_runtime_per_cell(common, label):
        for method, data in [("DAF", common[2]), ("FOF", common[0]), ("TLC", common[1])]:
            if not data or any(d is None for d in data):
                continue
            rates = [d["runtime_s"] / d["n_init_cells"] for d in data]
            runtimes = [d["runtime_s"] for d in data]
            print(f"  {method}: avg = {sum(rates)/len(rates):.6f} s/cell | "
                  f"max = {max(rates):.6f} s/cell || avg = {sum(runtimes)/(len(runtimes)*60):.6f} min | "
                  f"max = {max(runtimes)/60:.6f} min" )

    common_3 = compute_common_successes_spine(fof_data, tlc_data, ve_data, ignore_tlc=False)
    n_common_3 = len(common_3[0])
    print(f"\n--- Common success rate (DAF, FOF, TLC) ---")
    print(f"  {n_common_3} / {TOTAL_CASES} ({100 * n_common_3 / TOTAL_CASES:.1f}%)")
    print(f"\n--- Runtime per input cell (on {n_common_3} common successes) ---")
    print_runtime_per_cell(common_3, "all three")

    common_2 = compute_common_successes_spine(fof_data, tlc_data, ve_data, ignore_tlc=True)
    n_common_2 = len(common_2[0])
    print(f"\n--- Common success rate (DAF, FOF only) ---")
    print(f"  {n_common_2} / {TOTAL_CASES} ({100 * n_common_2 / TOTAL_CASES:.1f}%)")
    print(f"\n--- Runtime per input cell (on {n_common_2} common successes, excl. TLC) ---")
    print_runtime_per_cell(common_2, "excl. TLC")

    # 10 shortest-runtime TLC failures
    tlc_failures = sorted(
        (d for d in tlc_data if not is_success(d)),
        key=lambda d: d["runtime_s"]
    )[:10]
    print(f"\n--- 10 shortest-runtime TLC failures ---")
    if not tlc_failures:
        print("  (none)")
    else:
        print(f"  {'case':<25}  {'n_init_cells':>12}  {'runtime_s':>10}  "
              f"{'n_degen':>7}  {'n_flipped':>9}")
        for d in tlc_failures:
            stats = d.get("stats", {})
            print(f"  {d['subdir']:<25}  {d['n_init_cells']:>12}  "
                  f"{d['runtime_s']:>10.4f}  "
                  f"{stats.get('n_degenerate_cells', '?'):>7}  "
                  f"{stats.get('n_flipped_cells', '?'):>9}")

    def print_growth_ratio(common, label):
        valid = [d for d in common[2] if d is not None
                 and d.get("n_final_vertices") is not None
                 and d.get("n_init_vertices") not in (None, 0)]
        if valid:
            ratios = [d["n_final_vertices"] / d["n_init_vertices"] for d in valid]
            print(f"  DAF: avg = {sum(ratios)/len(ratios):.4f} | "
                  f"max = {max(ratios):.4f} | min = {min(ratios):.4f}")
        else:
            print("  DAF: no data")

    # Growth ratio (DAF only)
    print(f"\n--- DAF growth ratio (on {n_common_3} common successes) ---")
    print_growth_ratio(common_3, "all three")
    print(f"\n--- DAF growth ratio (on {n_common_2} common successes, excl. TLC) ---")
    print_growth_ratio(common_2, "excl. TLC")


# ---------------------------------------------------------------------------
# Volmap data collection (G1_star / G3_star / star, string IDs)
# ---------------------------------------------------------------------------

def collect_data_volmap(input_dir):
    """
    Walk input_dir -> G1_star and G3_star for FOF/TLC, star/remeshing-maps-data for DAF.

    Directory structure:
        <input_dir>/
        +-- G1_star/
        |   +-- <ID>_[star_]FOF.time
        |   +-- <ID>_[star_]FOF_stats.txt
        |   +-- <ID>_[star_]TLC.time
        |   +-- <ID>_[star_]TLC_stats.txt
        +-- G3_star/
        |   +-- ... (same pattern)
        +-- star/
            +-- remeshing-maps-data/
                +-- <ID>_star_boundary_deformation_map_output.json   (DAF)

    Cases are identified by (group, ID). G1 and G3 cases are merged.
    Returns fof_data, tlc_data, ve_data (SAE and FOL are empty lists).
    """
    fof_data = []
    tlc_data = []
    ve_data  = []

    input_path = Path(input_dir)
    if not input_path.is_dir():
        print(f"Error: {input_dir} is not a valid directory.")
        sys.exit(1)

    daf_dir = input_path / "star" / "remeshing-maps-data"
    if not daf_dir.is_dir():
        print(f"  Warning: star/remeshing-maps-data/ not found in {input_dir}. "
              f"DAF data will be missing.")

    group_dirs = [p for p in input_path.iterdir()
                  if p.is_dir() and p.name in ("G1_star",)]  # TODO: re-enable G3_star
    if not group_dirs:
        print(f"Error: no G1_star or G3_star directories found in {input_dir}.")
        sys.exit(1)

    print(f"Found group directories: {[p.name for p in sorted(group_dirs)]}")

    for group_dir in sorted(group_dirs):
        print(f"\nProcessing {group_dir.name}...")
        cases = {}  # ID (str) -> dict of file paths

        for f in group_dir.iterdir():
            if not f.is_file():
                continue
            name = f.name
            # Match <ID>_[star_]FOF.time  etc.
            # The ID is everything before the last _[star_]FOF/TLC suffix
            for suffix, key in [
                (r"^(.+?)_(?:star_)?FOF\.time$",       "FOF_time"),
                (r"^(.+?)_(?:star_)?FOF_stats\.txt$",  "FOF_stats"),
                (r"^(.+?)_(?:star_)?TLC\.time$",       "TLC_time"),
                (r"^(.+?)_(?:star_)?TLC_stats\.txt$",  "TLC_stats"),
            ]:
                m = re.match(suffix, name)
                if m:
                    idx = m.group(1)
                    cases.setdefault(idx, {})[key] = f
                    break

        # DAF lookup temporarily disabled
        #TODO: re-enable once DAF data is verified
        if daf_dir.is_dir():
            for idx in list(cases.keys()):
                daf_file = daf_dir / f"{idx}_star_boundary_deformation_map_output.json"
                if daf_file.is_file():
                    cases[idx]["VE_json"] = daf_file

        print(f"  Found {len(cases)} case(s).")

        for idx, file_map in cases.items():
            case_id = f"{group_dir.name}_{idx}"
            fof_stats = None
            tlc_stats = None
            ve_raw    = None

            if "FOF_stats" in file_map:
                fof_stats = parse_stats_txt(file_map["FOF_stats"])
            if "TLC_stats" in file_map:
                tlc_stats = parse_stats_txt(file_map["TLC_stats"])
            if "VE_json" in file_map:
                ve_raw = parse_json_file(file_map["VE_json"])

            n_init_vertices = None
            if fof_stats and "n_init_vertices" in fof_stats:
                n_init_vertices = fof_stats["n_init_vertices"]
            elif tlc_stats and "n_init_vertices" in tlc_stats:
                n_init_vertices = tlc_stats["n_init_vertices"]

            cell_sources = {}
            if fof_stats and "n_init_cells" in fof_stats:
                cell_sources["FOF"] = fof_stats["n_init_cells"]
            if tlc_stats and "n_init_cells" in tlc_stats:
                cell_sources["TLC"] = tlc_stats["n_init_cells"]
            unique_cell_values = set(cell_sources.values())
            if len(unique_cell_values) > 1:
                print(f"  Warning [{case_id}] n_init_cells mismatch: {cell_sources}")
            n_init_cells = next(iter(unique_cell_values)) if unique_cell_values else None
            if n_init_cells is None:
                print(f"  Warning [{case_id}] could not determine n_init_cells, skipping.")
                continue

            if "FOF_time" in file_map and fof_stats is not None:
                runtime = parse_time_file(file_map["FOF_time"])
                if runtime is not None:
                    fof_data.append({
                        "subdir": case_id,
                        "idx": idx,
                        "n_init_vertices": n_init_vertices,
                        "n_init_cells": n_init_cells,
                        "runtime_s": runtime,
                        "success": is_success_stats(fof_stats, n_init_cells),
                        "stats": fof_stats,
                    })
            else:
                missing = [k for k in ("FOF_time", "FOF_stats") if k not in file_map]
                if missing:
                    print(f"  [{case_id}] Missing FOF files: {missing}")

            if "TLC_time" in file_map and tlc_stats is not None:
                runtime = parse_time_file(file_map["TLC_time"])
                if runtime is not None:
                    tlc_data.append({
                        "subdir": case_id,
                        "idx": idx,
                        "n_init_vertices": n_init_vertices,
                        "n_init_cells": n_init_cells,
                        "runtime_s": runtime,
                        "success": is_success_stats(tlc_stats, n_init_cells),
                        "stats": tlc_stats,
                    })
            else:
                missing = [k for k in ("TLC_time", "TLC_stats") if k not in file_map]
                if missing:
                    print(f"  [{case_id}] Missing TLC files: {missing}")

            if ve_raw is not None:
                runtime = ve_raw.get("runtime_s")
                if runtime is not None:
                    ve_data.append({
                        "subdir": case_id,
                        "idx": idx,
                        "n_init_vertices": n_init_vertices,
                        "n_init_cells": n_init_cells,
                        "n_final_vertices": ve_raw.get("n_final_vertices"),
                        "runtime_s": runtime,
                        "success": ve_raw.get("opti_result") == 0,
                        "raw": ve_raw,
                    })
                else:
                    print(f"  [{case_id}] DAF JSON missing key: runtime_s")
            else:
                print(f"  [{case_id}] Missing DAF JSON file.")

    return fof_data, tlc_data, ve_data

# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

PLOT_WIDTH_PX = 1400
PLOT_HEIGHT_PX = 500
PLOT_DPI = 150
PLOT_XLIM_STARS = (1e1, 1e3)
PLOT_XLIM_SPINE = (1e5, 1e6)

def make_scatter_plot(fof_data, tlc_data, ve_data, sae_data, fol_data,
                      plot_output_path, show_plots, marker_size=25, ignore_tlc=False,
                      xlim=PLOT_XLIM_STARS, use_index=False):
    """Generate, save, and optionally show the two-panel scatter plot."""
    figsize = (PLOT_WIDTH_PX / PLOT_DPI, 2 * PLOT_HEIGHT_PX / PLOT_DPI)
    fig, (ax_runtime, ax_growth) = plt.subplots(2, 1, figsize=figsize)

    method_cfg = [
        ("DAF (ours)",  ve_data,  "#5CBF7A", "^"),
        ("FOF", fof_data, "#E05C5C", "o"),
        ("TLC", tlc_data, "#5C9BE0", "s"),
        ("FOL", fol_data, "#A05CE0", "P"),
        ("SAE", sae_data, "#E0A83C", "D"),
    ]
    if ignore_tlc:
        method_cfg = [(l, d, c, m) for l, d, c, m in method_cfg if l != "TLC"]

    # --- Top panel: runtime vs n_init_cells (or case index) ---
    for label, data, color, marker in method_cfg:
        if not data or any(d is None for d in data):
            continue
        xs = list(range(len(data))) if use_index else [d["n_init_cells"] for d in data]
        ys = [d["runtime_s"] for d in data]
        ax_runtime.scatter(
            xs, ys,
            label=label,
            color=color,
            marker=marker,
            alpha=0.99,
            edgecolors="white",
            linewidths=0.5,
            s=marker_size,
        )

    fontsize = 16
    ax_runtime.set_xscale("linear" if xlim == PLOT_XLIM_SPINE else "log")
    ax_runtime.set_yscale("log")
    ax_runtime.set_xlabel("case index" if use_index else "#cells", fontsize=fontsize)
    ax_runtime.set_ylabel("runtime (s)", fontsize=fontsize)
    ax_runtime.tick_params(axis="both", labelsize=fontsize)
    ax_runtime.set_ylim(top=ax_runtime.get_ylim()[1] * 10)
    if not use_index:
        ax_runtime.set_xlim(*xlim)
    all_handles, all_labels = ax_runtime.get_legend_handles_labels()
    if all_handles:
        legend = ax_runtime.legend(
            all_handles, all_labels,
            loc="upper center", ncol=len(all_handles),
            fontsize=fontsize - 2, framealpha=0.9,
        )
        for handle in legend.legend_handles:
            handle.set_sizes([80])
    ax_runtime.grid(True, linestyle="--", alpha=0.4, which="major")

    # --- Bottom panel: growth ratio vs n_init_cells (DAF, SAE, FOL only) ---
    growth_cfg = [(l, d, c, m) for l, d, c, m in method_cfg
                  if l in ("DAF (ours)", "SAE", "FOL")]
    for label, data, color, marker in growth_cfg:
        if not data or any(d is None for d in data):
            continue
        valid = [d for d in data
                 if d.get("n_final_vertices") is not None
                 and d.get("n_init_vertices") not in (None, 0)]
        points = [(i if use_index else d["n_init_cells"],
                   d["n_final_vertices"] / d["n_init_vertices"])
                  for i, d in enumerate(valid)]
        if not points:
            continue
        xs, ys = zip(*points)
        ax_growth.scatter(
            xs, ys,
            label=label,
            color=color,
            marker=marker,
            alpha=0.99,
            edgecolors="white",
            linewidths=0.5,
            s=marker_size,
        )

    ax_growth.set_xscale("linear" if xlim == PLOT_XLIM_SPINE else "log")
    ax_growth.set_yscale("log")
    ax_growth.set_xlabel("case index" if use_index else "#cells", fontsize=fontsize)
    ax_growth.set_ylabel("growth ratio", fontsize=fontsize)
    ax_growth.tick_params(axis="both", labelsize=fontsize)
    if not use_index:
        ax_growth.set_xlim(*xlim)
    ax_growth.grid(True, linestyle="--", alpha=0.4, which="major")

    fig.tight_layout()
    fig.savefig(plot_output_path, dpi=PLOT_DPI)
    print(f"\nPlot saved to: {plot_output_path}")

    if show_plots:
        plt.show()


# ---------------------------------------------------------------------------
# Remeshing-necessity plot
# ---------------------------------------------------------------------------

def make_remeshing_necessity_plot(fof_data, tlc_data, ve_data, plot_output_path, show_plots,
                                  xlim=PLOT_XLIM_STARS, use_index=False):
    """Plot DAF growth ratios split by whether TLC or FOF succeeded (no remeshing needed)
    vs. both failed (remeshing necessary). Single scatter plot with stats in legend."""

    fof_success_keys = {(d["subdir"], d["idx"]) for d in fof_data if is_success(d)}
    tlc_success_keys = {(d["subdir"], d["idx"]) for d in tlc_data if is_success(d)}
    no_remesh_keys   = fof_success_keys | tlc_success_keys

    no_remesh, remesh_needed = [], []
    for d in ve_data:
        if d.get("n_final_vertices") is None or d.get("n_init_vertices") in (None, 0):
            continue
        ratio = d["n_final_vertices"] / d["n_init_vertices"]
        entry = {"n_init_cells": d["n_init_cells"], "ratio": ratio}
        key = (d["subdir"], d["idx"])
        if key in no_remesh_keys:
            no_remesh.append(entry)
        else:
            remesh_needed.append(entry)

    ratios_no  = [e["ratio"] for e in no_remesh]
    ratios_yes = [e["ratio"] for e in remesh_needed]
    cells_no   = [e["n_init_cells"] for e in no_remesh]
    cells_yes  = [e["n_init_cells"] for e in remesh_needed]

    n_no_gt1    = sum(1 for r in ratios_no  if r > 1)
    n_yes_gt1   = sum(1 for r in ratios_yes if r > 1)
    pct_no_gt1  = 100 * n_no_gt1  / len(ratios_no)  if ratios_no  else float("nan")
    pct_yes_gt1 = 100 * n_yes_gt1 / len(ratios_yes) if ratios_yes else float("nan")

    print(f"\n--- Remeshing necessity split (DAF growth ratio) ---")
    print(f"  No remeshing needed (FOF or TLC succeeded): {len(no_remesh)} cases")
    print(f"  Remeshing needed   (both FOF and TLC failed): {len(remesh_needed)} cases")
    if ratios_no:
        print(f"  Growth ratio > 1 -- no remeshing needed: "
              f"{n_no_gt1} / {len(ratios_no)} ({pct_no_gt1:.1f}%)")
    if ratios_yes:
        print(f"  Growth ratio > 1 -- remeshing needed:    "
              f"{n_yes_gt1} / {len(ratios_yes)} ({pct_yes_gt1:.1f}%)")

    COLOR_NO_REMESH = "#5C9BE0"
    COLOR_REMESH    = "#E05C5C"
    LABEL_NO_REMESH = "no remeshing needed (FOF or TLC succeeded)"
    LABEL_REMESH    = "remeshing needed (both FOF and TLC failed)"

    def avg(vals):
        return sum(vals) / len(vals) if vals else float("nan")

    mu_no      = avg(ratios_no)
    mu_no_gt1  = avg([r for r in ratios_no  if r > 1])
    mu_yes     = avg(ratios_yes)
    mu_yes_gt1 = avg([r for r in ratios_yes if r > 1])

    label_no = (
        f"$\mu$ = {mu_no:.3f} \n"
        f"$\hat{{\mu}}$ = {mu_no_gt1:.3f}"
    )
    label_yes = (
        f"$\mu$ = {mu_yes:.3f}\n"
        f"$\hat{{\mu}}$ = {mu_yes_gt1:.3f}"
    )

    figsize = (PLOT_WIDTH_PX / PLOT_DPI, PLOT_HEIGHT_PX / PLOT_DPI)
    fig, ax = plt.subplots(1, 1, figsize=figsize)

    fontsize = 14
    xs_no  = list(range(len(no_remesh)))     if use_index else cells_no
    xs_yes = list(range(len(remesh_needed))) if use_index else cells_yes
    ax.scatter(xs_no,  ratios_no,  color=COLOR_NO_REMESH, alpha=0.6,
               s=15, label=label_no,  edgecolors="white", linewidths=0.5)
    ax.scatter(xs_yes, ratios_yes, color=COLOR_REMESH,    alpha=1.0,
               s=15, label=label_yes, edgecolors="white", linewidths=0.5)

    ax.set_xscale("linear" if xlim == PLOT_XLIM_SPINE else "log")
    ax.set_yscale("log")
    if not use_index:
        ax.set_xlim(*xlim)
    ax.set_xlabel("case index" if use_index else "#cells", fontsize=fontsize)
    ax.set_ylabel("growth ratio", fontsize=fontsize)
    ax.tick_params(labelsize=fontsize - 2)
    ax.legend(fontsize=fontsize - 3, loc="upper right")
    ax.grid(True, linestyle="--", alpha=0.4)

    if not ratios_no and not ratios_yes:
        print("  Skipping remeshing-necessity plot: no data.")
        plt.close(fig)
        return
    fig.tight_layout()
    fig.savefig(plot_output_path, dpi=PLOT_DPI)
    print(f"Remeshing-necessity plot saved to: {plot_output_path}")

    if show_plots:
        plt.show()


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Analyze mesh processing results and plot runtime vs. mesh size."
    )
    parser.add_argument(
        "--input-directory",
        required=True,
        help="Path to the main directory containing result subdirectories.",
    )
    parser.add_argument(
        "--plot-output-path",
        required=True,
        help="Path to save the plot image (e.g. plot.png).",
    )
    parser.add_argument(
        "--show-plots",
        action="store_true",
        help="If set, also display the plot interactively after saving.",
    )
    parser.add_argument(
        "--ignore-TLC",
        action="store_true",
        help="If set, skip plotting TLC data.",
    )
    parser.add_argument(
        "--remeshing-necessity-plot-path",
        default=None,
        help="If set, generate the remeshing-necessity plot and save it to this path.",
    )
    parser.add_argument(
        "--analyse-spine",
        action="store_true",
        help="If set, run spine (FSU/AIS) analysis instead of the default star analysis.",
    )
    parser.add_argument(
        "--analyse-volmap",
        action="store_true",
        help="If set, run volmap (G1_star/G3_star) analysis instead of the default star analysis.",
    )
    args = parser.parse_args()

    if args.analyse_volmap:
        fof_data, tlc_data, ve_data = collect_data_volmap(args.input_directory)
        print_summary_spine(fof_data, tlc_data, ve_data,
                            total_cases=TOTAL_CASES, title="volmap")
        common = compute_common_successes_spine(fof_data, tlc_data, ve_data,
                                                ignore_tlc=args.ignore_TLC)
        print(f"\nPlotting {len(common[0])} common-success cases "
              f"({'excl. TLC' if args.ignore_TLC else 'all three methods'}).")
        make_scatter_plot(
            *common,
            plot_output_path=args.plot_output_path,
            show_plots=args.show_plots,
            ignore_tlc=args.ignore_TLC,
        )
        if args.remeshing_necessity_plot_path:
            make_remeshing_necessity_plot(
                fof_data, tlc_data, ve_data,
                plot_output_path=args.remeshing_necessity_plot_path,
                show_plots=args.show_plots,
            )
    elif args.analyse_spine:
        fof_data, tlc_data, ve_data = collect_data_spine(args.input_directory)
        print_summary_spine(fof_data, tlc_data, ve_data, total_cases=TOTAL_CASES,
                            title="spine")
        common = compute_common_successes_spine(fof_data, tlc_data, ve_data,
                                                ignore_tlc=args.ignore_TLC)
        print(f"\nPlotting {len(common[0])} common-success cases "
              f"({'excl. TLC' if args.ignore_TLC else 'all three methods'}).")
        make_scatter_plot(
            *common,
            plot_output_path=args.plot_output_path,
            show_plots=args.show_plots,
            ignore_tlc=args.ignore_TLC,
            xlim=PLOT_XLIM_SPINE,
            use_index=True,
        )
        if args.remeshing_necessity_plot_path:
            make_remeshing_necessity_plot(
                fof_data, tlc_data, ve_data,
                plot_output_path=args.remeshing_necessity_plot_path,
                show_plots=args.show_plots,
                xlim=PLOT_XLIM_SPINE,
                use_index=True,
            )
    else:
        fof_data, tlc_data, ve_data, sae_data, fol_data = collect_data(args.input_directory)
        print_summary(fof_data, tlc_data, ve_data, sae_data, fol_data)
        print_daf_variants_summary(args.input_directory, fof_data, tlc_data, ve_data)
        common = compute_common_successes(fof_data, tlc_data, ve_data, sae_data, fol_data,
                                          ignore_tlc=args.ignore_TLC)
        print(f"\nPlotting {len(common[0])} common-success cases "
              f"({'excl. TLC' if args.ignore_TLC else 'all five methods'}).")
        make_scatter_plot(
            *common,
            plot_output_path=args.plot_output_path,
            show_plots=args.show_plots,
            ignore_tlc=args.ignore_TLC,
        )
        if args.remeshing_necessity_plot_path:
            make_remeshing_necessity_plot(
                fof_data, tlc_data, ve_data,
                plot_output_path=args.remeshing_necessity_plot_path,
                show_plots=args.show_plots,
            )


if __name__ == "__main__":
    main()