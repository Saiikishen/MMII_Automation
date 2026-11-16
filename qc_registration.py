#!/usr/bin/env python3
"""Independent, read-only registration checks and anatomical montages.

Example:
    python qc_registration.py --output PATH_TO_RUN

Reads the successful run's manifest.json and config_used.json. The manifest
selects the actual displayed main image, processed source images, and fixed-RAS
to moving-RAS (pull) matrices. These are Greedy/c3d RAS matrices, not ITK LPS
parameter files. Only groups with qc_anatomy=true receive intensity comparisons;
other processed groups receive geometry checks. Missing optional modalities
recorded as skipped by the builder are reported without being loaded.

Dependencies: numpy, scipy, matplotlib. No nibabel dependency is required;
the small NIfTI-1 reader accepts uncompressed .nii and NIfTI .hdr/.img pairs.
Only output/qc is written. Inputs and transforms are never changed.
Metrics and images support human review; they are not clinical validation.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import struct

import numpy as np
from scipy.ndimage import map_coordinates


class NiftiImage:
    """Read NIfTI-1 with read-only memory mapping and explicit physical axes."""

    DTYPES = {2: "u1", 4: "i2", 8: "i4", 16: "f4", 64: "f8",
              256: "i1", 512: "u2", 768: "u4", 1024: "i8", 1280: "u8"}

    def __init__(self, path: Path):
        self.path = path
        with path.open("rb") as stream:
            header = stream.read(352)
        if len(header) < 348:
            raise ValueError(f"Incomplete NIfTI header: {path}")
        endian = "<" if struct.unpack_from("<i", header)[0] == 348 else ">"
        unpack = lambda fmt, offset: struct.unpack_from(endian + fmt, header, offset)
        if unpack("i", 0)[0] != 348 or header[344:347] not in (b"n+1", b"ni1"):
            raise ValueError(f"Only NIfTI-1 is supported: {path}")
        dims = unpack("8h", 40)
        if dims[0] < 3 or any(n != 1 for n in dims[4:dims[0] + 1]):
            raise ValueError(f"Expected a single 3-D volume: {path}")
        self.shape = tuple(dims[1:4])
        if min(self.shape) <= 0:
            raise ValueError(f"Invalid volume shape: {path}")
        code = unpack("h", 70)[0]
        if code not in self.DTYPES:
            raise ValueError(f"Unsupported datatype {code}: {path}")
        self.slope, self.intercept = unpack("2f", 112)
        if self.slope == 0 or not math.isfinite(self.slope):
            self.slope, self.intercept = 1.0, 0.0
        if not math.isfinite(self.intercept):
            raise ValueError(f"Nonfinite intensity intercept: {path}")
        self.affine = np.eye(4)
        qcode, scode = unpack("2h", 252)
        if scode > 0:
            self.affine[:3] = [unpack("4f", 280 + i * 16) for i in range(3)]
            self.affine_source = "sform"
        elif qcode > 0:
            b, c, d = unpack("3f", 256)
            a2 = 1 - b*b - c*c - d*d
            if a2 < 1e-7:
                length = math.sqrt(b*b + c*c + d*d)
                b, c, d = b/length, c/length, d/length
                a = 0.0
            else:
                a = math.sqrt(a2)
            rotation = np.array([
                [a*a+b*b-c*c-d*d, 2*b*c-2*a*d, 2*b*d+2*a*c],
                [2*b*c+2*a*d, a*a+c*c-b*b-d*d, 2*c*d-2*a*b],
                [2*b*d-2*a*c, 2*c*d+2*a*b, a*a+d*d-c*c-b*b]])
            pixdim = unpack("8f", 76)
            self.affine[:3, :3] = rotation @ np.diag([
                pixdim[1], pixdim[2], pixdim[3] * (-1 if pixdim[0] < 0 else 1)])
            self.affine[:3, 3] = unpack("3f", 268)
            self.affine_source = "qform"
        else:
            raise ValueError(f"No physical qform or sform: {path}")
        units = header[123] & 7
        if units != 2:
            raise ValueError(f"Expected explicit millimetre units, found {units}: {path}")
        if not np.isfinite(self.affine).all() or abs(np.linalg.det(self.affine[:3, :3])) < 1e-12:
            raise ValueError(f"Invalid spatial affine: {path}")
        self.inverse_affine = np.linalg.inv(self.affine)
        pair = header[344:347] == b"ni1"
        data_path = path.with_suffix(".img") if pair else path
        offset = max(0, round(unpack("f", 108)[0]))
        if pair:
            offset = 0
        self.data = np.memmap(data_path, dtype=endian + self.DTYPES[code], mode="r",
                              offset=offset, shape=self.shape, order="F")

    def sample(self, ras: np.ndarray):
        vox = (self.inverse_affine @ ras)[:3]
        upper = (np.asarray(self.shape) - 1)[:, None]
        # Matrix round-off at exact image boundaries must not discard valid voxels.
        inside = np.all((vox >= -1e-5) & (vox <= upper + 1e-5), axis=0)
        values = map_coordinates(self.data, np.clip(vox, 0, upper), order=1, mode="constant", cval=np.nan,
                                 output=np.float64, prefilter=False)
        values[~inside] = np.nan
        return values * self.slope + self.intercept, inside

    def metadata(self):
        return {"filename": self.path.name, "shape": list(self.shape),
                "voxel_to_RAS_mm": self.affine.tolist(), "affine_source": self.affine_source}


def read_transform(path: Path):
    matrix = np.loadtxt(path)
    if matrix.shape != (4, 4) or not np.isfinite(matrix).all():
        raise ValueError(f"Expected a finite 4x4 RAS matrix: {path}")
    if not np.allclose(matrix[3], [0, 0, 0, 1], atol=1e-7):
        raise ValueError(f"Invalid homogeneous row: {path}")
    if abs(np.linalg.det(matrix[:3, :3])) < 1e-10:
        raise ValueError(f"Singular transform: {path}")
    return matrix


def scores(reference, moving, inside, mask):
    finite = np.isfinite(reference) & np.isfinite(moving)
    good = mask & inside & finite & (moving != 0)
    result = {"reference_foreground_samples": int(mask.sum()), "paired_samples": int(good.sum()),
              "source_FOV_coverage": float(np.mean(inside[mask])),
              "source_finite_coverage": float(np.mean(finite[mask])),
              "paired_nonzero_coverage": float(np.mean(good[mask])),
              "pearson_r": None, "normalized_mutual_information": None}
    if good.sum() < 100:
        return result
    x, y = reference[good], moving[good]
    if np.std(x) > 0 and np.std(y) > 0:
        result["pearson_r"] = float(np.corrcoef(x, y)[0, 1])
    # Entropy-ratio NMI, (H(reference)+H(moving))/H(joint), using 64 bins.
    hist = np.histogram2d(x, y, bins=64)[0]
    hist /= hist.sum()
    entropy = lambda p: -float(np.sum(p[p > 0] * np.log(p[p > 0])))
    joint = entropy(hist)
    if joint > 0:
        result["normalized_mutual_information"] = (entropy(hist.sum(0)) + entropy(hist.sum(1))) / joint
    return result


def normalize(values):
    finite = values[np.isfinite(values) & (values != 0)]
    if len(finite) < 10:
        return np.zeros_like(values)
    low, high = np.percentile(finite, [1, 99])
    return np.clip((np.nan_to_num(values, nan=low) - low) / max(high - low, 1e-12), 0, 1)


def montage(main, sources, ras_samples, mask, output, requested_center=None):
    import matplotlib
    matplotlib.use("Agg")
    from matplotlib import pyplot as plt

    if requested_center is None:
        # Whole-head scans can contain a long neck. Center the display in the
        # superior half of foreground samples; this is not a brain segmentation.
        superior = mask & (ras_samples[2] >= np.median(ras_samples[2, mask]))
        center = np.median(ras_samples[:3, superior], axis=1)
    else:
        center = np.asarray(requested_center, dtype=float)
    # Physical RAS planes: axial x/y, coronal x/z, sagittal y/z.
    specs = [(0, 1, 2, center[2] + dz, f"Axial z={center[2]+dz:.1f}") for dz in (-25, 0, 25)]
    specs += [(0, 2, 1, center[1], f"Coronal y={center[1]:.1f}"),
              (1, 2, 0, center[0], f"Sagittal x={center[0]:.1f}")]
    rows = [("Reference T1", main, np.eye(4))] + sources
    fig, axes = plt.subplots(len(rows), len(specs), figsize=(15, 3 * len(rows)), squeeze=False)
    side = 241
    extent = [-120, 120, -120, 120]
    grid = np.linspace(-120, 120, side)
    u, v = np.meshgrid(grid, grid)
    checker = ((np.indices((side, side))[0] // 24 + np.indices((side, side))[1] // 24) % 2).astype(bool)
    for col, (u_axis, v_axis, fixed_axis, level, title) in enumerate(specs):
        points = np.ones((4, side * side))
        points[u_axis] = center[u_axis] + u.ravel()
        points[v_axis] = center[v_axis] + v.ravel()
        points[fixed_axis] = level
        ref_values, _ = main.sample(points)
        ref_norm = normalize(ref_values).reshape(side, side)
        for row, (name, source, matrix) in enumerate(rows):
            ax = axes[row, col]
            if row == 0:
                display = ref_norm
            else:
                values, _ = source.sample(matrix @ points)
                display = np.where(checker, normalize(values).reshape(side, side), ref_norm)
            ax.imshow(display, cmap="gray", vmin=0, vmax=1, origin="lower", extent=extent)
            ax.axhline(0, color="#45c5e8", linewidth=.45, alpha=.7)
            ax.axvline(0, color="#45c5e8", linewidth=.45, alpha=.7)
            ax.set_xticks([])
            ax.set_yticks([])
            if row == 0:
                ax.set_title(title + " mm RAS", fontsize=9)
            if col == 0:
                ax.set_ylabel(name, fontsize=10)
            # Neurological display: the right side of axial/coronal panels is R.
            ax.text(.98, .03, "R" if u_axis == 0 else "A", color="white", fontsize=8,
                    transform=ax.transAxes, ha="right")
    if sources:
        fig.suptitle("Registration review: alternating T1 and transformed anatomy tiles\n"
                     "Crosshairs share a physical RAS center", fontsize=12)
    else:
        fig.suptitle("T1 reference review: no anatomical overlays were processed", fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, .88 if len(rows) == 1 else .96])
    target = output / "anatomy_checkerboard.png"
    fig.savefig(target, dpi=150, facecolor="white")
    plt.close(fig)
    return target, center.tolist()


def _recorded_path(value, base, label):
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"Run metadata has no valid {label} path")
    path = Path(value)
    return (path if path.is_absolute() else base / path).resolve()


def _required_file(path, label):
    if not path.is_file():
        raise FileNotFoundError(f"{label} is missing: {path}")
    if path.stat().st_size == 0:
        raise ValueError(f"{label} is empty: {path}")
    return path


def resolve_run_inputs(output, images=None, main_image=None):
    """Resolve a completed run without guessing filenames or optional groups.

    --images and --main are compatibility checks only. In particular, the
    original --main filename must never replace a prepared display copy recorded
    in the manifest. Relative manifest paths, when present in older runs, are
    interpreted relative to the run directory.
    """
    output = Path(output).resolve()
    if not output.is_dir():
        raise ValueError("--output must be an existing pipeline run directory")
    manifest_path = _required_file(output / "manifest.json", "Run manifest")
    config_path = _required_file(output / "config_used.json", "Recorded run configuration")
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    config = json.loads(config_path.read_text(encoding="utf-8"))
    success = {"built; anatomical review required", "built", "complete", "completed", "success", "successful"}
    if manifest.get("status") not in success or manifest.get("error"):
        raise ValueError(f"QC requires a successfully built run; manifest status is {manifest.get('status')!r}")
    portable = manifest.get("path_format") == "relative_to_run_directory"

    def run_path(value, label):
        path = _recorded_path(value, output, label)
        if portable and not path.is_relative_to(output):
            raise ValueError(f"Portable run {label} must be inside the run directory: {path}")
        return path

    images_dir = run_path(manifest.get("images_dir", config.get("images")), "input images directory")
    if images is not None and Path(images).resolve() != images_dir:
        raise ValueError(f"--images does not match the run's recorded input directory: {images_dir}")
    reference = _required_file(run_path(manifest.get("main"), "main image"), "Recorded main image")
    if main_image is not None:
        supplied = Path(main_image)
        supplied = (supplied if supplied.is_absolute() else images_dir / supplied).resolve()
        allowed = {reference}
        original_main = manifest.get("main_file", config.get("main"))
        if original_main:
            allowed.add(_recorded_path(original_main, images_dir, "original main image"))
        if supplied not in allowed:
            raise ValueError("--main does not match this run's original or displayed main image")

    modalities = manifest.get("modalities", {})
    skipped = {}
    for name, details in modalities.items():
        if details.get("status") == "skipped":
            reason = details.get("reason")
            if not reason:
                parts = []
                for field, label in [("missing_files", "Missing files"), ("empty_files", "Empty files")]:
                    if details.get(field):
                        parts.append(label + ": " + ", ".join(map(str, details[field])))
                reason = "; ".join(parts) or "Skipped by the workspace builder"
            skipped[name] = {**details, "reason": reason}

    recorded_groups = manifest.get("groups", {})
    configured_groups = config.get("groups", {})
    for name in configured_groups:
        if name not in recorded_groups and name not in skipped:
            raise ValueError(f"Configured active group {name!r} has no completed manifest entry")
    groups = {}
    for name, details in recorded_groups.items():
        if name in skipped:
            continue
        if name not in configured_groups:
            raise ValueError(f"Processed group {name!r} has no recorded configuration")
        spec = configured_groups[name]
        qc_anatomy = spec.get("qc_anatomy", False)
        if not isinstance(qc_anatomy, bool):
            raise ValueError(f"qc_anatomy must be true or false for group {name!r}")
        source = _required_file(run_path(details.get("source"), f"{name} source image"),
                                f"Active {name} source image")
        fallback = output / "transforms" / f"{name}_RAS.txt"
        transform = _required_file(run_path(details.get("RAS_transform", str(fallback)),
                                                f"{name} RAS transform"), f"Active {name} RAS transform")
        groups[name] = {"source": source, "transform": transform, "qc_anatomy": qc_anatomy}
    return {"output": output, "manifest_path": manifest_path, "config_path": config_path,
            "images_dir": images_dir, "reference": reference, "groups": groups, "skipped": skipped,
            "portable": portable,
            "workspace_kind": manifest.get("workspace_kind", "multimodal" if groups else "T1-only")}


def run_qc(output, images=None, main_image=None, stride=4, center_ras=None, make_montage=True):
    """Check the images and transforms actually recorded in a successful run."""
    if stride < 1:
        raise ValueError("Sampling stride must be positive")
    if center_ras is not None and (len(center_ras) != 3 or not np.isfinite(center_ras).all()):
        raise ValueError("Montage center must contain three finite RAS coordinates")
    inputs = resolve_run_inputs(output, images=images, main_image=main_image)

    def report_path(path):
        # Operational paths in a portable report are relative to the run root,
        # not the qc subdirectory and not the caller's working directory.
        path = Path(path).resolve()
        return path.relative_to(inputs["output"]).as_posix() if inputs["portable"] else str(path)

    qc_dir = inputs["output"] / "qc"
    reference = NiftiImage(inputs["reference"])
    voxel = np.array(np.meshgrid(*[np.arange(0, n, stride) for n in reference.shape], indexing="ij")).reshape(3, -1)
    points = reference.affine @ np.vstack([voxel, np.ones(voxel.shape[1])])
    values = reference.data[tuple(voxel)] * reference.slope + reference.intercept
    positive = values[np.isfinite(values) & (values > 0)]
    if len(positive) < 100:
        raise ValueError("Reference contains insufficient positive finite samples")
    threshold = float(np.percentile(positive, 70))
    mask = np.isfinite(values) & (values > threshold)
    if not mask.any():
        raise ValueError("Reference has no varying positive foreground for anatomical QC")
    report = {"reference": {**reference.metadata(), "path": report_path(inputs["reference"])},
              "manifest": report_path(inputs["manifest_path"]), "config_used": report_path(inputs["config_path"]),
              "images_dir": report_path(inputs["images_dir"]), "workspace_kind": inputs["workspace_kind"],
              "transform_convention": "fixed RAS mm to moving RAS mm",
              "sampling_stride_voxels": stride, "reference_foreground_threshold": threshold,
              "interpretation": "Descriptive metrics only. Human review of registration is required; correlation does not validate clinical use.",
              "checks": {}, "transform_geometry": {}, "skipped_modalities": inputs["skipped"], "warnings": []}
    if inputs["portable"]:
        report["path_format"] = "relative_to_run_directory"
    sources = []
    for name, group in inputs["groups"].items():
        source = NiftiImage(group["source"])
        matrix = read_transform(group["transform"])
        report["transform_geometry"][name] = {
            "transform_path": report_path(group["transform"]), "source": source.metadata(),
            "matrix_RAS": matrix.tolist(), "linear_determinant": float(np.linalg.det(matrix[:3, :3])),
            "linear_singular_values": np.linalg.svd(matrix[:3, :3], compute_uv=False).tolist()}
        source_metadata = {**source.metadata(), "path": report_path(group["source"])}
        if not group["qc_anatomy"]:
            report["checks"][name] = {"status": "geometry_only", "source": source_metadata,
                                     "reason": "Anatomical intensity comparison is not enabled for this group"}
            continue
        identity_values, identity_inside = source.sample(points)
        mapped_values, mapped_inside = source.sample(matrix @ points)
        report["checks"][name] = {"status": "measured", "source": source_metadata,
                                 "identity": scores(values, identity_values, identity_inside, mask),
                                 "registered": scores(values, mapped_values, mapped_inside, mask)}
        sources.append((name + " anatomy", source, matrix))
    report["anatomical_comparisons_performed"] = len(sources)
    qc_dir.mkdir(exist_ok=True)
    if make_montage:
        path, center = montage(reference, sources, points, mask, qc_dir, center_ras)
        report["montage"] = report_path(path)
        report["montage_center_RAS_mm"] = center
        report["montage_center_method"] = "user specified" if center_ras is not None else "median of superior half of foreground samples (display heuristic only)"
    destination = qc_dir / "registration_metrics.json"
    destination.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--images", type=Path, help="Optional compatibility check against the recorded input folder")
    parser.add_argument("--output", required=True, type=Path, help="Successful run directory containing manifest.json and config_used.json")
    parser.add_argument("--main", help="Optional compatibility check against the recorded original or displayed main image")
    parser.add_argument("--stride", type=int, default=4, help="Voxel sampling stride for numeric checks")
    parser.add_argument("--center-ras", type=float, nargs=3, metavar=("X", "Y", "Z"),
                        help="Optional montage crosshair position in reference RAS mm")
    parser.add_argument("--no-montage", action="store_true")
    args = parser.parse_args()
    try:
        report = run_qc(args.output, images=args.images, main_image=args.main, stride=args.stride,
                        center_ras=args.center_ras, make_montage=not args.no_montage)
    except (ValueError, FileNotFoundError) as exc:
        parser.error(str(exc))
    destination = args.output.resolve() / "qc" / "registration_metrics.json"
    print(destination)
    print(f"Workspace: {report['workspace_kind']}; anatomical comparisons: {report['anatomical_comparisons_performed']}")
    for name, check in report["checks"].items():
        if check["status"] == "measured":
            s = check["registered"]
            print(f"{name}: r={s['pearson_r']}, NMI={s['normalized_mutual_information']}, "
                  f"FOV={s['source_FOV_coverage']:.3f}")
        else:
            print(f"{name}: {check['status']}")
    for name, skipped in report["skipped_modalities"].items():
        print(f"{name}: skipped ({skipped['reason']})")


if __name__ == "__main__":
    main()
