"""Deterministic translation initialization for corresponding T1 anatomy.

This initializer uses image intensities and headers only. It never reads an
existing registration. Its output maps fixed RAS millimetres to moving RAS
millimetres (the pull convention used by Greedy). Use with corresponding T1
anatomical images, not PET, FLAIR, sparse statistical maps, or label images.
The numeric checks reject gross failures; they do not replace visual review.
"""

from __future__ import annotations

import itertools
from pathlib import Path

import numpy as np
from scipy.optimize import minimize

from qc_registration import NiftiImage


STRIDE = 8
FOREGROUND_PERCENTILE = 70
MIN_OVERLAP = 0.50
MIN_PAIRED_SAMPLES = 100
COARSE_HALF_WIDTHS = np.array([20.0, 40.0, 40.0])
REFINE_HALF_WIDTHS = np.array([30.0, 50.0, 50.0])


def _grid(image):
    # The fixed phase makes repeated runs deterministic and avoids sampling
    # exactly on every image boundary. All coordinate calculations remain RAS.
    indices = np.array(np.meshgrid(
        *[np.arange(min(2, n - 1), n, STRIDE) for n in image.shape],
        indexing="ij")).reshape(3, -1)
    points = image.affine @ np.vstack([indices, np.ones(indices.shape[1])])
    values = np.asarray(image.data[tuple(indices)], dtype=float) * image.slope + image.intercept
    return points, values


def _intensity_center(points, values):
    weights = np.where(np.isfinite(values) & (values > 0), values, 0.0)
    if np.count_nonzero(weights) < MIN_PAIRED_SAMPLES or weights.sum() <= 0:
        raise ValueError("T1 image has insufficient positive finite samples for center of mass")
    return np.sum(points[:3] * weights[None, :], axis=1) / weights.sum()


class _AlignmentContext:
    def __init__(self, fixed_path, moving_path):
        self.fixed = NiftiImage(Path(fixed_path))
        self.moving = NiftiImage(Path(moving_path))
        fixed_points, fixed_values = _grid(self.fixed)
        positive = fixed_values[np.isfinite(fixed_values) & (fixed_values > 0)]
        if len(positive) < MIN_PAIRED_SAMPLES:
            raise ValueError("Fixed T1 image has insufficient positive finite samples")
        self.threshold = float(np.percentile(positive, FOREGROUND_PERCENTILE))
        keep = np.isfinite(fixed_values) & (fixed_values > self.threshold)
        if keep.sum() < MIN_PAIRED_SAMPLES:
            raise ValueError("Fixed T1 foreground is too small for translation initialization")
        self.points, self.values = fixed_points[:, keep], fixed_values[keep]
        self.fixed_center = _intensity_center(fixed_points, fixed_values)
        moving_points, moving_values = _grid(self.moving)
        self.moving_center = _intensity_center(moving_points, moving_values)

    def score(self, matrix):
        matrix = np.asarray(matrix, dtype=float)
        if matrix.shape != (4, 4) or not np.isfinite(matrix).all():
            raise ValueError("Alignment must be a finite 4x4 RAS matrix")
        if not np.allclose(matrix[3], [0, 0, 0, 1], atol=1e-7):
            raise ValueError("Alignment has an invalid homogeneous row")
        moving, inside = self.moving.sample(matrix @ self.points)
        good = inside & np.isfinite(moving) & (moving > 0)
        overlap = float(np.mean(good))
        pearson = None
        if good.sum() >= MIN_PAIRED_SAMPLES:
            x, y = self.values[good], moving[good]
            if np.std(x) > 1e-10 and np.std(y) > 1e-10:
                pearson = float(np.corrcoef(x, y)[0, 1])
        return {
            "pearson_r": pearson,
            "overlap_fraction": overlap,
            "source_FOV_fraction": float(np.mean(inside)),
            "foreground_sample_count": int(len(self.values)),
            "paired_sample_count": int(good.sum()),
            "minimum_overlap_fraction": MIN_OVERLAP,
            "overlap_sufficient": bool(overlap >= MIN_OVERLAP and good.sum() >= MIN_PAIRED_SAMPLES),
            "foreground_threshold": self.threshold,
            "foreground_percentile_of_positive_samples": FOREGROUND_PERCENTILE,
            "sampling_stride_voxels": STRIDE,
        }

    def objective(self, translation):
        score = self.score(_translation_matrix(translation))
        if not score["overlap_sufficient"] or score["pearson_r"] is None:
            return 2.0 + max(0.0, MIN_OVERLAP - score["overlap_fraction"])
        return -score["pearson_r"]


def _translation_matrix(translation):
    matrix = np.eye(4)
    matrix[:3, 3] = translation
    return matrix


def score_t1_alignment(fixed_path, moving_path, matrix):
    """Return deterministic T1 Pearson and overlap measurements for a RAS map.

    overlap_fraction is the fraction of sampled fixed foreground with finite,
    positive moving intensity. It is not a brain-coverage or lesion measure.
    A small overlapping fragment must never be accepted based on Pearson alone.
    """
    return _AlignmentContext(fixed_path, moving_path).score(matrix)


def estimate_t1_translation(fixed_path, moving_path):
    """Estimate a translation without an existing transform; return matrix/report.

    Start at intensity-weighted centers, examine a deterministic 405-point
    translation grid, and refine the best candidate with bounded Powell search.
    This intentionally estimates translation only; a later rigid registration
    may refine rotation. The caller must inspect the returned quality scores.
    """
    context = _AlignmentContext(fixed_path, moving_path)
    start = context.moving_center - context.fixed_center
    grid_axes = [np.arange(-width, width + 1, 10) for width in COARSE_HALF_WIDTHS]
    best_translation = None
    best_objective = float("inf")
    grid_count = 0
    for offset in itertools.product(*grid_axes):
        translation = start + np.asarray(offset)
        value = context.objective(translation)
        grid_count += 1
        if value < best_objective:
            best_objective, best_translation = value, translation.copy()
    if best_objective > 1.0:
        raise ValueError("T1 translation search found no candidate with sufficient finite overlap")
    coarse_translation = best_translation.copy()
    coarse_score = context.score(_translation_matrix(coarse_translation))
    bounds = list(zip(start - REFINE_HALF_WIDTHS, start + REFINE_HALF_WIDTHS))
    refined = minimize(context.objective, coarse_translation, method="Powell", bounds=bounds,
                       options={"xtol": 0.001, "ftol": 1e-7, "maxiter": 100, "maxfev": 2000})
    # A numerical optimizer is allowed to stop unsuccessfully. Preserve the
    # measured best grid candidate if the refinement does not improve it.
    if np.isfinite(refined.fun) and refined.fun < best_objective:
        best_translation = np.asarray(refined.x).copy()
        selected = "bounded Powell refinement"
    else:
        selected = "coarse grid candidate"
    matrix = _translation_matrix(best_translation)
    selected_score = context.score(matrix)
    report = {
        "method": "intensity center of mass, deterministic translation grid, bounded Powell",
        "transform_convention": "fixed RAS mm to moving RAS mm",
        "uses_existing_transform": False,
        "fixed_image": Path(fixed_path).name,
        "moving_image": Path(moving_path).name,
        "fixed_intensity_center_RAS_mm": context.fixed_center.tolist(),
        "moving_intensity_center_RAS_mm": context.moving_center.tolist(),
        "center_of_mass_translation_RAS_mm": start.tolist(),
        "center_of_mass_score": context.score(_translation_matrix(start)),
        "coarse_grid_half_widths_mm": COARSE_HALF_WIDTHS.tolist(),
        "coarse_grid_step_mm": 10,
        "coarse_grid_evaluations": grid_count,
        "coarse_translation_RAS_mm": coarse_translation.tolist(),
        "coarse_score": coarse_score,
        "refinement_bounds_RAS_mm": [[float(lo), float(hi)] for lo, hi in bounds],
        "refinement_optimizer_success": bool(refined.success),
        "refinement_optimizer_message": str(refined.message),
        "refinement_evaluations": int(refined.nfev),
        "selected_candidate": selected,
        "translation_RAS_mm": best_translation.tolist(),
        "score": selected_score,
        "pearson_r": selected_score["pearson_r"],
        "overlap_fraction": selected_score["overlap_fraction"],
        "overlap_sufficient": selected_score["overlap_sufficient"],
        "interpretation": "Initialization and gross-failure measurements only; visual anatomical review is required.",
    }
    return matrix, report
