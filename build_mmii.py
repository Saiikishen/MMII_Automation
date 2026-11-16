"""Generate fresh registration transforms and an ITK-SNAP MMII workspace.

Requires Python 3.10+, numpy, scipy and installed ITK-SNAP command-line tools.
Input files are read only. Every invocation requires a new output directory.
This assembles existing modality outputs; it does not rerun their analyses.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import struct
import subprocess
import sys
import xml.etree.ElementTree as ET

import numpy as np

from modality_inputs import plan_modalities, print_input_report, report_lines
from portable_run import make_portable


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(4 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def read_volume(path):
    """Read uncompressed scalar NIfTI-1, including .hdr/.img NIfTI pairs.

    Refuse legacy Analyze (ambiguous orientation) instead of guessing it.
    Affines map zero-based voxels to RAS millimetres.
    """
    path = Path(path)
    if path.suffix.lower() == '.img':
        path = path.with_suffix('.hdr')
    with path.open("rb") as f:
        h = f.read(348)
    if len(h) != 348:
        raise ValueError(f"Truncated header: {path}")
    endian = "<" if struct.unpack_from("<i", h)[0] == 348 else ">"
    if struct.unpack_from(endian + "i", h)[0] != 348 or h[344:347] not in (b"n+1", b"ni1"):
        raise ValueError(f"Expected an uncompressed NIfTI-1 file, not legacy Analyze: {path}")
    def unpack(fmt, offset):
        return struct.unpack_from(endian + fmt, h, offset)
    dims = unpack("8h", 40)
    if dims[0] < 3 or any(v != 1 for v in dims[4:dims[0] + 1]):
        raise ValueError(f"Expected a single 3-D scalar volume: {path}")
    shape = tuple(dims[1:4])
    if any(v <= 0 for v in shape):
        raise ValueError(f"Invalid dimensions: {path}")
    types = {2:"u1", 4:"i2", 8:"i4", 16:"f4", 64:"f8", 256:"i1", 512:"u2", 768:"u4", 1024:"i8", 1280:"u8"}
    code = unpack("h", 70)[0]
    if code not in types:
        raise ValueError(f"Unsupported scalar datatype {code}: {path}")
    dtype = np.dtype(endian + types[code])
    affine = np.eye(4)
    qcode, scode = unpack("2h", 252)
    if scode:
        affine[:3] = np.array(unpack("12f", 280)).reshape(3, 4)
    elif qcode:
        b, c, d = unpack("3f", 256)
        norm = b*b+c*c+d*d
        if norm > 1:
            b, c, d = np.array([b,c,d]) / np.sqrt(norm)
        a = np.sqrt(max(0, 1-b*b-c*c-d*d))
        rot = np.array([[a*a+b*b-c*c-d*d, 2*(b*c-a*d), 2*(b*d+a*c)],
                        [2*(b*c+a*d), a*a+c*c-b*b-d*d, 2*(c*d-a*b)],
                        [2*(b*d-a*c), 2*(c*d+a*b), a*a+d*d-c*c-b*b]])
        pixdim = unpack("8f", 76)
        affine[:3,:3] = rot @ np.diag([pixdim[1],pixdim[2],pixdim[3]*(-1 if pixdim[0]<0 else 1)])
        affine[:3,3] = unpack("3f", 268)
    else:
        raise ValueError(f"No coded physical-space affine: {path}")
    if h[123] & 7 != 2:
        raise ValueError(f"Only explicit millimetre units are supported: {path}")
    if not np.isfinite(affine).all() or abs(np.linalg.det(affine[:3,:3])) < 1e-9:
        raise ValueError(f"Invalid affine: {path}")
    data_path = path.with_suffix(".img") if h[344:347] == b"ni1" else path
    offset = int(unpack("f", 108)[0])
    if data_path.stat().st_size < offset + int(np.prod(shape)) * dtype.itemsize:
        raise ValueError(f"Truncated image data: {data_path}")
    slope, intercept = unpack("2f",112)
    if not np.isfinite(slope) or slope == 0:
        slope, intercept = 1., 0.
    if not np.isfinite(intercept):
        raise ValueError(f"Invalid scaling intercept: {path}")
    data = np.memmap(data_path, dtype=dtype, mode="r", offset=offset, shape=shape, order="F")
    finite = np.isfinite(data)
    if not finite.any():
        raise ValueError(f"No finite voxels: {path}")
    extremes = np.array([float(data[finite].min()), float(data[finite].max())]) * slope + intercept
    return {"path": path, "data_path": data_path, "shape": shape, "affine": affine,
            "min": float(extremes.min()), "max": float(extremes.max()),
            "nonfinite": int(data.size - np.count_nonzero(finite)), "data": data,
            "slope": slope, "intercept": intercept}


def same_grid(a, b):
    return a["shape"] == b["shape"] and np.allclose(a["affine"], b["affine"], atol=1e-4, rtol=0)


def center_transform(fixed, moving):
    transform = np.eye(4)
    fc = fixed["affine"] @ np.r_[(np.array(fixed["shape"])-1)/2, 1]
    mc = moving["affine"] @ np.r_[(np.array(moving["shape"])-1)/2, 1]
    transform[:3,3] = mc[:3] - fc[:3]
    return transform


def check_rigid(matrix):
    if matrix.shape != (4,4) or not np.isfinite(matrix).all():
        raise ValueError("Invalid 4x4 transform")
    if not np.allclose(matrix[3], [0,0,0,1], atol=1e-7):
        raise ValueError("Invalid homogeneous row")
    r = matrix[:3,:3]
    if not np.allclose(r.T @ r, np.eye(3), atol=1e-3) or not np.isclose(np.linalg.det(r),1,atol=1e-3):
        raise ValueError("Expected a proper rigid transform with no reflection")


def foreground_fov_fraction(fixed, moving, matrix):
    """A gross-failure guard, not a registration-quality acceptance test."""
    ijk = np.indices(tuple((np.array(fixed["shape"])+7)//8)).reshape(3,-1) * 8
    values = np.asarray(fixed["data"][tuple(ijk)]) * fixed["slope"] + fixed["intercept"]
    nonzero = values[np.isfinite(values) & (values > 0)]
    if len(nonzero) < 100:
        raise ValueError("Too little positive reference anatomy for a coverage check")
    foreground = np.isfinite(values) & (values >= np.percentile(nonzero,25))
    if not np.any(foreground):
        raise ValueError("No reference foreground available for a coverage check")
    points = np.vstack((ijk[:,foreground], np.ones(np.count_nonzero(foreground))))
    moving_voxels = (np.linalg.inv(moving["affine"]) @ matrix @ fixed["affine"] @ points)[:3]
    inside = np.all((moving_voxels >= -1e-4) & (moving_voxels <= (np.array(moving["shape"])-1)[:,None]+1e-4),axis=0)
    return float(np.mean(inside))


def display_settings(layer, volume):
    """Keep valid empty result maps loadable without colouring the whole T1."""
    low, high = volume['min'], volume['max']
    lo = layer.get('minimum', low)
    hi = layer.get('maximum', high)
    if not np.isfinite([lo, hi, low, high]).all():
        raise ValueError(f"Nonfinite contrast limits: {layer['file']}")
    if 'maximum' in layer and hi <= lo:
        raise ValueError(f"Configured maximum must exceed minimum: {layer['file']}")
    if high == low:
        return {'range': [lo, max(hi, lo+1.)], 'normalized': [0., 1.], 'alpha': 0.,
                'warning': f"{layer['file']}: constant value {low:g}; layer retained with opacity zero"}
    warning = None
    if hi <= lo:
        hi = lo + max(high-low, 1.)
        warning = f"{layer['file']}: no values exceed display minimum {lo:g}; overlay remains transparent"
    return {'range': [lo, hi], 'normalized': [(lo-low)/(high-low), (hi-low)/(high-low)],
            'alpha': layer['alpha'], 'warning': warning}


def default_output(source_dir):
    name = source_dir.name
    patient = name[:-7] if name.lower().endswith('_images') else name
    return source_dir.parent / ((patient or name) + '_MMII_' + datetime.now().strftime('%Y%m%d_%H%M%S_%f'))


def build(args):
    config_path = Path(args.config).resolve()
    config = json.loads(config_path.read_text(encoding="utf-8-sig"))
    input_location = args.images or config.get('images')
    if not input_location:
        raise ValueError('Supply --images with the patientName_Images folder containing t1.nii')
    source_dir = Path(input_location).expanduser()
    if not source_dir.is_absolute() and not args.images:
        source_dir = config_path.parent / source_dir
    source_dir = source_dir.resolve()
    plan = plan_modalities(config, source_dir)
    print_input_report(plan)
    if plan['main']['status'] != 'ready':
        raise ValueError(f"Required main image {config['main']} is {plan['main']['status']}; cannot continue")
    if getattr(args, 'check_inputs', False):
        print('File availability check only; image contents and registration were not checked.', flush=True)
        return plan
    output = Path(args.output).resolve() if args.output else default_output(source_dir)
    if output == source_dir or source_dir in output.parents:
        raise ValueError("Choose an output directory outside the original image folder")
    if output.exists():
        raise FileExistsError(f"Use a new output directory; will not overwrite: {output}")
    bindir = Path(args.itksnap_bin or config["itksnap_bin"])
    exe = {n: bindir / (n + (".exe" if sys.platform == "win32" else "")) for n in ("greedy", "c3d", "c3d_affine_tool", "itksnap-wt")}
    for p in exe.values():
        if not p.is_file():
            raise FileNotFoundError(p)
    # Remove every reference to an incomplete modality BEFORE reading images,
    # registering, writing layers, or asking QC to inspect it.
    config['groups'] = {n:g for n,g in config['groups'].items() if n in plan['active_groups']}
    config['layers'] = [layer for layer in config['layers'] if layer['group'] in plan['active_groups']]
    names = {config["main"], *(g["source"] for g in config["groups"].values()), *(l["file"] for l in config["layers"])}
    output.mkdir(parents=True)
    for sub in ("transforms", "logs", "prepared"):
        (output/sub).mkdir()
    config["images"] = str(source_dir)
    (output/"config_used.json").write_text(json.dumps(config, indent=2), encoding="utf-8")
    manifest = {"created_utc": datetime.now(timezone.utc).isoformat(), "main": str(source_dir/config['main']),
                "main_file": config['main'], "images_dir": str(source_dir),
                "workspace_kind": 'multimodal' if plan['active_groups'] else 'T1-only',
                "status": "incomplete", "clinical_qc": "not performed", "inputs_sha256": {},
                "builder_sha256": sha256(Path(__file__)),
                "modalities": plan['modalities'],
                "groups": {}, "layers": [], "commands": [], "warnings": []}
    def save_manifest():
        (output/"manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
        (output/"modality_report.json").write_text(json.dumps(plan, indent=2), encoding='utf-8')
        (output/"modality_report.txt").write_text('\n'.join(report_lines(plan)) + '\n', encoding='utf-8')
    def run(name, argv, label):
        cmd = [str(exe[name]), *map(str, argv)]
        manifest["commands"].append(cmd)
        save_manifest()
        with (output/"logs"/(label+".log")).open("w", encoding="utf-8") as log:
            result = subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT)
        if result.returncode:
            raise RuntimeError(f"{name} failed ({result.returncode}); see logs/{label}.log")
    processing_group = None
    save_manifest()
    try:
        print("Checking complete modalities' image contents and geometry...", flush=True)
        volumes = {}
        for name in sorted(names):
            try:
                volumes[name] = read_volume(source_dir/name)
            except Exception as exc:
                if name == config['main']:
                    plan['main'].update(status='failed', error=str(exc))
                else:
                    for group, record in plan['modalities'].items():
                        if group in plan['active_groups'] and name in record['required_files']:
                            record.update(status='failed', error=str(exc))
                raise
        main = volumes[config['main']]
        if main['max'] <= main['min']:
            plan['main'].update(status='failed', error='Main image has no intensity variation')
            raise ValueError('Main image has no intensity variation; cannot register anatomy')
        for layer in config['layers']:
            processing_group = layer['group']
            source = volumes[config['groups'][processing_group]['source']]
            if not same_grid(source, volumes[layer['file']]):
                raise ValueError(f"{processing_group}: transform propagation requires a verified shared grid: {layer['file']}")
        processing_group = None
        inputs = {str(p): sha256(p) for v in volumes.values() for p in {v['path'], v['data_path']}}
        manifest['inputs_sha256'] = inputs
        # Replace nonfinite values only in a new display copy, never in inputs.
        display_paths = {}
        for name, v in volumes.items():
            display_paths[name] = v["path"]
            if v["nonfinite"]:
                cleaned = output/"prepared"/(Path(name).stem+"_finite.nii")
                run("c3d", [v["path"], "-replace", "nan", "0", "inf", "0", "-inf", "0", "-o", cleaned], "clean_"+Path(name).stem)
                clean = read_volume(cleaned)
                if clean["nonfinite"] or not same_grid(v, clean):
                    raise ValueError(f"Failed finite display preparation: {name}")
                display_paths[name] = cleaned
                manifest["warnings"].append(f"{name}: {v['nonfinite']} nonfinite voxels replaced with zero in display copy only")
                volumes[name] = clean
        main = volumes[config["main"]]
        manifest['main'] = str(main['path'])
        transforms = {}
        for group, spec in config["groups"].items():
            processing_group = group
            plan['modalities'][group]['status'] = 'processing'
            save_manifest()
            method = spec["method"]
            source = volumes[spec["source"]]
            matrix_path = output/"transforms"/(group+"_RAS.txt")
            initialization_report = None
            final_anatomy_score = None
            print(f"Generating {group} transform ({method})...", flush=True)
            if method == "centers":
                matrix = center_transform(main, source)
                np.savetxt(matrix_path, matrix, fmt="%.17g")
                manifest["warnings"].append(f"{group}: geometric center alignment only; manual anatomical review still required")
            elif (method == "rigid_or_identical" and same_grid(main,source)
                  and sha256(main["path"]) == sha256(source["path"])
                  and sha256(main["data_path"]) == sha256(source["data_path"])):
                method = "identical_images_identity"
                np.savetxt(matrix_path, np.eye(4), fmt="%.17g")
            elif method in ("rigid", "rigid_or_identical"):
                reg = config["registration"]
                init_args = ["-ia-image-centers"]
                if spec.get("initialize_t1_translation"):
                    from anatomical_init import estimate_t1_translation, score_t1_alignment
                    initial, initialization_report = estimate_t1_translation(main["path"],display_paths[spec["source"]])
                    if (initialization_report["pearson_r"] is None or not initialization_report["overlap_sufficient"]
                        or initialization_report["pearson_r"] < 0.5):
                        raise ValueError(f"{group}: anatomical initialization correlation is too low; review source anatomy")
                    init_path = output/"transforms"/(group+"_initial_RAS.txt")
                    np.savetxt(init_path,initial,fmt="%.17g")
                    init_args = ["-ia",init_path]
                extra = ["-jitter",spec["jitter"]] if "jitter" in spec else []
                run("greedy", ["-d",3,"-a","-dof",6,"-m",reg["metric"],"-i",main["path"],display_paths[spec["source"]],
                               *init_args,"-n",spec.get("iterations",reg["iterations"]),"-threads",reg["threads"],*extra,"-o",matrix_path], group+"_registration")
                if initialization_report:
                    final_anatomy_score = score_t1_alignment(main["path"],display_paths[spec["source"]],np.loadtxt(matrix_path))
                    if (final_anatomy_score["pearson_r"] is None or not final_anatomy_score["overlap_sufficient"]
                        or final_anatomy_score["pearson_r"] < max(0.5,initialization_report["pearson_r"]-0.03)):
                        raise ValueError(f"{group}: rigid refinement worsened anatomical alignment; inspect initialization and fit")
            elif method == "external_ras":
                supplied = Path(spec["transform"])
                np.savetxt(matrix_path, np.loadtxt(supplied), fmt="%.17g")
            else:
                raise ValueError(f"Unsupported transform method {method}")
            matrix = np.loadtxt(matrix_path)
            check_rigid(matrix)
            coverage = foreground_fov_fraction(main, source, matrix)
            if not np.isfinite(coverage):
                raise ValueError(f"Invalid field-of-view coverage for {group}")
            if coverage < spec.get("minimum_fov_fraction",0.3):
                raise ValueError(f"Gross alignment failure for {group}: only {coverage:.1%} reference foreground inside moving FOV")
            reg_log = output/"logs"/(group+"_registration.log")
            if reg_log.exists() and "LINE SEARCH FAILED" in reg_log.read_text(errors="replace"):
                manifest["warnings"].append(f"{group}: optimizer reported a line-search warning; inspect registration log and anatomical QC")
            itk_path = output/"transforms"/("trans_"+group+".txt")
            run("c3d_affine_tool", [matrix_path,"-oitk",itk_path], group+"_export_itk")
            # Round-trip the GUI-loadable transform to catch inversion/convention mistakes.
            roundtrip = output/"transforms"/(group+"_roundtrip_RAS.txt")
            run("c3d_affine_tool", ["-itk",itk_path,"-o",roundtrip], group+"_roundtrip")
            # c3d_affine_tool prints about six significant digits in RAS text.
            if not np.allclose(matrix, np.loadtxt(roundtrip), atol=1e-4, rtol=0):
                raise ValueError(f"Transform round-trip failed for {group}")
            transforms[group] = matrix_path
            manifest["groups"][group] = {"method": method, "source": str(source["path"]),
                "RAS_transform": str(matrix_path),
                "fixed_to_moving_RAS": matrix.tolist(), "itk_transform": str(itk_path),
                "reference_foreground_in_source_fov": coverage,
                "initialization": initialization_report, "final_anatomy_score": final_anatomy_score}
            plan['modalities'][group]['status'] = 'processed'
            save_manifest()
        processing_group = None
        print("Building workspace and display settings...", flush=True)
        workspace = output/"MMII.itksnap"
        command = ["-layers-set-main",main["path"],"-props-set-nickname","T1 reference"]
        for layer in config["layers"]:
            processing_group = layer['group']
            vol = volumes[layer["file"]]
            settings = display_settings(layer, vol)
            lo, hi = settings['range']
            norm_lo, norm_hi = settings['normalized']
            if settings['warning']:
                manifest['warnings'].append(settings['warning'])
            command += ["-layers-add-anat", display_paths[layer["file"]],"-props-set-nickname",layer["name"],
                        "-props-set-colormap",layer["color"],"-props-set-contrast","LINEAR",format(norm_lo,".17g"),format(norm_hi,".17g"),
                        "-props-registry-set","LayerMetaData.Sticky",1,
                        "-props-registry-set","LayerMetaData.Alpha",settings['alpha'],
                        "-props-set-transform",transforms[layer["group"]]]
            if "midpoint_y" in layer:
                command += ["-props-registry-set","LayerMetaData.DisplayMapping.Curve.ControlPoint[1].xValue",layer["midpoint_y"]]
            manifest["layers"].append({**layer, "alpha":settings['alpha'], "display_file":str(display_paths[layer["file"]]),
                "data_range":[vol["min"],vol["max"]], "display_range":[lo,hi], "normalized_range":[norm_lo,norm_hi]})
        processing_group = None
        command += ["-o",workspace]
        run("itksnap-wt", command, "workspace_build")
        run("itksnap-wt", ["-i",workspace,"-layers-list"], "workspace_validate")
        ET.parse(workspace)
        for index, layer in enumerate(config["layers"],1):
            run("itksnap-wt", ["-i",workspace,"-layers-pick",index,"-P","-props-get-transform"], f"validate_transform_{index:02}")
            lines = (output/"logs"/f"validate_transform_{index:02}.log").read_text().splitlines()
            rows = []
            for line in lines:
                try:
                    row = [float(value) for value in line.split()]
                except ValueError:
                    continue
                if len(row) == 4:
                    rows.append(row)
            actual = np.array(rows)
            if actual.shape != (4,4):
                raise ValueError(f"Could not read stored transform: {layer['name']}")
            if not np.allclose(actual, np.loadtxt(transforms[layer["group"]]), atol=1e-4, rtol=0):
                raise ValueError(f"Stored workspace transform mismatch: {layer['name']}")
        if any(sha256(p) != digest for p,digest in inputs.items()):
            raise ValueError("An input file changed during the run; investigate before use")
        manifest["status"] = "built; anatomical review required"
        save_manifest()
        print('Copying required images and packaging a self-contained run...', flush=True)
        make_portable(output, backup=False)
        print_input_report(plan)
        if not plan['active_groups']:
            print('No complete optional modalities: created a T1-only workspace.', flush=True)
        print(f"Created {workspace}", flush=True)
        print('Portable output: move or copy this entire output folder, including its subfolders.', flush=True)
        if transforms:
            print("Fresh ITK transforms: " + str(output/"transforms"), flush=True)
        print('Modality report: ' + str(output/'modality_report.txt'), flush=True)
        return output
    except Exception as exc:
        manifest['status'] = 'failed'
        manifest["error"] = str(exc)
        if processing_group:
            plan['modalities'][processing_group].update(status='failed', error=str(exc))
        for record in plan['modalities'].values():
            if record['status'] in ('ready','processing'):
                record.update(status='not_run', detail='Run stopped after an error; see manifest.json')
        save_manifest()
        print_input_report(plan)
        raise


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", default=str(Path(__file__).with_name("mmii_config.json")))
    parser.add_argument("--images", help="Override input folder for the explicit filenames in the config")
    parser.add_argument("--output", help="New output directory (must not already exist)")
    parser.add_argument("--itksnap-bin", help="Override directory containing ITK-SNAP command-line executables")
    parser.add_argument('--check-inputs', action='store_true', help='Report missing/empty modality files without processing or writing output')
    try:
        build(parser.parse_args())
    except (OSError, ValueError, RuntimeError) as exc:
        print(f'ERROR: {exc}', file=sys.stderr)
        sys.exit(1)
