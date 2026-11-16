"""Plan complete MMII modality groups before opening any image data.

This module checks file availability only. It does not infer anatomical
compatibility, read voxel data, substitute files, or write patient files.
Geometry and registration validation remain the builder's responsibility.
"""
from __future__ import annotations

from collections.abc import Mapping
from pathlib import Path, PurePosixPath, PureWindowsPath


def _component(value, location):
    """Require a single filename/path component on both Windows and POSIX."""
    if not isinstance(value, str) or not value or value.strip() != value:
        raise ValueError(f"{location} must be a nonempty name without surrounding whitespace")
    if value in (".", "..") or any(c in value for c in '<>:"/\\|?*') or any(ord(c) < 32 for c in value):
        raise ValueError(f"{location} must be a simple relative name without path components: {value!r}")
    if PureWindowsPath(value).name != value or PurePosixPath(value).name != value or value.endswith("."):
        raise ValueError(f"{location} must be a simple relative name: {value!r}")
    return value


def _image_filename(value, location):
    name = _component(value, location)
    if Path(name).suffix.lower() not in (".nii", ".hdr", ".img"):
        raise ValueError(f"{location} must name an uncompressed .nii or .hdr/.img image: {name!r}")
    return name


def _add_image(required, seen, name):
    """Add a volume's configured name and its paired sidecar if applicable."""
    pair = {".hdr": ".img", ".img": ".hdr"}
    names = [name]
    suffix = Path(name).suffix.lower()
    if suffix in pair:
        names.append(str(Path(name).with_suffix(pair[suffix])))
    for candidate in names:
        # Windows filenames are case insensitive. Preserve the first spelling
        # in the report while avoiding duplicate dependency entries.
        key = candidate.casefold()
        if key not in seen:
            required.append(candidate)
            seen.add(key)


def _availability(images, required):
    missing, empty = [], []
    for name in required:
        path = images / name
        if not path.is_file():
            missing.append(name)
        elif path.stat().st_size == 0:
            empty.append(name)
    return {"required_files": list(required), "missing_files": missing, "empty_files": empty}


def plan_modalities(config, images: Path):
    """Return a JSON-serializable plan for the mandatory main and all groups.

    A modality is ready only if every source/layer dependency exists as a
    nonempty regular file. NIfTI .hdr/.img pairs always require both members.
    A missing main is reported here; the caller decides when to stop the run.
    Invalid configuration raises ValueError instead of being called missing.
    """
    if not isinstance(config, Mapping):
        raise ValueError("Configuration must be an object")
    main_name = _image_filename(config.get("main"), "main")
    groups = config.get("groups")
    layers = config.get("layers")
    if not isinstance(groups, Mapping):
        raise ValueError("groups must be an object mapping modality names to definitions")
    if not isinstance(layers, list):
        raise ValueError("layers must be a list of layer objects")

    dependencies = {}
    seen = {}
    group_names = set()
    for group, spec in groups.items():
        _component(group, "group name")
        if group.casefold() in group_names:
            raise ValueError(f"Modality names must not differ only in case: {group!r}")
        group_names.add(group.casefold())
        if not isinstance(spec, Mapping):
            raise ValueError(f"groups.{group} must be an object")
        source = _image_filename(spec.get("source"), f"groups.{group}.source")
        dependencies[group], seen[group] = [], set()
        _add_image(dependencies[group], seen[group], source)

    for index, layer in enumerate(layers):
        if not isinstance(layer, Mapping):
            raise ValueError(f"layers[{index}] must be an object")
        group = layer.get("group")
        if not isinstance(group, str) or group not in groups:
            raise ValueError(f"layers[{index}].group refers to an unknown modality: {group!r}")
        filename = _image_filename(layer.get("file"), f"layers[{index}].file")
        _add_image(dependencies[group], seen[group], filename)

    images = Path(images).expanduser().resolve()
    if not images.is_dir():
        raise FileNotFoundError(f"Patient image folder does not exist or is not a directory: {images}")
    main_required = []
    _add_image(main_required, set(), main_name)
    main = {"file": main_name, **_availability(images, main_required)}
    main["status"] = "missing" if main["missing_files"] else "empty" if main["empty_files"] else "ready"
    modalities, active, skipped = {}, [], []
    for group, required in dependencies.items():
        result = _availability(images, required)
        unavailable = bool(result["missing_files"] or result["empty_files"])
        result["status"] = "skipped" if unavailable else "ready"
        modalities[group] = result
        (skipped if unavailable else active).append(group)
    return {"images": str(images), "main": main, "modalities": modalities,
            "active_groups": active, "skipped_groups": skipped}


def _unavailable_reason(record):
    parts = []
    if record.get("missing_files"):
        parts.append("missing: " + ", ".join(record["missing_files"]))
    if record.get("empty_files"):
        parts.append("empty: " + ", ".join(record["empty_files"]))
    for key in ("error", "detail"):
        if record.get(key):
            value = str(record[key]).strip()
            if value and value not in parts:
                parts.append(value)
    return "; ".join(parts)


def report_lines(plan):
    """Return a readable input report with explicit reasons for each skip."""
    main = plan["main"]
    lines = ["Patient images: " + plan["images"],
             f"Main reference ({main['file']}): {main['status'].upper()}"]
    reason = _unavailable_reason(main)
    if reason:
        lines[-1] += " - " + reason
    for group, record in plan["modalities"].items():
        line = f"{group}: {record['status'].upper()}"
        reason = _unavailable_reason(record)
        if reason:
            line += " - " + reason
        lines.append(line)
    lines.append("Active modalities: " + (", ".join(plan["active_groups"]) or "none"))
    if not plan["active_groups"]:
        lines.append("No complete modalities are available.")
    return lines


def print_input_report(plan):
    for line in report_lines(plan):
        print(line, flush=True)
