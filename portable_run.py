"""Copy a run's required images locally without recomputing registration.

ITK-SNAP relocates its absolute child image paths using SaveLocation. JSON paths
used by Python QC are relative to the run root. Historical commands and original
input hashes remain untouched as provenance, not runtime dependencies.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import shutil
import xml.etree.ElementTree as ET

from modality_inputs import report_lines


def digest(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(4 * 1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def _write_json(path, value):
    temporary = path.with_name(path.name + '.portable.tmp')
    temporary.write_text(json.dumps(value, indent=2, allow_nan=False) + '\n', encoding='utf-8')
    temporary.replace(path)


def make_portable(output, images=None, backup=True):
    """Package a completed builder run; optionally locate moved original inputs.

    No image contents, spatial transforms or display settings are changed.
    Existing copies are accepted only when their SHA-256 digest matches.
    A failure before metadata publication leaves the original workspace intact.
    """
    output = Path(output).resolve()
    manifest_path = output/'manifest.json'
    manifest = json.loads(manifest_path.read_text(encoding='utf-8'))
    config_path = output/'config_used.json'
    config = json.loads(config_path.read_text(encoding='utf-8'))
    if manifest.get('status') != 'built; anatomical review required' or manifest.get('error'):
        raise ValueError('Packaging requires a successfully completed MMII run')
    workspace = output/'MMII.itksnap'
    tree = ET.parse(workspace)
    save_location = tree.getroot().find("entry[@key='SaveLocation']")
    if save_location is None:
        raise ValueError('Workspace has no SaveLocation for resolving its recorded paths')
    old_root = Path(save_location.attrib['value'])
    was_portable = manifest.get('path_format') == 'relative_to_run_directory'
    if was_portable:
        prior_hashes = manifest.get('packaged_files_sha256')
        if not prior_hashes:
            raise ValueError('Portable run has no packaged-file hashes')
        for name, expected in prior_hashes.items():
            path = (output/name).resolve()
            if not path.is_relative_to(output) or not path.is_file() or digest(path) != expected:
                raise ValueError(f'Packaged file is missing or changed: {name}')
    qc_path = output/'qc'/'registration_metrics.json'
    qc = json.loads(qc_path.read_text(encoding='utf-8')) if qc_path.exists() else None
    # A workspace can have been resaved after a move while its Python run
    # metadata still records the earlier root. Retain those explicit roots too.
    historical_roots = {old_root}
    for group in manifest.get('groups', {}).values():
        for key in ('RAS_transform', 'itk_transform'):
            path = Path(group.get(key, ''))
            if path.is_absolute() and path.parent.name == 'transforms':
                historical_roots.add(path.parent.parent)
    main_path = Path(manifest['main'])
    if main_path.is_absolute() and main_path.parent.name == 'prepared':
        historical_roots.add(main_path.parent.parent)
    if qc and qc.get('manifest'):
        path = Path(qc['manifest'])
        if path.is_absolute() and path.name == 'manifest.json':
            historical_roots.add(path.parent)
    historical_roots = sorted(historical_roots, key=lambda p:len(p.parts), reverse=True)
    originals = manifest.get('inputs_sha256', {})
    if not originals:
        raise ValueError('Run has no original input hashes; cannot verify packaged inputs')
    input_dir = Path(images).resolve() if images else (
        output/'images' if was_portable else Path(config['images']))
    if not input_dir.is_absolute():
        input_dir = (output/input_dir).resolve()

    # Resolve every raw dependency and verify provenance before any copying.
    raw_map = {}
    dependencies = []
    destinations = set()
    for recorded, expected in originals.items():
        old_path = Path(recorded)
        target = output/'images'/old_path.name
        if target.name.casefold() in destinations:
            raise ValueError(f'Ambiguous input filename: {target.name}')
        destinations.add(target.name.casefold())
        candidate = input_dir/old_path.name
        if not candidate.is_file() and not images and not was_portable:
            candidate = old_path
        if not candidate.is_file():
            raise FileNotFoundError(f'Required original input missing: {candidate}; supply --images if moved')
        if digest(candidate) != expected:
            raise ValueError(f'Input differs from the completed registration run: {candidate}')
        raw_map[str(old_path.resolve()).casefold()] = target
        dependencies.append((candidate, target, expected))

    def resolve_recorded(value):
        path = Path(value)
        if not path.is_absolute():
            result = (output/path).resolve()
        elif str(path.resolve()).casefold() in raw_map:
            result = raw_map[str(path.resolve()).casefold()]
        else:
            result = path.resolve()
            for root in historical_roots:
                if path.is_relative_to(root):
                    result = (output/path.relative_to(root)).resolve()
                    break
        if not result.is_relative_to(output):
            raise ValueError(f'Untracked dependency outside the run: {value}')
        return result

    def relative(value):
        return resolve_recorded(value).relative_to(output).as_posix()

    # Prepare all metadata in memory first. Unknown dependencies fail closed.
    manifest.setdefault('original_images_dir', manifest['images_dir'])
    manifest['main'] = relative(manifest['main'])
    manifest['images_dir'] = 'images'
    required = {resolve_recorded(manifest['main'])}
    for group in manifest['groups'].values():
        for key in ['source', 'RAS_transform', 'itk_transform']:
            if key in group:
                group[key] = relative(group[key])
                required.add(resolve_recorded(group[key]))
    for layer in manifest['layers']:
        layer['display_file'] = relative(layer['display_file'])
        required.add(resolve_recorded(layer['display_file']))
    entries = tree.getroot().findall(".//entry[@key='AbsolutePath']")
    if not entries:
        raise ValueError('Workspace has no image layers')
    for entry in entries:
        local = resolve_recorded(entry.attrib['value'])
        required.add(local)
        # This is the native workspace convention, not a dependency on this
        # directory forever: SNAP rebases children when SaveLocation changes.
        entry.set('value', local.as_posix())
    save_location.set('value', output.as_posix())
    config['images'] = 'images'
    config['path_format'] = 'relative_to_run_directory'
    for name, spec in config['groups'].items():
        if spec.get('method') == 'external_ras':
            spec['transform'] = manifest['groups'][name]['RAS_transform']

    (output/'images').mkdir(exist_ok=True)
    for source, target, expected in dependencies:
        if target.exists():
            if not target.is_file() or digest(target) != expected:
                raise ValueError(f'Conflicting packaged image: {target}')
        else:
            temporary = target.with_name(target.name + '.portable.tmp')
            shutil.copy2(source, temporary)
            if digest(temporary) != expected:
                raise ValueError(f'Copied image verification failed: {source}')
            temporary.replace(target)
        required.add(target)

    # A NIfTI pair is one image, but both files must travel with the bundle.
    for path in list(required):
        if path.suffix.lower() in ('.hdr', '.img'):
            required.add(path.with_suffix('.img' if path.suffix.lower() == '.hdr' else '.hdr'))
    for path in required:
        if not path.is_file() or path.stat().st_size == 0:
            raise FileNotFoundError(f'Incomplete portable bundle: {path}')
    manifest['packaged_files_sha256'] = {
        path.relative_to(output).as_posix(): digest(path) for path in sorted(required)}
    manifest['path_format'] = 'relative_to_run_directory'
    manifest['self_contained'] = True
    manifest['packaged_utc'] = datetime.now(timezone.utc).isoformat()
    manifest['packaging_module_sha256'] = digest(Path(__file__))
    manifest['portability_note'] = ('Move or copy the entire run folder. ITK-SNAP rebases workspace child paths using '
                                  'SaveLocation; Python operational paths are relative to this directory. '
                                  'Original paths in hashes and commands are historical provenance only.')

    report_path = output/'modality_report.json'
    plan = json.loads(report_path.read_text(encoding='utf-8')) if report_path.exists() else None
    if plan:
        plan['images'] = 'images'
    if qc:
        def rewrite_qc(record):
            if isinstance(record, dict):
                for key, value in record.items():
                    if key == 'images_dir':
                        record[key] = 'images'
                    elif key in ('path', 'manifest', 'config_used', 'transform_path', 'montage') and isinstance(value, str):
                        record[key] = relative(value)
                    else:
                        rewrite_qc(value)
            elif isinstance(record, list):
                for value in record:
                    rewrite_qc(value)
        rewrite_qc(qc)
        qc['path_format'] = 'relative_to_run_directory'

    if backup:
        backup_dir = output/('portability_backup_' + datetime.now().strftime('%Y%m%d_%H%M%S_%f'))
        backup_dir.mkdir()
        for path in [workspace, manifest_path, config_path, report_path, output/'modality_report.txt', qc_path]:
            if path.exists():
                dest = backup_dir/path.relative_to(output)
                dest.parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(path, dest)
    tmp_workspace = workspace.with_suffix('.portable.tmp')
    tree.write(tmp_workspace, encoding='utf-8', xml_declaration=True)
    tmp_workspace.replace(workspace)
    _write_json(config_path, config)
    if plan:
        _write_json(report_path, plan)
        (output/'modality_report.txt').write_text('\n'.join(report_lines(plan))+'\n', encoding='utf-8')
    if qc:
        _write_json(qc_path, qc)
    _write_json(manifest_path, manifest)
    return manifest


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', required=True, help='Existing successfully completed run directory')
    parser.add_argument('--images', help='Current original image directory, if it has moved')
    args = parser.parse_args()
    try:
        result = make_portable(args.output, args.images)
        print('Self-contained workspace: ' + str(Path(args.output).resolve()/'MMII.itksnap'))
        print(f"Verified bundled files: {len(result['packaged_files_sha256'])}")
    except (OSError, ValueError, KeyError) as exc:
        parser.exit(1, f'ERROR: {exc}\n')
