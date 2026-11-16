# MMII workspace automation

Create an ITK-SNAP workspace from one patient's existing imaging outputs. Supply a `patientName_Images` folder containing `t1.nii` and whichever optional modalities are available. The pipeline creates fresh transforms, applies them to the corresponding result maps, and saves `MMII.itksnap` in a self-contained output folder. It does not rerun VBM, MELD, PASCOM or SISCOM analysis.

## Quick start

Use these steps on Windows with Python 3.10 or newer and ITK-SNAP installed. Copy the **complete** `mmii_automation` folder to the PC first. The examples assume it is in Downloads; change the paths if you put it elsewhere.

1. Open PowerShell and install the Python packages once:

   ```powershell
   Set-Location "$env:USERPROFILE\Downloads\mmii_automation"
   python -m pip install -r requirements.txt
   ```

2. Open `mmii_config.json` and set `itksnap_bin` to your ITK-SNAP tools folder. For example, the tested setup uses `C:/Program Files/ITK-SNAP 4.4/bin`. That folder must contain `greedy.exe`, `c3d.exe`, `c3d_affine_tool.exe` and `itksnap-wt.exe`. See **Setup on another PC or ITK-SNAP version** below if your installation differs.

3. Put `t1.nii` and your available modality files directly in `patientName_Images`, using the filenames in the table below. Run:

   ```powershell
   .\run_mmii.ps1 -Images "$env:USERPROFILE\Downloads\patientName_Images" -WithQC
   ```

   Replace `patientName_Images` with your actual folder name or full path. An incomplete optional modality is skipped as a whole and reported; the other complete modalities continue.

4. Open `MMII.itksnap` from the new `patientName_MMII_<timestamp>` folder beside your input folder. Review the alignment and the images in its `qc` folder. To move the result, copy the **whole output folder**; its required images are included automatically.

## Run options and input checks

With the scripts in your Downloads folder, this command works from any PowerShell working directory:

```powershell
& "$env:USERPROFILE\Downloads\mmii_automation\run_mmii.ps1" -Images "$env:USERPROFILE\Downloads\patientName_Images" -WithQC
```

Replace `patientName_Images` with the actual patient's folder. The launcher uses `python` from PATH. To choose another interpreter, add `-PythonExe 'C:\path with spaces\python.exe'`.

The default output is a new sibling folder such as `patientName_MMII_20261114_143012_123`, beside `patientName_Images`. The input suffix `_Images` is removed without regard to case. Use `-Output 'D:\MMII results\patientName_run1'` to choose a different new folder. Existing output folders are refused rather than overwritten.

Check which modalities can be included before running registration:

```powershell
& "$env:USERPROFILE\Downloads\mmii_automation\run_mmii.ps1" -Images "$env:USERPROFILE\Downloads\patientName_Images" -CheckInputs
```

This checks whether required input files are present and nonempty, without creating an output folder or running registration or QC. It does not validate image contents or native-tool compatibility. `-CheckInputs` takes precedence if `-WithQC` is also supplied.

## Required filenames

All files belong directly inside the folder passed to `-Images`. Filenames are not inferred from the patient name. In particular, the main image is always `t1.nii`, and the VBM pair is `probFCD.hdr` plus `probFCD.img`.

| Modality | Files needed to include the entire group | Alignment method |
|---|---|---|
| Main T1 — required | `t1.nii` | Fixed reference for every overlay |
| FLAIR | `flair.nii` | Automatic rigid NMI |
| PASCOM | `wwt1_u_rc1t1_Template.nii`, `wwpet_u_rc1t1_Template.nii`, `Hypo_PASCOM_without_C.nii`, `Hyper_PASCOM_without_C.nii` | Register the accompanying T1, then propagate to all four layers |
| SISCOM | `mr.nii`, `SISCOM_without_C.nii` | Register `mr.nii`, then propagate to the result map |
| MELD | `T1w.nii`, `prediction.nii` | Identity if the source T1 is byte-identical to the main T1; otherwise rigid NMI |
| VBM | `probFCD.hdr`, `probFCD.img` | Physical image-center alignment |
| Texture | `Texture_MRI.nii` | Physical image-center alignment |

Missing or empty `t1.nii` stops the run. If any required file is missing or empty, that optional group is skipped in full, including its registration, layers and QC: for example, `prediction.nii` without `T1w.nii` does not create a MELD layer. This rule applies to FLAIR, PASCOM, SISCOM, MELD, VBM and Texture. Other complete groups still run. A folder containing only `t1.nii` is allowed and produces a main-image-only workspace with a warning. Extra files are ignored.

Missing optional files and invalid data are different conditions. Invalid images, unsupported geometry, mismatched grids within a group, or failed registration stop the build instead of silently discarding a problematic available modality. Supported inputs are uncompressed, scalar 3-D NIfTI-1 volumes with explicit millimetre geometry, including NIfTI `.hdr`/`.img` pairs. Legacy Analyze pairs without reliable NIfTI spatial information are rejected.

## Requirements and configuration

Use Python 3.10 or newer with NumPy and SciPy. Optional QC also requires Matplotlib. The package list is in `requirements.txt`. Install these packages in the Python environment chosen with `-PythonExe`.

Install ITK-SNAP with its `greedy.exe`, `c3d.exe`, `c3d_affine_tool.exe` and `itksnap-wt.exe` command-line executables. Set `itksnap_bin` in `mmii_config.json` to their folder; the supplied setting points to `C:/Program Files/ITK-SNAP 4.4/bin`. MATLAB and SPM are not required for this pipeline.

The supplied configuration has `images: null`; it does not select a previous patient's data. The PowerShell launcher requires `-Images` for each invocation. The Python builder requires `--images` unless a custom configuration supplies a nonempty `images` value. `-Config` selects an alternative JSON configuration; layer names, display settings and registration methods can be changed there.

### Setup on another PC or ITK-SNAP version

Copy all five Python modules, `run_mmii.ps1`, `mmii_config.json` and `requirements.txt`. In particular, `portable_run.py` is required by the builder. Install the packages with the same Python interpreter that will run the pipeline; the quick start uses `python` from PATH. For a different interpreter, use its full path for installation and pass it to the launcher:

```powershell
& 'C:\path with spaces\python.exe' -m pip install -r requirements.txt
.\run_mmii.ps1 -Images 'D:\Patients\patientName_Images' -PythonExe 'C:\path with spaces\python.exe' -WithQC
```

Run these examples from the pipeline folder and replace the example paths. In **`mmii_config.json`** (the supplied configuration filename), update the existing `itksnap_bin` field to the actual tools directory on the destination PC. For example, if that installation is in an ITK-SNAP 4.2 folder, the field would be:

```json
"itksnap_bin": "C:/Program Files/ITK-SNAP 4.2/bin"
```

This is one field within the existing JSON object, not a replacement for the entire file. Use forward slashes as shown, or escape Windows backslashes as `\\`. Preserve valid JSON, including commas between fields, and do not add comments. The directory must contain **all four executables** listed above; the builder checks their presence before a full build. It does not automatically search for another installation.

**Changing the path selects an installation; it does not establish compatibility with that version.** ITK-SNAP 4.4 is the tested setup for this pipeline. Other versions can differ in command-line options or workspace behavior. Review a representative full run, workspace and QC on the new setup, including reopening a copy of the complete output folder from another location. `-CheckInputs` alone cannot validate this. If a required executable or option is unavailable, use a compatible installation rather than changing registration or workspace commands blindly.

If the four-module `tests` folder is supplied with your copy, transfer it too. From the pipeline folder, run `python -B -m unittest discover -s tests -v`. Check failures, zero discovered tests and skipped native-tool tests; an empty or partially skipped run is not evidence of compatibility. These development checks supplement anatomical review. They are not required to open a completed workspace.

### Direct Python commands and later QC

To run Python directly from any working directory:

```powershell
python "$env:USERPROFILE\Downloads\mmii_automation\build_mmii.py" --images "$env:USERPROFILE\Downloads\patientName_Images"
```

The builder also accepts `--check-inputs`, `--config`, `--output` and `--itksnap-bin`. The last option overrides the configured tools folder for that invocation; it is a **Python builder** option, not a `run_mmii.ps1` parameter. The builder's automatic output name uses a microsecond timestamp. To create QC later, point it at a completed run; the run manifest supplies the images and reference:

```powershell
python "$env:USERPROFILE\Downloads\mmii_automation\qc_registration.py" --output 'D:\MMII results\patientName_run1'
```

## Outputs and review

Open the generated `MMII.itksnap` using ITK-SNAP's Workspace menu. Every new successful run is self-contained: the builder calls `portable_run.make_portable()` automatically, with no extra launcher switch. Its `images` directory contains verified copies of all required inputs for included groups (including both `.hdr` and `.img`), while `prepared` contains any cleaned display copies. Files from skipped groups and unrelated input files are not copied. This uses additional disk space approximately equal to the included original inputs, on top of any prepared copies.

Move, rename, copy or zip the **entire output folder**, keeping its internal files and subfolders together. After extracting or moving it, open `MMII.itksnap` from the new location. The original patient image folder and pipeline installation are not needed to view it; ITK-SNAP must be installed. Copying only `MMII.itksnap` is insufficient.

The workspace uses ITK-SNAP's native `SaveLocation` relocation mechanism: its stored absolute image paths point inside the saved run, so ITK-SNAP resolves the same child path beneath the moved run. Keep the internal filenames and folder arrangement unchanged, and do not manually edit `SaveLocation`. Python operational paths in `manifest.json`, `config_used.json` and new QC reports are relative to the run directory, not the `qc` subfolder. QC can be regenerated after moving using the same `--output` interface, with the pipeline code and Python dependencies available. Original absolute paths in historical command lists and input hashes remain as provenance only.

### Convert an older run to a portable output

Use `portable_run.py` to package a successfully completed run without repeating registration. If the original images have moved, supply their current folder with `--images`; omit that option if they are still available at the recorded location:

```powershell
python "$env:USERPROFILE\Downloads\mmii_automation\portable_run.py" --output 'D:\Results\patientName_MMII_run' --images 'D:\Patients\patientName_Images'
```

The helper checks original input SHA-256 hashes before copying, verifies the copies, and preserves image contents, transforms and display settings. Before publishing updated metadata, standalone conversion backs up the previous workspace and metadata to `portability_backup_<timestamp>` inside the run. Automatic packaging during a new build does not create that backup. A missing or changed input stops conversion. Repackaging an already portable run also checks its previous packaged hashes and rejects missing or changed prepared images and referenced transforms. Do not open the backup workspace for normal use. Existing unrelated files are preserved.

### Run records and review

Each run records:

- `modality_report.txt` and `modality_report.json`: available and skipped groups, including missing dependencies.
- `manifest.json`: input paths and hashes, selected groups and layers, transforms, exact commands, warnings and build status. Portable runs also record `self_contained`, `path_format`, `packaged_files_sha256`, `packaged_utc`, `packaging_module_sha256`, `original_images_dir` and `portability_note`.
- `images`: byte-for-byte copies of required included input files; `packaged_files_sha256` in the manifest records the packaged images and referenced transforms.
- `config_used.json`: the run configuration, with `images` set to `images` and operational paths relative to the run root after packaging.
- `transforms`: generated ITK LPS transforms and corresponding `_RAS.txt` matrices for included groups.
- `logs`: registration output and workspace verification.
- With `-WithQC`, `qc/anatomy_checkerboard.png` and `qc/registration_metrics.json`: anatomical checkerboards and descriptive registration measurements.

The original images and supplied transform files are not changed. Input SHA-256 hashes are checked before and after a successful build. Registration uses newly calculated transforms, not previously supplied `trans_*.txt` files.

Check the overall manifest status before relying on an output: a processed group alone does not mean the build succeeded. Invalid data, registration errors, failed packaging or failed verification stop the build. If optional QC fails after a successful build, the workspace remains available and the launcher reports the QC failure separately.

Review the midline, cortical outlines, ventricles, cerebellum and full anatomical coverage across slices. The numerical checks detect some gross failures but do not establish anatomical or clinical acceptance. Use each group's anatomical reference when judging alignment of its statistical result map. QC draws on the included groups recorded in the manifest, so it does not require absent modalities.

## Python files

| File | Responsibility |
|---|---|
| `build_mmii.py` | Validate images, estimate transforms, configure layers, verify the workspace and package a new run. |
| `modality_inputs.py` | Determine each group's required files, include complete groups and report missing or empty dependencies. |
| `anatomical_init.py` | Estimate PASCOM's initial translation and score anatomical alignment. |
| `qc_registration.py` | Read completed-run metadata, calculate registration measurements and create checkerboards. |
| `portable_run.py` | Copy and verify included inputs, update workspace and metadata paths, and package older runs without repeating registration. |

See [MMII_Pipeline_Technical_Guide.docx](MMII_Pipeline_Technical_Guide.docx) for the 15-page explanation of configuration, algorithms, transform conventions, portability and validation.


## References

- [ITK-SNAP command-line tools](https://itksnap.org/pmwiki/pmwiki.php?n=Documentation.CommandLine)
- [Official workspace creation example](https://ashs.itksnap.org/ashs/cloud/command-line/)
- [Greedy registration options](https://greedy.readthedocs.io/en/latest/reference.html)

The supplied MMII tutorial informed the overlay order, colors, contrast settings and registration workflow. Its instructions to install or run upstream analysis software are separate from this workspace-assembly pipeline.
