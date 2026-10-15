# MAS glass simulation data and analysis code

This repository contains molecular dynamics simulation inputs and outputs, selected glass configurations, a structural-analysis script, and experimental Fourier transform infrared (FTIR) spectra for MgO–Al₂O₃–SiO₂ (MAS) glasses. The files support the associated manuscript:

**Effect of the Al₂O₃/SiO₂ ratio on the structure and Young’s modulus of MgO–Al₂O₃–SiO₂ glasses: Experiments and molecular dynamics simulations**

The shared package covers compositions B1–B5. MgO is fixed at 23.4 mol%, while Al₂O₃ replaces SiO₂.

## Repository structure

```text
.
├── README.md
├── Analysis/
│   └── multi_xyzFiles_process.py
├── Data/
│   └── FTIR/
│       ├── B1.CSV
│       ├── B2.CSV
│       ├── B3.CSV
│       ├── B4.CSV
│       └── B5.CSV
└── Simulation/
    ├── B1/
    ├── B2/
    ├── B3/
    ├── B4/
    └── B5/
```

Directory and file names are case-sensitive on Linux. `Analysis/` contains the post-processing code; `Simulation/` contains the simulation files; `Data/FTIR/` contains the experimental spectra.

## Compositions and atom counts

Oxide compositions are expressed in mol%. The atom counts below correspond to the supplied initial LAMMPS data files.

| Sample | SiO₂ | Al₂O₃ | MgO | Si atoms | Al atoms | Mg atoms | O atoms | Total atoms |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| B1 | 57.2 | 19.4 | 23.4 | 1716 | 1164 | 702 | 5880 | 9462 |
| B2 | 55.2 | 21.4 | 23.4 | 1656 | 1284 | 702 | 5940 | 9582 |
| B3 | 53.2 | 23.4 | 23.4 | 1596 | 1404 | 702 | 6000 | 9702 |
| B4 | 51.2 | 25.4 | 23.4 | 1536 | 1524 | 702 | 6060 | 9822 |
| B5 | 49.2 | 27.4 | 23.4 | 1476 | 1644 | 702 | 6120 | 9942 |

## Simulation files

Each `Simulation/Bn/` directory contains the following files, where `n` is the sample number.

| File | Description |
|---|---|
| `n.data` | Initially constructed atomic model read by the corresponding `in.Du` file. For example, B1 uses `1.data`. |
| `400.xyz` | Single glass configuration supplied for structural analysis. It contains element labels, Cartesian coordinates, atom IDs, and simulation-cell information in an extended-XYZ comment line. The filename identifies the snapshot; it does not mean that the file contains 400 frames. |
| `in.Du` | LAMMPS input script containing the charges, Buckingham potential parameters, electrostatic settings, melt–quench and relaxation stages, and subsequent tensile deformation along the z direction. |
| `log.lammps` | LAMMPS run log, including echoed commands, thermodynamic quantities, timing, and runtime messages. This is a text log rather than an atomic-coordinate trajectory. |
| `stretch.xyz` | Multi-frame atomic-coordinate trajectory recorded during tensile deformation. |
| `stretch.data` | LAMMPS data file written at the end of the tensile run; it represents the deformed configuration. |
| `stressz-strainz.form` | Two-column tensile output: engineering strain in the z direction, followed by the corresponding stress in GPa. The supplied tables have no column headers. |
| `n.xlsx` | Supplied structural-analysis workbook, where present. The current package includes `2.xlsx`, `3.xlsx`, `4.xlsx`, and `5.xlsx`; no B1 workbook is included. |

The Excel workbooks contain structural statistics such as coordination numbers, oxygen types, Qⁿ statistics, network linkages, and oxygen-tricluster connections. Some workbooks also include Mg-centered distribution statistics. Sheet availability follows the individual supplied workbook. The analysis script writes a new workbook named `multifile_cn_and_oxygen_statistics.xlsx` by default.

### Simulation settings and units

The supplied logs identify **LAMMPS 22 Dec 2022**. The inputs use `units metal`, `atom_style charge`, periodic boundaries in all three directions, and a timestep of 0.001 ps (1 fs). Compatibility with other LAMMPS releases has not been assessed here.

| Atom type | Element | Assigned charge (elementary charge units) |
|---|---|---:|
| 1 | Si | +2.4 |
| 2 | Al | +1.8 |
| 3 | Mg | +1.2 |
| 4 | O | −1.2 |

Potential coefficients are defined directly in each `in.Du`; a separate potential file is not required by these inputs. The B1 input uses `pppm`, while B2–B5 use `pppm/stagger`, as recorded in the supplied files. Use a LAMMPS build that supports the specified pair style and long-range solver.

The scripts contain 300 K relaxation, a 4500 K melt stage, cooling from 4500 to 300 K at 2.5 K/ps, further relaxation, and tensile deformation. The tensile engineering strain rate specified by `fix deform` is 0.05 ps⁻¹. The exact commands and stage durations are retained in the individual input files.

Relevant units are:

- Coordinates and box lengths: Å.
- Time: ps; temperature: K; energy: eV.
- Pressure and stress in the thermodynamic log: bar, following LAMMPS `metal` units.
- Tensile stress in `stressz-strainz.form`: GPa, calculated in the script as `-pzz/10000`.
- Engineering strain: dimensionless, calculated as `(Lz − Lz₀)/Lz₀`.

The `.form` file records the supplied deformation response; it does not itself contain a fitted Young’s modulus. An elastic fit requires an appropriate, explicitly documented strain interval and analysis method.

### Running a supplied simulation

Run from the sample's working directory so that `read_data n.data` resolves correctly. To preserve the archived output files, copy the input and initial model into a separate working directory first. For B1, the following example assumes a POSIX shell and starts at the repository root:

```bash
mkdir -p work/lammps_B1
cp Simulation/B1/in.Du Simulation/B1/1.data work/lammps_B1/
cd work/lammps_B1
lmp -in in.Du
```

Replace `lmp` with the executable or parallel launch command for the local installation. Repeat with the matching directory and numbered data file for B2–B5. The simulation inputs also request coordinate outputs such as `quench.coordinate` and `NVE.coordinate`; these are generated on execution and are not included in the current shared package.

## Structural analysis

### Dependencies

The Python script imports MDAnalysis, NumPy, and pandas. An Excel writer such as openpyxl is needed to save `.xlsx` output. Install compatible versions in a Python environment, for example:

```bash
python -m pip install MDAnalysis numpy pandas openpyxl
```

The original Python package versions were not recorded in the supplied files. The script is interactive and displays prompts in Chinese.

### Coordination settings

Use the following settings for the supplied MAS glass configurations:

| Cation–oxygen pair | Distance cutoff (Å) |
|---|---:|
| Si–O | 2.19 |
| Al–O | 2.53 |
| Mg–O | 2.70 |

Neighbors are selected using distances **less than or equal to** the specified cutoff. Enable periodic boundary conditions and use Mg as the modifier species. The script reads the orthorhombic box lengths from the `Lattice` field in the supplied extended-XYZ files when necessary.

### Running the analysis on one supplied configuration

The script scans **all `.xyz` files in the selected directory** and analyzes the initially loaded configuration of each file. It does not iterate over every frame of a multi-frame trajectory. Therefore, place only the intended structural-analysis snapshots in the selected input directory. In particular, keep `stretch.xyz` separate from the analysis snapshots.

For B1, start at the repository root and use a dedicated working directory:

```bash
mkdir -p work/analysis_B1/snapshots
cp Simulation/B1/400.xyz work/analysis_B1/snapshots/
cd work/analysis_B1
python ../../Analysis/multi_xyzFiles_process.py
```

Enter the following responses, without quotation marks, when prompted:

| Prompt meaning | Response |
|---|---|
| Directory containing the `.xyz` files | `snapshots` |
| Cation species to analyze | `Si Al Mg` |
| Corresponding cation–oxygen cutoffs, in the same order | `2.19 2.53 2.70` |
| Periodic boundary conditions | `Yes` |
| Modifier species for distribution analysis | `Mg` |

No boron-specific question is expected for these MAS inputs. The output workbook is saved in the current working directory as `multifile_cn_and_oxygen_statistics.xlsx`. Use separate working directories for different compositions to keep their outputs distinct.

The workbook includes coordination counts, oxygen classifications, Qⁿ and detailed Qⁿ statistics, bridging-oxygen linkages, tricluster connection statistics, non-bridging-oxygen-associated modifier counts, and tricluster counts within Al coordination units. Mg-centered distribution statistics are also requested by the example above. Many sheets contain counts or per-file averages; only columns explicitly labeled as percentages should be read directly as percentages.

For oxygen classification, the script treats Si and four-coordinated Al as network formers in its bridging/non-bridging oxygen rules. Its oxygen-tricluster rule also includes five- and six-coordinated Al. Consult `classify_oxygen` and the other analysis functions for the implemented definitions. Output categories should be interpreted using those definitions.

## Experimental FTIR spectra

`Data/FTIR/B1.CSV` through `B5.CSV` contain the spectra for the corresponding glass compositions. Each file is a comma-separated, two-column numerical table without a header:

| Column | Meaning |
|---|---|
| 1 | Wavenumber, in cm⁻¹; approximately 399–4001 cm⁻¹ in the supplied exports. |
| 2 | FTIR transmittance signal on the supplied numerical scale. |

The manuscript describes transmittance measurements and labels the plotted FTIR ordinate **“Transmittance (a.u.)”**. The CSV files do not specify a percent unit, so the second column is preserved as supplied and is not assumed or converted to percent transmittance. Plotting in descending wavenumber order does not require changing the stored data order.

## Scope of this shared package

One initial model and one selected analysis configuration are included for each of B1–B5, together with the supplied simulation logs and tensile outputs. `400.xyz` is a single snapshot in each directory. Reproducing a multi-frame or independent-run average requires the actual frames and runs used for that average; the single-snapshot example above demonstrates the analysis workflow rather than reproducing an entire sampling ensemble. B6 and the separate elastic-tensor calculation workflow are not part of this package.

## Citation and contact

When using these files, cite the associated manuscript and the LAMMPS and analysis tools as appropriate. The manuscript title is given above; no publication DOI is assigned in this README.

Corresponding author: **Junfeng Kang**, School of Materials Science and Engineering, University of Jinan, Jinan 250002, China. Email: **mse_kangjf@ujn.edu.cn**.
