# UV-Vis Simulator API

A Python API for simulating UV-Vis absorption spectra from quantum-chemistry output files.

Supported features:

- ORCA, OpenMolcas, and cclib-supported formats
- Gaussian, Lorentzian, and Pseudo-Voigt broadening
- Individual and ensemble spectrum analysis
- Average spectra and standard deviation envelopes
- Excel export with raw transitions and simulated spectra
- Publication-quality PNG generation
- Fully scriptable Python interface

---

# Installation

```bash
pip install numpy pandas matplotlib xlsxwriter cclib
```

---

# Quick Start

```python
from uvvis_api import simulate_uvvis

simulate_uvvis(
    ["mol1.out", "mol2.out"],
    png_out="uvvis.png",
    excel_out="uvvis.xlsx"
)
```

This will:

- Parse the output files
- Simulate UV-Vis spectra using the default Pseudo-Voigt model
- Plot individual spectra (default mode)
- Export results to Excel

---

# Supported File Formats

The parser automatically detects the source software.

| Software | Parser |
|-----------|----------|
| ORCA | Native parser |
| OpenMolcas | Native parser |
| Gaussian | cclib |
| Q-Chem | cclib |
| NWChem | cclib |
| Other cclib-supported formats | cclib |

---

# High-Level API

## simulate_uvvis

```python
simulate_uvvis(
    files,
    *,
    sigma=0.33,
    n_states=3,
    weight=0.85,
    wavelength_range=(200, 800),
    n_components=3,
    plot_mode="individual",
    excel_out=None,
    png_out=None,
)
```

### Parameters

#### files

List of output files to process.

```python
files=["conf1.out", "conf2.out"]
```

---

#### sigma

Spectral broadening parameter (eV).

Default:

```python
sigma=0.33
```

Example:

```python
sigma=0.20
```

---

#### n_states

Number of excited states extracted from each file.

Default:

```python
n_states=3
```

Example:

```python
n_states=10
```

---

#### weight

Gaussian fraction in the Pseudo-Voigt profile.

| Value | Profile |
|---------|---------|
| 1.0 | Pure Gaussian |
| 0.0 | Pure Lorentzian |
| 0.85 | Default Pseudo-Voigt |

Examples:

```python
weight=1.0
```

Pure Gaussian broadening.

```python
weight=0.0
```

Pure Lorentzian broadening.

```python
weight=0.7
```

70% Gaussian, 30% Lorentzian.

---

#### wavelength_range

Wavelength range used for spectrum generation.

Default:

```python
(200, 800)
```

Example:

```python
wavelength_range=(180, 700)
```

---

#### n_components

Number of individual transition contributions displayed and exported.

Default:

```python
n_components=3
```

Example:

```python
n_components=5
```

---

#### plot_mode

Simulation and export mode.

Available values:

```python
"individual"
"average"
"combined"
```

Default:

```python
plot_mode="individual"
```

---

#### excel_out

Path to Excel output file.

Example:

```python
excel_out="spectrum.xlsx"
```

---

#### png_out

Path to PNG output file.

Example:

```python
png_out="spectrum.png"
```

---

# Plot Modes

## Individual (Default)

```python
simulate_uvvis(
    files,
    plot_mode="individual"
)
```

Features:

- One spectrum per file
- Individual transition components
- No averaging
- No standard deviation
- Filenames preserved in exported Excel columns

### Plot

```text
conf1
conf2
conf3
```

### Excel Export

```text
Wavelength (nm)
conf1_TOTAL
conf1_E1
conf1_E2

conf2_TOTAL
conf2_E1
conf2_E2
...
```

---

## Average

```python
simulate_uvvis(
    files,
    plot_mode="average"
)
```

Features:

- Average spectrum only
- No individual spectra
- No standard deviation envelope

### Summary Sheet

```text
Wavelength (nm)
Average
```

---

## Combined

```python
simulate_uvvis(
    files,
    plot_mode="combined"
)
```

Features:

- Ensemble average
- Standard deviation envelope
- Suitable for conformational or MD ensemble analysis

### Summary Sheet

```text
Wavelength (nm)
Average
SD
```

---

# Output Files

## PNG Export

```python
simulate_uvvis(
    files,
    png_out="spectrum.png"
)
```

Generated figure:

```text
spectrum.png
```

Characteristics:

- 300 dpi
- Transparent background
- Publication-ready output

---

## Excel Export

```python
simulate_uvvis(
    files,
    excel_out="spectrum.xlsx"
)
```

Generated workbook:

```text
spectrum.xlsx
```

### Sheet: Configuration

Contains simulation metadata:

```text
Date
Sigma (eV)
Model
Files
Plot mode
```

---

### Sheet: Summary

Generated only for:

```python
plot_mode="average"
```

or

```python
plot_mode="combined"
```

Contents:

**Average mode**

```text
Wavelength (nm)
Average
```

**Combined mode**

```text
Wavelength (nm)
Average
SD
```

---

### Sheet: Series_and_Components

Contains all simulated spectra and their individual transitions.

Example:

```text
Wavelength (nm)
conf1_TOTAL
conf1_E1
conf1_E2
conf1_E3

conf2_TOTAL
conf2_E1
conf2_E2
conf2_E3
```

---

### Sheet: Raw_Data

Contains the raw transition information extracted from the QM outputs.

Example:

```text
conf1_E_nm
conf1_f_osc

conf2_E_nm
conf2_f_osc
```

---

# Low-Level API

The simulation engine can be used directly.

```python
from uvvis_api import UVVisSimulator

sim = UVVisSimulator(
    sigma=0.33,
    n_states=10,
    weight=0.85,
    x_range=(200, 800),
)

sim.process_files(
    ["conf1.out", "conf2.out"]
)
```

---

# Accessing Results

## Simulated Spectra

```python
for result in sim.data["results"]:
    print(result["name"])
    print(result["total"])
```

Each entry contains:

```python
{
    "name": str,
    "total": ndarray,
    "peaks": list,
    "en_raw": list,
    "f_raw": list,
}
```

---

## Wavelength Grid

```python
x = sim.x
```

Returns the wavelength axis used for all calculations.

---

## Average and Standard Deviation

```python
avg, std = sim.calculate_statistics()
```

or

```python
avg = sim.data["avg"]
std = sim.data["std"]
```

---

# Custom Plotting

You may access the calculated spectra and create your own plots.

```python
import matplotlib.pyplot as plt

sim = simulate_uvvis(
    ["conf1.out", "conf2.out"]
)

for result in sim.data["results"]:
    plt.plot(
        sim.x,
        result["total"],
        label=result["name"]
    )

plt.xlabel("Wavelength (nm)")
plt.ylabel("Extinction Coefficient")
plt.legend()
plt.show()
```

---

# Examples

## Default Workflow

```python
simulate_uvvis(
    ["conf1.out", "conf2.out", "conf3.out"],
    png_out="individual.png",
    excel_out="individual.xlsx",
)
```

---

## Average Spectrum

```python
simulate_uvvis(
    ["conf1.out", "conf2.out", "conf3.out"],
    plot_mode="average",
    png_out="average.png",
    excel_out="average.xlsx",
)
```

---

## Average + Standard Deviation

```python
simulate_uvvis(
    ["conf1.out", "conf2.out", "conf3.out"],
    plot_mode="combined",
    png_out="combined.png",
    excel_out="combined.xlsx",
)
```

---

## Pure Gaussian Broadening

```python
simulate_uvvis(
    ["molecule.out"],
    weight=1.0,
)
```

---

## Pure Lorentzian Broadening

```python
simulate_uvvis(
    ["molecule.out"],
    weight=0.0,
)
```

---

# Error Handling

```python
try:
    simulate_uvvis(["missing.out"])
except ValueError as exc:
    print(exc)
```

Typical exceptions:

```text
At least one input file is required
```

```text
No spectra could be parsed from the supplied files
```

```text
plot_mode must be 'combined', 'average', or 'individual'
```

```text
No spectra have been processed
```

---

# License

Specify your project's license (MIT, BSD-3-Clause, GPL, etc.) here.
