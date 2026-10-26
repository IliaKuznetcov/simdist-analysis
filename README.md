# SimDist — v1.0.0

Python/PySide6 desktop application for processing high-temperature gas chromatography data into n-paraffin-equivalent carbon distributions and boiling-point–recovery curves. Developed by Ilia Kuznetcov for heavy-hydrocarbon research.

## From chromatogram to engineering result

The workflow turns detector measurements into an interpretable heavy-hydrocarbon distribution: review standard peaks, map retention time to carbon number and boiling point, inspect a gravimetric-blend reference comparison, and process the unknown sample.

![Bitumen 1 carbon distribution and cumulative recovery](docs/images/bitumen_1_distribution.png)

*Bitumen 1: the captured run reports 60.282% recovery within C5–C112. Blue bars show recovery per carbon interval; the red curve accumulates the same bin contributions. This is selected-range recovery, not recovery normalized to 100%.*

The example demonstrates experimental data processing, numerical integration, interactive scientific software, and explicit handling of calibration limits. See [EXAMPLE.md](EXAMPLE.md) for settings and image captions, and [PRESENTATION.md](PRESENTATION.md) for a short technical walkthrough.

## Features

- Interactive retention-standard and Polywax calibration; JSON import/export.
- Blank correction and gravimetric-standard response scaling.
- Carbon distributions, cumulative recovery, BP curves, and optional QC comparison.
- Optional high-carbon extension, BP extrapolation, and sample solvent-spike removal.
- CSV exports and Matplotlib plot controls.

## Install and run

Use your existing virtual environment, or create one:

```powershell
python -m venv .venv
.\.venv\Scripts\Activate.ps1
python -m pip install -r requirements.txt
python simdist_app.py
```

Dependencies are listed without version pins; use your tested environment for the release build.

## Example

See [EXAMPLE.md](EXAMPLE.md) for the supplied three-bitumen dataset, masses, settings, and expected results. Load `calibration.json` before analysis. For calibration editing, use the supplied paraffin BP spreadsheet; startup defaults are not the example's calibration.

Optional peak editing, plot navigation, and a gravimetric-blend reference check are also described in [EXAMPLE.md](EXAMPLE.md).

## Calculation and interpretation

Signals are aligned to the standard time grid, clipped at zero, blank-subtracted, and clipped again. The standard area and material/solvent mass ratios set the response scale. One piecewise-linear integral supplies carbon-bin areas and the BP curve, including partial acquisition intervals at bin boundaries.

Cumulative recovery starts at zero at the lower selected carbon boundary and ends at the selected-range total. Results remain on a sample-mass basis; they are not normalized to 100%. BP mapping uses shape-preserving interpolation of the calibration points. Extension beyond the measured range is an estimate, marked on the plots.

The quantity `100 − selected-range recovery` is material not accounted for within that range, not an independently measured non-eluting residue. Carbon numbers are n-paraffin equivalents, not molecular identification. Blank selection, calibration applicability, baseline treatment, and solvent exclusion affect interpretation. The optional QC comparison does not automatically enforce the reference's `Deviation` tolerances.

Research software supporting an ASTM D7169-related workflow; no ASTM compliance certification is claimed. Numerical integration and endpoint checks passed on the supplied examples; the author also checked the corrected app locally on Windows.

