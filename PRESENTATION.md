# SimDist: processing heavy-hydrocarbon GC data

## Repository description

A Python/PySide6 app for processing high-temperature GC data, reviewing calibration peaks, and plotting carbon distributions and boiling-point recovery curves.

## About the project

I developed SimDist to simplify the analysis of high-temperature GC measurements from heavy-oil and bitumen samples. It brings calibration, blank correction, mass inputs, plotting, and CSV export into one desktop application.

## Calibration peaks

![Standard traces and selected peaks](docs/images/calibration_peaks.png)

The standard and Polywax chromatograms provide peaks for calibration. Users can detect peaks automatically, add or remove them, and review the carbon-number assignments before saving the calibration.

## Retention-time calibration

![Retention-time calibration curves](docs/images/calibration_mapping.png)

The selected peaks link retention time to equivalent carbon number and boiling point. The supplied calibration covers C5–C112, corresponding to 36–744 °C.

## From signal area to recovery

The standard's material fraction and blank-corrected area give the response factor:

```math
RF=\frac{M_{\mathrm{std}}}{M_{\mathrm{std}}+M_{\mathrm{sol,std}}}\frac{1}{A_{\mathrm{std}}}
```

The sample's material/solvent masses then convert its selected-range area into recovery:

```math
R_{\mathrm{range}}=100\,RF\,\frac{M_{\mathrm{smp}}+M_{\mathrm{sol,smp}}}{M_{\mathrm{smp}}}\,A_{\mathrm{smp,range}}
```

Carbon boundaries use linearly interpolated retention times. Integrating between them gives each histogram bar. The BP curve uses the same cumulative integral, with retention times converted to temperature by PCHIP interpolation. See [README.md](README.md#calculation-and-interpretation) for the variable definitions.

## Optional extrapolation

To extend carbon bins beyond the calibration, the app fits:

```math
t_r=a+b\ln C+c(\ln C)^2
```

It uses high-carbon calibration points, shifts the fit to meet the final calibrated point, and keeps predicted times increasing. The fit estimates additional bin boundaries within the recorded run.

Boiling-point extrapolation is separate: beyond the calibration, temperature continues linearly using the PCHIP curve's endpoint slope. For the supplied calibration, results beyond C112 are extrapolated. The logarithmic fit does not itself calculate boiling point.

## Gravimetric-blend comparison

![Gravimetric blend compared with reference values](docs/images/gravimetric_reference_comparison.png)

The measured blend curve is plotted against the supplied reference values. Its two main rises and plateau near 50% recovery reflect the blend's two fractions. The comparison also reports temperature differences for review.

This run reports 99.875% recovery with extension to C118. It uses the same chromatogram as standard and sample, so it checks the distribution shape without providing an independent recovery check. Reference tolerances still need to be assessed separately.

## Bitumen 1

![Bitumen 1 carbon distribution and cumulative recovery](docs/images/bitumen_1_distribution.png)

Blue bars show recovery in each carbon interval, and the red curve shows cumulative recovery. This run reports 60.282% within C5–C112, with a broad distribution and a tail toward higher carbon numbers.

The histogram and cumulative curve use the same integral. Signal outside the selected range is excluded from the reported recovery. The completed-analysis summary shows that no extension was applied in this screenshot.

## Development

The application replaces a script-based analysis with an interface for reviewing calibration, entering sample masses, comparing results, and exporting data. Recent fixes made recovery calculations consistent across plots and ensured imported calibration peaks refresh on screen.
