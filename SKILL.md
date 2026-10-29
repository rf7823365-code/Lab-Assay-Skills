# Lab Assay Calculator Skill

Calculate enzyme activity rates (nmol/min) from kinetic spectrophotometry data.

Supports the three assay types used in the G3DH / 3-ketoglucose / 3-ketosucrose work.

## When to use this skill

- Rob provides absorbance vs. time data from a DU800 or similar spectrophotometer
- Calculating enzyme activity from ferricyanide, phosphate, or DCIP kinetics
- Comparing rates between fractions, dilutions, or substrates
- Checking whether a reaction is linear (depletion warning)

## Tool location

`assay_calculator.py` — in this repo (`Lab-Assay-Skills/assay_calculator.py`)

## Quick reference

```bash
# Ferricyanide kinetic assay at 455 nm (default)
python3 assay_calculator.py --assay ferricyanide --file data.csv

# Ferricyanide at 420 nm (stopped endpoint)
python3 assay_calculator.py --assay ferricyanide --file data.csv --wavelength 420

# Phosphate assay (Hayano-Fukui protocol)
# Auto-skips first 180 s; requires two-column absorbance file (310 nm, 800 nm)
python3 assay_calculator.py --assay phosphate --file data.csv

# DCIP assay
python3 assay_calculator.py --assay dcip --file data.csv

# Full example with all corrections
python3 assay_calculator.py --assay ferricyanide --file data.csv \
    --wavelength 455 --path-length 1.0 --dilution 50 \
    --enzyme-volume 0.02 --start 60 --end 300 \
    --starting-conc 0.5

# Print built-in extinction coefficients
python3 assay_calculator.py --assay ferricyanide --file data.csv --list-defaults
```

## CSV format

**Single-wavelength (ferricyanide, DCIP):**
```
time_s,absorbance
0,1.234
15,1.198
30,1.162
```

**Two-wavelength (phosphate — 310 nm analytical + 800 nm background):**
```
time_s,absorbance_310,absorbance_800
0,0.823,0.012
15,0.841,0.011
30,0.860,0.011
```

Header row is auto-detected. Column names are flexible (also accepts `time`, `abs`, `Abs`). No header: columns assumed to be [time_s, absorbance] or [time_s, abs_310, abs_800].

## Assay-specific parameters

### Ferricyanide

| Parameter | Default | Notes |
|-----------|---------|-------|
| `--wavelength` | 455 nm | Use 455 for kinetic; 420 for stopped assay |
| `--path-length` | 1.0 cm | DU800 standard cuvette |
| `--dilution` | 1 | e.g. `50` for 1-in-50 dilution |
| `--starting-conc` | 0.5 mM | For depletion check |
| `--enzyme-volume` | — | mL of enzyme in cuvette; enables per-mL-enzyme rate |

**Built-in extinction coefficients:**
- 455 nm: **329 M⁻¹cm⁻¹**
- 420 nm: **1040 M⁻¹cm⁻¹**

Absorbance **decreases** as ferricyanide is reduced. Tool takes the absolute value of slope.

**Turbidity:** constant absorbance offset at any wavelength; cancels out of the slope. No turbidity correction needed. (If turbidity changes over time due to precipitation, flag the run as unreliable.)

**Depletion warning:** triggered when >50% of starting ferricyanide is consumed in the analysis window. Reaction kinetics become non-linear as substrate depletes. Shorten the window, dilute the enzyme, or use lower ferricyanide if this fires.

**Concentration range note:** At 25 mM ferricyanide (stock), absorbance at 455 nm would be 8.2 AU — off scale. Typical working concentration: 0.5–2 mM.

### Phosphate (Hayano-Fukui protocol)

| Parameter | Default | Notes |
|-----------|---------|-------|
| `--wavelength` | 310 nm | Analytical; 800 nm background built-in |
| `--path-length` | 1.0 cm | DU800 standard cuvette |
| `--start` | 180 s | Skip first 3 min post-phosphate addition |
| `--end` | last data point | Typically 360 s (6 min run) |
| `--extinction` | required for nmol | Provide from standard curve |

**Dilution factor:** 1.25 (1.6 mL sample + 0.4 mL 1M KₓPO₄ pH 7.0 = 2.0 mL total) — applied automatically.

**Background correction:** corrected A = A_310 − A_800. Tool does this automatically from the two-column CSV.

**Read window:** Start at 3–4 min (180–240 s) after phosphate addition. Skip early points where mixing artifacts appear.

To get nmol/min rates you need the molar absorptivity from a phosphate standard curve. Pass it with `--extinction`. Without it, the tool reports ΔAbs/min only.

### DCIP

| Parameter | Default | Notes |
|-----------|---------|-------|
| `--wavelength` | 600 nm | |
| `--path-length` | 1.0 cm | |
| `--extinction` | 16,100 M⁻¹cm⁻¹ | Built-in for reduced DCIP at 600 nm |

Absorbance **decreases** as DCIP is reduced (accepts electrons from substrate oxidation). Measures enolization of 3-ketoglucose or 3-ketosucrose.

## Interpreting the output

```
============================================================
  FERRICYANIDE ASSAY RESULT  (λ = 455 nm)
============================================================
  Window:          60–300 s  (17 points)
  R²:              0.9987
  Slope:           -0.001234 abs/s
  Intercept:       1.2104
  ΔAbs/min:        0.0740
  Rate (rxn):      136.8 nmol/mL_rxn/min
  Rate (enzyme):   6840.0 nmol/mL_enzyme/min
============================================================
```

- **R² < 0.95** → linear fit is poor; check for a lag phase, burst, or depletion
- **Rate (rxn)** → nmol per mL of reaction mixture per minute
- **Rate (enzyme)** → nmol per mL of enzyme fraction per minute (requires `--enzyme-volume`)
- **Depletion warning** → reaction is non-linear in this window; results unreliable

## Common workflows

### Fraction screen (ferricyanide, 455 nm, 1-in-50 dilution)
```bash
for f in fraction_*.csv; do
    echo "=== $f ==="; 
    python3 assay_calculator.py --assay ferricyanide --file "$f" \
        --dilution 50 --enzyme-volume 0.02 --start 30 --end 180;
done
```

### DEAE vs. raw comparison
```bash
python3 assay_calculator.py --assay ferricyanide --file raw.csv \
    --dilution 100 --enzyme-volume 0.01
python3 assay_calculator.py --assay ferricyanide --file deae_peak.csv \
    --dilution 100 --enzyme-volume 0.01
```

### Phosphate assay with known absorptivity
```bash
python3 assay_calculator.py --assay phosphate \
    --file phos_310_800.csv \
    --extinction 3200 \
    --enzyme-volume 0.05
```

## Dependencies

Standard library only: `csv`, `math`, `argparse`, `pathlib`. No pip installs needed. Works with Python 3.7+.

## Known issues and gotchas

- **Column naming:** The CSV reader tries common column name variants but is not magic. If it loads 0 rows, check that the header uses `time_s` (or `time`, `Time`) and `absorbance` (or `abs`, `Abs`).
- **Phosphate one-wavelength data:** If you only have A_310 (no A_800), pass the file as single-column and the tool will use 0 for background. Results will include any turbidity contribution.
- **Non-zero start:** Some DU800 exports begin at t=0 but the reaction isn't started until t=30. Use `--start 30` (or whatever the dead-time is) to exclude the pre-reaction baseline.
- **Negative rates on phosphate:** Phosphate A_310 should rise over time. A negative rate means background drift or wrong file orientation (time reversed). Check raw data.
