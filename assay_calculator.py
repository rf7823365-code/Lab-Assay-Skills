#!/usr/bin/env python3
"""
assay_calculator.py — Kinetic enzyme assay rate calculator.

Supports three assay types used in the G3DH / 3-ketoglucose / 3-ketosucrose work:

1. FERRICYANIDE (oxidation assay)
   - Wavelength: 455 nm (kinetic) or 420 nm (stopped endpoint)
   - Extinction coefficient: 329 M⁻¹cm⁻¹ at 455 nm; ~1040 M⁻¹cm⁻¹ at 420 nm
   - Measures: decrease in absorbance as ferricyanide is reduced
   - Turbidity: constant offset, cancels from slope (no correction needed)
   - Depletion warning: if >~50% ferricyanide consumed, linearity is lost

2. PHOSPHATE (Hayano-Fukui protocol)
   - Analytical wavelength: 310 nm; background: 800 nm
   - Background-corrected: A_corrected = A_310 - A_800
   - Read window: 3–4 min after phosphate addition (skip first 3 min)
   - Sample ratio: 1.6 mL sample + 0.4 mL 1M KxPO4 pH 7.0 = 2.0 mL total
   - Dilution factor: 2.0 / 1.6 = 1.25
   - Molar absorptivity: derived from standard curve (provide as argument)

3. DCIP (2,6-dichlorophenolindophenol — electron acceptor assay)
   - Measures enolization of 3-ketoglucose / 3-ketosucrose
   - Wavelength: 600 nm (or user-specified)
   - Extinction coefficient: ~16,100 M⁻¹cm⁻¹ at 600 nm (reduced DCIP form)
   - Decrease in absorbance = electron transfer = substrate oxidation

Common to all:
   - Input: CSV with columns [time_s, absorbance] (and optionally [absorbance_800])
   - Linear regression over specified time window → slope → rate in nmol/min
   - Path-length correction using user-supplied path length (cm)
   - All rates reported in nmol enzyme⁻¹ min⁻¹ if enzyme volume/amount given

Usage:
    python3 assay_calculator.py --assay ferricyanide --file data.csv
    python3 assay_calculator.py --assay phosphate --file data.csv --start 180 --end 240
    python3 assay_calculator.py --assay dcip --file data.csv --wavelength 600
    python3 assay_calculator.py --assay ferricyanide --file data.csv \\
        --wavelength 455 --path-length 1.0 --dilution 50 \\
        --enzyme-volume 0.02 --start 60 --end 300

CSV format (with header row):
    time_s,absorbance
    0,1.234
    15,1.198
    ...

For phosphate (two-wavelength):
    time_s,absorbance_310,absorbance_800
    0,0.823,0.012
    15,0.841,0.011
    ...
"""

import argparse
import csv
import sys
import math
from pathlib import Path

# ---------------------------------------------------------------------------
# Extinction coefficients (M⁻¹ cm⁻¹)
# ---------------------------------------------------------------------------
EXTINCTION = {
    'ferricyanide_455': 329.0,    # mM⁻¹cm⁻¹ reported as 0.329; using M⁻¹cm⁻¹ = 329
    'ferricyanide_420': 1040.0,
    'dcip_600': 16100.0,
}

# Default starting ferricyanide concentration for depletion check
DEFAULT_FERRICYANIDE_CONC_mM = 0.5   # mM; typical working range 0.5–2 mM

# Phosphate assay dilution factor: 1.6 mL sample in 2.0 mL total
PHOSPHATE_DILUTION = 2.0 / 1.6  # = 1.25


# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------

def load_csv(filepath, assay_type):
    """
    Load a CSV file. Returns list of dicts with keys:
      time_s, absorbance (and absorbance_800 for phosphate)
    """
    rows = []
    with open(filepath, newline='', encoding='utf-8') as f:
        # Try to detect header
        sample = f.read(1024)
        f.seek(0)
        has_header = not sample.strip()[0].isdigit()
        reader = csv.DictReader(f) if has_header else csv.reader(f)

        if has_header:
            for row in reader:
                try:
                    entry = {
                        'time_s': float(row.get('time_s', row.get('time', row.get('Time', 0)))),
                        'absorbance': float(row.get('absorbance', row.get('abs', row.get('Abs', 0)))),
                    }
                    if assay_type == 'phosphate':
                        entry['absorbance_800'] = float(
                            row.get('absorbance_800', row.get('abs_800', row.get('Abs_800', 0)))
                        )
                    rows.append(entry)
                except (ValueError, KeyError) as e:
                    print(f"  [skip] row {row}: {e}", file=sys.stderr)
        else:
            for row in reader:
                if len(row) < 2:
                    continue
                try:
                    entry = {'time_s': float(row[0]), 'absorbance': float(row[1])}
                    if assay_type == 'phosphate' and len(row) >= 3:
                        entry['absorbance_800'] = float(row[2])
                    rows.append(entry)
                except ValueError as e:
                    print(f"  [skip] row {row}: {e}", file=sys.stderr)

    return sorted(rows, key=lambda r: r['time_s'])


# ---------------------------------------------------------------------------
# Linear regression (least squares, no scipy dependency)
# ---------------------------------------------------------------------------

def linear_regression(x_vals, y_vals):
    """
    Returns (slope, intercept, r_squared).
    slope is in y-units per x-unit.
    """
    n = len(x_vals)
    if n < 2:
        raise ValueError(f"Need at least 2 data points for regression; got {n}")

    sum_x = sum(x_vals)
    sum_y = sum(y_vals)
    sum_xx = sum(xi ** 2 for xi in x_vals)
    sum_xy = sum(xi * yi for xi, yi in zip(x_vals, y_vals))

    denom = n * sum_xx - sum_x ** 2
    if denom == 0:
        raise ValueError("All time values are identical — cannot compute slope")

    slope = (n * sum_xy - sum_x * sum_y) / denom
    intercept = (sum_y - slope * sum_x) / n

    # R²
    mean_y = sum_y / n
    ss_tot = sum((yi - mean_y) ** 2 for yi in y_vals)
    if ss_tot == 0:
        r_squared = 1.0 if all(abs(yi - mean_y) < 1e-12 for yi in y_vals) else 0.0
    else:
        ss_res = sum((yi - (slope * xi + intercept)) ** 2 for xi, yi in zip(x_vals, y_vals))
        r_squared = 1.0 - ss_res / ss_tot

    return slope, intercept, r_squared


# ---------------------------------------------------------------------------
# Window filtering
# ---------------------------------------------------------------------------

def filter_window(rows, start_s, end_s):
    """Return only rows with start_s <= time_s <= end_s."""
    filtered = [r for r in rows if start_s <= r['time_s'] <= end_s]
    if not filtered:
        times = [r['time_s'] for r in rows]
        raise ValueError(
            f"No data in window [{start_s}–{end_s}s]. "
            f"Data range: [{min(times):.0f}–{max(times):.0f}s]."
        )
    return filtered


# ---------------------------------------------------------------------------
# Assay-specific rate calculations
# ---------------------------------------------------------------------------

def calc_ferricyanide(rows, wavelength, path_length, dilution_factor,
                      start_s, end_s, enzyme_volume_mL,
                      starting_conc_mM=DEFAULT_FERRICYANIDE_CONC_mM,
                      sample_volume_mL=None):
    """
    Ferricyanide kinetic assay.

    Absorbance DECREASES as ferricyanide is reduced.
    slope (abs/s) is negative → rate is positive.

    rate (nmol/min) = |slope| / (epsilon * path_length) * dilution_factor
                       * 1e9 [mol→nmol] / 1000 [mL→L] * sample_volume_mL * 60 [s→min]

    But often reported per mL enzyme:
    rate = |slope| * 60 / (epsilon * path_length) * dilution_factor * 1e6  (nmol/mL/min)
    """
    eps = EXTINCTION.get(f'ferricyanide_{wavelength}')
    if eps is None:
        # Allow user to override by wavelength not in table — warn
        print(f"  [warn] No built-in extinction coefficient for {wavelength} nm. "
              f"Using 455 nm value ({EXTINCTION['ferricyanide_455']}). "
              f"Pass --extinction to override.", file=sys.stderr)
        eps = EXTINCTION['ferricyanide_455']

    window = filter_window(rows, start_s, end_s)
    x = [r['time_s'] for r in window]
    y = [r['absorbance'] for r in window]

    slope, intercept, r2 = linear_regression(x, y)

    # Depletion check
    a_start = intercept + slope * start_s
    a_end = intercept + slope * end_s
    delta_a = a_start - a_end   # positive if absorbance fell (ferricyanide consumed)
    # Convert delta_A to fraction consumed
    # Starting A ≈ eps * path_length * conc_M
    # If starting_conc_mM given, use that; else estimate from first absorbance
    start_conc_M = starting_conc_mM / 1000.0
    start_abs_expected = eps * path_length * start_conc_M
    if start_abs_expected > 0:
        fraction_consumed = delta_a / start_abs_expected
    else:
        fraction_consumed = None

    # Rate in nmol/min per mL reaction (before enzyme normalisation)
    # |slope| abs/s → abs/min; then Beer-Lambert → ΔM/min; then nmol; then dilution
    rate_abs_per_min = abs(slope) * 60.0
    rate_M_per_min = rate_abs_per_min / (eps * path_length)
    rate_nmol_per_mL_per_min = rate_M_per_min * 1e9 / 1000.0  # M→nmol/L→nmol/mL

    # Apply sample dilution
    rate_nmol_per_mL_per_min *= dilution_factor

    # If enzyme volume given, normalise to per-mL-enzyme
    if enzyme_volume_mL and enzyme_volume_mL > 0:
        rate_per_mL_enzyme = rate_nmol_per_mL_per_min / enzyme_volume_mL
    else:
        rate_per_mL_enzyme = None

    return {
        'assay': 'ferricyanide',
        'wavelength_nm': wavelength,
        'window_s': (start_s, end_s),
        'n_points': len(window),
        'slope_abs_per_s': slope,
        'intercept': intercept,
        'r_squared': r2,
        'rate_abs_per_min': rate_abs_per_min,
        'rate_nmol_per_mL_rxn_per_min': rate_nmol_per_mL_per_min,
        'rate_nmol_per_mL_enzyme_per_min': rate_per_mL_enzyme,
        'fraction_ferricyanide_consumed': fraction_consumed,
        'depletion_warning': (fraction_consumed is not None and fraction_consumed > 0.5),
        'path_length_cm': path_length,
        'dilution_factor': dilution_factor,
        'extinction_M_per_cm': eps,
    }


def calc_phosphate(rows, path_length, start_s, end_s,
                   enzyme_volume_mL, molar_absorptivity=None):
    """
    Phosphate assay (Hayano-Fukui protocol).

    Corrected absorbance = A_310 - A_800 (background at 800 nm)
    Read window starts 3 min (180 s) after phosphate addition.

    Increasing absorbance = phosphate release.

    Dilution factor built-in: 1.6 mL sample + 0.4 mL KxPO4 = 2.0 mL; factor = 1.25
    """
    # Apply background correction
    corrected = []
    for r in rows:
        bg = r.get('absorbance_800', 0.0)
        corrected.append({
            'time_s': r['time_s'],
            'absorbance': r['absorbance'] - bg,
        })

    window = filter_window(corrected, start_s, end_s)
    x = [r['time_s'] for r in window]
    y = [r['absorbance'] for r in window]

    slope, intercept, r2 = linear_regression(x, y)

    rate_abs_per_min = slope * 60.0  # positive = increasing A = phosphate release

    if molar_absorptivity:
        eps = molar_absorptivity
        rate_M_per_min = rate_abs_per_min / (eps * path_length)
        rate_nmol_per_mL_per_min = rate_M_per_min * 1e9 / 1000.0 * PHOSPHATE_DILUTION
        rate_per_mL_enzyme = (rate_nmol_per_mL_per_min / enzyme_volume_mL
                               if enzyme_volume_mL and enzyme_volume_mL > 0 else None)
    else:
        rate_nmol_per_mL_per_min = None
        rate_per_mL_enzyme = None
        eps = None

    return {
        'assay': 'phosphate',
        'analytical_wavelength_nm': 310,
        'background_wavelength_nm': 800,
        'window_s': (start_s, end_s),
        'n_points': len(window),
        'slope_corrected_abs_per_s': slope,
        'intercept': intercept,
        'r_squared': r2,
        'rate_corrected_abs_per_min': rate_abs_per_min,
        'rate_nmol_per_mL_rxn_per_min': rate_nmol_per_mL_per_min,
        'rate_nmol_per_mL_enzyme_per_min': rate_per_mL_enzyme,
        'path_length_cm': path_length,
        'dilution_factor': PHOSPHATE_DILUTION,
        'molar_absorptivity_used': eps,
        'note': ('Corrected A = A_310 − A_800. Dilution factor 1.25 applied '
                 '(1.6 mL sample + 0.4 mL KxPO4 = 2.0 mL).'),
    }


def calc_dcip(rows, wavelength, path_length, dilution_factor,
              start_s, end_s, enzyme_volume_mL, extinction=None):
    """
    DCIP assay (2,6-dichlorophenolindophenol).

    Absorbance DECREASES as DCIP is reduced (electron acceptor consumed).
    Rate = substrate oxidation rate.
    """
    eps = extinction or EXTINCTION.get(f'dcip_{wavelength}')
    if eps is None:
        eps = EXTINCTION['dcip_600']
        print(f"  [warn] Using default DCIP extinction at 600 nm ({eps} M⁻¹cm⁻¹). "
              f"Pass --extinction for wavelength {wavelength} nm.", file=sys.stderr)

    window = filter_window(rows, start_s, end_s)
    x = [r['time_s'] for r in window]
    y = [r['absorbance'] for r in window]

    slope, intercept, r2 = linear_regression(x, y)

    rate_abs_per_min = abs(slope) * 60.0
    rate_M_per_min = rate_abs_per_min / (eps * path_length)
    rate_nmol_per_mL_per_min = rate_M_per_min * 1e9 / 1000.0 * dilution_factor

    rate_per_mL_enzyme = (rate_nmol_per_mL_per_min / enzyme_volume_mL
                           if enzyme_volume_mL and enzyme_volume_mL > 0 else None)

    return {
        'assay': 'dcip',
        'wavelength_nm': wavelength,
        'window_s': (start_s, end_s),
        'n_points': len(window),
        'slope_abs_per_s': slope,
        'intercept': intercept,
        'r_squared': r2,
        'rate_abs_per_min': rate_abs_per_min,
        'rate_nmol_per_mL_rxn_per_min': rate_nmol_per_mL_per_min,
        'rate_nmol_per_mL_enzyme_per_min': rate_per_mL_enzyme,
        'path_length_cm': path_length,
        'dilution_factor': dilution_factor,
        'extinction_M_per_cm': eps,
    }


# ---------------------------------------------------------------------------
# Output formatting
# ---------------------------------------------------------------------------

def print_result(result):
    """Pretty-print the result dict."""
    assay = result['assay'].upper()
    wl = result.get('wavelength_nm') or result.get('analytical_wavelength_nm', '?')
    w_start, w_end = result['window_s']
    n = result['n_points']
    r2 = result['r_squared']

    print(f"\n{'='*60}")
    print(f"  {assay} ASSAY RESULT  (λ = {wl} nm)")
    print(f"{'='*60}")
    print(f"  Window:          {w_start}–{w_end} s  ({n} points)")
    print(f"  R²:              {r2:.4f}", end='')
    if r2 < 0.95:
        print("  ← LOW — check linearity", end='')
    print()

    slope_key = 'slope_abs_per_s' if 'slope_abs_per_s' in result else 'slope_corrected_abs_per_s'
    print(f"  Slope:           {result[slope_key]:+.6f} abs/s")
    print(f"  Intercept:       {result['intercept']:.4f}")

    rate_key = 'rate_abs_per_min' if 'rate_abs_per_min' in result else 'rate_corrected_abs_per_min'
    print(f"  ΔAbs/min:        {result[rate_key]:.4f}")

    if result.get('rate_nmol_per_mL_rxn_per_min') is not None:
        print(f"  Rate (rxn):      {result['rate_nmol_per_mL_rxn_per_min']:.2f} nmol/mL_rxn/min")
    else:
        print(f"  Rate (nmol):     [requires --extinction or --molar-absorptivity]")

    if result.get('rate_nmol_per_mL_enzyme_per_min') is not None:
        print(f"  Rate (enzyme):   {result['rate_nmol_per_mL_enzyme_per_min']:.2f} nmol/mL_enzyme/min")

    if result.get('depletion_warning'):
        frac = result.get('fraction_ferricyanide_consumed', 0)
        print(f"\n  ⚠ DEPLETION WARNING: ~{frac*100:.0f}% ferricyanide consumed in window.")
        print(f"    Reaction may not be linear. Consider shorter window or lower enzyme.")

    if result.get('note'):
        print(f"\n  Note: {result['note']}")

    print(f"{'='*60}\n")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description='Kinetic enzyme assay rate calculator (ferricyanide, phosphate, DCIP).',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument('--assay', required=True,
                        choices=['ferricyanide', 'phosphate', 'dcip'],
                        help='Assay type')
    parser.add_argument('--file', required=True,
                        help='CSV data file (time_s, absorbance [, absorbance_800])')
    parser.add_argument('--wavelength', type=int, default=None,
                        help='Measurement wavelength in nm (default: 455 for ferricyanide, '
                             '310 for phosphate, 600 for DCIP)')
    parser.add_argument('--path-length', type=float, default=1.0,
                        help='Cuvette path length in cm (default: 1.0)')
    parser.add_argument('--dilution', type=float, default=1.0,
                        help='Sample dilution factor (e.g. 50 for 1-in-50). '
                             'Phosphate protocol applies its own 1.25 internally.')
    parser.add_argument('--enzyme-volume', type=float, default=None,
                        help='Volume of enzyme added to cuvette in mL (for per-mL-enzyme rates)')
    parser.add_argument('--start', type=float, default=None,
                        help='Start of linear window in seconds (default: auto)')
    parser.add_argument('--end', type=float, default=None,
                        help='End of linear window in seconds (default: auto)')
    parser.add_argument('--extinction', type=float, default=None,
                        help='Molar extinction coefficient (M⁻¹ cm⁻¹) — overrides built-in values')
    parser.add_argument('--starting-conc', type=float, default=DEFAULT_FERRICYANIDE_CONC_mM,
                        help=f'Starting ferricyanide concentration in mM '
                             f'(for depletion check; default {DEFAULT_FERRICYANIDE_CONC_mM})')
    parser.add_argument('--list-defaults', action='store_true',
                        help='Print built-in extinction coefficients and exit')
    args = parser.parse_args()

    if args.list_defaults:
        print("\nBuilt-in extinction coefficients:")
        for k, v in sorted(EXTINCTION.items()):
            print(f"  {k:30s}  {v:8.1f} M⁻¹cm⁻¹")
        print(f"\nPhosphate dilution factor: {PHOSPHATE_DILUTION:.4f} "
              f"(1.6 mL sample + 0.4 mL KxPO4 = 2.0 mL total)")
        return

    # Defaults by assay type
    if args.wavelength is None:
        args.wavelength = {'ferricyanide': 455, 'phosphate': 310, 'dcip': 600}[args.assay]

    print(f"\nLoading: {args.file}")
    rows = load_csv(args.file, args.assay)
    if not rows:
        print("ERROR: No data loaded. Check file format.", file=sys.stderr)
        sys.exit(1)

    times = [r['time_s'] for r in rows]
    print(f"  {len(rows)} rows loaded, time range {min(times):.0f}–{max(times):.0f} s")

    # Default windows
    t_min, t_max = min(times), max(times)
    if args.start is None:
        if args.assay == 'phosphate':
            args.start = 180.0  # skip first 3 min (Hayano protocol)
        else:
            args.start = t_min
    if args.end is None:
        args.end = t_max

    print(f"  Window: {args.start:.0f}–{args.end:.0f} s")

    try:
        if args.assay == 'ferricyanide':
            result = calc_ferricyanide(
                rows,
                wavelength=args.wavelength,
                path_length=args.path_length,
                dilution_factor=args.dilution,
                start_s=args.start,
                end_s=args.end,
                enzyme_volume_mL=args.enzyme_volume,
                starting_conc_mM=args.starting_conc,
            )
            if args.extinction:
                result['extinction_M_per_cm'] = args.extinction
                # Recalculate rate with user extinction
                rate_M = result['rate_abs_per_min'] / (args.extinction * args.path_length)
                result['rate_nmol_per_mL_rxn_per_min'] = (
                    rate_M * 1e9 / 1000.0 * args.dilution
                )
                if args.enzyme_volume:
                    result['rate_nmol_per_mL_enzyme_per_min'] = (
                        result['rate_nmol_per_mL_rxn_per_min'] / args.enzyme_volume
                    )

        elif args.assay == 'phosphate':
            result = calc_phosphate(
                rows,
                path_length=args.path_length,
                start_s=args.start,
                end_s=args.end,
                enzyme_volume_mL=args.enzyme_volume,
                molar_absorptivity=args.extinction,
            )

        elif args.assay == 'dcip':
            result = calc_dcip(
                rows,
                wavelength=args.wavelength,
                path_length=args.path_length,
                dilution_factor=args.dilution,
                start_s=args.start,
                end_s=args.end,
                enzyme_volume_mL=args.enzyme_volume,
                extinction=args.extinction,
            )

        print_result(result)

    except ValueError as e:
        print(f"\nERROR: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == '__main__':
    main()
