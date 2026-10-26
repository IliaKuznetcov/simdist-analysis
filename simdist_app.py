import sys
import json
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

from scipy import interpolate
from scipy.interpolate import PchipInterpolator
from scipy.integrate import cumulative_trapezoid, simpson, trapezoid
from scipy.signal import find_peaks

from PySide6 import QtWidgets, QtCore
from PySide6.QtWidgets import (
    QApplication, QMainWindow, QWidget, QFileDialog, QMessageBox,
    QVBoxLayout, QHBoxLayout, QFormLayout, QTabWidget, QPushButton,
    QLabel, QLineEdit, QTableWidget, QTableWidgetItem, QSplitter,
    QGroupBox, QDoubleSpinBox, QHeaderView, QAbstractItemView,
    QCheckBox, QSpinBox
)

from matplotlib.figure import Figure
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qtagg import NavigationToolbar2QT as NavigationToolbar
from PySide6.QtGui import QIcon

# -----------------------------
# Data structures
# -----------------------------
@dataclass
class AnalysisInputs:
    blank_file: str
    standard_file: str
    sample_file: str
    sample_blank_file: str = ""
    m_std: float = 0.0
    m_slstd: float = 0.0
    m_smp: float = 0.0
    m_slsmp: float = 0.0
    recovery_name: str = ""
    enable_extension: bool = False
    carbon_max: int = 118
    extrapolate_bp: bool = False
    fit_min_carbon: float = 44.0
    remove_sample_solvent_spike: bool = False
    sample_spike_cutoff_time: float = 0.0
    sample_spike_search_max_time: float = 5.0


@dataclass
class AnalysisResult:
    time: np.ndarray
    blank: np.ndarray
    sample_blank: np.ndarray
    standard: np.ndarray
    sample: np.ndarray
    standard_sub: np.ndarray
    sample_sub: np.ndarray
    areas: np.ndarray
    carbon_edges: np.ndarray
    cumulative_edges: np.ndarray
    carbon_centers: np.ndarray
    bp_recovery: np.ndarray
    half_percent_slices: np.ndarray
    rc_total: float
    rf: float
    me: float
    residual: float
    requested_carbon_max: int
    used_carbon_max: int
    measured_carbon_max: int
    used_extension: bool
    bp_extrapolated: bool
    sample_spike_removed: bool
    sample_spike_cutoff_time: float
    assigned_recovery: float
    unassigned_recovery: float


# -----------------------------
# Helpers
# -----------------------------
def show_error(parent, text: str):
    QMessageBox.critical(parent, "Error", text)


def safe_float(item_text: str, default=np.nan):
    try:
        return float(item_text)
    except Exception:
        return default


def load_chromatogram_excel(file_path: str) -> tuple[np.ndarray, np.ndarray]:
    df = pd.read_excel(file_path)
    cols = {c.lower().strip(): c for c in df.columns}

    if "time" not in cols or "value" not in cols:
        raise ValueError(
            f"{Path(file_path).name} must contain columns named 'Time' and 'Value'. "
            f"Found: {list(df.columns)}"
        )

    time = pd.to_numeric(df[cols["time"]], errors="coerce").to_numpy()
    value = pd.to_numeric(df[cols["value"]], errors="coerce").to_numpy()

    mask = np.isfinite(time) & np.isfinite(value)
    time = time[mask]
    value = value[mask]

    if len(time) < 2:
        raise ValueError(f"{Path(file_path).name} does not contain enough valid data.")

    return time, value


def ensure_same_timebase(t_ref: np.ndarray, t_other: np.ndarray, y_other: np.ndarray) -> np.ndarray:
    if len(t_ref) == len(t_other) and np.allclose(t_ref, t_other):
        return y_other

    f = interpolate.interp1d(t_other, y_other, bounds_error=False, fill_value=0.0)
    return f(t_ref)


def subtract_blank(
    ref_time: np.ndarray,
    signal: np.ndarray,
    blank_time: np.ndarray | None,
    blank_signal: np.ndarray | None,
):
    if blank_time is None or blank_signal is None:
        return np.maximum(signal, 0.0)

    blank_interp = interpolate.interp1d(blank_time, blank_signal, bounds_error=False, fill_value=0.0)
    corrected = signal - blank_interp(ref_time)
    return np.maximum(corrected, 0.0)


def _moving_average(y: np.ndarray, window: int = 7) -> np.ndarray:
    y = np.asarray(y, dtype=float)
    if window <= 1 or len(y) < 3:
        return y.copy()

    window = max(3, int(window))
    if window % 2 == 0:
        window += 1
    if len(y) < window:
        window = len(y) if len(y) % 2 == 1 else max(1, len(y) - 1)
    if window <= 1:
        return y.copy()

    kernel = np.ones(window, dtype=float) / float(window)
    return np.convolve(y, kernel, mode="same")


def detect_initial_solvent_spike_cutoff(
    time: np.ndarray,
    signal: np.ndarray,
    search_max_time: float = 5.0,
) -> float:
    """
    Detect the end of a dominant early solvent spike.
    Returns a cutoff time; data at times <= cutoff should be ignored/zeroed.
    """
    time = np.asarray(time, dtype=float)
    signal = np.asarray(signal, dtype=float)

    if len(time) < 5 or len(signal) != len(time):
        return 0.0

    search_max_time = float(max(search_max_time, time[0]))
    mask = (time >= time[0]) & (time <= search_max_time)
    idx = np.flatnonzero(mask)
    if len(idx) < 5:
        return 0.0

    t = time[idx]
    y = np.maximum(signal[idx], 0.0)
    ys = _moving_average(y, window=min(9, max(3, len(y) // 15 * 2 + 1)))

    peak_rel = int(np.argmax(ys))
    peak_height = float(ys[peak_rel])
    if not np.isfinite(peak_height) or peak_height <= 0:
        return 0.0

    baseline = float(np.median(ys[:max(1, min(peak_rel + 1, 5))]))
    threshold = max(baseline + 0.01 * peak_height, 0.02 * peak_height)

    for j in range(peak_rel + 1, len(ys) - 1):
        if ys[j] <= threshold and ys[j] <= ys[j - 1] and ys[j] <= ys[j + 1]:
            return float(t[j])

    for j in range(peak_rel + 1, len(ys)):
        if ys[j] <= threshold:
            return float(t[j])

    return float(t[min(len(t) - 1, peak_rel)])


def remove_initial_region(time: np.ndarray, signal: np.ndarray, cutoff_time: float) -> np.ndarray:
    out = np.asarray(signal, dtype=float).copy()
    if cutoff_time > float(time[0]):
        out[np.asarray(time, dtype=float) <= cutoff_time] = 0.0
    return out


def default_bp_dataframe() -> pd.DataFrame:
    carbon_numbers = np.array(
        [5, 6, 7, 8, 9, 10, 12, 14, 16, 18, 20, 24, 30, 40, 50, 60, 80, 98, 118],
        dtype=float
    )
    boiling_points = np.array(
        [36, 69, 98, 126, 151, 174, 216, 254, 287, 316, 343, 391, 449, 525, 600, 675, 810, 930, 1080],
        dtype=float
    )
    return pd.DataFrame({"CarbonNumber": carbon_numbers, "BoilingPoint": boiling_points})


def load_bp_table(file_path: str | None) -> pd.DataFrame:
    if not file_path:
        return default_bp_dataframe()

    df = pd.read_excel(file_path)
    cols = {c.lower().strip(): c for c in df.columns}

    cn_col = cols.get("cn") or cols.get("carbon number") or cols.get("carbonnumber")
    bp_col = cols.get("bp") or cols.get("boiling point") or cols.get("boilingpoint")

    if not cn_col or not bp_col:
        raise ValueError("Boiling point file must contain CN / CarbonNumber and BP / BoilingPoint columns.")

    out = pd.DataFrame({
        "CarbonNumber": pd.to_numeric(df[cn_col], errors="coerce"),
        "BoilingPoint": pd.to_numeric(df[bp_col], errors="coerce"),
    }).dropna()

    out = out.sort_values("CarbonNumber").reset_index(drop=True)
    if len(out) < 2:
        raise ValueError("Boiling point table must contain at least 2 valid rows.")
    return out


def load_qc_reference(file_path: str) -> pd.DataFrame:
    df = pd.read_excel(file_path)
    cols = {c.lower().strip(): c for c in df.columns}

    rec_col = (
        cols.get("recovery")
        or cols.get("% recovery")
        or cols.get("recovery %")
        or cols.get("percent recovery")
    )
    bp_col = (
        cols.get("bp")
        or cols.get("boiling point")
        or cols.get("boilingpoint")
    )

    if not rec_col or not bp_col:
        raise ValueError(
            "QC file must contain columns like 'Recovery' and 'BP' / 'BoilingPoint'."
        )

    out = pd.DataFrame({
        "Recovery": pd.to_numeric(df[rec_col], errors="coerce"),
        "BP": pd.to_numeric(df[bp_col], errors="coerce"),
    }).dropna()

    out = out.sort_values("Recovery").reset_index(drop=True)

    if len(out) < 2:
        raise ValueError("QC reference file must contain at least 2 valid rows.")

    return out


def default_expected_carbon_sequence() -> np.ndarray:
    poly1 = np.append(np.arange(5, 19), 20).astype(int)
    poly2 = np.arange(22, 121, 2).astype(int)
    return np.append(poly1, poly2).astype(float)


def default_calibration_dataframe() -> pd.DataFrame:
    expected = default_expected_carbon_sequence()
    retention_times = np.linspace(0.6, 30.0, len(expected))
    bp_df = default_bp_dataframe()

    cn = bp_df["CarbonNumber"].to_numpy(dtype=float)
    bp = bp_df["BoilingPoint"].to_numpy(dtype=float)

    interp_bp = interpolate.interp1d(
        cn, bp,
        kind="linear",
        bounds_error=False,
        fill_value="extrapolate",
        assume_sorted=True
    )
    bps = interp_bp(expected)

    return pd.DataFrame({
        "CarbonNumber": expected,
        "RetentionTime": retention_times,
        "BoilingPoint": bps
    })


def build_rt_bp_model(calib_df: pd.DataFrame):
    """
    Build a shape-preserving monotonic RT -> BP model directly from calibration points.
    """
    df = calib_df.dropna(subset=["RetentionTime", "BoilingPoint"]).copy()
    df = df.sort_values("RetentionTime").reset_index(drop=True)

    if len(df) < 2:
        raise ValueError("Calibration table must contain at least 2 valid RT/BP points.")

    rt = df["RetentionTime"].to_numpy(dtype=float)
    bp = df["BoilingPoint"].to_numpy(dtype=float)

    if np.any(np.diff(rt) <= 0):
        raise ValueError("Calibration retention times must be strictly increasing.")
    if np.any(np.diff(bp) <= 0):
        raise ValueError("Calibration boiling points must be strictly increasing.")

    pchip = PchipInterpolator(rt, bp, extrapolate=False)
    deriv = pchip.derivative()

    slope_left = float(deriv(rt[0]))
    slope_right = float(deriv(rt[-1]))

    return {
        "pchip": pchip,
        "rt_min": float(rt[0]),
        "rt_max": float(rt[-1]),
        "bp_min": float(bp[0]),
        "bp_max": float(bp[-1]),
        "slope_left": slope_left,
        "slope_right": slope_right,
    }


def evaluate_rt_bp_model(rt_values: np.ndarray, rt_bp_model: dict, extrapolate_bp: bool = False):
    """
    Evaluate RT -> BP mapping.
    Inside range: PCHIP.
    Outside range:
      - if extrapolate_bp=False: clamp
      - if extrapolate_bp=True: linear extrapolation using end slope
    """
    rt_values = np.asarray(rt_values, dtype=float)

    pchip = rt_bp_model["pchip"]
    rt_min = rt_bp_model["rt_min"]
    rt_max = rt_bp_model["rt_max"]
    bp_min = rt_bp_model["bp_min"]
    bp_max = rt_bp_model["bp_max"]
    slope_left = rt_bp_model["slope_left"]
    slope_right = rt_bp_model["slope_right"]

    bp_out = np.empty_like(rt_values, dtype=float)

    inside = (rt_values >= rt_min) & (rt_values <= rt_max)
    left = rt_values < rt_min
    right = rt_values > rt_max

    if np.any(inside):
        bp_out[inside] = pchip(rt_values[inside])

    if extrapolate_bp:
        if np.any(left):
            bp_out[left] = bp_min + slope_left * (rt_values[left] - rt_min)
        if np.any(right):
            bp_out[right] = bp_max + slope_right * (rt_values[right] - rt_max)
    else:
        if np.any(left):
            bp_out[left] = bp_min
        if np.any(right):
            bp_out[right] = bp_max

    return bp_out


def fit_log_poly_retention_model(
    carbon_numbers: np.ndarray,
    retention_times: np.ndarray,
    fit_min_carbon: float = 44.0
) -> np.ndarray:
    carbon_numbers = np.asarray(carbon_numbers, dtype=float)
    retention_times = np.asarray(retention_times, dtype=float)

    mask = carbon_numbers >= fit_min_carbon
    if np.count_nonzero(mask) < 3:
        n_tail = min(max(3, len(carbon_numbers) // 3), len(carbon_numbers))
        mask = np.zeros(len(carbon_numbers), dtype=bool)
        mask[-n_tail:] = True

    cn_fit = carbon_numbers[mask]
    rt_fit = retention_times[mask]

    if len(cn_fit) < 3:
        raise ValueError("Not enough calibration points to fit high-carbon logarithmic extension.")

    coeffs = np.polyfit(np.log(cn_fit), rt_fit, deg=2)
    return coeffs


def eval_log_poly_retention(carbon_numbers: np.ndarray, coeffs: np.ndarray) -> np.ndarray:
    carbon_numbers = np.asarray(carbon_numbers, dtype=float)
    ln_c = np.log(carbon_numbers)
    return coeffs[0] * ln_c ** 2 + coeffs[1] * ln_c + coeffs[2]


def build_extended_calibration_arrays(
    calib_df: pd.DataFrame,
    requested_carbon_max: int,
    enable_extension: bool,
    extrapolate_bp: bool,
    fit_min_carbon: float = 44.0,
):
    calib_df = calib_df.dropna().sort_values("CarbonNumber").copy()

    if len(calib_df) < 2:
        raise ValueError("Calibration table must have at least 2 valid rows.")

    cn_meas = calib_df["CarbonNumber"].to_numpy(dtype=float)
    rt_meas = calib_df["RetentionTime"].to_numpy(dtype=float)

    if np.any(np.diff(cn_meas) <= 0):
        raise ValueError("Carbon numbers in calibration must be strictly increasing.")
    if np.any(np.diff(rt_meas) <= 0):
        raise ValueError("Retention times in calibration must be strictly increasing.")

    measured_carbon_max = int(np.floor(cn_meas.max()))
    carbon_min = max(5, int(np.floor(cn_meas.min())))

    used_extension = False
    used_carbon_max = min(int(requested_carbon_max), measured_carbon_max)

    cn_all = cn_meas.copy()
    rt_all = rt_meas.copy()

    if enable_extension and int(requested_carbon_max) > measured_carbon_max:
        coeffs = fit_log_poly_retention_model(cn_meas, rt_meas, fit_min_carbon=fit_min_carbon)

        cn_ext = np.arange(measured_carbon_max + 1, int(requested_carbon_max) + 1, dtype=float)
        if len(cn_ext) > 0:
            rt_ext = eval_log_poly_retention(cn_ext, coeffs)

            anchor_cn = float(measured_carbon_max)
            anchor_rt = float(np.interp(anchor_cn, cn_meas, rt_meas))
            anchor_model = float(eval_log_poly_retention(np.array([anchor_cn]), coeffs)[0])
            rt_ext = rt_ext + (anchor_rt - anchor_model)

            prev_rt = anchor_rt
            fixed_rt_ext = []
            min_step = 1e-6
            for val in rt_ext:
                prev_rt = max(float(val), prev_rt + min_step)
                fixed_rt_ext.append(prev_rt)
            rt_ext = np.array(fixed_rt_ext, dtype=float)

            cn_all = np.concatenate([cn_meas, cn_ext])
            rt_all = np.concatenate([rt_meas, rt_ext])
            used_extension = True
            used_carbon_max = int(requested_carbon_max)

    carbon_edges = np.arange(carbon_min, used_carbon_max + 1, dtype=float)
    rt_interp = interpolate.interp1d(
        cn_all,
        rt_all,
        kind="linear",
        bounds_error=False,
        fill_value="extrapolate",
        assume_sorted=True
    )
    retention_time_edges = rt_interp(carbon_edges)

    rt_bp_model = build_rt_bp_model(calib_df)

    bp_all = evaluate_rt_bp_model(
        rt_all,
        rt_bp_model=rt_bp_model,
        extrapolate_bp=extrapolate_bp
    )

    return {
        "carbon_edges": carbon_edges,
        "retention_time_edges": retention_time_edges,
        "cn_all": cn_all,
        "rt_all": rt_all,
        "bp_all": bp_all,
        "rt_bp_model": rt_bp_model,
        "measured_carbon_max": measured_carbon_max,
        "used_carbon_max": used_carbon_max,
        "used_extension": used_extension,
    }


def integral_at(time, signal, query):
    """Exact integral of the linearly interpolated signal, from time[0] to query."""
    t, y, q = np.asarray(time, float), np.asarray(signal, float), np.asarray(query, float)
    if t.ndim != 1 or y.shape != t.shape or len(t) < 2:
        raise ValueError('Time and signal must be matching one-dimensional arrays.')
    if not np.isfinite(t).all() or not np.isfinite(y).all() or np.any(np.diff(t) <= 0):
        raise ValueError('Time and signal must be finite; time must increase strictly.')
    if not np.isfinite(q).all() or np.any(q < t[0]) or np.any(q > t[-1]):
        raise ValueError('Integration limits must lie within the measured time interval.')
    cumulative = cumulative_trapezoid(y, t, initial=0.)
    k = np.clip(np.searchsorted(t, q, side='right') - 1, 0, len(t)-2)
    delta = q - t[k]
    slope = (y[k+1] - y[k]) / (t[k+1] - t[k])
    return cumulative[k] + y[k]*delta + .5*slope*delta**2


def inverse_integral(time, signal, areas):
    """Invert the same integral for a nonnegative signal, including zero plateaus."""
    t, y, a = np.asarray(time, float), np.asarray(signal, float), np.asarray(areas, float)
    total = float(integral_at(t, y, t[-1]))
    tolerance = 64 * np.finfo(float).eps * max(total, 1.)
    if np.any(y < 0) or total <= 0 or not np.isfinite(a).all() or np.any(a < -tolerance) or np.any(a > total + tolerance):
        raise ValueError('Inverse areas must lie in the integral range of a nonnegative signal.')
    a = np.clip(a, 0., total)
    c = cumulative_trapezoid(y, t, initial=0.)
    k = np.clip(np.searchsorted(c, a, side='right') - 1, 0, len(t)-2)
    delta_area = np.maximum(a-c[k], 0.)
    slope = (y[k+1] - y[k]) / (t[k+1]-t[k])
    discriminant = np.maximum(y[k]**2 + 2*slope*delta_area, 0.)
    denominator = y[k] + np.sqrt(discriminant)
    delta_time = np.divide(2*delta_area, denominator, out=np.zeros_like(delta_area), where=denominator > 0)
    out = np.clip(t[k] + delta_time, t[0], t[-1])
    return np.where(a == 0, t[0], np.where(a == total, t[-1], out))


def run_d7169_analysis(inputs: AnalysisInputs, calib_df: pd.DataFrame) -> AnalysisResult:
    t_blank, y_blank = load_chromatogram_excel(inputs.blank_file)
    t_std, y_std = load_chromatogram_excel(inputs.standard_file)
    t_smp, y_smp = load_chromatogram_excel(inputs.sample_file)

    if inputs.sample_blank_file and inputs.sample_blank_file.strip():
        t_smp_blank, y_smp_blank = load_chromatogram_excel(inputs.sample_blank_file)
    else:
        t_smp_blank, y_smp_blank = t_blank, y_blank

    time = t_std
    blank = ensure_same_timebase(time, t_blank, y_blank)
    sample_blank = ensure_same_timebase(time, t_smp_blank, y_smp_blank)
    standard = ensure_same_timebase(time, t_std, y_std)
    sample = ensure_same_timebase(time, t_smp, y_smp)

    blank = np.maximum(blank, 0)
    sample_blank = np.maximum(sample_blank, 0)
    standard = np.maximum(standard, 0)
    sample = np.maximum(sample, 0)

    sample_spike_cutoff_time = 0.0
    if inputs.remove_sample_solvent_spike:
        if inputs.sample_spike_cutoff_time > 0:
            sample_spike_cutoff_time = float(inputs.sample_spike_cutoff_time)
        else:
            sample_spike_cutoff_time = max(
                detect_initial_solvent_spike_cutoff(
                    time,
                    sample_blank,
                    search_max_time=inputs.sample_spike_search_max_time,
                ),
                detect_initial_solvent_spike_cutoff(
                    time,
                    sample,
                    search_max_time=inputs.sample_spike_search_max_time,
                ),
            )

        if sample_spike_cutoff_time > float(time[0]):
            sample_blank = remove_initial_region(time, sample_blank, sample_spike_cutoff_time)
            sample = remove_initial_region(time, sample, sample_spike_cutoff_time)
        else:
            sample_spike_cutoff_time = 0.0

    standard_sub = np.maximum(standard - blank, 0)
    sample_sub = np.maximum(sample - sample_blank, 0)

    a_std = trapezoid(standard_sub, time)
    a_smp = trapezoid(sample_sub, time)

    if a_std <= 0:
        raise ValueError("Standard corrected area is zero or negative.")
    if a_smp <= 0:
        raise ValueError("Sample corrected area is zero or negative.")

    rf = inputs.m_std / (inputs.m_std + inputs.m_slstd) / a_std
    me = a_smp * rf
    full_rc_total = a_smp * rf * (inputs.m_smp + inputs.m_slsmp) / inputs.m_smp * 100.0


    ext = build_extended_calibration_arrays(
        calib_df=calib_df,
        requested_carbon_max=inputs.carbon_max,
        enable_extension=inputs.enable_extension,
        extrapolate_bp=inputs.extrapolate_bp,
        fit_min_carbon=inputs.fit_min_carbon,
    )

    carbon_edges = ext["carbon_edges"]
    retention_time_edges = ext["retention_time_edges"]

    # One piecewise-linear integral defines both carbon bins and the BP curve.
    # Recovery is zero at the first selected carbon boundary, not time[0].
    integral_edges = integral_at(time, sample_sub, retention_time_edges)
    scale = full_rc_total / a_smp
    cumulative_edges = (integral_edges - integral_edges[0]) * scale
    areas = np.diff(cumulative_edges)
    final_recovery = float(cumulative_edges[-1])
    if final_recovery <= 0:
        raise ValueError("No corrected sample signal lies within the selected carbon range.")

    half_percent_slices = np.r_[0.0, np.arange(0.5, final_recovery, 0.5), final_recovery]
    target_areas = integral_edges[0] + half_percent_slices / scale
    sample_times = inverse_integral(
        time, sample_sub, np.clip(target_areas, 0.0, a_smp)
    )
    sample_times = np.clip(
        sample_times, retention_time_edges[0], retention_time_edges[-1]
    )
    # Explicit endpoints retain the selected calibration/extension boundaries,
    # including when a flat cumulative segment makes inversion non-unique.
    sample_times[0] = retention_time_edges[0]
    sample_times[-1] = retention_time_edges[-1]
    bp_recovery = evaluate_rt_bp_model(
        sample_times, rt_bp_model=ext["rt_bp_model"],
        extrapolate_bp=inputs.extrapolate_bp
    )

    carbon_centers = (carbon_edges[:-1] + carbon_edges[1:]) / 2.0
    assigned_recovery = float(np.sum(areas))
    rc_total = assigned_recovery
    residual = 100.0 - rc_total
    unassigned_recovery = 0.0

    return AnalysisResult(
        time=time,
        blank=blank,
        sample_blank=sample_blank,
        standard=standard,
        sample=sample,
        standard_sub=standard_sub,
        sample_sub=sample_sub,
        areas=areas,
        carbon_edges=carbon_edges,
        cumulative_edges=cumulative_edges,
        carbon_centers=carbon_centers,
        bp_recovery=bp_recovery,
        half_percent_slices=half_percent_slices,
        rc_total=rc_total,
        rf=rf,
        me=me,
        residual=residual,
        requested_carbon_max=int(inputs.carbon_max),
        used_carbon_max=int(ext["used_carbon_max"]),
        measured_carbon_max=int(ext["measured_carbon_max"]),
        used_extension=bool(ext["used_extension"]),
        bp_extrapolated=bool(inputs.extrapolate_bp),
        sample_spike_removed=bool(inputs.remove_sample_solvent_spike and sample_spike_cutoff_time > 0),
        sample_spike_cutoff_time=float(sample_spike_cutoff_time),
        assigned_recovery=assigned_recovery,
        unassigned_recovery=unassigned_recovery,
    )


# -----------------------------
# Plot widget with toolbar + hover annotation
# -----------------------------
class MplCanvas(FigureCanvas):
    clicked = QtCore.Signal(object)
    moved = QtCore.Signal(object)
    released = QtCore.Signal(object)
    scrolled = QtCore.Signal(object)

    def __init__(self):
        self.figure = Figure(figsize=(8, 5), tight_layout=True)
        super().__init__(self.figure)

        self.setFocusPolicy(QtCore.Qt.StrongFocus)
        self.setFocus()

        self._cid_click = self.mpl_connect("button_press_event", self._on_click)
        self._cid_move = self.mpl_connect("motion_notify_event", self._on_move)
        self._cid_release = self.mpl_connect("button_release_event", self._on_release)
        self._cid_scroll = self.mpl_connect("scroll_event", self._on_scroll)

    def _on_click(self, event):
        self.clicked.emit(event)

    def _on_move(self, event):
        self.moved.emit(event)

    def _on_release(self, event):
        self.released.emit(event)

    def _on_scroll(self, event):
        self.scrolled.emit(event)


class PlotPanel(QWidget):
    def __init__(self):
        super().__init__()
        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        self.canvas = MplCanvas()
        self.toolbar = NavigationToolbar(self.canvas, self)
        layout.addWidget(self.toolbar)
        layout.addWidget(self.canvas)


# -----------------------------
# Calibration tab
# -----------------------------
class CalibrationTab(QWidget):
    calibration_changed = QtCore.Signal()

    def __init__(self):
        super().__init__()
        self.calibration_path = None
        self._table_updating = False
        # display precision for table only
        self.cn_decimals = 0
        self.rt_decimals = 4
        self.bp_decimals = 2
        self.ph_decimals = 2

        self.bp_df = default_bp_dataframe()

        self.std_time = None
        self.std_signal = None
        self.poly_time = None
        self.poly_signal = None
        self.blank_time = None
        self.blank_signal = None

        self.peak_full_xlim = None
        self.peak_full_ylim = None
        self.peak_current_xlim = None
        self.peak_current_ylim = None

        self.std_corr = None
        self.poly_corr = None

        self.candidate_peaks_df = pd.DataFrame(columns=["Source", "RetentionTime", "PeakHeight"])
        self.selected_peaks_df = pd.DataFrame(columns=["CarbonNumber", "RetentionTime", "BoilingPoint", "PeakHeight", "Source"])

        self.pick_mode = "none"
        self.expected_cn = default_expected_carbon_sequence()

        self._peak_hover_annotation = None

        self._middle_pan_active = False
        self._middle_pan_last_px = None
        self._middle_pan_ax = None

        main_layout = QVBoxLayout(self)

        top_splitter = QSplitter()
        top_splitter.setOrientation(QtCore.Qt.Horizontal)

        left = QWidget()
        left_layout = QVBoxLayout(left)

        source_group = QGroupBox("Calibration Source Files")
        source_form = QFormLayout(source_group)

        self.std_file_edit = QLineEdit()
        self.poly_file_edit = QLineEdit()
        self.blank_file_edit = QLineEdit()
        self.bp_file_edit = QLineEdit()

        self.std_browse_btn = QPushButton("Browse")
        self.poly_browse_btn = QPushButton("Browse")
        self.blank_browse_btn = QPushButton("Browse")
        self.bp_browse_btn = QPushButton("Browse")

        source_form.addRow("Retention standard:", self._file_row(self.std_file_edit, self.std_browse_btn))
        source_form.addRow("Polywax:", self._file_row(self.poly_file_edit, self.poly_browse_btn))
        source_form.addRow("Blank (optional):", self._file_row(self.blank_file_edit, self.blank_browse_btn))
        source_form.addRow("BP table (optional):", self._file_row(self.bp_file_edit, self.bp_browse_btn))

        action_group = QGroupBox("Calibration Actions")
        action_layout = QVBoxLayout(action_group)

        self.btn_auto_detect = QPushButton("Auto Detect")
        self.btn_add_mode = QPushButton("Add Peak Mode")
        self.btn_remove_mode = QPushButton("Remove Peak Mode")
        self.btn_reset_view = QPushButton("Reset Peak View")
        self.btn_clear_picks = QPushButton("Reset Picks")
        self.btn_import_excel = QPushButton("Import Calibration from Excel")
        self.btn_load = QPushButton("Load Calibration JSON")
        self.btn_save = QPushButton("Save Calibration JSON")

        self.chk_prefer_visible_trace = QCheckBox("Prefer nearest visible trace")
        self.chk_prefer_visible_trace.setChecked(True)

        self.chk_middle_mouse_pan = QCheckBox("Enable middle-mouse pan")
        self.chk_middle_mouse_pan.setChecked(True)

        self.mode_label = QLabel("Mode: None")
        self.mode_label.setWordWrap(True)

        action_layout.addWidget(self.btn_auto_detect)
        action_layout.addWidget(self.btn_add_mode)
        action_layout.addWidget(self.btn_remove_mode)
        action_layout.addWidget(self.btn_reset_view)
        action_layout.addWidget(self.btn_clear_picks)
        action_layout.addSpacing(8)
        action_layout.addWidget(self.chk_prefer_visible_trace)
        action_layout.addWidget(self.chk_middle_mouse_pan)
        action_layout.addSpacing(8)
        action_layout.addWidget(self.btn_import_excel)
        action_layout.addWidget(self.btn_load)
        action_layout.addWidget(self.btn_save)
        action_layout.addSpacing(8)
        action_layout.addWidget(self.mode_label)

        info = QLabel(
            "Workflow:\n"
            "1) Load retention standard and Polywax.\n"
            "2) Click Auto Detect.\n"
            "3) Zoom if needed, then use Add Peak / Remove Peak.\n"
            "4) Edit the calibration table if needed.\n"
            "5) Results tab will use this calibration.\n\n"
            "Controls:\n"
            "• Mouse wheel = zoom\n"
            "• Middle mouse drag = pan (optional)\n"
            "• Left click in Add/Remove mode = interactively pick peaks"
        )
        info.setWordWrap(True)

        left_layout.addWidget(source_group)
        left_layout.addWidget(action_group)
        left_layout.addWidget(info)
        left_layout.addStretch(1)

        right = QTabWidget()
        self.peak_panel = PlotPanel()
        self.calib_panel = PlotPanel()
        right.addTab(self.peak_panel, "Peak Picking")
        right.addTab(self.calib_panel, "Calibration Curves")

        top_splitter.addWidget(left)
        top_splitter.addWidget(right)
        top_splitter.setSizes([360, 980])

        table_group = QGroupBox("Calibration Table")
        table_layout = QVBoxLayout(table_group)

        table_button_row = QHBoxLayout()
        self.btn_add_row = QPushButton("Add Row")
        self.btn_remove_row = QPushButton("Remove Selected Row")
        self.btn_reassign = QPushButton("Reassign Carbon Numbers by Time")
        table_button_row.addWidget(self.btn_add_row)
        table_button_row.addWidget(self.btn_remove_row)
        table_button_row.addWidget(self.btn_reassign)
        table_button_row.addStretch(1)

        self.table = QTableWidget(0, 5)
        self.table.setHorizontalHeaderLabels([
            "Carbon Number", "Retention Time (min)", "Boiling Point (°C)", "Peak Height", "Source"
        ])
        self.table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        self.table.setSelectionBehavior(QAbstractItemView.SelectRows)

        table_layout.addLayout(table_button_row)
        table_layout.addWidget(self.table)

        main_layout.addWidget(top_splitter, stretch=3)
        main_layout.addWidget(table_group, stretch=2)

        self.std_browse_btn.clicked.connect(lambda: self.browse_file(self.std_file_edit))
        self.poly_browse_btn.clicked.connect(lambda: self.browse_file(self.poly_file_edit))
        self.blank_browse_btn.clicked.connect(lambda: self.browse_file(self.blank_file_edit))
        self.bp_browse_btn.clicked.connect(lambda: self.browse_file(self.bp_file_edit))

        self.btn_auto_detect.clicked.connect(self.auto_detect_peaks)
        self.btn_add_mode.clicked.connect(lambda: self.set_pick_mode("add"))
        self.btn_remove_mode.clicked.connect(lambda: self.set_pick_mode("remove"))
        self.btn_clear_picks.clicked.connect(self.reset_picks)

        self.btn_load.clicked.connect(self.load_json)
        self.btn_save.clicked.connect(self.save_json)
        self.btn_import_excel.clicked.connect(self.import_excel)

        self.btn_add_row.clicked.connect(self.add_row)
        self.btn_remove_row.clicked.connect(self.remove_selected_rows)
        self.btn_reassign.clicked.connect(self.reassign_carbon_numbers_from_table)

        self.table.itemChanged.connect(self.on_table_item_changed)
        self.peak_panel.canvas.clicked.connect(self.on_peak_plot_clicked)
        self.peak_panel.canvas.moved.connect(self.on_peak_plot_hover)
        self.peak_panel.canvas.released.connect(self.on_peak_plot_released)
        self.peak_panel.canvas.scrolled.connect(self.on_peak_plot_scrolled)

        self.btn_reset_view.clicked.connect(self.reset_peak_view)

        self.set_dataframe(default_calibration_dataframe())

    
    def _format_number(self, value, decimals=4):
        try:
            v = float(value)
        except Exception:
            return ""
        if not np.isfinite(v):
            return ""
        return f"{v:.{decimals}f}"

    def _make_table_item(self, value, decimals=4):
        return QTableWidgetItem(self._format_number(value, decimals))
    
    
    def _file_row(self, line_edit, button):
        w = QWidget()
        l = QHBoxLayout(w)
        l.setContentsMargins(0, 0, 0, 0)
        l.addWidget(line_edit)
        l.addWidget(button)
        return w

    def browse_file(self, line_edit):
        path, _ = QFileDialog.getOpenFileName(self, "Select file", "", "Excel Files (*.xls *.xlsx)")
        if path:
            line_edit.setText(path)

    def peak_toolbar_mode_active(self) -> bool:
        mode = getattr(self.peak_panel.toolbar, "mode", "")
        return bool(mode)

    def store_peak_view(self):
        if self.peak_panel.canvas.figure.axes:
            ax = self.peak_panel.canvas.figure.axes[0]
            self.peak_current_xlim = ax.get_xlim()
            self.peak_current_ylim = ax.get_ylim()

    def on_peak_plot_released(self, event):
        if event is not None and event.button == 2:
            self._middle_pan_active = False
            self._middle_pan_last_px = None
            self._middle_pan_ax = None

            if self.pick_mode in {"add", "remove"}:
                self.peak_panel.canvas.setCursor(QtCore.Qt.CrossCursor)
            else:
                self.peak_panel.canvas.setCursor(QtCore.Qt.ArrowCursor)

        self.store_peak_view()

    def on_peak_plot_scrolled(self, event):
        self.zoom_peak_plot_with_mouse(event)

    def zoom_peak_plot_with_mouse(self, event, scale_step=1.2):
        if event is None or event.inaxes is None:
            return

        ax = event.inaxes
        xdata = event.xdata
        ydata = event.ydata
        if xdata is None or ydata is None:
            return

        cur_xlim = ax.get_xlim()
        cur_ylim = ax.get_ylim()

        x_left, x_right = cur_xlim
        y_bottom, y_top = cur_ylim

        x_range = x_right - x_left
        y_range = y_top - y_bottom
        if x_range <= 0 or y_range <= 0:
            return

        if event.button == "up":
            scale_factor = 1.0 / scale_step
        elif event.button == "down":
            scale_factor = scale_step
        else:
            step = getattr(event, "step", 0)
            scale_factor = 1.0 / scale_step if step > 0 else scale_step

        new_x_range = x_range * scale_factor
        new_y_range = y_range * scale_factor

        x_frac = (xdata - x_left) / x_range
        y_frac = (ydata - y_bottom) / y_range

        new_xlim = (
            xdata - new_x_range * x_frac,
            xdata + new_x_range * (1.0 - x_frac)
        )
        new_ylim = (
            ydata - new_y_range * y_frac,
            ydata + new_y_range * (1.0 - y_frac)
        )

        ax.set_xlim(new_xlim)
        ax.set_ylim(new_ylim)

        self.peak_current_xlim = new_xlim
        self.peak_current_ylim = new_ylim
        self.peak_panel.canvas.draw_idle()

    def reset_peak_view(self):
        if self.peak_full_xlim is None or self.peak_full_ylim is None:
            return
        self.update_peak_plot(force_view=(tuple(self.peak_full_xlim), tuple(self.peak_full_ylim)))

    def set_pick_mode(self, mode: str):
        if self.pick_mode == mode:
            self.pick_mode = "none"
        else:
            self.pick_mode = mode

        self.mode_label.setText(f"Mode: {self.pick_mode.capitalize()}")

        if self.pick_mode in {"add", "remove"}:
            self.peak_panel.canvas.setCursor(QtCore.Qt.CrossCursor)
        else:
            self.peak_panel.canvas.setCursor(QtCore.Qt.ArrowCursor)

        self.update_peak_plot()

    def load_calibration_sources(self):
        std_path = self.std_file_edit.text().strip()
        poly_path = self.poly_file_edit.text().strip()
        blank_path = self.blank_file_edit.text().strip()
        bp_path = self.bp_file_edit.text().strip()

        if not std_path:
            raise ValueError("Retention standard file is required.")
        if not poly_path:
            raise ValueError("Polywax file is required.")

        self.std_time, self.std_signal = load_chromatogram_excel(std_path)
        self.poly_time, self.poly_signal = load_chromatogram_excel(poly_path)

        self.blank_time = None
        self.blank_signal = None
        if blank_path:
            self.blank_time, self.blank_signal = load_chromatogram_excel(blank_path)

        self.bp_df = load_bp_table(bp_path if bp_path else None)

        self.std_corr = subtract_blank(self.std_time, self.std_signal, self.blank_time, self.blank_signal)
        self.poly_corr = subtract_blank(self.poly_time, self.poly_signal, self.blank_time, self.blank_signal)

    def detect_candidate_peaks(self) -> pd.DataFrame:
        if self.std_time is None or self.poly_time is None:
            raise ValueError("Calibration source files are not loaded.")

        mask_std_early = (self.std_time >= 0.0) & (self.std_time <= 2.0)
        peaks_std_early, _ = find_peaks(
            self.std_corr * mask_std_early,
            prominence=2000.0,
            height=1000.0,
            distance=5
        )

        mask_std_mid = (self.std_time > 2.0) & (self.std_time <= 5.0)
        peaks_std_mid, _ = find_peaks(
            self.std_corr * mask_std_mid,
            prominence=500.0,
            height=400.0,
            distance=10
        )

        std_idx = np.unique(np.concatenate([peaks_std_early, peaks_std_mid]))
        std_df = pd.DataFrame({
            "Source": "Retention Standard",
            "RetentionTime": self.std_time[std_idx],
            "PeakHeight": self.std_corr[std_idx]
        })

        last_std_time = std_df["RetentionTime"].max() if not std_df.empty else 5.0

        mask_poly = (
            (self.poly_time > last_std_time + 0.2) &
            (self.poly_time > 5.0) &
            (self.poly_time < 25.4)
        )
        peaks_poly, _ = find_peaks(
            self.poly_corr * mask_poly,
            prominence=200.0,
            distance=10
        )

        poly_df = pd.DataFrame({
            "Source": "Polywax",
            "RetentionTime": self.poly_time[peaks_poly],
            "PeakHeight": self.poly_corr[peaks_poly]
        })

        out = pd.concat([std_df, poly_df], ignore_index=True)
        out = out.sort_values("RetentionTime").reset_index(drop=True)
        return out

    def auto_detect_peaks(self):
        try:
            self.load_calibration_sources()
            self.candidate_peaks_df = self.detect_candidate_peaks()

            if self.candidate_peaks_df.empty:
                raise ValueError("No candidate peaks were detected.")

            self.selected_peaks_df = self.candidate_peaks_df.copy()
            self.reassign_selected_peaks()
            self.update_all_views()
            
            self.selected_peaks_df = self.candidate_peaks_df.copy()
            self.reassign_selected_peaks()
            self.update_all_views()
            self.reset_peak_view()


        except Exception as e:
            show_error(self, str(e))

    def reset_picks(self):
        self.selected_peaks_df = pd.DataFrame(
            columns=["CarbonNumber", "RetentionTime", "BoilingPoint", "PeakHeight", "Source"]
        )
        self.update_all_views()

    def on_peak_plot_clicked(self, event):
        if event is None or event.inaxes is None:
            return

        if event.button == 2 and self.chk_middle_mouse_pan.isChecked():
            self._middle_pan_active = True
            self._middle_pan_last_px = (event.x, event.y)
            self._middle_pan_ax = event.inaxes
            self.peak_panel.canvas.setCursor(QtCore.Qt.ClosedHandCursor)
            return

        if event.xdata is None or event.ydata is None:
            return

        if event.button == 1:
            if self.peak_toolbar_mode_active():
                return

            self.store_peak_view()

            if self.pick_mode == "add":
                self.add_peak_near(event.xdata, event.ydata)
            elif self.pick_mode == "remove":
                self.remove_peak_near(event.xdata, event.ydata)

    def get_peak_axes_ranges(self):
        ax = self.peak_panel.canvas.figure.axes[0] if self.peak_panel.canvas.figure.axes else None
        if ax is None:
            return None, None, None, None, None

        xlim = ax.get_xlim()
        ylim = ax.get_ylim()
        x_range = max(xlim[1] - xlim[0], 1e-12)
        y_range = max(ylim[1] - ylim[0], 1e-12)
        return ax, xlim, ylim, x_range, y_range

    def find_closest_candidate_peak(self, x_click: float, y_click: float):
        if self.candidate_peaks_df.empty:
            return None

        _, _, _, x_range, y_range = self.get_peak_axes_ranges()
        if x_range is None or y_range is None:
            return None

        rt = self.candidate_peaks_df["RetentionTime"].to_numpy(dtype=float)
        ph = self.candidate_peaks_df["PeakHeight"].to_numpy(dtype=float)

        if self.chk_prefer_visible_trace.isChecked():
            dxn = (rt - x_click) / x_range
            dyn = (ph - y_click) / y_range
            dist2 = dxn ** 2 + dyn ** 2
            idx = int(np.argmin(dist2))

            if abs(rt[idx] - x_click) <= max(0.02 * x_range, 0.08) and abs(ph[idx] - y_click) <= max(0.08 * y_range, 3000):
                return idx
            return None

        dx = np.abs(rt - x_click)
        idx = int(np.argmin(dx))
        if dx[idx] <= max(0.02 * x_range, 0.08):
            return idx
        return None

    def add_peak_near(self, x_click: float, y_click: float):
        try:
            if self.std_time is None or self.poly_time is None:
                self.load_calibration_sources()

            if self.candidate_peaks_df.empty:
                self.candidate_peaks_df = self.detect_candidate_peaks()

            if not self.candidate_peaks_df.empty:
                idx = self.find_closest_candidate_peak(x_click, y_click)
                if idx is not None:
                    row = self.candidate_peaks_df.iloc[idx]
                    rt = float(row["RetentionTime"])
                    ph = float(row["PeakHeight"])
                    src = str(row["Source"])
                    self._append_selected_if_new(rt, ph, src)
                    return

            rt, ph, src = self.snap_to_local_max(x_click, y_click, window=0.10)
            if rt is not None:
                self._append_selected_if_new(rt, ph, src)

        except Exception as e:
            show_error(self, str(e))

    def snap_to_local_max(self, x_click: float, y_click: float, window: float = 0.10):
        _, _, _, x_range, y_range = self.get_peak_axes_ranges()
        if x_range is None or y_range is None:
            return None, None, None

        candidates = []

        for time_arr, sig_arr, src in [
            (self.std_time, self.std_corr, "Retention Standard"),
            (self.poly_time, self.poly_corr, "Polywax")
        ]:
            if time_arr is None or sig_arr is None:
                continue

            mask = (time_arr >= x_click - window) & (time_arr <= x_click + window)
            if not np.any(mask):
                continue

            local_t = time_arr[mask]
            local_y = sig_arr[mask]

            peak_idx, _ = find_peaks(local_y, distance=2)
            if len(peak_idx) == 0:
                peak_idx = np.array([int(np.argmax(local_y))])

            for i in peak_idx:
                px = float(local_t[i])
                py = float(local_y[i])

                if self.chk_prefer_visible_trace.isChecked():
                    dxn = (px - x_click) / x_range
                    dyn = (py - y_click) / y_range
                    dist2 = dxn ** 2 + dyn ** 2
                else:
                    dxn = (px - x_click) / x_range
                    dist2 = dxn ** 2

                candidates.append((dist2, px, py, src))

        if not candidates:
            return None, None, None

        candidates.sort(key=lambda t: t[0])
        _, rt, ph, src = candidates[0]
        return rt, ph, src

    def _append_selected_if_new(self, rt: float, ph: float, src: str):
        if not self.selected_peaks_df.empty:
            dx = np.abs(self.selected_peaks_df["RetentionTime"].to_numpy(dtype=float) - rt)
            if np.min(dx) < 1e-3:
                return

        new_row = pd.DataFrame([{
            "CarbonNumber": np.nan,
            "RetentionTime": rt,
            "BoilingPoint": np.nan,
            "PeakHeight": ph,
            "Source": src
        }])

        self.selected_peaks_df = pd.concat([self.selected_peaks_df, new_row], ignore_index=True)
        self.reassign_selected_peaks()
        self.update_all_views()

    def remove_peak_near(self, x_click: float, y_click: float):
        if self.selected_peaks_df.empty:
            return

        _, _, _, x_range, y_range = self.get_peak_axes_ranges()
        if x_range is None or y_range is None:
            return

        rt = self.selected_peaks_df["RetentionTime"].to_numpy(dtype=float)
        ph = self.selected_peaks_df["PeakHeight"].to_numpy(dtype=float)

        if self.chk_prefer_visible_trace.isChecked():
            dxn = (rt - x_click) / x_range
            dyn = (ph - y_click) / y_range
            dist2 = dxn ** 2 + dyn ** 2
            idx = int(np.argmin(dist2))
            accept = (
                abs(rt[idx] - x_click) <= max(0.03 * x_range, 0.10) and
                abs(ph[idx] - y_click) <= max(0.08 * y_range, 3000)
            )
        else:
            dx = np.abs(rt - x_click)
            idx = int(np.argmin(dx))
            accept = dx[idx] <= max(0.03 * x_range, 0.10)

        if accept:
            self.selected_peaks_df = self.selected_peaks_df.drop(
                self.selected_peaks_df.index[idx]
            ).reset_index(drop=True)
            self.reassign_selected_peaks()
            self.update_all_views()

    def _interp_bp_from_cn(self, cn_vals: np.ndarray) -> np.ndarray:
        cn_vals = np.asarray(cn_vals, dtype=float)
        cn = self.bp_df["CarbonNumber"].to_numpy(dtype=float)
        bp = self.bp_df["BoilingPoint"].to_numpy(dtype=float)

        interp_bp = interpolate.interp1d(
            cn, bp,
            kind="linear",
            bounds_error=False,
            fill_value="extrapolate",
            assume_sorted=True
        )
        return interp_bp(cn_vals)

    def reassign_selected_peaks(self):
        if self.selected_peaks_df.empty:
            return

        self.selected_peaks_df = self.selected_peaks_df.sort_values("RetentionTime").reset_index(drop=True)

        n = len(self.selected_peaks_df)
        if n > len(self.expected_cn):
            raise ValueError(
                f"Detected/selected {n} peaks, but only {len(self.expected_cn)} expected carbon-number slots are defined."
            )

        self.selected_peaks_df["CarbonNumber"] = self.expected_cn[:n]
        self.selected_peaks_df["BoilingPoint"] = self._interp_bp_from_cn(
            self.selected_peaks_df["CarbonNumber"].to_numpy(dtype=float)
        )

    def set_dataframe(self, df: pd.DataFrame):
        self.selected_peaks_df = df.copy()
        self.selected_peaks_df = self.selected_peaks_df.sort_values("CarbonNumber").reset_index(drop=True)

        if "PeakHeight" not in self.selected_peaks_df.columns:
            self.selected_peaks_df["PeakHeight"] = np.nan
        if "Source" not in self.selected_peaks_df.columns:
            self.selected_peaks_df["Source"] = "Imported"

        self.update_all_views()

    def get_dataframe(self) -> pd.DataFrame:
        if self.selected_peaks_df.empty:
            return pd.DataFrame(columns=["CarbonNumber", "RetentionTime", "BoilingPoint"])
        cols = ["CarbonNumber", "RetentionTime", "BoilingPoint"]
        return self.selected_peaks_df[cols].dropna().sort_values("CarbonNumber").reset_index(drop=True)

    def get_bp_dataframe(self) -> pd.DataFrame:
        return self.bp_df.copy()

    def update_table_from_selected(self):
        self._table_updating = True
        self.table.blockSignals(True)
        self.table.setRowCount(0)

        for _, row in self.selected_peaks_df.iterrows():
            r = self.table.rowCount()
            self.table.insertRow(r)

            self.table.setItem(r, 0, self._make_table_item(row.get("CarbonNumber", np.nan), self.cn_decimals))
            self.table.setItem(r, 1, self._make_table_item(row.get("RetentionTime", np.nan), self.rt_decimals))
            self.table.setItem(r, 2, self._make_table_item(row.get("BoilingPoint", np.nan), self.bp_decimals))
            self.table.setItem(r, 3, self._make_table_item(row.get("PeakHeight", np.nan), self.ph_decimals))

            src = str(row.get("Source", "")) if pd.notna(row.get("Source", "")) else ""
            self.table.setItem(r, 4, QTableWidgetItem(src))

        self.table.blockSignals(False)
        self._table_updating = False

    def rebuild_selected_from_table(self):
        rows = []
        for r in range(self.table.rowCount()):
            cn_item = self.table.item(r, 0)
            rt_item = self.table.item(r, 1)
            bp_item = self.table.item(r, 2)
            ph_item = self.table.item(r, 3)
            src_item = self.table.item(r, 4)

            cn = safe_float(cn_item.text() if cn_item else "")
            rt = safe_float(rt_item.text() if rt_item else "")
            bp = safe_float(bp_item.text() if bp_item else "")
            ph = safe_float(ph_item.text() if ph_item else "")
            src = src_item.text().strip() if src_item else ""

            if np.isfinite(cn) and np.isfinite(rt):
                if not np.isfinite(bp):
                    bp = float(self._interp_bp_from_cn(np.array([cn], dtype=float))[0])
                rows.append({
                    "CarbonNumber": cn,
                    "RetentionTime": rt,
                    "BoilingPoint": bp,
                    "PeakHeight": ph if np.isfinite(ph) else np.nan,
                    "Source": src or "Manual"
                })

        self.selected_peaks_df = pd.DataFrame(rows)
        if not self.selected_peaks_df.empty:
            self.selected_peaks_df = self.selected_peaks_df.sort_values("CarbonNumber").reset_index(drop=True)

    def on_table_item_changed(self, _item):
        if self._table_updating:
            return
        try:
            self.rebuild_selected_from_table()
            self.update_peak_plot()
            self.update_calibration_plot()
            self.calibration_changed.emit()
        except Exception as e:
            show_error(self, str(e))

    def add_row(self):
        r = self.table.rowCount()
        self.table.insertRow(r)
        self.table.setItem(r, 0, QTableWidgetItem(""))
        self.table.setItem(r, 1, QTableWidgetItem(""))
        self.table.setItem(r, 2, QTableWidgetItem(""))
        self.table.setItem(r, 3, QTableWidgetItem(""))
        self.table.setItem(r, 4, QTableWidgetItem("Manual"))

    def remove_selected_rows(self):
        rows = sorted({idx.row() for idx in self.table.selectionModel().selectedRows()}, reverse=True)
        for row in rows:
            self.table.removeRow(row)
        self.rebuild_selected_from_table()
        self.update_peak_plot()
        self.update_calibration_plot()
        self.calibration_changed.emit()

    def reassign_carbon_numbers_from_table(self):
        try:
            self.rebuild_selected_from_table()
            if self.selected_peaks_df.empty:
                return

            self.selected_peaks_df = self.selected_peaks_df.sort_values("RetentionTime").reset_index(drop=True)
            n = len(self.selected_peaks_df)
            if n > len(self.expected_cn):
                raise ValueError("More selected rows than expected carbon-number sequence supports.")

            self.selected_peaks_df["CarbonNumber"] = self.expected_cn[:n]
            self.selected_peaks_df["BoilingPoint"] = self._interp_bp_from_cn(
                self.selected_peaks_df["CarbonNumber"].to_numpy(dtype=float)
            )
            self.update_all_views()
        except Exception as e:
            show_error(self, str(e))

    def update_all_views(self):
        self.update_table_from_selected()
        self.update_peak_plot()
        self.update_calibration_plot()
        self.calibration_changed.emit()

    def update_peak_plot(self, force_view=None):
        fig = self.peak_panel.canvas.figure

        if force_view is None and fig.axes:
            old_ax = fig.axes[0]
            self.peak_current_xlim = old_ax.get_xlim()
            self.peak_current_ylim = old_ax.get_ylim()

        fig.clear()
        ax = fig.add_subplot(111)

        if self.std_time is not None and self.std_corr is not None:
            ax.plot(self.std_time, self.std_corr, label="Retention Standard")
        if self.poly_time is not None and self.poly_corr is not None:
            ax.plot(self.poly_time, self.poly_corr, label="Polywax")

        if not self.candidate_peaks_df.empty:
            ax.scatter(
                self.candidate_peaks_df["RetentionTime"],
                self.candidate_peaks_df["PeakHeight"],
                marker="o",
                s=20,
                alpha=0.35,
                label="Candidate Peaks"
            )

        if not self.selected_peaks_df.empty:
            ax.scatter(
                self.selected_peaks_df["RetentionTime"],
                self.selected_peaks_df["PeakHeight"],
                marker="x",
                s=50,
                label="Selected Peaks"
            )
            for _, row in self.selected_peaks_df.iterrows():
                try:
                    cn = int(float(row["CarbonNumber"]))
                    lbl = f"C{cn}"
                except Exception:
                    lbl = "?"
                ax.annotate(
                    lbl,
                    (row["RetentionTime"], row["PeakHeight"]),
                    xytext=(0, 18),
                    textcoords="offset points",
                    ha="center",
                    rotation=90
                )

        ax.set_xlabel("Retention time, min")
        ax.set_ylabel("Signal amplitude")
        ax.set_title(f"Peak Picking ({self.pick_mode.capitalize()} mode)")
        ax.grid(True)
        ax.legend(loc="best")

        ax.relim()
        ax.autoscale_view()
        self.peak_full_xlim = tuple(ax.get_xlim())
        self.peak_full_ylim = tuple(ax.get_ylim())

        if force_view is not None:
            xlim, ylim = force_view
            ax.set_xlim(xlim)
            ax.set_ylim(ylim)
            self.peak_current_xlim = tuple(xlim)
            self.peak_current_ylim = tuple(ylim)
        elif self.peak_current_xlim is not None and self.peak_current_ylim is not None:
            try:
                ax.set_xlim(self.peak_current_xlim)
                ax.set_ylim(self.peak_current_ylim)
            except Exception:
                ax.set_xlim(self.peak_full_xlim)
                ax.set_ylim(self.peak_full_ylim)
                self.peak_current_xlim = tuple(self.peak_full_xlim)
                self.peak_current_ylim = tuple(self.peak_full_ylim)
        else:
            self.peak_current_xlim = tuple(self.peak_full_xlim)
            self.peak_current_ylim = tuple(self.peak_full_ylim)

        self._peak_hover_annotation = ax.annotate(
            "",
            xy=(0, 0),
            xytext=(15, 15),
            textcoords="offset points",
            bbox=dict(boxstyle="round", fc="w", alpha=0.85),
            arrowprops=dict(arrowstyle="->")
        )
        self._peak_hover_annotation.set_visible(False)

        self.peak_panel.canvas.draw_idle()

    def on_peak_plot_hover(self, event):
        if self._middle_pan_active and self.chk_middle_mouse_pan.isChecked():
            ax = self._middle_pan_ax
            if ax is None or self._middle_pan_last_px is None:
                return

            if event.x is None or event.y is None:
                return

            last_px = self._middle_pan_last_px
            cur_px = (event.x, event.y)

            if cur_px == last_px:
                return

            inv = ax.transData.inverted()
            last_data = inv.transform(last_px)
            cur_data = inv.transform(cur_px)

            dx = cur_data[0] - last_data[0]
            dy = cur_data[1] - last_data[1]

            cur_xlim = ax.get_xlim()
            cur_ylim = ax.get_ylim()

            new_xlim = (cur_xlim[0] - dx, cur_xlim[1] - dx)
            new_ylim = (cur_ylim[0] - dy, cur_ylim[1] - dy)

            ax.set_xlim(new_xlim)
            ax.set_ylim(new_ylim)

            self.peak_current_xlim = new_xlim
            self.peak_current_ylim = new_ylim
            self._middle_pan_last_px = cur_px

            if self._peak_hover_annotation is not None:
                self._peak_hover_annotation.set_visible(False)

            self.peak_panel.canvas.draw_idle()
            return

        if event.inaxes is None or self.selected_peaks_df.empty or self._peak_hover_annotation is None:
            if self._peak_hover_annotation is not None and self._peak_hover_annotation.get_visible():
                self._peak_hover_annotation.set_visible(False)
                self.peak_panel.canvas.draw_idle()
            return

        x = event.xdata
        y = event.ydata
        if x is None or y is None:
            if self._peak_hover_annotation.get_visible():
                self._peak_hover_annotation.set_visible(False)
                self.peak_panel.canvas.draw_idle()
            return

        rt = self.selected_peaks_df["RetentionTime"].to_numpy(dtype=float)
        ph = self.selected_peaks_df["PeakHeight"].to_numpy(dtype=float)

        if self.chk_prefer_visible_trace.isChecked():
            x_range = event.inaxes.get_xlim()[1] - event.inaxes.get_xlim()[0]
            y_range = event.inaxes.get_ylim()[1] - event.inaxes.get_ylim()[0]
            dxn = (rt - x) / max(x_range, 1e-12)
            dyn = (ph - y) / max(y_range, 1e-12)
            idx = int(np.argmin(dxn ** 2 + dyn ** 2))
            close_x = abs(rt[idx] - x) <= max(0.03 * x_range, 0.08)
            close_y = abs(ph[idx] - y) <= max(0.06 * y_range, 1e-12)
            close_enough = close_x and close_y
        else:
            dx = np.abs(rt - x)
            idx = int(np.argmin(dx))
            x_range = event.inaxes.get_xlim()[1] - event.inaxes.get_xlim()[0]
            close_enough = dx[idx] <= max(0.03 * x_range, 0.08)

        if close_enough:
            row = self.selected_peaks_df.iloc[idx]
            self._peak_hover_annotation.xy = (row["RetentionTime"], row["PeakHeight"])
            self._peak_hover_annotation.set_text(
                f"C{int(row['CarbonNumber'])}\n"
                f"RT = {row['RetentionTime']:.3f} min\n"
                f"Height = {row['PeakHeight']:.3f}\n"
                f"{row['Source']}"
            )
            if not self._peak_hover_annotation.get_visible():
                self._peak_hover_annotation.set_visible(True)
                self.peak_panel.canvas.draw_idle()
            else:
                self.peak_panel.canvas.draw_idle()
        else:
            if self._peak_hover_annotation.get_visible():
                self._peak_hover_annotation.set_visible(False)
                self.peak_panel.canvas.draw_idle()

    def update_calibration_plot(self):
        fig = self.calib_panel.canvas.figure
        fig.clear()

        ax1 = fig.add_subplot(211)
        ax2 = fig.add_subplot(212)

        calib_df = self.get_dataframe()
        if not calib_df.empty and len(calib_df) >= 2:
            ax1.plot(calib_df["RetentionTime"], calib_df["CarbonNumber"], "o-")
            ax1.set_xlabel("Retention time, min")
            ax1.set_ylabel("Carbon number")
            ax1.set_title("Carbon number vs retention time")
            ax1.grid(True)

            ax2.plot(calib_df["RetentionTime"], calib_df["BoilingPoint"], "o-")
            ax2.set_xlabel("Retention time, min")
            ax2.set_ylabel("Boiling point, °C")
            ax2.set_title("Boiling point vs retention time")
            ax2.grid(True)
        else:
            ax1.text(0.5, 0.5, "Not enough calibration points", ha="center", va="center", transform=ax1.transAxes)
            ax2.text(0.5, 0.5, "Not enough calibration points", ha="center", va="center", transform=ax2.transAxes)

        fig.tight_layout()
        self.calib_panel.canvas.draw()

    def save_json(self):
        path, _ = QFileDialog.getSaveFileName(self, "Save Calibration", "calibration.json", "JSON Files (*.json)")
        if not path:
            return

        df = self.selected_peaks_df.copy()
        if df.empty:
            show_error(self, "Calibration table is empty.")
            return

        with open(path, "w", encoding="utf-8") as f:
            json.dump(df.to_dict(orient="records"), f, indent=2)

        self.calibration_path = path
        QMessageBox.information(self, "Saved", f"Calibration saved to:\n{path}")

    def load_json(self):
        path, _ = QFileDialog.getOpenFileName(self, "Load Calibration", "", "JSON Files (*.json)")
        if not path:
            return

        with open(path, "r", encoding="utf-8") as f:
            records = json.load(f)

        df = pd.DataFrame(records)
        required = {"CarbonNumber", "RetentionTime", "BoilingPoint"}
        if not required.issubset(df.columns):
            show_error(self, "JSON must contain CarbonNumber, RetentionTime, and BoilingPoint.")
            return

        if "PeakHeight" not in df.columns:
            df["PeakHeight"] = np.nan
        if "Source" not in df.columns:
            df["Source"] = "Imported"

        self.set_dataframe(df)
        self.calibration_path = path

    def import_excel(self):
        path, _ = QFileDialog.getOpenFileName(self, "Import Calibration Excel", "", "Excel Files (*.xls *.xlsx)")
        if not path:
            return

        df = pd.read_excel(path)
        normalized = {c.lower().strip(): c for c in df.columns}

        possible_cn = normalized.get("cn") or normalized.get("carbonnumber") or normalized.get("carbon number")
        possible_rt = normalized.get("retentiontime") or normalized.get("retention time") or normalized.get("rt")
        possible_bp = normalized.get("bp") or normalized.get("boilingpoint") or normalized.get("boiling point")
        possible_ph = normalized.get("peakheight") or normalized.get("peak height")
        possible_src = normalized.get("source")

        if not (possible_cn and possible_rt and possible_bp):
            show_error(
                self,
                "Excel must contain columns similar to:\nCN / CarbonNumber\nRetentionTime / RT\nBP / BoilingPoint"
            )
            return

        out = pd.DataFrame({
            "CarbonNumber": pd.to_numeric(df[possible_cn], errors="coerce"),
            "RetentionTime": pd.to_numeric(df[possible_rt], errors="coerce"),
            "BoilingPoint": pd.to_numeric(df[possible_bp], errors="coerce"),
            "PeakHeight": pd.to_numeric(df[possible_ph], errors="coerce") if possible_ph else np.nan,
            "Source": df[possible_src].astype(str) if possible_src else "Imported",
        }).dropna(subset=["CarbonNumber", "RetentionTime", "BoilingPoint"])

        self.set_dataframe(out)


# -----------------------------
# Results tab
# -----------------------------
class ResultsTab(QWidget):
    def __init__(self, calibration_tab: CalibrationTab):
        super().__init__()
        self.calibration_tab = calibration_tab
        self.result: AnalysisResult | None = None

        self._dist_hover_annotation = None
        self._bp_hover_annotation = None

        main_layout = QVBoxLayout(self)

        splitter = QSplitter()
        splitter.setOrientation(QtCore.Qt.Horizontal)

        left = QWidget()
        left_layout = QVBoxLayout(left)

        file_group = QGroupBox("Input Files")
        file_form = QFormLayout(file_group)

        self.blank_edit = QLineEdit()
        self.std_edit = QLineEdit()
        self.sample_edit = QLineEdit()
        self.sample_blank_edit = QLineEdit()

        self.blank_btn = QPushButton("Browse")
        self.std_btn = QPushButton("Browse")
        self.sample_btn = QPushButton("Browse")
        self.sample_blank_btn = QPushButton("Browse")

        file_form.addRow("Blank for standard / CS2:", self._file_row(self.blank_edit, self.blank_btn))
        file_form.addRow("Gravimetric standard:", self._file_row(self.std_edit, self.std_btn))
        file_form.addRow("Sample:", self._file_row(self.sample_edit, self.sample_btn))
        file_form.addRow("Sample blank (optional):", self._file_row(self.sample_blank_edit, self.sample_blank_btn))

        mass_group = QGroupBox("Mass Inputs")
        mass_form = QFormLayout(mass_group)

        self.m_std = self._spin(0.0321)
        self.m_slstd = self._spin(1.7723)
        self.m_smp = self._spin(0.0475)
        self.m_slsmp = self._spin(4.2730)

        mass_form.addRow("M_std:", self.m_std)
        mass_form.addRow("M_slstd:", self.m_slstd)
        mass_form.addRow("M_smp:", self.m_smp)
        mass_form.addRow("M_slsmp:", self.m_slsmp)

        meta_group = QGroupBox("Metadata")
        meta_form = QFormLayout(meta_group)
        self.recovery_name = QLineEdit("Sample")
        meta_form.addRow("Recovery name:", self.recovery_name)

        extension_group = QGroupBox("Extension Options")
        extension_form = QFormLayout(extension_group)

        self.chk_enable_extension = QCheckBox("Enable high-carbon logarithmic extension")
        self.chk_enable_extension.setChecked(False)

        self.spin_carbon_max = QSpinBox()
        self.spin_carbon_max.setRange(20, 400)
        self.spin_carbon_max.setValue(118)

        self.chk_extrapolate_bp = QCheckBox("Extrapolate boiling point beyond last RT/BP point")
        self.chk_extrapolate_bp.setChecked(False)

        self.fit_min_carbon_spin = QDoubleSpinBox()
        self.fit_min_carbon_spin.setDecimals(1)
        self.fit_min_carbon_spin.setRange(5.0, 300.0)
        self.fit_min_carbon_spin.setValue(44.0)

        extension_form.addRow("", self.chk_enable_extension)
        extension_form.addRow("Extend to carbon number:", self.spin_carbon_max)
        extension_form.addRow("Log-fit minimum carbon:", self.fit_min_carbon_spin)
        extension_form.addRow("", self.chk_extrapolate_bp)

        spike_group = QGroupBox("Sample Solvent Spike Handling")
        spike_form = QFormLayout(spike_group)

        self.chk_remove_sample_spike = QCheckBox("Remove initial sample solvent spike")
        self.chk_remove_sample_spike.setChecked(False)

        self.sample_spike_cutoff_spin = QDoubleSpinBox()
        self.sample_spike_cutoff_spin.setDecimals(3)
        self.sample_spike_cutoff_spin.setRange(0.0, 20.0)
        self.sample_spike_cutoff_spin.setValue(0.0)
        self.sample_spike_cutoff_spin.setSingleStep(0.05)
        self.sample_spike_cutoff_spin.setToolTip("0 = auto-detect using the early sample/sample-blank spike")

        self.sample_spike_search_max_spin = QDoubleSpinBox()
        self.sample_spike_search_max_spin.setDecimals(3)
        self.sample_spike_search_max_spin.setRange(0.1, 20.0)
        self.sample_spike_search_max_spin.setValue(5.0)
        self.sample_spike_search_max_spin.setSingleStep(0.1)

        spike_form.addRow("", self.chk_remove_sample_spike)
        spike_form.addRow("Cutoff time, min (0 = auto):", self.sample_spike_cutoff_spin)
        spike_form.addRow("Auto-detect search window, min:", self.sample_spike_search_max_spin)

        qc_group = QGroupBox("Quality Control (optional)")
        qc_form = QFormLayout(qc_group)
        self.qc_ref_edit = QLineEdit()
        self.qc_ref_btn = QPushButton("Browse")
        qc_form.addRow("QC reference (Recovery/BP Excel):", self._file_row(self.qc_ref_edit, self.qc_ref_btn))

        self.qc_summary_label = QLabel("No QC comparison run.")
        self.qc_summary_label.setWordWrap(True)

        self.run_btn = QPushButton("Run Analysis")
        self.export_btn = QPushButton("Export Results CSV")
        self.export_btn.setEnabled(False)

        self.summary_label = QLabel("No analysis run yet.")
        self.summary_label.setWordWrap(True)

        left_layout.addWidget(file_group)
        left_layout.addWidget(mass_group)
        left_layout.addWidget(meta_group)
        left_layout.addWidget(extension_group)
        left_layout.addWidget(spike_group)
        left_layout.addWidget(qc_group)
        left_layout.addWidget(self.run_btn)
        left_layout.addWidget(self.export_btn)
        left_layout.addWidget(self.summary_label)
        left_layout.addWidget(self.qc_summary_label)
        left_layout.addStretch(1)

        right = QTabWidget()

        self.dist_panel = PlotPanel()
        self.bp_panel = PlotPanel()
        self.raw_panel = PlotPanel()
        self.qc_panel = PlotPanel()

        right.addTab(self.dist_panel, "Carbon Distribution")
        right.addTab(self.bp_panel, "Boiling Point vs Recovery")
        right.addTab(self.raw_panel, "Corrected Chromatograms")
        right.addTab(self.qc_panel, "QC Comparison")

        splitter.addWidget(left)
        splitter.addWidget(right)
        splitter.setSizes([390, 980])

        main_layout.addWidget(splitter)

        self.blank_btn.clicked.connect(lambda: self.browse_file(self.blank_edit))
        self.std_btn.clicked.connect(lambda: self.browse_file(self.std_edit))
        self.sample_btn.clicked.connect(lambda: self.browse_file(self.sample_edit))
        self.sample_blank_btn.clicked.connect(lambda: self.browse_file(self.sample_blank_edit))
        self.qc_ref_btn.clicked.connect(lambda: self.browse_file(self.qc_ref_edit))
        self.run_btn.clicked.connect(self.run_analysis)
        self.export_btn.clicked.connect(self.export_results)

        self.dist_panel.canvas.moved.connect(self.on_dist_hover)
        self.bp_panel.canvas.moved.connect(self.on_bp_hover)

    def _file_row(self, line_edit, button):
        w = QWidget()
        l = QHBoxLayout(w)
        l.setContentsMargins(0, 0, 0, 0)
        l.addWidget(line_edit)
        l.addWidget(button)
        return w

    def _spin(self, value):
        box = QDoubleSpinBox()
        box.setDecimals(6)
        box.setRange(0.0, 1e9)
        box.setValue(value)
        return box

    def _position_annotation_horizontally(self, annotation, ax, canvas):
        renderer = canvas.get_renderer()
        if renderer is None:
            return

        annotation.set_ha("left")
        annotation.xyann = (15, 15)
        bbox = annotation.get_window_extent(renderer=renderer)
        ax_bbox = ax.get_window_extent(renderer=renderer)

        if bbox.x1 > ax_bbox.x1:
            annotation.set_ha("right")
            annotation.xyann = (-15, 15)
            bbox = annotation.get_window_extent(renderer=renderer)

        if bbox.x0 < ax_bbox.x0:
            annotation.set_ha("left")
            annotation.xyann = (15, 15)

    def browse_file(self, line_edit):
        path, _ = QFileDialog.getOpenFileName(self, "Select chromatogram file", "", "Excel Files (*.xls *.xlsx)")
        if path:
            line_edit.setText(path)

    def run_analysis(self):
        try:
            calib_df = self.calibration_tab.get_dataframe()

            if calib_df.empty:
                raise ValueError("Calibration table is empty.")

            inputs = AnalysisInputs(
                blank_file=self.blank_edit.text().strip(),
                standard_file=self.std_edit.text().strip(),
                sample_file=self.sample_edit.text().strip(),
                sample_blank_file=self.sample_blank_edit.text().strip(),
                m_std=self.m_std.value(),
                m_slstd=self.m_slstd.value(),
                m_smp=self.m_smp.value(),
                m_slsmp=self.m_slsmp.value(),
                recovery_name=self.recovery_name.text().strip() or "Sample",
                enable_extension=self.chk_enable_extension.isChecked(),
                carbon_max=self.spin_carbon_max.value(),
                extrapolate_bp=self.chk_extrapolate_bp.isChecked(),
                fit_min_carbon=self.fit_min_carbon_spin.value(),
                remove_sample_solvent_spike=self.chk_remove_sample_spike.isChecked(),
                sample_spike_cutoff_time=self.sample_spike_cutoff_spin.value(),
                sample_spike_search_max_time=self.sample_spike_search_max_spin.value(),
            )

            for label, path in [
                ("Blank for standard / CS2", inputs.blank_file),
                ("Standard", inputs.standard_file),
                ("Sample", inputs.sample_file),
            ]:
                if not path:
                    raise ValueError(f"{label} file is not selected.")

            self.result = run_d7169_analysis(inputs, calib_df)
            self.plot_distribution(inputs.recovery_name)
            self.plot_bp_curve(inputs.recovery_name)
            self.plot_raw(inputs.recovery_name)

            qc_text = ""
            qc_path = self.qc_ref_edit.text().strip()
            if qc_path:
                qc_text = self.plot_qc_comparison(inputs.recovery_name, qc_path)
            else:
                self.clear_qc_plot()
                self.qc_summary_label.setText("No QC comparison run.")

            extension_text = "Yes" if self.result.used_extension else "No"
            bp_extrap_text = "Yes" if self.result.bp_extrapolated else "No"

            if (
                not self.result.used_extension
                and inputs.enable_extension
                and self.result.used_carbon_max < inputs.carbon_max
            ):
                extension_note = (
                    f"\nRequested extension to C{inputs.carbon_max}, but only used up to C{self.result.used_carbon_max}."
                )
            else:
                extension_note = ""

            sample_blank_text = (
                Path(inputs.sample_blank_file).name
                if inputs.sample_blank_file.strip()
                else "Same as standard / CS2 blank"
            )

            spike_text = (
                f"Removed to {self.result.sample_spike_cutoff_time:.3f} min"
                if self.result.sample_spike_removed
                else "No"
            )

            self.summary_label.setText(
                f"Analysis complete.\n"
                f"RF = {self.result.rf:.6e}\n"
                f"ME = {self.result.me:.6f}\n"
                f"Recovery RC (selected carbon range) = {self.result.rc_total:.3f} %\n"
                f"Not accounted for within selected range = {self.result.residual:.3f} %\n"
                f"Measured calibration max = C{self.result.measured_carbon_max}\n"
                f"Used carbon max = C{self.result.used_carbon_max}\n"
                f"High-carbon extension used = {extension_text}\n"
                f"Boiling-point extrapolation = {bp_extrap_text}\n"
                f"Sample blank used = {sample_blank_text}\n"
                f"Initial sample solvent spike removal = {spike_text}"
                f"{extension_note}"
            )

            if qc_text:
                self.qc_summary_label.setText(qc_text)
            self.export_btn.setEnabled(True)

        except Exception as e:
            show_error(self, str(e))

    def plot_distribution(self, title):
        if self.result is None:
            return

        fig = self.dist_panel.canvas.figure
        fig.clear()

        ax1 = fig.add_subplot(111)
        ax2 = ax1.twinx()

        widths = np.diff(self.result.carbon_edges)
        right_edges = self.result.carbon_edges[1:]
        cumulative_recovery = self.result.cumulative_edges[1:]

        bar_centers = self.result.carbon_centers.astype(float)
        bar_heights = self.result.areas.astype(float)
        bar_widths = widths.astype(float)

        ax1.bar(
            bar_centers,
            bar_heights,
            width=bar_widths,
            align="center",
            edgecolor="black"
        )
        ax1.set_xlabel("Carbon number")
        ax1.set_ylabel("Fractional recovery, %")

        tick_positions = list(self.result.carbon_edges[::max(1, int(np.ceil(len(self.result.carbon_edges) / 20)))])
        tick_labels = [f"C{int(x)}" for x in tick_positions]
        ax1.set_xticks(tick_positions)
        ax1.set_xticklabels(tick_labels, rotation=90)

        ax2.plot(
            right_edges,
            cumulative_recovery,
            marker="o",
            linewidth=2.0,
            color="tab:red",
            label="Cumulative recovery"
        )
        ax2.set_ylabel("Cumulative recovery, %", color="tab:red")
        ax2.tick_params(axis="y", colors="tab:red")

        ax2.axhline(
            self.result.rc_total,
            linestyle=":",
            linewidth=1.2,
            color="gray",
            label="Recovery RC"
        )

        ax1.set_title(title)

        if self.result.used_extension:
            ax1.axvline(
                self.result.measured_carbon_max,
                linestyle="--",
                linewidth=1.2,
                color="gray"
            )

        self._dist_hover_annotation = ax1.annotate(
            "",
            xy=(0, 0),
            xytext=(15, 15),
            textcoords="offset points",
            bbox=dict(boxstyle="round", fc="w", alpha=0.85),
            arrowprops=dict(arrowstyle="->")
        )
        self._dist_hover_annotation.set_visible(False)

        fig.tight_layout()
        self.dist_panel.canvas.draw()

    def plot_bp_curve(self, title):
        if self.result is None:
            return

        fig = self.bp_panel.canvas.figure
        fig.clear()
        ax = fig.add_subplot(111)

        ax.plot(self.result.bp_recovery, self.result.half_percent_slices, marker="o")
        ax.set_xlabel("Boiling point, °C")
        ax.set_ylabel("Recovery, %")
        ax.set_title(f"{title} - Boiling Point vs Recovery")
        ax.grid(True)

        self._bp_hover_annotation = ax.annotate(
            "",
            xy=(0, 0),
            xytext=(15, 15),
            textcoords="offset points",
            bbox=dict(boxstyle="round", fc="w", alpha=0.85),
            arrowprops=dict(arrowstyle="->")
        )
        self._bp_hover_annotation.set_visible(False)

        fig.tight_layout()
        self.bp_panel.canvas.draw()

    def plot_raw(self, title):
        if self.result is None:
            return

        fig = self.raw_panel.canvas.figure
        fig.clear()
        ax = fig.add_subplot(111)

        same_blank = np.allclose(self.result.blank, self.result.sample_blank)

        ax.plot(self.result.time, self.result.blank, label="Standard Blank (CS2)", linestyle="--")
        if not same_blank:
            ax.plot(self.result.time, self.result.sample_blank, label="Sample Blank", linestyle=":")

        ax.plot(self.result.time, self.result.standard_sub, label="Standard - Standard Blank")
        ax.plot(self.result.time, self.result.sample_sub, label="Sample - Sample Blank")

        if self.result.sample_spike_removed and self.result.sample_spike_cutoff_time > 0:
            ax.axvline(
                self.result.sample_spike_cutoff_time,
                linestyle="-.",
                linewidth=1.2,
                color="gray",
                label=f"Sample spike cutoff = {self.result.sample_spike_cutoff_time:.3f} min"
            )

        ax.set_xlabel("Retention time, min")
        ax.set_ylabel("Signal amplitude, pA")
        ax.set_title(f"{title} - Corrected Chromatograms")
        ax.legend()
        ax.grid(True)

        fig.tight_layout()
        self.raw_panel.canvas.draw()

    def clear_qc_plot(self):
        fig = self.qc_panel.canvas.figure
        fig.clear()
        ax = fig.add_subplot(111)
        ax.text(0.5, 0.5, "No QC reference loaded.", ha="center", va="center", transform=ax.transAxes)
        ax.set_xlabel("Boiling point, °C")
        ax.set_ylabel("Recovery, %")
        ax.grid(False)
        self.qc_panel.canvas.draw()

    def plot_qc_comparison(self, title: str, qc_file_path: str) -> str:
        if self.result is None:
            self.clear_qc_plot()
            return ""

        ref_df = load_qc_reference(qc_file_path)

        exp_rec = self.result.half_percent_slices.astype(float)
        exp_bp = self.result.bp_recovery.astype(float)

        ref_rec = ref_df["Recovery"].to_numpy(dtype=float)
        ref_bp = ref_df["BP"].to_numpy(dtype=float)

        exp_interp = interpolate.interp1d(
            exp_rec,
            exp_bp,
            bounds_error=False,
            fill_value=np.nan,
            assume_sorted=True
        )
        exp_bp_at_ref = exp_interp(ref_rec)

        valid = np.isfinite(exp_bp_at_ref) & np.isfinite(ref_bp)

        if np.any(valid):
            diff = exp_bp_at_ref[valid] - ref_bp[valid]
            mae = float(np.mean(np.abs(diff)))
            rmse = float(np.sqrt(np.mean(diff ** 2)))
            max_abs = float(np.max(np.abs(diff)))
        else:
            mae = np.nan
            rmse = np.nan
            max_abs = np.nan

        fig = self.qc_panel.canvas.figure
        fig.clear()

        ax = fig.add_subplot(111)
        ax.plot(ref_bp, ref_rec, "o", label="QC Reference")
        ax.plot(exp_bp, exp_rec, "-o", label="Experiment")
        ax.set_xlabel("Boiling point, °C")
        ax.set_ylabel("Recovery, %")
        ax.set_title(f"{title} - QC Boiling Point vs Recovery")
        ax.grid(True)
        ax.legend()

        fig.tight_layout()
        self.qc_panel.canvas.draw()

        if np.any(valid):
            return (
                f"QC comparison complete.\n"
                f"MAE = {mae:.2f} °C\n"
                f"RMSE = {rmse:.2f} °C\n"
                f"Max |ΔBP| = {max_abs:.2f} °C"
            )

        return "QC comparison could not be computed because the recovery ranges do not overlap enough."

    def on_dist_hover(self, event):
        if self.result is None or self._dist_hover_annotation is None or event.inaxes is None:
            if self._dist_hover_annotation is not None and self._dist_hover_annotation.get_visible():
                self._dist_hover_annotation.set_visible(False)
                self.dist_panel.canvas.draw_idle()
            return

        x = event.xdata
        y = event.ydata
        if x is None or y is None:
            if self._dist_hover_annotation.get_visible():
                self._dist_hover_annotation.set_visible(False)
                self.dist_panel.canvas.draw_idle()
            return

        centers = self.result.carbon_centers.astype(float)
        areas = self.result.areas.astype(float)
        edges = self.result.carbon_edges.astype(float)
        right_edges = edges[1:]
        cumulative = self.result.cumulative_edges[1:].astype(float)

        x_range = event.inaxes.get_xlim()[1] - event.inaxes.get_xlim()[0]

        candidates = []
        for idx, cx in enumerate(centers):
            candidates.append((abs(cx - x), "assigned", idx))

        if self.result.unassigned_recovery > 1e-9:
            extra_center = float(self.result.used_carbon_max + 1.0)
            candidates.append((abs(extra_center - x), "unassigned", 0))

        best_dx, kind, idx = min(candidates, key=lambda t: t[0])

        if best_dx > max(0.03 * x_range, 0.8):
            if self._dist_hover_annotation.get_visible():
                self._dist_hover_annotation.set_visible(False)
                self.dist_panel.canvas.draw_idle()
            return

        if kind == "assigned":
            left_cn = int(edges[idx])
            right_cn = int(edges[idx + 1])
            self._dist_hover_annotation.xy = (centers[idx], areas[idx])
            self._dist_hover_annotation.set_text(
                f"Bin: C{left_cn}–C{right_cn}\n"
                f"Fractional assigned = {areas[idx]:.3f} %\n"
                f"Cumulative = {cumulative[idx]:.3f} %"
            )
        else:
            extra_center = float(self.result.used_carbon_max + 1.0)
            self._dist_hover_annotation.xy = (extra_center, self.result.unassigned_recovery)
            self._dist_hover_annotation.set_text(
                f"Bin: >C{self.result.used_carbon_max}\n"
                f"Unassigned recovery = {self.result.unassigned_recovery:.3f} %\n"
                f"Total recovery RC = {self.result.rc_total:.3f} %"
            )

        self._dist_hover_annotation.set_visible(True)
        self._position_annotation_horizontally(
            self._dist_hover_annotation,
            event.inaxes,
            self.dist_panel.canvas,
        )
        self.dist_panel.canvas.draw_idle()

    def on_bp_hover(self, event):
        if self.result is None or self._bp_hover_annotation is None or event.inaxes is None:
            if self._bp_hover_annotation is not None and self._bp_hover_annotation.get_visible():
                self._bp_hover_annotation.set_visible(False)
                self.bp_panel.canvas.draw_idle()
            return

        x = event.xdata
        y = event.ydata
        if x is None or y is None:
            if self._bp_hover_annotation.get_visible():
                self._bp_hover_annotation.set_visible(False)
                self.bp_panel.canvas.draw_idle()
            return

        bp = self.result.bp_recovery
        rec = self.result.half_percent_slices

        idx = int(np.argmin((bp - x) ** 2 + (rec - y) ** 2))

        x_range = event.inaxes.get_xlim()[1] - event.inaxes.get_xlim()[0]
        y_range = event.inaxes.get_ylim()[1] - event.inaxes.get_ylim()[0]

        close_x = abs(bp[idx] - x) <= max(0.03 * x_range, 3.0)
        close_y = abs(rec[idx] - y) <= max(0.05 * y_range, 0.6)

        if close_x and close_y:
            self._bp_hover_annotation.xy = (bp[idx], rec[idx])
            self._bp_hover_annotation.set_text(
                f"Boiling point = {bp[idx]:.2f} °C\n"
                f"Recovery = {rec[idx]:.2f} %"
            )
            if not self._bp_hover_annotation.get_visible():
                self._bp_hover_annotation.set_visible(True)
                self.bp_panel.canvas.draw_idle()
            else:
                self.bp_panel.canvas.draw_idle()
        else:
            if self._bp_hover_annotation.get_visible():
                self._bp_hover_annotation.set_visible(False)
                self.bp_panel.canvas.draw_idle()

    def export_results(self):
        if self.result is None:
            return

        path, _ = QFileDialog.getSaveFileName(self, "Export Results", "d7169_results.csv", "CSV Files (*.csv)")
        if not path:
            return

        right_edges = self.result.carbon_edges[1:]
        n_main = min(len(right_edges), len(self.result.carbon_centers), len(self.result.areas), len(self.result.cumulative_edges[1:]))

        df_main = pd.DataFrame({
            "BinLabel": [f"C{int(self.result.carbon_edges[i])}-C{int(self.result.carbon_edges[i + 1])}" for i in range(n_main)],
            "CarbonBinRightEdge": right_edges[:n_main],
            "CarbonBinCenter": self.result.carbon_centers[:n_main],
            "FractionalRecoveryPercent": self.result.areas[:n_main],
            "CumulativeRecoveryPercent": self.result.cumulative_edges[1:][:n_main]
        })

        df_bp = pd.DataFrame({
            "RecoveryPercent_BPcurve": self.result.half_percent_slices,
            "BoilingPointC_BPcurve": self.result.bp_recovery
        })

        base = Path(path)
        df_main.to_csv(base, index=False)
        df_bp.to_csv(base.with_name(base.stem + "_bp_curve.csv"), index=False)

        QMessageBox.information(
            self,
            "Exported",
            f"Saved:\n{base}\n{base.with_name(base.stem + '_bp_curve.csv')}"
        )


# -----------------------------
# Main window
# -----------------------------
class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("D7169 Recovery Calculator")
        self.resize(1550, 920)

        tabs = QTabWidget()
        self.calibration_tab = CalibrationTab()
        self.results_tab = ResultsTab(self.calibration_tab)

        tabs.addTab(self.calibration_tab, "Calibration")
        tabs.addTab(self.results_tab, "Results")

        self.setCentralWidget(tabs)


def main():
    app = QApplication(sys.argv)
    app.setWindowIcon(QIcon(str(Path(__file__).resolve().parent / "SIMDIST.png")))
    win = MainWindow()
    win.show()
    sys.exit(app.exec())


if __name__ == "__main__":
    main()