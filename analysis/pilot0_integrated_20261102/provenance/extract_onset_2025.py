#!/usr/bin/env python3
"""
Extract the instability onset time from the per-step `*_sigmatop.csv` written by the
instrumented DuMuX model.

Same diagnostic as `code/paper_comparison/extract_onset.py` -- onset is the time at which
the spread of the salt mole fraction along the top boundary is minimal, the Phase I -> II
transition (Bringedal et al. 2025, section 5.3; eqs. 5.1-5.2). The CSV route was validated
against the VTU route to 3e-12 relative, so the two are interchangeable; the CSV is used
here because it survives VTU output throttling and has full temporal resolution.

Note the paper's sigma^top is the unnormalised line integral sqrt(int (x - xbar)^2 dx)
while the CSV holds the population standard deviation. For the uniform lateral spacing
used here they differ by the constant sqrt(domain width) and share the same argmin.
"""
import glob
import json
import os
import sys

import numpy as np


def _boxcar(values, width):
    """Boxcar smoothing that does not invent a minimum at the ends.

    `np.convolve(..., mode="same")` zero-pads, so the last half-window of the smoothed
    curve is scaled down -- by 0.8 and 0.6 for a five-point window -- which plants a
    spurious global minimum at the final sample. It stays hidden while the post-onset
    rise exceeds 1/0.6, which is why the saturated golden set never showed it, and bites
    exactly when a run stops shortly after onset. Replicating the edge keeps a monotone
    tail monotone and leaves every interior value unchanged.
    """
    if values.size < width:
        return values
    half = width // 2
    padded = np.pad(values, half, mode="edge")
    return np.convolve(padded, np.ones(width) / width, mode="valid")


# Set from the measured cross-build variation of the sigma trace: below this the two
# builds' traces are indistinguishable, so samples this close to the minimum cannot be
# ordered by the data.
PLATEAU_TOLERANCE = 1.0e-3
# A plateau spanning this much of the record means the minimum is not localised at all.
PLATEAU_MAX_FRACTION = 0.25
# Growth after the minimum, below which the Phase-I -> II transition is not confirmed.
RISE_MIN = 1.002


def extract_onset(csv_path, *, edge_guard=2):
    """Return a dict with the onset time and enough context to judge it.

    `edge_guard`: how many samples from each end are disqualified as the minimum. A
    minimum sitting at the very first or last sample means the run either started
    supercritical or has not yet reached Phase II, and must not be reported as an onset.
    """
    a = np.loadtxt(csv_path)
    if a.ndim != 2 or a.shape[0] < 5:
        return dict(onset_s=None, status="too few samples", csv=csv_path)
    t_all, mean_all, sigma_all = a[:, 0], a[:, 1], a[:, 2]

    # Same cascade as code/paper_comparison/extract_onset.py, so the saturated and
    # partially saturated legs use an identical diagnostic. A plain global argmin is not
    # enough: when the initial state carries a transient (or the mode is the longest
    # wavelength) sigma can rise first and never drop below its t=0 value, which parks the
    # global minimum at the first sample even though a genuine Phase I -> II dip exists.
    m = t_all > 0
    t, mean, sigma = t_all[m], mean_all[m], sigma_all[m]
    ss = _boxcar(sigma, 5)

    n = len(ss)
    i = int(np.argmin(ss))
    if 0 < i < n - 1 - edge_guard:
        status = "ok (interior minimum)"
    elif i >= n - 1 - edge_guard:
        status = "minimum still at the end (run not past onset yet)"
    else:
        loc = [k for k in range(1, n - 1) if ss[k] < ss[k - 1] and ss[k] <= ss[k + 1]]
        if loc:
            i = loc[0]
            status = "ok (first interior local minimum after initial transient)"
        else:
            i = 0
            status = "immediate onset (monotonic growth, no Phase-I decay)"

    # Depth of the minimum relative to the later rise: a shallow minimum is
    # sampling-sensitive (the saturated suite documents this for high-K cases).
    rise = float(sigma[-1] / sigma[i]) if sigma[i] > 0 else float("nan")

    # An argmin reports a sample index with full apparent precision whether or not the
    # trace resolves it. The plateau is the contiguous run of samples within
    # PLATEAU_TOLERANCE of the minimum: on a well-resolved dip it is a few samples wide,
    # on a flat trace it can span the whole record. Its width is the honest uncertainty of
    # the onset, and the tolerance is set by the measured cross-build variation of the
    # trace itself (2.4e-4 strong to 7.4e-4 weak), not chosen for convenience.
    lo = hi = i
    while lo > 0 and ss[lo - 1] <= (1.0 + PLATEAU_TOLERANCE) * ss[i]:
        lo -= 1
    while hi < n - 1 and ss[hi + 1] <= (1.0 + PLATEAU_TOLERANCE) * ss[i]:
        hi += 1

    # Sub-sample vertex from a quadratic fit of log sigma about the minimum. With the
    # record sampled uniformly this removes the quantisation that otherwise makes every
    # reported ratio a multiple of one sample.
    onset_refined = float(t[i])
    lo_fit, hi_fit = max(0, i - 3), min(n, i + 4)
    if hi_fit - lo_fit >= 3 and np.all(ss[lo_fit:hi_fit] > 0):
        coeff = np.polyfit(t[lo_fit:hi_fit], np.log(ss[lo_fit:hi_fit]), 2)
        if coeff[0] > 0:
            vertex = -coeff[1] / (2.0 * coeff[0])
            if t[lo] <= vertex <= t[hi]:
                onset_refined = float(vertex)

    pre = float(np.polyfit(t[max(0, i - 10):i + 1],
                           np.log(ss[max(0, i - 10):i + 1]), 1)[0]) if i >= 2 else float("nan")
    post = float(np.polyfit(t[i:min(n, i + 11)],
                            np.log(ss[i:min(n, i + 11)]), 1)[0]) if n - i >= 3 else float("nan")

    # A four-way classification the downstream fit can act on, kept beside `status` rather
    # than replacing it so existing consumers are unaffected.
    if status.startswith("minimum still at the end"):
        onset_class = "right_censored"
    elif status.startswith("immediate"):
        onset_class = "immediate_growth"
    elif (hi - lo + 1) > PLATEAU_MAX_FRACTION * n or rise < RISE_MIN:
        onset_class = "weak_interval_onset"
    else:
        onset_class = "resolved_onset"

    return dict(
        # `immediate_growth` has no located minimum: the reported time would be the first
        # sample, an artefact of the sampling grid rather than a measurement. Returning
        # None keeps it from being read as an onset by anything downstream.
        onset_s=(float(t[i]) if status.startswith("ok") else None),
        onset_refined_s=(onset_refined if status.startswith("ok") else None),
        onset_class=onset_class,
        onset_low_s=float(t[lo]),
        onset_high_s=float(t[hi]),
        plateau_width_s=float(t[hi] - t[lo]),
        plateau_samples=int(hi - lo + 1),
        pre_slope_per_s=pre,
        post_slope_per_s=post,
        minimum_position_fraction=float(i / max(1, n - 1)),
        status=status,
        argmin_index=i,
        n_samples=n,
        t_last_s=float(t[-1]),
        sigma_min=float(sigma[i]),
        sigma_first=float(sigma[0]),
        sigma_last=float(sigma[-1]),
        amplification_since_min=rise,
        x_top_mean_at_min=float(mean[i]),
        csv=csv_path,
    )


def extract_run(workdir):
    """Extract from a run directory and write result.json next to the CSV."""
    hits = glob.glob(os.path.join(workdir, "*_sigmatop.csv"))
    if not hits:
        return dict(onset_s=None, status="no *_sigmatop.csv", workdir=workdir)
    res = extract_onset(hits[0])
    meta_path = os.path.join(workdir, "meta.json")
    if os.path.exists(meta_path):
        with open(meta_path) as f:
            res["meta"] = json.load(f)
    with open(os.path.join(workdir, "result.json"), "w") as f:
        json.dump(res, f, indent=2)
    return res


if __name__ == "__main__":
    targets = sys.argv[1:] or sorted(
        d for d in glob.glob(os.path.join(
            os.path.dirname(os.path.abspath(__file__)), "run_*"))
        if os.path.isdir(d))
    for wd in targets:
        r = extract_run(wd)
        o = r.get("onset_s")
        print(f"{os.path.basename(wd):32s} onset = "
              f"{('%9.1f s' % o) if o else '     ---  '}   "
              f"[{r['status']}]  t_last = {r.get('t_last_s', float('nan')):.0f} s")
