# Reproduce WP14 Phase A

From `paper/2`:

```powershell
python scripts/wp14_phase_a.py
python scripts/wp14_phase_a_figures.py
python scripts/verify_wp14_phase_a.py
```

The analysis reads only the frozen WP12/WP13 reduced products. It does not connect to hpcc, open raw production directories or execute DuMuX.

Figure QA uses the locally installed `nature-figure` audit scripts for source validation, PDF glyph size, rendered panel alignment and PDF collision checks. The plotting script is the authoritative source for regenerating TIFF files when needed.
