
### code/lsa_unsaturated.py @ 158395f3d25a29680dea5bcd0d3a24e52a28b4cd
```text
571: def critical_rayleigh(t, n, m, L, Swr, delta, pbot, beta, *, nx=None, nx_max=30,
572:                       ell_hat=3.0, **kw):
573:     """R_c(t) = min over the quantized wavenumbers, and the minimizing n_x.
574:
575:     With `nx` given (int or sequence) the minimization is restricted to it, which is how
576:     the paper's Table 5.9 column is defined.
577:     """
578:     nx_arr = np.arange(1, nx_max + 1) if nx is None else np.atleast_1d(np.asarray(nx, int))
579:     a = nx_arr * 2.0 * math.pi / ell_hat
580:     R = eig_val_problem(a, t, n, m, L, Swr, delta, pbot, beta, **kw)
581:     if not np.any(np.isfinite(R)):
582:         return np.nan, -1
583:     j = int(np.nanargmin(R))
584:     return float(R[j]), int(nx_arr[j])
585:
586:
587: def onset_time_hat(R_E_target, n, m, L, Swr, delta, pbot, beta, *, nx=None, nx_max=30,
588:                    ell_hat=3.0, bracket=(1e-6, 0.1), tol=1e-6, itmax=60, **kw):
589:     """Dimensionless onset time: the first t_hat where R_c(t_hat) crosses R_E_target.
590:
591:     R_c decreases with time over the range of interest, so this bisects on log t. The
592:     bracket is deliberately wider than the reference implementation's [1e-4, 1e-2], which
593:     excludes Table 5.6 Case 6 (t_hat = 0.0254).
594:
595:     Returns (t_hat, status). status is "ok", or "no crossing" when R_c never reaches the
596:     target inside the bracket -- the censored / never-unstable case, which the paper
597:     discusses in section 3.5 and which must be propagated rather than clipped.
598:     """
599:     lo, hi = bracket
600:
601:     def f(t):
602:         Rc, _ = critical_rayleigh(t, n, m, L, Swr, delta, pbot, beta,
603:                                   nx=nx, nx_max=nx_max, ell_hat=ell_hat, **kw)
604:         return Rc - R_E_target
605:
...
```

### code/parameter_database/solvers.py @ 158395f3d25a29680dea5bcd0d3a24e52a28b4cd
```text
286: def run_lsa(point: DimensionlessPoint, realization: PhysicalRealization, *, nx: int | None = None,
287:             nx_max: int = 30, n_chebyshev: int = 96, t_hat_min: float = 1e-7,
288:             t_hat_max: float = 20.0, scan_points: int = 33,
289:             tolerance: float = 1e-6,
290:             chi_coeff_mode: str = "reference") -> dict[str, Any]:
291:     """Find the first LSA crossing on a declared dimensionless time horizon.
292:
293:     A log-time scan precedes bisection.  This is intentionally more expensive than
294:     checking two bracket endpoints: the unsaturated neutral curve need not be monotone,
295:     and equal endpoint signs do not prove that there was no crossing in between.
296:     """
297:     point.validate()
298:     if not (0.0 < t_hat_min < t_hat_max):
299:         raise ValueError("LSA horizon must satisfy 0 < t_hat_min < t_hat_max")
300:     if scan_points < 3:
301:         raise ValueError("LSA scan needs at least three log-time points")
302:
303:     times = np.geomspace(t_hat_min, t_hat_max, scan_points)
304:     critical = np.empty(scan_points)
305:     modes = np.empty(scan_points, dtype=int)
306:     for idx, t_hat in enumerate(times):
307:         critical[idx], modes[idx] = unsaturated_lsa.critical_rayleigh(
308:             float(t_hat), point.vg_n, point.vg_m, point.mualem_l,
309:             point.residual_saturation, 1.0 / point.gamma, point.p_hat_b, point.beta,
310:             nx=nx, nx_max=nx_max, ell_hat=point.width_over_height, N=n_chebyshev,
311:             chi_coeff_mode=chi_coeff_mode,
312:         )
313:
314:     finite = np.isfinite(critical) & (critical > 0.0)
315:     ratios = np.full(scan_points, np.nan)
316:     ratios[finite] = np.log10(critical[finite] / point.r_e)
317:     result: dict[str, Any] = {
318:         "onset_status": "solver_failure",
319:         "diagnostic_status": "solver_failure_during_time_scan",
320:         "t_hat_onset_lsa": None,
321:         "t_onset_lsa_s": None,
322:         "critical_nx": None,
323:         "critical_rayleigh": None,
324:         "reference_time_s": realization.reference_time_s,
325:         "requested_nx": nx,
326:         "nx_max": nx_max,
327:         "n_chebyshev": n_chebyshev,
328:         "search_t_hat_min": t_hat_min,
329:         "search_t_hat_max": t_hat_max,
330:         "search_scan_points": scan_points,
331:         "bisection_log_tolerance": tolerance,
332:         "chi_coeff_mode": chi_coeff_mode,
333:         "crossing_count": 0,
334:         "bracket_sufficient": False,
335:         "stability_margin_log10": (float(np.nanmin(ratios))
...
402:         if log_hi - log_lo < tolerance:
403:             break
404:     t_hat = math.exp(0.5 * (log_lo + log_hi))
405:     critical_nx = None
406:     critical_r = None
407:     critical_r, critical_nx = unsaturated_lsa.critical_rayleigh(
408:         t_hat, point.vg_n, point.vg_m, point.mualem_l,
409:         point.residual_saturation, 1.0 / point.gamma, point.p_hat_b, point.beta,
410:         nx=nx, nx_max=nx_max, ell_hat=point.width_over_height, N=n_chebyshev,
411:         chi_coeff_mode=chi_coeff_mode,
412:     )
413:     bracket_sufficient = 1 < crossing_index < scan_points - 1
414:     diagnostic = "interior_first_crossing"
415:     if not bracket_sufficient:
416:         diagnostic = ("crossing_near_lower_search_bound" if crossing_index <= 1
417:                       else "crossing_near_upper_search_bound")
418:     if len(transitions) > 1:
419:         diagnostic += ";multiple_crossings_in_scan"
420:     result.update({
421:         "onset_status": "onset_found",
422:         "diagnostic_status": diagnostic,
423:         "t_hat_onset_lsa": float(t_hat),
424:         "t_onset_lsa_s": float(t_hat * realization.reference_time_s),
425:         "critical_nx": critical_nx,
426:         "critical_rayleigh": critical_r,
427:         "onset_bracket_t_hat": [lo, hi],
428:         "eigenprofile_t_hat": float(t_hat),
429:         "bracket_sufficient": bracket_sufficient,
430:     })
431:     return result
432:
433:
434: def run_precipitation_clock(
...
730:
731: def select_dominant_mode(
732:     point: DimensionlessPoint, lsa_result: dict[str, Any], *,
733:     n_chebyshev: int = 96, chi_coeff_mode: str = "paper_eq367",
734:     candidates: Sequence[int] = DOMINANT_MODE_CANDIDATES, sharpness_min: float = SHARPNESS_MIN,
735:     sharpness_floor: float = SIGMA_SHARP_MIN, t_hat_stop: float | None = None,
736:     t_hat_min: float = 1e-7, t_hat_max: float = 0.1,
737: ) -> dict[str, Any]:
738:     """Choose the lateral mode to seed. One rule for every case, no class distinction.
739:
740:     For each candidate n_x >= n_c (in increasing order, starting at the critical mode
741:     itself) the mode's own neutral time t_c(n_x) is found by
742:     `lsa_unsaturated.onset_time_hat` restricted to that n_x, and its sharpness
743:     sigma(n_x, 2 t_c(n_x)) * t_c(n_x) computed. A candidate whose crossing falls after
744:     `t_hat_stop` -- the competing-risk stop time min(precipitation, horizon) -- is not
745:     an event and is rejected: without that ceiling the search happily returns a mode
746:     crossing at 711 t_c, 268 times later than the case's own precipitation.
747:
748:     Of the surviving candidates the *first* to reach `sharpness_min` is returned
749:     (`selection = "preferred"`): above the threshold the earliest event is the one the
750:     sigma_top diagnostic sees. If none reaches it, the candidate with the *highest*
751:     sharpness is returned provided it clears `sharpness_floor`
752:     (`selection = "fallback_best"`): below the preferred threshold there is nothing to
753:     be early about, so the best available seed wins. Only when even that fails is the
754:     status "no_sharp_mode".
755:
756:     When the critical mode already reaches `sharpness_min` -- the "sharp" class -- the
757:     rule returns n_c and `t_ratio_to_critical` is exactly 1, i.e. the pre-existing
758:     protocol is a special case of this one, not a separate branch. The fallback exists
759:     for the same reason: on a 180-case sample four grey-band cases (critical sharpness
760:     0.05-0.08, no higher mode crossing at all) would otherwise be dropped although they
761:     resolve in DuMuX today.
762:
763:     The DuMuX onset is always compared against the returned `t_hat_onset_dominant`, the
764:     seeded mode's own crossing time, never against t_c(n_c). Which mode is seeded sets
765:     that reference time but not the LSA-to-DuMuX discrepancy: on the paper's four cases
766:     the ratio is 1.03-1.11 whichever of two modes is seeded, and across 109 HF cases
767:     the rank correlation between critical-mode sharpness and |ln(ratio) - median| is
768:     -0.04. Record the seeded n_x and its sharpness with every row so that this stays
769:     checkable.
770:     """
771:     if lsa_result.get("onset_status") != "onset_found":
772:         raise ValueError("dominant-mode selection needs an onset_found LSA result")
773:     n_c = int(lsa_result["critical_nx"])
774:     t_c = float(lsa_result["t_hat_onset_lsa"])
775:     args = (point.vg_n, point.vg_m, point.mualem_l, point.residual_saturation,
776:             1.0 / point.gamma, point.p_hat_b, point.beta)
777:     ell = point.width_over_height
778:     tried: list[dict[str, Any]] = []
779:     viable: list[tuple[int, float, float]] = []          # (nx, t_nx, sharpness)
780:
781:     def _hit(nx: int, t_nx: float, sharpness: float, selection: str) -> dict[str, Any]:
782:         return {"status": "ok", "selection": selection, "nx_dominant": nx,
783:                 "t_hat_onset_dominant": float(t_nx),
784:                 "t_ratio_to_critical": float(t_nx / t_c), "sharpness": float(sharpness),
785:                 "cells_x": cells_x_for_mode(nx, ell), "critical_nx": n_c,
786:                 "sharpness_min": sharpness_min, "sharpness_floor": sharpness_floor,
787:                 "t_hat_stop": t_hat_stop, "tried": tried}
788:
789:     for nx in sorted(set(int(v) for v in candidates) | {n_c}):
790:         if nx < n_c:
791:             continue
792:         if nx == n_c:
793:             # The crossing the LSA result already found; re-solving for it would only
794:             # add root-finder noise to a time we hold exactly.
795:             t_nx, status = t_c, "ok"
796:         else:
797:             t_nx, status = unsaturated_lsa.onset_time_hat(
798:                 point.r_e, *args, nx=nx, ell_hat=ell, bracket=(t_hat_min, t_hat_max),
799:                 N=n_chebyshev, chi_coeff_mode=chi_coeff_mode)
800:         if status != "ok":
801:             tried.append({"nx": nx, "status": status})
802:             continue
803:         if t_hat_stop is not None and t_nx > t_hat_stop:
804:             tried.append({"nx": nx, "t_hat_onset": float(t_nx),
805:                           "status": "crossing_after_stop"})
806:             continue
807:         sigma = unsaturated_lsa.growth_rate(
808:             nx * 2.0 * math.pi / ell, 2.0 * t_nx, point.r_e, *args,
809:             N=n_chebyshev, chi_coeff_mode=chi_coeff_mode)[0]
810:         sharpness = float(sigma * t_nx)
811:         tried.append({"nx": nx, "t_hat_onset": float(t_nx), "sharpness": sharpness})
812:         if sharpness >= sharpness_min:
813:             return _hit(nx, t_nx, sharpness, "preferred")
814:         viable.append((nx, float(t_nx), sharpness))
815:
816:     if viable:
817:         nx, t_nx, sharpness = max(viable, key=lambda v: v[2])
818:         if sharpness >= sharpness_floor:
819:             return _hit(nx, t_nx, sharpness, "fallback_best")
820:     return {"status": "no_sharp_mode", "nx_dominant": None, "critical_nx": n_c,
821:             "sharpness_min": sharpness_min, "sharpness_floor": sharpness_floor,
822:             "t_hat_stop": t_hat_stop, "tried": tried}
823:
...
890: def run_lsa_resolution_audit(
891:     point: DimensionlessPoint, realization: PhysicalRealization, *,
892:     n_chebyshev_values: Sequence[int], relative_tolerance: float = 0.05,
893:     **kwargs: Any,
894: ) -> dict[str, Any]:
895:     """Escalate an ordered resolution ladder until two adjacent levels agree."""
896:     levels = [int(value) for value in n_chebyshev_values]
897:     if len(levels) < 2 or levels != sorted(set(levels)) or levels[0] < 8:
898:         raise ValueError("resolution ladder requires at least two unique increasing N >= 8")
899:
900:     history: list[dict[str, Any]] = []
901:     previous: dict[str, Any] | None = None
902:     current: dict[str, Any] | None = None
903:     relative_error = None
904:     accepted_pair = None
905:     for level in levels:
906:         current = run_lsa(point, realization, n_chebyshev=level, **kwargs)
907:         history.append({
908:             "n_chebyshev": level, "onset_status": current["onset_status"],
909:             "t_hat_onset_lsa": current.get("t_hat_onset_lsa"),
910:             "critical_nx": current.get("critical_nx"),
911:             "stability_margin_log10": current.get("stability_margin_log10"),
912:         })
913:         if previous is None:
914:             previous = current
915:             continue
916:         same_status = previous["onset_status"] == current["onset_status"]
917:         same_mode = previous.get("critical_nx") == current.get("critical_nx")
918:         relative_error = None
919:         if (same_status and current["onset_status"] == "onset_found" and
920:                 previous.get("t_hat_onset_lsa") is not None):
921:             relative_error = abs(
922:                 current["t_hat_onset_lsa"] / previous["t_hat_onset_lsa"] - 1.0)
923:             converged = relative_error <= relative_tolerance and same_mode
924:         elif same_status and current["onset_status"] == "right_censored_at_horizon":
925:             lower_margin = previous.get("stability_margin_log10")
926:             upper_margin = current.get("stability_margin_log10")
927:             relative_error = (abs(upper_margin - lower_margin)
928:                               if lower_margin is not None and upper_margin is not None else None)
...
950:     if accepted_pair is None:
951:         current["diagnostic_status"] += ";resolution_unresolved_ladder_exhausted"
952:     return current
953:
954:
955: def write_lsa_eigenprofile(path: str | Path, point: DimensionlessPoint,
956:                            realization: PhysicalRealization, lsa_result: dict[str, Any],
957:                            *, n_chebyshev: int = 96, n_output: int = 400) -> Path:
958:     if lsa_result.get("onset_status") not in {"onset_found", "right_censored_at_horizon"}:
959:         raise ValueError("an LSA event or finite-horizon classification is required")
960:     if not lsa_result.get("critical_nx") or lsa_result.get("eigenprofile_t_hat") is None:
961:         raise ValueError("critical_nx and eigenprofile_t_hat are required for an eigenprofile")
962:     nx = int(lsa_result["critical_nx"])
963:     profile_t_hat = float(lsa_result["eigenprofile_t_hat"])
964:     a_hat = nx * 2.0 * math.pi / point.width_over_height
965:     critical_r, z_cheb, chi = unsaturated_lsa.critical_eigenmode(
966:         a_hat, profile_t_hat, point.vg_n, point.vg_m,
967:         point.mualem_l, point.residual_saturation, 1.0 / point.gamma,
968:         point.p_hat_b, point.beta, N=n_chebyshev,
969:         chi_coeff_mode=lsa_result.get("chi_coeff_mode", "reference"),
970:     )
971:     y = np.linspace(0.0, realization.column_height_m, n_output)
972:     amplitude = np.interp(y, realization.column_height_m * z_cheb, chi)
973:     peak = np.max(np.abs(amplitude))
974:     if peak > 0.0:
975:         amplitude /= peak
...
```

### code/parameter_database/pipeline.py @ 158395f3d25a29680dea5bcd0d3a24e52a28b4cd
```text
245:         # retention curve. The published work's two legs agreed on the first and differed
246:         # on the second, which is what put its DuMuX/LSA onset ratio between 0.668 and
247:         # 1.465.
248:         raise ValueError(
249:             f"DuMuX and LSA disagree on the relative permeability closure "
250:             f"(max relative error {krw_error:.3g}): {audit}")
251:
252:     resolved_t_hat_max = spec.censor_time_hat if lsa_t_hat_max is None else lsa_t_hat_max
253:     lsa_kwargs = dict(nx=nx, nx_max=nx_max, t_hat_min=lsa_t_hat_min,
254:                       t_hat_max=resolved_t_hat_max, scan_points=lsa_scan_points,
255:                       tolerance=lsa_tolerance, chi_coeff_mode=lsa_chi_coeff_mode)
256:     if lsa_resolution_ladder is not None:
257:         lsa = run_lsa_resolution_audit(
258:             spec.point, realization, n_chebyshev_values=lsa_resolution_ladder,
259:             relative_tolerance=lsa_resolution_tolerance, **lsa_kwargs)
260:     elif lsa_check_chebyshev is not None:
261:         lsa = run_lsa_resolution_audit(
262:             spec.point, realization,
263:             n_chebyshev_values=(lsa_check_chebyshev, n_chebyshev),
264:             relative_tolerance=lsa_resolution_tolerance, **lsa_kwargs)
265:     else:
266:         lsa = run_lsa(spec.point, realization, n_chebyshev=n_chebyshev, **lsa_kwargs)
267:     precipitation = run_precipitation_clock(
268:         spec.point, realization,
269:         salt_solubility_mole_fraction=spec.salt_solubility_mole_fraction,
270:         n_chebyshev=int(lsa["n_chebyshev"]), t_hat_min=lsa_t_hat_min,
271:         t_hat_max=resolved_t_hat_max, scan_points=lsa_scan_points,
272:         tolerance=lsa_tolerance,
273:     )
274:     lsa.update(precipitation)
275:     # A censored LSA leg is not yet a censored observation: R_c may settle above R_E, in
...
280:                     n_chebyshev=int(lsa["n_chebyshev"]),
281:                     chi_coeff_mode=lsa_chi_coeff_mode)
282:                 if lsa["onset_status"] == "right_censored_at_horizon" else None)
283:     if eventual is not None:
284:         lsa.update(eventual)
285:     lsa.update(classify_competing_events(lsa, precipitation, eventual))
286:     case_numerics = numerics.with_lsa_result(lsa)
287:     workdir = Path(workdir)
288:     workdir.mkdir(parents=True, exist_ok=True)
289:     safe_key = "".join(ch if ch.isalnum() or ch in "-_" else "_" for ch in spec.case_key)
290:     name = f"{safe_key}--{policy}"
291:     lsa_protocol = "adaptive_log_time_competing_risk_v2"
292:     dumux_protocol = None
293:     profile_path = None
294:     input_path = None
295:     if (lsa["first_event_type"] != "precipitation_first" and
296:             lsa["onset_status"] in {"onset_found", "right_censored_at_horizon"}):
297:         if dumux_perturbation == "Random":
298:             # The published protocol: Gaussian noise of standard deviation `amplitude`
299:             # on the top row of vertices, every lateral mode seeded alike. The seed is
300:             # part of the protocol string so that realizations stay distinct rows.
301:             dumux_protocol = f"top_row_random_v1:seed={int(dumux_random_seed)}"
302:             input_text = render_dumux_input(
303:                 realization, name=name, numerics=case_numerics,
304:                 constitutive=constitutive, perturbation="Random",
305:                 critical_nx=int(lsa["critical_nx"]), random_seed=int(dumux_random_seed),
306:             )
307:         else:
308:             dumux_protocol = (
309:                 "lsa_eigenprofile_at_onset_v1" if lsa["onset_status"] == "onset_found"
310:                 else "lsa_least_stable_eigenprofile_for_horizon_classification_v1"
311:             )
312:             # 400 samples over the column is 0.5 mm, coarser than the sub-millimetre
313:             # eigenmodes of thin boundary layers; the count is a protocol choice.
314:             profile_path = write_lsa_eigenprofile(
315:                 workdir / f"{name}_eigenprofile.dat", spec.point, realization, lsa,
316:                 n_chebyshev=int(lsa["n_chebyshev"]), n_output=int(eigenprofile_samples))
317:             input_text = render_dumux_input(
318:                 realization, name=name, numerics=case_numerics,
319:                 constitutive=constitutive, perturbation="Eigenprofile",
320:                 critical_nx=int(lsa["critical_nx"]), eigenprofile_path=profile_path,
321:             )
322:         input_path = workdir / f"{name}.input"
323:         input_path.write_text(input_text, encoding="utf-8")
324:     store_lsa_result(con, spec, realization, lsa, audit, solver_version=solver_version,
325:                      constitutive=constitutive, perturbation_protocol=lsa_protocol)
326:
327:     manifest: dict[str, Any] = {
328:         "case": {"case_id": spec.case_id, "case_key": spec.case_key},
329:         "policy": policy, "dimensionless_point": asdict(spec.point),
330:         "realization": asdict(realization), "constitutive_audit": audit,
331:         "dumux_numerics": {
332:             "requested": asdict(numerics),
333:             "case_resolved": asdict(case_numerics),
334:             "physical": case_numerics.resolved(realization),
335:         },
336:         "lsa": lsa, "input": str(input_path) if input_path else None,
337:         "eigenprofile": str(profile_path) if profile_path else None,
338:         "dumux": None,
339:     }
340:     if execute_dumux:
341:         if not lsa.get("bracket_sufficient"):
342:             raise ValueError(
343:                 "DuMuX execution requires an LSA crossing away from the search-horizon "
344:                 "endpoints; enlarge the LSA time horizon"
345:             )
...
424:     nx_by_case: dict[int, int] | None = None,
425:     eigenprofile_samples: int = 400,
426: ) -> dict[str, Any]:
427:     con.row_factory = sqlite3.Row
428:     rows = con.execute(
429:         """SELECT cdp.selection_rank,pc.case_id,pc.case_key
430:         FROM campaign_design_points cdp JOIN parameter_cases pc USING(case_id)
431:         WHERE cdp.design_name=? ORDER BY cdp.selection_rank""", (design_name,)).fetchall()
432:     if not rows:
433:         raise ValueError(f"campaign design not found or empty: {design_name}")
434:     if limit is not None:
435:         rows = rows[:limit]
436:     root = Path(workdir)
437:     results: list[dict[str, Any]] = []
438:     for row in rows:
439:         case_dir = root / f"{row['selection_rank']:04d}--case-{row['case_id']}--{policy}"
440:         try:
441:             manifest = prepare_and_run_case(
442:                 con, row["case_id"], policy=policy, workdir=case_dir,
443:                 numerics=numerics, constitutive=constitutive,
444:                 execute_dumux=execute_dumux, executable=executable,
445:                 nx=(nx_by_case or {}).get(int(row["case_id"]), nx),
446:                 nx_max=nx_max, n_chebyshev=n_chebyshev,
447:                 lsa_check_chebyshev=lsa_check_chebyshev,
448:                 lsa_resolution_ladder=lsa_resolution_ladder,
449:                 lsa_resolution_tolerance=lsa_resolution_tolerance,
450:                 lsa_t_hat_min=lsa_t_hat_min, lsa_t_hat_max=lsa_t_hat_max,
451:                 lsa_scan_points=lsa_scan_points,
...
```

### code/parameter_database/run_design_campaign_parallel.py @ 158395f3d25a29680dea5bcd0d3a24e52a28b4cd
```text
163:     ladder = ([int(v) for v in args.lsa_resolution_ladder.split(",")]
164:               if args.lsa_resolution_ladder else None)
165:     kwargs = dict(nx_max=args.nx_max, n_chebyshev=args.n_chebyshev,
166:                   lsa_resolution_ladder=ladder,
167:                   lsa_resolution_tolerance=args.lsa_resolution_tolerance,
168:                   lsa_t_hat_min=args.lsa_t_hat_min,
169:                   lsa_scan_points=args.lsa_scan_points,
170:                   lsa_chi_coeff_mode=args.lsa_chi_coeff_mode,
171:                   dumux_perturbation=args.perturbation,
172:                   dumux_random_seed=args.random_seed,
173:                   nx=args.nx, eigenprofile_samples=args.eigenprofile_samples,
174:                   nx_by_case=({int(c): int(v["nx"]) for c, v in
175:                                json.load(open(args.nx_map))["nx_map"].items()}
176:                               if args.nx_map else None))
177:     # The frozen protocol, not the published one it was measured against. Its home is
178:     # solvers.canonical_numerics(); restating the numbers here is how they drift.
...
```

### code/parameter_database/dominant_mode_cli.py @ 158395f3d25a29680dea5bcd0d3a24e52a28b4cd
```text
90:                              "limit that passing it on the command line is not safe")
91:     parser.add_argument("--out", required=True)
92:     parser.add_argument("--workers", type=int, default=8)
93:     parser.add_argument("--sharpness-min", type=float, default=None,
94:                         help="override solvers.SHARPNESS_MIN; the threshold sensitivity "
95:                              "study needs the same selection at several values")
96:     args = parser.parse_args(argv)
97:     if (args.case_ids is None) == (args.case_ids_file is None):
98:         parser.error("give exactly one of --case-ids or --case-ids-file")
99:     raw_ids = (args.case_ids if args.case_ids is not None
100:                else open(args.case_ids_file, encoding="utf-8").read())
101:
102:     con = sqlite3.connect(args.db)
103:     con.row_factory = sqlite3.Row
104:     rows = []
105:     for cid in (int(v) for v in raw_ids.replace("\n", ",").split(",") if v.strip()):
106:         row = _campaign_row(con, cid)
107:         if row is None:
108:             raise SystemExit(f"case {cid}: no searched onset_found LSA result in {args.db}")
109:         rows.append(row)
110:     con.close()
111:
112:     with ProcessPoolExecutor(max_workers=args.workers) as pool:
113:         results = dict(pool.map(_job, rows, repeat(args.sharpness_min)))
114:     nx_map = {str(c): {"nx": r["nx_dominant"], "cells_x": r["cells_x"],
115:                        "t_hat_onset_dominant": r["t_hat_onset_dominant"],
116:                        "t_ratio_to_critical": r["t_ratio_to_critical"],
117:                        "sharpness": r["sharpness"], "selection": r["selection"],
118:                        "critical_nx": r["critical_nx"]}
119:               for c, r in results.items() if r["status"] == "ok"}
120:     no_mode = {str(c): r["tried"] for c, r in results.items() if r["status"] != "ok"}
121:     payload = {"sharpness_min": next(iter(results.values()))["sharpness_min"] if results else None,
122:                "nx_map": nx_map, "no_sharp_mode": no_mode}
123:     with open(args.out, "w", encoding="utf-8") as fh:
124:         json.dump(payload, fh, indent=2, sort_keys=True)
125:     selection: dict[str, int] = {}
126:     for v in nx_map.values():
127:         selection[v["selection"]] = selection.get(v["selection"], 0) + 1
128:     same_as_critical = sum(1 for v in nx_map.values() if v["nx"] == v["critical_nx"])
129:     print(json.dumps({"selected": len(nx_map), "no_sharp_mode": len(no_mode),
130:                       "by_selection": selection, "same_as_critical_mode": same_as_critical,
131:                       "cells_x": sorted({v["cells_x"] for v in nx_map.values()})}, indent=2))
132:     return 0
133:
134:
135: if __name__ == "__main__":
136:     sys.exit(main())
...
```

### code/parameter_database/seeded_mode_cli.py @ 158395f3d25a29680dea5bcd0d3a24e52a28b4cd
```text
70:                "candidates": list(DOMINANT_MODE_CANDIDATES),
71:                "solvers_sha256": solvers_sha256}
72:     encoded = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
73:     return hashlib.sha256(encoded).hexdigest()[:16]
74:
75:
76: def _reject_reason(tried: list[dict[str, Any]]) -> str:
77:     """Why no mode was seeded: the swept threshold, or something independent of it.
78:
79:     The distinction matters for the data contract. A case rejected because every
80:     candidate crosses after its own competing risk, or because no candidate clears the
81:     resolvability floor, is rejected identically at every `sharpness_min`; only a case
82:     with viable candidates could ever have been rescued by a lower threshold.
83:     """
84:     if not tried:
85:         return "no_candidate_evaluated"
86:     if any("sharpness" in t for t in tried):
87:         return "below_resolvability_floor"
88:     if all(t.get("status") == "crossing_after_stop" for t in tried):
89:         return "all_crossings_after_competing_risk"
90:     return "no_candidate_solved"
91:
92:
93: def _rows_from_file(path: Path, solvers_sha256: str, stamp: str) -> list[tuple[Any, ...]]:
94:     payload = json.loads(path.read_text(encoding="utf-8"))
95:     theta = float(payload["sharpness_min"])
96:     phash = _protocol_hash(theta, solvers_sha256)
97:     common = dict(sharpness_min=theta, sharpness_floor=SIGMA_SHARP_MIN,
98:                   seeding_protocol_hash=phash, solvers_sha256=solvers_sha256,
99:                   source_file=path.name, created_utc=stamp)
100:     rows = []
101:     for case_id, v in payload["nx_map"].items():
102:         rows.append({**common, "case_id": int(case_id), "status": "ok",
103:                      "selection": v["selection"], "nx_seeded": v["nx"],
104:                      "critical_nx": v["critical_nx"], "cells_x": v["cells_x"],
105:                      "t_hat_onset_seeded": v["t_hat_onset_dominant"],
106:                      "t_ratio_to_critical": v["t_ratio_to_critical"],
107:                      "sharpness": v["sharpness"], "reject_reason": None})
108:     for case_id, tried in payload["no_sharp_mode"].items():
109:         rows.append({**common, "case_id": int(case_id), "status": "no_sharp_mode",
110:                      "selection": None, "nx_seeded": None, "critical_nx": None,
111:                      "cells_x": None, "t_hat_onset_seeded": None,
112:                      "t_ratio_to_critical": None, "sharpness": None,
113:                      "reject_reason": _reject_reason(tried)})
114:     return [tuple(r[c] for c in COLUMNS) for r in rows]
115:
116:
117: def main(argv: list[str] | None = None) -> int:
118:     parser = argparse.ArgumentParser(description=(__doc__ or "").split("\n\n")[0])
119:     parser.add_argument("--db", required=True)
120:     parser.add_argument("--json-dir", required=True)
121:     parser.add_argument("--glob", default="nx_*.json")
122:     parser.add_argument("--solvers-sha256", required=True,
123:                         help="sha256 of the solvers.py the sweep ran with")
124:     args = parser.parse_args(argv)
125:
126:     files = sorted(Path(args.json_dir).glob(args.glob))
127:     if not files:
128:         raise SystemExit(f"no {args.glob} under {args.json_dir}")
129:     stamp = datetime.now(timezone.utc).isoformat(timespec="seconds")
130:
131:     con = sqlite3.connect(args.db, timeout=120.0)
132:     con.execute(TABLE_DDL)
133:     known = {int(r[0]) for r in con.execute("SELECT case_id FROM parameter_cases")}
134:     insert = (f"INSERT OR REPLACE INTO lsa_seeded_modes({','.join(COLUMNS)}) "
135:               f"VALUES ({','.join('?' for _ in COLUMNS)})")
136:     total = 0
137:     for path in files:
138:         rows = _rows_from_file(path, args.solvers_sha256, stamp)
139:         unknown = {r[0] for r in rows} - known
140:         if unknown:
141:             raise SystemExit(f"{path.name}: {len(unknown)} case ids absent from the "
142:                              f"database, e.g. {sorted(unknown)[:5]}")
143:         con.executemany(insert, rows)
144:         total += len(rows)
145:     con.commit()
...
```
