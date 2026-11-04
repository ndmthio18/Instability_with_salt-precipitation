"""Final read-only delivery gate for WP10."""
import hashlib,json,sys
from pathlib import Path

root=Path(sys.argv[1]);run=root/"run"
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
scalar=json.loads((root/"analysis/scalar_flux_metrics.json").read_text())
field=json.loads((root/"analysis/field_metrics.json").read_text())
checks=json.loads((root/"local_post/independent_checks.json").read_text())
bench=json.loads((root/"benchmark_result.json").read_text())
completion=json.loads((root/"completion.json").read_text())
tests={
 "production_input_hash":sha(root/"params.input")=="c1441be8b69a11a23e9c8ce36839c127866d1d5567207fb0c63b6645c33d1bd1",
 "fields_archive_hash":sha(root/"fields.tar.gz")=="9863089b39ca6cc5b81558086f9a0b20554e1589507242ff8790bfbba21419f2",
 "products_archive_hash":sha(root/"products.tar.gz")=="4c05c3a219ad1c9d78d03c43cea9f79aabd86bbae9314b4140ffd64f599af3e3",
 "selected_archive_hash":sha(root/"selected_fields.tar.gz")=="2ce1729a8d159974e38b6525959ee0cb8db84a10771eefd12b5f881f42bd4b0a",
 "corrected_post_archive_hash":sha(root/"post_corrected.tar.gz")=="190bd4fcd0e1b974c8446f6ea8e7a134e35fb9def0f75071347eccf6ab015c27",
 "same_node_benchmark_identical":bench["all_bitwise_identical"],
 "solver_success":completion["scientific_run"]["solver_returncode"]==0,
 "threshold_stop":completion["scientific_run"]["stop_reason"]=="surface_solubility_threshold",
 "no_retry":not completion["scientific_run"]["automatic_restart_or_retry"],
 "old_4T_prefix_identical":all(v["bitwise_prefix"] for v in scalar["through_4T"].values()),
 "terminal_field_preserved":checks["terminal_irregular_stride"] and checks["frame_step_intervals"]==[2,50],
 "independent_conservation_passed":checks["passed"],
 "all_sparse_fields_present":len(list(run.glob("*.vtu")))==checks["frames"]==512,
 "velocity_checkpoints_complete":len(field["velocity"])==8,
 "full_mode_timeseries_complete":len(field["mode_timeseries"])==512,
 "report_present":(root/"REPORT_ZH.md").stat().st_size>5000,
}
result={"schema_version":1,"tests":tests,"passed":all(tests.values()),
        "scope":"delivery integrity and registered numerical checks; not physical validation"}
(root/"delivery_verification.json").write_text(json.dumps(result,indent=2))
print(json.dumps(result,indent=2));assert result["passed"]
