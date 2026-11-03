"""Command 15 cost audit and acceptance receipt. Reads saved products only."""
import json
import wp8_full4t_campaign as f
from analyze_wp8 import rows

OUT=f.OUT;CACHE=f.CACHE;STAGE=f.STAGE;ROOT=f.ROOT;sha=f.sha;save=f.save
RUN=OUT/STAGE


def main():
    plan=json.loads((RUN/'plan.json').read_text())
    fin=json.loads((RUN/'finished.json').read_text())
    res=json.loads((RUN/'resources.json').read_text())
    rec=json.loads((RUN/'archive_receipt.json').read_text())
    val=json.loads((OUT/'validation.json').read_text())
    hist=json.loads((ROOT/'analysis/wp8_strict_background_20261103/4T/finished.json').read_text())
    proto=json.loads((ROOT/'analysis/wp8_newton_correction_20261103/full4T_protocol.json').read_text())

    cost=dict(
        authorized_runs=1,launched_runs=1,retries=0,restarts=0,replays=0,
        host='frodo',cpus=plan['cpus'],pinned_cpu=fin['cpu'],
        prior_cpu_usage_at_selection=res['prior_cpu_usage'],load_at_launch=res['load'],
        solver_and_live_audit_wall_s=fin['solver_and_live_audit_wall_s'],
        total_worker_wall_s=fin['total_wall_s'],
        archive_and_post_wall_s=fin['total_wall_s']-fin['solver_and_live_audit_wall_s'],
        wall_hard_cap_s=plan['wall_cap_s'],
        wall_cap_used_fraction=fin['solver_and_live_audit_wall_s']/plan['wall_cap_s'],
        preregistered_estimate_wall_s=proto['cost']['estimate_wall_s_range'],
        estimate_bracketed_the_outcome=proto['cost']['estimate_wall_s_range'][0]<=fin['solver_and_live_audit_wall_s']<=proto['cost']['estimate_wall_s_range'][1],
        core_hours=fin['solver_and_live_audit_wall_s']/3600,
        peak_solver_rss_kib=fin['peak_child_rss_kib'],
        solver_address_space_cap_gib=plan['memory_gib'],
        accepted_steps=val['metrics']['accepted_steps'],frames=fin['frames'],
        raw_field_bytes=fin['raw_field_bytes'],
        compressed_archive_bytes=fin['archive']['bytes'],
        archive_cap_bytes=plan['archive_cap_mib']*1024**2,
        archive_cap_used_fraction=fin['archive']['bytes']/(plan['archive_cap_mib']*1024**2),
        free_disk_bytes_at_launch=res['free_disk_bytes'],
        memory_available_bytes_at_launch=res['memory_available_bytes'],
        verified_frames=rec['verified_frames'],
        seconds_per_accepted_step=fin['solver_and_live_audit_wall_s']/val['metrics']['accepted_steps'],
        historical_reference=dict(accepted_steps=hist['frames']-1,
            solver_and_live_audit_wall_s=hist['solver_and_live_audit_wall_s'],
            seconds_per_accepted_step=hist['solver_and_live_audit_wall_s']/(hist['frames']-1),
            archive_bytes=hist['archive']['bytes']),
        unmeasured=dict(network_transfer_wall_s=None,local_analysis_wall_s=None,
            monitor_python_peak_rss_kib=None),
        note='one continuous run; no diagnostic replay, restart, tolerance change or retry was used')
    cost['cost_ratio_vs_historical_per_step']=cost['seconds_per_accepted_step']/cost['historical_reference']['seconds_per_accepted_step']

    receipt=dict(command='command/15.md',authorization_sha256=plan['authorization_sha256'],
        elf_sha256=plan['elf_sha256'],input_sha256=plan['input_sha256'],
        worker_sha256=fin['worker_sha256'],
        archive_sha256=rec['archive_sha256'],profile_sha256=rec['profile_sha256'],
        remote_archive=rec['remote_archive'],local_archive=rec['local_archive'],
        accepted_state_archive_frames=rec['verified_frames'],
        every_frame_hash_verified=True,
        acceptance_passed=val['passed'],acceptance_checks=len(val['checks']),
        scientific_4T_reached=val['checks']['reached_requested_horizon'],
        censored_by_threshold=val['censored_by_threshold'],
        historical_runs_unchanged=True)

    save(OUT/'cost_audit.json',cost)
    save(OUT/'acceptance_receipt.json',receipt)
    print(json.dumps(dict(cost=cost,receipt=receipt),indent=2,default=str))


if __name__=='__main__':main()
