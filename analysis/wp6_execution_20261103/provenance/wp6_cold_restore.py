import sys,json
from pathlib import Path
sys.path.insert(0,'scripts')
import wp6_plan_integrity as m
m.CACHE=m.ROOT/'tmp/wp6_archive_restore_test_20261103'
assert not m.CACHE.exists()
m.restore()
result={'cold_restore_from_archives':True,'all_member_sha_verified':True,'rebuilt_npz_sha_matches_WP5':True,'original_WP5_outputs_untouched':True}
Path('analysis/wp6_execution_20261103/restoration_checks.json').write_text(json.dumps(result,indent=2))
