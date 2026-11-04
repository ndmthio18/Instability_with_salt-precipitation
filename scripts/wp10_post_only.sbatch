#!/bin/bash
#SBATCH --job-name=WP10-post-only
#SBATCH --partition=cnmix
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=00:30:00
#SBATCH --output=post_only.stdout.%j
#SBATCH --error=post_only.stderr.%j
set -euo pipefail
export PATH=/usr/local/bin:/usr/bin:/bin
case "$(hostname)" in ibcln*|test0*|gln0*) echo "refusing login/test node"; exit 2;; esac
root="$1"
work="${TMPDIR:-/tmp}/wp10_post_${SLURM_JOB_ID}"
mkdir -p "$work/run"
trap 'rm -rf "$work"' EXIT
tar -xzf "$root/products.tar.gz" -C "$work/run"
tar -xzf "$root/fields.tar.gz" -C "$work/run"
rm -rf "$work/run/post"
/WORK/hulm_work/ys/miniconda3/bin/python3 "$root/wp9_post.py" "$work/run" "$work/run/post" > "$work/post.log" 2>&1
tar -czf "$root/post_corrected.tar.gz" -C "$work/run" post post.log
sha256sum "$root/post_corrected.tar.gz" > "$root/post_corrected.sha256"
