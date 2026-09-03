#!/usr/bin/env bash
set -euo pipefail

usage() {
  echo "usage: $0 PLAN.tsv OUTPUT_DIR" >&2
}

if [[ $# -ne 2 ]]; then
  usage
  exit 2
fi

plan_path=$1
result_root=$2
if [[ ! -f "$plan_path" ]]; then
  echo "benchmark plan not found: $plan_path" >&2
  exit 2
fi
if rg -n 'OWNER_TBD' "$plan_path" >/dev/null; then
  echo "benchmark plan still contains OWNER_TBD fields; refusing to run" >&2
  exit 2
fi
if [[ ! -x /usr/bin/time ]]; then
  echo "/usr/bin/time is required" >&2
  exit 2
fi

declare -A plan
while IFS=$'\t' read -r key value; do
  [[ "$key" == "key" || -z "$key" ]] && continue
  plan["$key"]=$value
done < "$plan_path"

required_keys=(reference_path query_path threads repetitions ramag_binary
  ramag_seed_mode min_match headline_endpoint mummer_command_template)
for key in "${required_keys[@]}"; do
  if [[ -z "${plan[$key]:-}" ]]; then
    echo "missing benchmark-plan value: $key" >&2
    exit 2
  fi
done

mkdir -p "$result_root"
cp "$plan_path" "$result_root/accepted-benchmark-plan.tsv"
uname -a > "$result_root/host.uname.txt"
lscpu > "$result_root/host.lscpu.txt"
"${plan[ramag_binary]}" --version > "$result_root/ramag.version.txt"

run_ramag() {
  local repetition=$1
  local run_dir="$result_root/ramag-delta-r${repetition}"
  local prefix="$run_dir/result"
  mkdir -p "$run_dir/work"
  local -a command=("${plan[ramag_binary]}" align
    --reference "${plan[reference_path]}"
    --query "${plan[query_path]}"
    --output-prefix "$prefix"
    --work-dir "$run_dir/work"
    --threads "${plan[threads]}"
    --formats delta
    --seed-mode "${plan[ramag_seed_mode]}"
    --min-match "${plan[min_match]}")
  printf '%q ' "${command[@]}" > "$run_dir/command.sh.txt"
  printf '\n' >> "$run_dir/command.sh.txt"
  set +e
  /usr/bin/time -v -o "$run_dir/time.raw.txt" \
    "${command[@]}" > "$run_dir/stdout.txt" 2> "$run_dir/stderr.txt"
  local exit_code=$?
  set -e
  printf '%s\n' "$exit_code" > "$run_dir/exit-code.txt"
  if [[ $exit_code -ne 0 || ! -s "$prefix.delta" || \
        ! -s "$prefix.manifest.json" || ! -s "$prefix.complete" ]]; then
    touch "$run_dir.INELIGIBLE"
    return 1
  fi
  if [[ -n "${plan[mummer_show_coords_binary]:-}" ]]; then
    "${plan[mummer_show_coords_binary]}" -rcl "$prefix.delta" \
      > "$run_dir/show-coords.txt"
  fi
  wc -c "$prefix.delta" "$prefix.manifest.json" > "$run_dir/artifact-bytes.txt"
  touch "$run_dir.ELIGIBLE"
}

run_mummer() {
  local repetition=$1
  local run_dir="$result_root/mummer-delta-r${repetition}"
  mkdir -p "$run_dir"
  # The accepted command is intentionally explicit: the harness does not guess
  # nucmer flags, direction, seed semantics, or endpoint on the owner's behalf.
  local command=${plan[mummer_command_template]}
  command=${command//\{reference\}/${plan[reference_path]}}
  command=${command//\{query\}/${plan[query_path]}}
  command=${command//\{threads\}/${plan[threads]}}
  command=${command//\{min_match\}/${plan[min_match]}}
  command=${command//\{output_dir\}/$run_dir}
  printf '%s\n' "$command" > "$run_dir/command.sh.txt"
  set +e
  /usr/bin/time -v -o "$run_dir/time.raw.txt" \
    bash -c "$command" > "$run_dir/stdout.txt" 2> "$run_dir/stderr.txt"
  local exit_code=$?
  set -e
  printf '%s\n' "$exit_code" > "$run_dir/exit-code.txt"
  if [[ $exit_code -ne 0 ]]; then
    touch "$run_dir.INELIGIBLE"
    return 1
  fi
  touch "$run_dir.ELIGIBLE"
}

failed=0
for ((repetition = 1; repetition <= plan[repetitions]; ++repetition)); do
  run_ramag "$repetition" || failed=1
  run_mummer "$repetition" || failed=1
done

if [[ $failed -ne 0 ]]; then
  echo "one or more runs were ineligible; inspect *.INELIGIBLE" >&2
  exit 1
fi
echo "all raw runs completed; aggregation still requires correctness gates"
