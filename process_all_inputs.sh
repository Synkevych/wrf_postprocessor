#!/usr/bin/env bash
set -euo pipefail

input_dir="${1:-input}"
output_dir="${2:-output}"
ntimes1="${3:-28}"

if [[ ! -d "$input_dir" ]]; then
  echo "Input directory not found: $input_dir" >&2
  exit 1
fi

if [[ ! -x ./extract_wrf_fields ]]; then
  echo "Missing executable ./extract_wrf_fields. Run ./compile.sh first." >&2
  exit 1
fi

mkdir -p "$output_dir"

backup_config=""
if [[ -f config.nml ]]; then
  backup_config="$(mktemp)"
  cp config.nml "$backup_config"
fi

restore_config() {
  if [[ -n "$backup_config" && -f "$backup_config" ]]; then
    cp "$backup_config" config.nml
    rm -f "$backup_config"
  else
    rm -f config.nml
  fi
}
trap restore_config EXIT

shopt -s nullglob
files=("$input_dir"/*)

processed=0
for fpath in "${files[@]}"; do
  [[ -f "$fpath" ]] || continue

  fname="$(basename "$fpath")"
  stem="${fname%.*}"
  out_subdir="$output_dir/$stem"
  mkdir -p "$out_subdir"

  cat > config.nml <<NML
&io_nml
  infile        = "$fpath"
  grid_outfile  = "$out_subdir/grid.dat"
  pmsl_outfile  = "$out_subdir/pmsl_"
  ntimes1       = $ntimes1
/
NML

  echo "Processing: $fpath -> $out_subdir"
  ./extract_wrf_fields
  processed=$((processed + 1))
done

if [[ "$processed" -eq 0 ]]; then
  echo "No input files found in: $input_dir" >&2
  exit 1
fi

echo "Processed $processed file(s). Outputs are in: $output_dir"
