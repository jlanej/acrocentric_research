#!/usr/bin/env bash
# Fetch the public annotation tracks the analyses run on. No genome assemblies
# are needed. Total ~105 MB. Usage: ./fetch_data.sh ../work
set -euo pipefail
W="${1:-../work}"
mkdir -p "$W/hprc" "$W/aux"
S3="https://s3-us-west-2.amazonaws.com/human-pangenomics"
GH="https://raw.githubusercontent.com/Platinum-Pedigree-Consortium/AcroMutRecomb/main"
T2T="$S3/T2T/CHM13/assemblies/annotation"

echo "== Platinum Pedigree acrocentric annotation =="
for f in annotation/all_samples_censat.bed annotation/flagger_nucfreq_merged.bed; do
  curl -sSL --max-time 300 -o "$W/$(basename "$f")" "$GH/$f"
done

echo "== T2T-CHM13v2.0 reference annotation (used by the report figures) =="
for f in chm13v2.0_censat_v2.1.bed chm13v2.0_SD.full.bed chm13v2.0_cytobands_allchrs.bed \
         accessibility/combined_mask.bed.gz accessibility/reference_accessibility_comparison.txt \
         chm13v1.1.rdna_model.bed; do
  curl -sSL --max-time 600 -o "$W/$(basename "$f")" "$T2T/$f"
done

echo "== HPRC release 2 per-haplotype CenSat annotations (paginated listing) =="
tok=""; : > "$W/keys.txt"
while :; do
  if [ -z "$tok" ]; then url="$S3?list-type=2&prefix=working/HPRC/&max-keys=1000"
  else enc=$(python3 -c "import urllib.parse,sys;print(urllib.parse.quote(sys.argv[1],safe=''))" "$tok")
       url="$S3?list-type=2&prefix=working/HPRC/&max-keys=1000&continuation-token=$enc"; fi
  curl -sSL --max-time 180 "$url" -o "$W/page.xml"
  tr '<' '\n' < "$W/page.xml" | grep -oE '^Key>.*' | sed 's/^Key>//' >> "$W/keys.txt"
  tok=$(tr '<' '\n' < "$W/page.xml" | grep -oE '^NextContinuationToken>.*' | sed 's/^NextContinuationToken>//' | head -1)
  [ -z "$tok" ] && break
done
grep -E 'cenSat\.bed$' "$W/keys.txt" > "$W/dl.txt"
grep -E '\.fa\.gz\.fai$' "$W/keys.txt" >> "$W/dl.txt"
grep -E 'active\.centromeres\.bed$' "$W/keys.txt" >> "$W/dl.txt"
echo "   $(wc -l < "$W/dl.txt") objects"

cat > "$W/_one.sh" <<'SH'
#!/bin/bash
k="$1"; W="$2"
case "$k" in *active.centromeres.bed|*.fai) d="$W/aux";; *) d="$W/hprc";; esac
o="$d/$(basename "$k")"
[ -s "$o" ] || curl -sSL --max-time 300 -o "$o" "https://s3-us-west-2.amazonaws.com/human-pangenomics/$k"
SH
chmod +x "$W/_one.sh"
xargs -P 12 -n 1 -I{} "$W/_one.sh" {} "$W" < "$W/dl.txt"
rm -f "$W/page.xml" "$W/_one.sh"
echo "== done: $(ls "$W"/hprc/*.cenSat.bed | wc -l) haplotype annotations in $W/hprc =="
