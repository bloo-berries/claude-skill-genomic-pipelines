#!/usr/bin/env bash
# Reference sanity gate: verify that every input agrees with the reference on contig naming and length.
# Usage: ref_sanity.sh <reference.fa> <file...>     files: BAM/CRAM, VCF/BCF, GTF/GFF, BED-like, .fai, .dict
# Exit 0 = consistent; exit 1 = mismatch found (stop and reconcile before running anything).
set -uo pipefail

ref=${1:?usage: ref_sanity.sh <reference.fa> <file...>}
shift
[ $# -ge 1 ] || { echo "ref_sanity: give at least one input file" >&2; exit 2; }

fai="$ref.fai"
if [ ! -e "$fai" ]; then
  if [ -e "${ref%.fa}.fai" ]; then fai="${ref%.fa}.fai"
  elif command -v samtools >/dev/null 2>&1; then
    echo "ref_sanity: no .fai next to $ref; running samtools faidx" >&2
    samtools faidx "$ref" || { echo "ref_sanity: faidx failed" >&2; exit 2; }
  else
    echo "ref_sanity: need $fai (samtools faidx $ref)" >&2; exit 2
  fi
fi

tmp=$(mktemp -d); trap 'rm -rf "$tmp"' EXIT
cut -f1,2 "$fai" | sort -k1,1 > "$tmp/ref.tsv"
ref_n=$(wc -l < "$tmp/ref.tsv" | tr -d ' ')
ref_style=$(head -n1 "$tmp/ref.tsv" | cut -f1 | grep -q '^chr' && echo chr || echo nochr)

guess_build() {   # arg: TSV name<TAB>len ; print assembly guess from chr1/1 length
  local l
  l=$(awk '$1=="chr1"||$1=="1"{print $2; exit}' "$1")
  case "${l:-}" in
    248956422) echo GRCh38 ;;
    249250621) echo GRCh37/hg19 ;;
    248387328) echo T2T-CHM13v2.0 ;;
    195154279) echo GRCm39 ;;
    195471971) echo GRCm38/mm10 ;;
    "")        echo "unknown(no chr1)" ;;
    *)         echo "unknown(chr1=$l)" ;;
  esac
}
echo "reference: $ref  contigs=$ref_n  style=$ref_style  build=$(guess_build "$tmp/ref.tsv")"

reader() { case "$1" in *.gz|*.bgz) gzip -cd "$1" ;; *) cat "$1" ;; esac; }
status=0

contigs_of() {    # emit name<TAB>len (len may be empty) for a file
  local f=$1
  case "$f" in
    *.bam|*.cram|*.sam)
      samtools view -H "$f" | awk -F'\t' '$1=="@SQ"{n="";l="";for(i=2;i<=NF;i++){if($i~/^SN:/)n=substr($i,4);if($i~/^LN:/)l=substr($i,4)} print n"\t"l}' ;;
    *.vcf|*.vcf.gz|*.vcf.bgz|*.bcf)
      if command -v bcftools >/dev/null 2>&1; then
        bcftools view -h "$f" | awk -F'[<>,=]' '/^##contig/{n="";l="";for(i=1;i<=NF;i++){if($i=="ID")n=$(i+1);if($i=="length")l=$(i+1)} print n"\t"l}'
      else
        reader "$f" | awk -F'[<>,=]' '/^##contig/{n="";l="";for(i=1;i<=NF;i++){if($i=="ID")n=$(i+1);if($i=="length")l=$(i+1)} print n"\t"l} !/^#/{exit}'
      fi
      # fall back to body contigs if header lacks ##contig lines
      ;;
    *.fai)  cut -f1,2 "$f" ;;
    *.dict) awk -F'\t' '$1=="@SQ"{n="";l="";for(i=2;i<=NF;i++){if($i~/^SN:/)n=substr($i,4);if($i~/^LN:/)l=substr($i,4)} print n"\t"l}' "$f" ;;
    *)      reader "$f" | grep -v '^#' | head -n 2000000 | cut -f1 | sort -u | awk '{print $1"\t"}' ;;
  esac
}

for f in "$@"; do
  [ -e "$f" ] || { echo "MISSING  $f"; status=1; continue; }
  contigs_of "$f" | awk -F'\t' '$1!=""' | sort -k1,1 -u > "$tmp/in.tsv"
  n=$(wc -l < "$tmp/in.tsv" | tr -d ' ')
  if [ "$n" -eq 0 ] && [[ $f == *.vcf* || $f == *.bcf ]]; then
    (command -v bcftools >/dev/null 2>&1 && bcftools view -H "$f" | head -n 2000000 | cut -f1 || reader "$f" | grep -v '^#' | head -n 2000000 | cut -f1) \
      | sort -u | awk '{print $1"\t"}' > "$tmp/in.tsv"
    n=$(wc -l < "$tmp/in.tsv" | tr -d ' ')
  fi
  style=$(head -n1 "$tmp/in.tsv" | cut -f1 | grep -q '^chr' && echo chr || echo nochr)
  missing=$(join -t $'\t' -v1 "$tmp/in.tsv" "$tmp/ref.tsv" | cut -f1)
  missing_n=$(printf '%s' "$missing" | grep -c . || true)
  lenmis=$(join -t $'\t' "$tmp/in.tsv" "$tmp/ref.tsv" | awk -F'\t' '$2!="" && $2!=$3 {print $1"("$2" vs "$3")"}')
  lenmis_n=$(printf '%s' "$lenmis" | grep -c . || true)
  build=$(guess_build "$tmp/in.tsv")

  verdict=OK
  notes=""
  [ "$style" != "$ref_style" ] && { verdict=MISMATCH; notes+=" naming:$style-vs-$ref_style"; }
  [ "$missing_n" -gt 0 ] && { verdict=MISMATCH; notes+=" ${missing_n}_contigs_not_in_ref($(printf '%s' "$missing" | head -n 3 | tr '\n' ',')…)"; }
  [ "$lenmis_n" -gt 0 ] && { verdict=MISMATCH; notes+=" ${lenmis_n}_length_mismatches($(printf '%s' "$lenmis" | head -n 2 | tr '\n' ',')…)"; }
  case "$build" in unknown*) ;; *) [ "$build" != "$(guess_build "$tmp/ref.tsv")" ] && { verdict=MISMATCH; notes+=" build:$build"; } ;; esac
  [ "$verdict" = MISMATCH ] && status=1
  printf '%-8s %s  contigs=%s style=%s build=%s%s\n' "$verdict" "$f" "$n" "$style" "$build" "$notes"
done

if [ $status -ne 0 ]; then
  cat >&2 <<'EOF'
ref_sanity: at least one input disagrees with the reference.
  naming  -> bcftools annotate --rename-chrs map.txt (VCF) | samtools reheader (BAM) | sed on BED/GTF; never mix styles
  build   -> stop: realign/lift over (CrossMap, gatk LiftoverVcf) rather than renaming
  missing -> alt/decoy/HLA contigs: harmless for BAM if reads there are ignored; fatal for GATK -L/VCF contig order
EOF
fi
exit $status
