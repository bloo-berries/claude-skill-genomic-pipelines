#!/usr/bin/env bash
# Smoke test for .claude/skills/genomics/scripts/*.sh using synthetic files (no bioinformatics tools required).
set -euo pipefail
S=$(cd "$(dirname "$0")/.." && pwd)/.claude/skills/genomics/scripts
T=$(mktemp -d); trap 'rm -rf "$T"' EXIT
cd "$T"

printf '>chr1\nACGTACGTACGT\n>chr2\nACGTAC\n' > ref.fa
printf 'chr1\t248956422\t6\t12\t13\nchr2\t6\t19\t6\t7\n' > ref.fa.fai
printf '##fileformat=VCFv4.2\n##contig=<ID=chr1,length=248956422>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchr1\t5\t.\tA\tC\t50\tPASS\t.\n' > good.vcf
printf '##fileformat=VCFv4.2\n##contig=<ID=1,length=249250621>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n1\t5\t.\tA\tC\t50\tPASS\t.\n' > bad.vcf
printf 'chr1\t10\t20\n' > ok.bed
printf '@r1 1:N:0:ACGT\nACGTACGTAC\n+\nIIIIIIIIII\n' | gzip > r.fq.gz

fail=0
# ref_sanity: consistent inputs pass
bash "$S/ref_sanity.sh" ref.fa good.vcf ok.bed >/dev/null || { echo "FAIL ref_sanity should pass on consistent inputs"; fail=1; }
# ref_sanity: build + naming mismatch fails with exit 1
if bash "$S/ref_sanity.sh" ref.fa bad.vcf >/dev/null 2>&1; then echo "FAIL ref_sanity should reject GRCh37/nochr VCF"; fail=1; fi
out=$(bash "$S/ref_sanity.sh" ref.fa bad.vcf 2>/dev/null || true)
grep -q 'MISMATCH' <<<"$out" || { echo "FAIL ref_sanity verdict missing"; fail=1; }
grep -q 'GRCh37' <<<"$out" || { echo "FAIL ref_sanity build guess missing"; fail=1; }

# peek: bounded output, detects Illumina header, never fails on odd inputs
out=$(bash "$S/peek.sh" r.fq.gz 1)
grep -q 'Illumina' <<<"$out" || { echo "FAIL peek fastq header detection"; fail=1; }
[ "$(wc -l <<<"$out")" -le 12 ] || { echo "FAIL peek fastq output too long"; fail=1; }
bash "$S/peek.sh" good.vcf 1 >/dev/null || { echo "FAIL peek vcf"; fail=1; }
bash "$S/peek.sh" ok.bed >/dev/null || { echo "FAIL peek bed"; fail=1; }
bash "$S/peek.sh" ref.fa >/dev/null || { echo "FAIL peek fasta"; fail=1; }
if bash "$S/peek.sh" does-not-exist 2>/dev/null; then echo "FAIL peek should error on missing file"; fail=1; fi

# check_env: always exit 0, one row per tool, MISSING for unknown tool
out=$(bash "$S/check_env.sh" bash definitely-not-a-tool)
grep -q $'^bash\t' <<<"$out" || { echo "FAIL check_env bash row"; fail=1; }
grep -q $'^definitely-not-a-tool\tMISSING' <<<"$out" || { echo "FAIL check_env MISSING row"; fail=1; }
grep -q $'^cpus\t' <<<"$out" || { echo "FAIL check_env cpus row"; fail=1; }

[ $fail -eq 0 ] && echo "scripts OK"
exit $fail
