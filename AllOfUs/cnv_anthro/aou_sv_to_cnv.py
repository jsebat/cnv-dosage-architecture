"""All of Us adapter (1/2): GATK-SV short-read WGS SV VCF -> per-person CNV table in PLINK .cnv layout
(FID IID CHR BP1 BP2 TYPE SCORE SITES), the same table every downstream module reads.

AoU provides (support-hub docs, CDRv8/v9): a joint-called GATK-SV srWGS SV VCF (~97k participants), a sites-only VCF,
and a maximal-unrelated sample list. There is NO curated array CNV callset (only GDA array genotypes/IDATs), so the
srWGS SV callset is the CNV source.

Mapping rules (written down because they change what a 'probe' means; see README):
  * DEL / DUP records (biallelic, SVTYPE=DEL|DUP): carrier if GT contains a non-reference allele (0/1, 1/1, 1|0 ...).
  * multiallelic CNV records (SVTYPE=CNV, GATK-SV FILTER=MULTIALLELIC): copy number from FORMAT CN (or RD_CN);
    CN < 2 -> DEL, CN > 2 -> DUP (autosomes only, diploid reference).
  * FILTER: PASS (plus MULTIALLELIC for SVTYPE=CNV) by default.
  * size >= 10 kb (Omar's 10 kb rule). The array '10 probes' rule has no WGS analogue and is not applied.
  * frequency: the PGC Freeze_2024 rule (data_10perc, filter TAG=FREQUENCY): a call is removed when > 10% of the
    cohort's samples carry a same-type call with >= 50% reciprocal overlap (call-level, as in the PGC calling pipeline's
    PLTF_*_50_reciprocal frequencies). UKBB calls (from the freeze) already passed exactly this rule.
  * CHR coded numerically (chr1 -> 1); chroms 1-22 by default (male X dosage is not comparable between GT and CN records).
Consequence for Omar's probe framework: breakpoints of joint-called WGS sites are shared exactly by all carriers of a
site, so each site contributes two probes with identical carrier sets (the Jaccard merge collapses them); probe density
follows SV sites, not array markers. Use --bp_dir with the PGC master breakpoints for array-comparable probes, or the
cohort breakpoints (default) for WGS-native probes.
Reads plain .vcf/.vcf.gz or streams gs:// via `gsutil cat` (Workbench). Only lines that pass the INFO prefilter
have their genotype columns split, so the 97k-sample VCF streams in one pass.
"""
import argparse
import bisect
import gzip
import io
import os
import re
import subprocess
import sys


def open_vcf(path):
    if path.startswith("gs://"):
        proj = os.environ.get("GOOGLE_PROJECT")                      # requester-pays billing project in the Workbench
        cmd = ["gsutil"] + (["-u", proj] if proj else []) + ["cat", path]
        raw = subprocess.Popen(cmd, stdout=subprocess.PIPE).stdout
        return io.TextIOWrapper(gzip.GzipFile(fileobj=raw) if path.endswith(".gz") else raw)
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def info_get(info, key):
    m = re.search(r"(?:^|;)" + key + r"=([^;]+)", info)
    return m.group(1) if m else None


def read_sites(a, chroms, okf, keep):
    """Returns (records, n_samples): records = (chrom, bp1, bp2, typecode, [carrier IDs])."""
    recs = []; use = []; n_sites = 0
    with open_vcf(a.vcf) as f:
        samples = None
        for line in f:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                samples = line.rstrip("\n").split("\t")[9:]
                use = [i for i, s in enumerate(samples) if keep is None or s in keep]
                if a.write_samples:
                    with open(a.write_samples, "w") as fs:
                        fs.write("\n".join(samples[i] for i in use) + "\n")
                continue
            head = line.split("\t", 8)                    # cheap prefilter before touching genotypes
            chrom = head[0].replace("chr", "")
            if chrom not in chroms or head[6] not in okf:
                continue
            info = head[7]
            svt = info_get(info, "SVTYPE")
            if svt not in ("DEL", "DUP", "CNV"):
                continue
            if svt != "CNV" and head[6] != "PASS":
                continue
            bp1 = int(head[1]); end = info_get(info, "END"); bp2 = int(end) if end else bp1
            if bp2 - bp1 < a.min_size:
                continue
            fields = line.rstrip("\n").split("\t")
            fmt = fields[8].split(":"); gts = fields[9:]
            dels, dups = [], []
            if svt in ("DEL", "DUP"):
                gi = fmt.index("GT"); dst = dels if svt == "DEL" else dups
                for i in use:
                    if "1" in gts[i].split(":")[gi].replace("|", "/").split("/"):
                        dst.append(samples[i])
            else:
                ci = fmt.index("CN") if "CN" in fmt else (fmt.index("RD_CN") if "RD_CN" in fmt else None)
                if ci is None:
                    continue
                for i in use:
                    v = gts[i].split(":")
                    if ci < len(v) and v[ci] not in (".", ""):
                        cn = int(float(v[ci]))
                        if cn < 2:
                            dels.append(samples[i])
                        elif cn > 2:
                            dups.append(samples[i])
            n_sites += 1
            if dels:
                recs.append((chrom, bp1, bp2, 1, dels))
            if dups:
                recs.append((chrom, bp1, bp2, 3, dups))
    return recs, max(len(use), 1), n_sites


def pgc_frequency_drop(recs, n, max_freq):
    """Indices of records failing the PGC rule: > max_freq of samples carry a same-type call with >= 50% reciprocal overlap."""
    by = {}
    for j, r in enumerate(recs):
        by.setdefault((r[0], r[3]), []).append(j)
    drop = set()
    for idx in by.values():
        idx.sort(key=lambda j: recs[j][1])
        starts = [recs[j][1] for j in idx]
        maxlen = max(recs[j][2] - recs[j][1] for j in idx)
        for j in idx:
            _, s1, e1, _, _ = recs[j]
            carriers = set()
            for m in range(bisect.bisect_left(starts, s1 - maxlen), bisect.bisect_right(starts, e1)):   # possible overlaps only
                _, s2, e2, _, c2 = recs[idx[m]]
                ov = min(e1, e2) - max(s1, s2)
                if ov > 0 and ov >= 0.5 * (e1 - s1) and ov >= 0.5 * (e2 - s2):
                    carriers.update(c2)
            if len(carriers) / n > max_freq:
                drop.add(j)
    return drop


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--vcf", required=True); ap.add_argument("--out", required=True)
    ap.add_argument("--min_size", type=int, default=10000); ap.add_argument("--max_freq", type=float, default=0.10)
    ap.add_argument("--chroms", default="1-22"); ap.add_argument("--keep_samples", default=None, help="file of sample IDs to keep")
    ap.add_argument("--filters", default="PASS,MULTIALLELIC")
    ap.add_argument("--write_samples", default=None, help="write the VCF sample list (people in the CNV source)")
    a = ap.parse_args()
    lo, hi = map(int, a.chroms.split("-")); chroms = {str(c) for c in range(lo, hi + 1)}
    okf = set(a.filters.split(","))
    keep = set(l.strip() for l in open(a.keep_samples)) if a.keep_samples else None
    recs, n, n_sites = read_sites(a, chroms, okf, keep)
    drop = pgc_frequency_drop(recs, n, a.max_freq)
    n_calls = 0
    with open(a.out, "w") as o:
        o.write("FID\tIID\tCHR\tBP1\tBP2\tTYPE\tSCORE\tSITES\n")
        for j, (ch, s1, e1, tc, cs) in enumerate(recs):
            if j in drop:
                continue
            for s in cs:
                o.write(f"0\t{s}\t{ch}\t{s1}\t{e1}\t{tc}\t0\t0\n")
                n_calls += 1
    print(f"sites used {n_sites}; calls written {n_calls}; site-dosages dropped by the PGC rule "
          f"(> {a.max_freq:.0%} of samples with a same-type call at >= 50% reciprocal overlap) {len(drop)}", file=sys.stderr)


if __name__ == "__main__":
    main()
