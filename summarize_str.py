"""Combine per-sample STR BEDs into one table and flag alleles against known_loci.loci.tsv."""
import argparse, csv, glob, os, re

ap = argparse.ArgumentParser()
ap.add_argument("--results", required=True)           # folder with *_genotypes.bed
ap.add_argument("--loci", default="catalogs/known_loci.loci.tsv")
ap.add_argument("--out", default="cohort_str_calls.tsv")
args = ap.parse_args()

limits = {r["LocusId"]: r for r in csv.DictReader(open(args.loci), delimiter="\t")}
num = lambda s: int(s) if s not in ("", "None") else None

with open(args.out, "w") as out:
    out.write("sample\tlocus\tchrom\tpos\tgenotype\tsupport\tspanning_reads\tnormal_max\tpathogenic_min\tflag\n")
    for bed in sorted(glob.glob(os.path.join(args.results, "*_genotypes.bed"))):
        sample = os.path.basename(bed).replace("_genotypes.bed", "")
        for line in open(bed):
            chrom, start, end, info = line.rstrip("\n").split("\t")[:4]
            f = dict(kv.split("=", 1) for kv in info.split(";") if "=" in kv)
            lim = limits.get(f["VARID"], {})        # thresholds apply to the main variant of a locus
            nmax, pmin = num(lim.get("NormalMax", "")), num(lim.get("PathogenicMin", ""))
            alleles = [int(a) for a in re.findall(r"\d+", f.get("REPCN", ""))]
            top = max(alleles) if alleles else None
            flag = ""
            if top is not None and pmin is not None and top >= pmin:
                flag = "PATHOGENIC_RANGE"
            elif top is not None and nmax is not None and top > nmax:
                flag = "INTERMEDIATE"
            if "SPANNING" not in f.get("SO", ""):
                flag = (flag + "," if flag else "") + "FLANKING_ONLY"
            out.write("\t".join(map(str, [sample, f["VARID"], chrom, start, f.get("GT", ""),
                      f.get("SO", ""), f.get("ADSP", ""), nmax or "", pmin or "", flag])) + "\n")
print("wrote", args.out)
