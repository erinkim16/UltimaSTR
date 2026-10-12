"""Build sex-specific STR catalogs: <prefix>.female.json, <prefix>.male.json, <prefix>.loci.tsv"""
import argparse, json
KEEP = ["LocusId", "LocusStructure", "ReferenceRegion", "VariantType", "VariantId"]

def chroms(locus):
    regions = locus["ReferenceRegion"]
    regions = regions if isinstance(regions, list) else [regions]
    return {r.split(":")[0] for r in regions}

ap = argparse.ArgumentParser()
ap.add_argument("--catalog", required=True); ap.add_argument("--prefix", required=True)
args = ap.parse_args()
loci = json.load(open(args.catalog))
female, male = [], []
with open(f"{args.prefix}.loci.tsv", "w") as meta:
    meta.write("LocusId\tchrom\tGene\tInheritance\tNormalMax\tPathogenicMin\tDisease\n")
    for locus in loci:
        slim = {k: locus[k] for k in KEEP if k in locus}
        female.append(dict(slim))
        male.append(dict(slim, Haploid=True) if chroms(locus) & {"chrX", "chrY"} else dict(slim))
        d = (locus.get("Diseases") or [{}])[0]
        meta.write("\t".join(str(x) for x in [locus["LocusId"], ",".join(sorted(chroms(locus))),
            locus.get("Gene", ""), locus.get("Inheritance", d.get("Inheritance", "")),
            d.get("NormalMax", ""), d.get("PathogenicMin", ""), d.get("Symbol", "")]) + "\n")
json.dump(female, open(f"{args.prefix}.female.json", "w"), indent=1)
json.dump(male, open(f"{args.prefix}.male.json", "w"), indent=1)
print(f"{len(loci)} loci; {sum('Haploid' in l for l in male)} haploid in the male catalog")