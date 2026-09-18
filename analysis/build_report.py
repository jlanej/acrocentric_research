import json, re, os, datetime, importlib
import render
importlib.reload(render)

REFS = json.load(open("handoff/refs.json"))
PARTS = ["report_01.md","report_02.md","report_03.md","report_04.md",
         "report_05.md","report_06.md","report_07.md"]
md = "\n\n".join(open(p).read() for p in PARTS)

# ---- citation numbering in order of first appearance -------------------
order, seen = [], {}
missing = set()
def repl(m):
    keys = [k.strip().lstrip("@") for k in m.group(1).split(",")]
    nums = []
    for k in keys:
        if k not in REFS:
            missing.add(k); continue
        if k not in seen:
            order.append(k); seen[k] = len(order)
        nums.append(seen[k])
    if not nums: return ""
    nums = sorted(set(nums))
    # collapse runs
    out, i = [], 0
    while i < len(nums):
        j = i
        while j+1 < len(nums) and nums[j+1] == nums[j]+1: j += 1
        out.append(str(nums[i]) if j == i else "%d\u2013%d" % (nums[i], nums[j]))
        i = j+1
    return "<super>%s</super>" % ",".join(out)

md = re.sub(r"\[@([^\]]+)\]", repl, md)
if missing: raise SystemExit("MISSING REFERENCE KEYS: %s" % sorted(missing))

# ---- bibliography ------------------------------------------------------
bib = ["## 29. References", "",
       "References are numbered in order of first appearance. Every entry was verified "
       "against bibliographic metadata retrieved programmatically; data resources are "
       "cited by file name and repository.", ""]
for i, k in enumerate(order, 1):
    au, ti, yr, ven, doi = REFS[k]
    doi_s = ("doi:" + doi) if doi and not doi.startswith(("KY9", "chm13")) else doi
    bib.append("%d. %s. *%s*. %s **%s**. %s" % (i, au, ti, ven, yr, doi_s))
md += "\n\n" + "\n".join(bib) + "\n"

TITLE = "The Acrocentric Short Arms"
SUB = ("Sequence, function, assembly, alignment and variant analysis of the last 66 Mb "
       "of the human genome \u2014 a consolidated review, with an agenda")
BLURB = ("The five acrocentric short arms comprise 66.07 Mb of T2T-CHM13v2.0 and 2.1% of the human "
         "genome. GRCh38 represents 83% of that territory as N; only 4.3% of it is short-read "
         "accessible, against 86.9% of the autosomal genome; between 22% and 61% of each arm has a "
         "\u226599%-identical copy on a different acrocentric chromosome; and the de novo "
         "single-nucleotide mutation rate within it is roughly tenfold the autosomal euchromatic rate. "
         "This report consolidates what is known about these regions, sets out how to assemble, align "
         "and analyse them from single nucleotides to whole-arm rearrangements, and proposes specific "
         "solutions to the coordinate, assignability, representation and benchmarking problems that "
         "currently prevent them from being studied at scale.")
META = [
 "Prepared %s." % datetime.date.today().strftime("%-d %B %Y"),
 "",
 "All quantitative statements about reference content were computed for this report from the "
 "CenSat v2.1 satellite annotation, the SEDEF segmental-duplication call set, the short-read "
 "accessibility mask and the minimum-unique-<i>k</i>-mer-length tracks distributed for "
 "T2T-CHM13v2.0. Methods in section 25; data resources in section 27.",
 "",
 "%d figures, %d tables, %d references." % (11, 8, len(order)),
]
out = render.build("acrocentric_short_arms_report.pdf", TITLE, SUB, META, BLURB, md, figdir=".")
print("built:", out, os.path.getsize(out), "bytes | refs cited:", len(order))
unused = sorted(set(REFS) - set(order))
print("uncited entries dropped from bibliography:", len(unused))
