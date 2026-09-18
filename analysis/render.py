"""Minimal markup -> PDF renderer (reportlab Platypus) for long technical reports."""
import re, csv, os
from reportlab.lib import colors
from reportlab.lib.pagesizes import letter
from reportlab.lib.styles import ParagraphStyle
from reportlab.lib.units import inch
from reportlab.lib.enums import TA_JUSTIFY, TA_LEFT, TA_CENTER
from reportlab.platypus import (BaseDocTemplate, PageTemplate, Frame, Paragraph, Spacer,
                                Image, Table, TableStyle, KeepTogether, PageBreak, Flowable)
from reportlab.platypus.tableofcontents import TableOfContents
from reportlab.pdfbase.pdfmetrics import stringWidth

PW, PH = letter
LM = RM = 0.82 * inch
TM = 0.78 * inch
BM = 0.85 * inch
AVAIL = PW - LM - RM

INK = colors.HexColor("#111318")
MUTED = colors.HexColor("#5c6470")
RULE = colors.HexColor("#c9cfd8")
ACC = colors.HexColor("#1f4e79")
BOXBG = colors.HexColor("#f2f5f9")
BOXED = colors.HexColor("#b9c6d6")
HYPBG = colors.HexColor("#fbf4e8")
HYPED = colors.HexColor("#e0c99a")

BODY = ParagraphStyle("body", fontName="Times-Roman", fontSize=9.7, leading=13.3,
                      alignment=TA_JUSTIFY, textColor=INK, spaceAfter=5.2)
BODY1 = ParagraphStyle("body1", parent=BODY, firstLineIndent=0)
H1 = ParagraphStyle("h1", fontName="Helvetica-Bold", fontSize=17, leading=20.5,
                    textColor=ACC, spaceBefore=0, spaceAfter=3, alignment=TA_LEFT)
H1SUB = ParagraphStyle("h1sub", fontName="Helvetica-Oblique", fontSize=9.6, leading=12.5,
                       textColor=MUTED, spaceAfter=13)
H2 = ParagraphStyle("h2", fontName="Helvetica-Bold", fontSize=12.4, leading=15,
                    textColor=ACC, spaceBefore=14, spaceAfter=4.5)
H3 = ParagraphStyle("h3", fontName="Helvetica-Bold", fontSize=10.3, leading=13,
                    textColor=INK, spaceBefore=10, spaceAfter=3)
H4 = ParagraphStyle("h4", fontName="Helvetica-BoldOblique", fontSize=9.6, leading=12.4,
                    textColor=colors.HexColor("#2c3542"), spaceBefore=7, spaceAfter=2)
BULL = ParagraphStyle("bull", parent=BODY, leftIndent=13, bulletIndent=3, spaceAfter=3.0)
NUMB = ParagraphStyle("numb", parent=BODY, leftIndent=17, bulletIndent=2, spaceAfter=3.4)
CAP = ParagraphStyle("cap", fontName="Helvetica", fontSize=8.1, leading=10.6,
                     textColor=colors.HexColor("#3a424f"), alignment=TA_LEFT, spaceBefore=4)
TCELL = ParagraphStyle("tcell", fontName="Helvetica", fontSize=7.5, leading=9.3, textColor=INK)
THEAD = ParagraphStyle("thead", fontName="Helvetica-Bold", fontSize=7.5, leading=9.3,
                       textColor=colors.white)
BOXBODY = ParagraphStyle("boxbody", parent=BODY, fontSize=9.2, leading=12.5, spaceAfter=4)
BOXTITLE = ParagraphStyle("boxtitle", fontName="Helvetica-Bold", fontSize=9.0, leading=11.5,
                          textColor=ACC, spaceAfter=3)
REFST = ParagraphStyle("ref", fontName="Times-Roman", fontSize=8.2, leading=10.5,
                       leftIndent=15, firstLineIndent=-15, spaceAfter=2.6, textColor=INK)
TITLE = ParagraphStyle("title", fontName="Helvetica-Bold", fontSize=25, leading=29,
                       textColor=ACC, alignment=TA_LEFT)
STITLE = ParagraphStyle("stitle", fontName="Helvetica", fontSize=13, leading=17,
                        textColor=colors.HexColor("#33404f"), alignment=TA_LEFT)
META = ParagraphStyle("meta", fontName="Helvetica", fontSize=9, leading=13, textColor=MUTED)

TAGRE = re.compile(r"</?(?:b|i|u|sup|super|sub|br|font|link|para)\b[^>]*/?>", re.I)


def esc(t):
    out, last = [], 0
    for m in TAGRE.finditer(t):
        out.append(t[last:m.start()].replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;"))
        out.append(m.group(0)); last = m.end()
    out.append(t[last:].replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;"))
    return "".join(out)


def inline(t):
    t = esc(t)
    t = re.sub(r"\*\*(.+?)\*\*", r"<b>\1</b>", t)
    t = re.sub(r"(?<![\w*])\*([^*\n]+?)\*(?![\w*])", r"<i>\1</i>", t)
    t = re.sub(r"`([^`]+?)`", r'<font face="Courier" size="8.6">\1</font>', t)
    return t


class HRule(Flowable):
    def __init__(self, w=AVAIL, thick=0.6, color=RULE, pad=3):
        Flowable.__init__(self); self.w = w; self.t = thick; self.c = color; self.pad = pad
    def wrap(self, aw, ah): return (self.w, self.t + 2 * self.pad)
    def draw(self):
        self.canv.setStrokeColor(self.c); self.canv.setLineWidth(self.t)
        self.canv.line(0, self.pad, self.w, self.pad)


SEARCH = [".", "figs", "handoff", "tables"]


def resolve(p):
    if os.path.isabs(p) and os.path.exists(p): return p
    for d in SEARCH:
        c = os.path.join(d, os.path.basename(p))
        if os.path.exists(c): return c
    if os.path.exists(p): return p
    raise IOError("cannot find %r in %s" % (p, SEARCH))


def img(path, width=None, maxh=7.1 * inch):
    from reportlab.lib.utils import ImageReader
    path = resolve(path)
    iw, ih = ImageReader(path).getSize()
    w = width or AVAIL
    h = w * ih / iw
    if h > maxh:
        h = maxh; w = h * iw / ih
    return Image(path, width=w, height=h)


def csv_table(path, widths=None, fontsize=7.5, align=None):
    rows = list(csv.reader(open(resolve(path))))
    return make_table(rows, widths, fontsize, align)


def make_table(rows, widths=None, fontsize=7.5, align=None):
    tc = ParagraphStyle("tc", parent=TCELL, fontSize=fontsize, leading=fontsize * 1.24)
    th = ParagraphStyle("th", parent=THEAD, fontSize=fontsize, leading=fontsize * 1.24)
    tcr = ParagraphStyle("tcr", parent=tc, alignment=2)
    data = []
    ncol = max(len(r) for r in rows)
    for ri, r in enumerate(rows):
        r = list(r) + [""] * (ncol - len(r))
        out = []
        for ci, cell in enumerate(r):
            st = th if ri == 0 else (tcr if (align and ci < len(align) and align[ci] == "r") else tc)
            out.append(Paragraph(inline(str(cell)), st))
        data.append(out)
    if widths is None:
        raw = [max(stringWidth(str(rows[ri][ci]) if ci < len(rows[ri]) else "", "Helvetica", fontsize)
                   for ri in range(len(rows))) for ci in range(ncol)]
        # no column may be narrower than its longest unbreakable word
        word = [max([stringWidth(w, "Helvetica-Bold", fontsize)
                     for ri in range(len(rows))
                     for w in str(rows[ri][ci] if ci < len(rows[ri]) else "").split()] or [0])
                for ci in range(ncol)]
        tot = sum(raw) or 1
        widths = [max(0.42 * inch, word[ci] + 9.0, AVAIL * raw[ci] / tot)
                  for ci in range(ncol)]
        s = sum(widths)
        if s > AVAIL:  # shrink only the columns that have slack above their word floor
            slack = [w - max(0.42 * inch, word[ci] + 9.0) for ci, w in enumerate(widths)]
            over, ts = s - AVAIL, sum(slack) or 1
            widths = [w - over * slack[ci] / ts for ci, w in enumerate(widths)]
        s = sum(widths)
        widths = [w * AVAIL / s for w in widths]
    t = Table(data, colWidths=widths, repeatRows=1, hAlign="LEFT")
    t.setStyle(TableStyle([
        ("BACKGROUND", (0, 0), (-1, 0), ACC),
        ("TEXTCOLOR", (0, 0), (-1, 0), colors.white),
        ("VALIGN", (0, 0), (-1, -1), "TOP"),
        ("TOPPADDING", (0, 0), (-1, -1), 2.4), ("BOTTOMPADDING", (0, 0), (-1, -1), 2.4),
        ("LEFTPADDING", (0, 0), (-1, -1), 3.6), ("RIGHTPADDING", (0, 0), (-1, -1), 3.6),
        ("LINEBELOW", (0, 0), (-1, -2), 0.3, RULE),
        ("LINEBELOW", (0, -1), (-1, -1), 0.7, colors.HexColor("#8b97a6")),
        ("ROWBACKGROUNDS", (0, 1), (-1, -1), [colors.white, colors.HexColor("#f6f8fa")]),
    ]))
    return t


def boxed(flows, bg=BOXBG, ed=BOXED):
    t = Table([[flows]], colWidths=[AVAIL], hAlign="LEFT")
    t.setStyle(TableStyle([("BACKGROUND", (0, 0), (-1, -1), bg),
                           ("BOX", (0, 0), (-1, -1), 0.7, ed),
                           ("LEFTPADDING", (0, 0), (-1, -1), 8), ("RIGHTPADDING", (0, 0), (-1, -1), 8),
                           ("TOPPADDING", (0, 0), (-1, -1), 7), ("BOTTOMPADDING", (0, 0), (-1, -1), 6)]))
    return t


class Doc(BaseDocTemplate):
    def __init__(self, fn, title):
        BaseDocTemplate.__init__(self, fn, pagesize=letter, leftMargin=LM, rightMargin=RM,
                                 topMargin=TM, bottomMargin=BM, title=title,
                                 author="Claude Science", subject=title)
        fr = Frame(LM, BM, AVAIL, PH - TM - BM, id="n", leftPadding=0, rightPadding=0,
                   topPadding=0, bottomPadding=0)
        self.addPageTemplates([PageTemplate(id="plain", frames=[fr], onPage=self.deco)])
        self.shorttitle = title
        self._first = True

    def deco(self, c, d):
        c.saveState()
        pn = c.getPageNumber()
        if pn > 1:
            c.setFont("Helvetica", 7.2); c.setFillColor(MUTED)
            c.drawString(LM, PH - TM + 17, self.shorttitle)
            c.setStrokeColor(RULE); c.setLineWidth(0.4)
            c.line(LM, PH - TM + 12, PW - RM, PH - TM + 12)
            c.setFont("Helvetica", 8.2)
            c.drawRightString(PW - RM, BM - 22, "%d" % pn)
        c.restoreState()

    def afterFlowable(self, flow):
        if hasattr(flow, "style") and getattr(flow, "_tocLevel", None) is not None:
            self.notify("TOCEntry", (flow._tocLevel, flow.getPlainText(), self.page))


def H(text, style, level=None):
    p = Paragraph(inline(text), style)
    p._tocLevel = level
    return p


def parse(md, figdir="."):
    """Parse markup text into a flowable list."""
    story = []
    lines = md.split("\n")
    i = 0
    para, bullets, nums, table, boxbuf, boxkind, boxtitle = [], [], [], [], [], None, None

    def flush_para():
        nonlocal para
        if para:
            tgt = boxbuf if boxkind else story
            tgt.append(Paragraph(inline(" ".join(para)), BOXBODY if boxkind else BODY))
            para = []

    def flush_bul():
        nonlocal bullets
        if bullets:
            tgt = boxbuf if boxkind else story
            for b in bullets:
                tgt.append(Paragraph(inline(b), BULL, bulletText="\u2022"))
            bullets = []

    def flush_num():
        nonlocal nums
        if nums:
            tgt = boxbuf if boxkind else story
            for n, b in nums:
                tgt.append(Paragraph(inline(b), NUMB, bulletText=n + "."))
            nums = []

    def flush_tbl():
        nonlocal table
        if table:
            hdr = table[0]
            body = [r for r in table[1:] if not re.match(r"^[-: ]+$", "".join(r))]
            tgt = boxbuf if boxkind else story
            tgt.append(Spacer(1, 3))
            tgt.append(make_table([hdr] + body))
            tgt.append(Spacer(1, 2))
            table = []

    def flush_all():
        flush_para(); flush_bul(); flush_num(); flush_tbl()

    while i < len(lines):
        ln = lines[i].rstrip()
        s = ln.strip()
        if s.startswith(":::"):
            flush_all()
            if boxkind is None:
                parts = s[3:].strip().split(" ", 1)
                boxkind = parts[0] or "note"
                boxtitle = parts[1] if len(parts) > 1 else None
                boxbuf = []
                if boxtitle:
                    boxbuf.append(Paragraph(inline(boxtitle), BOXTITLE))
            else:
                bg, ed = (HYPBG, HYPED) if boxkind in ("idea", "hypothesis", "proposal") else (BOXBG, BOXED)
                story.append(Spacer(1, 3)); story.append(boxed(boxbuf, bg, ed)); story.append(Spacer(1, 6))
                boxkind, boxbuf, boxtitle = None, [], None
            i += 1; continue
        if not s:
            flush_all(); i += 1; continue
        if s == "!PAGEBREAK":
            flush_all(); story.append(PageBreak()); i += 1; continue
        if s == "!HR":
            flush_all(); (boxbuf if boxkind else story).append(HRule()); i += 1; continue
        if s.startswith("!FIG"):
            flush_all()
            body = s[4:].strip()
            parts = [x.strip() for x in body.split("|")]
            path = os.path.join(figdir, parts[0])
            capt = parts[1] if len(parts) > 1 else ""
            wfrac = float(parts[2]) if len(parts) > 2 else 1.0
            blk = [img(path, width=AVAIL * wfrac)]
            if capt: blk.append(Paragraph(inline(capt), CAP))
            story.append(Spacer(1, 4)); story.append(KeepTogether(blk)); story.append(Spacer(1, 7))
            i += 1; continue
        if s.startswith("!TBL"):
            flush_all()
            parts = [x.strip() for x in s[4:].strip().split("|")]
            path = os.path.join(figdir, parts[0])
            capt = parts[1] if len(parts) > 1 else ""
            blk = []
            if capt: blk.append(Paragraph(inline(capt), CAP))
            blk.append(Spacer(1, 2)); blk.append(csv_table(path))
            story.append(Spacer(1, 4)); story.extend(blk); story.append(Spacer(1, 7))
            i += 1; continue
        if s.startswith("#### "):
            flush_all(); (boxbuf if boxkind else story).append(Paragraph(inline(s[5:]), H4)); i += 1; continue
        if s.startswith("### "):
            flush_all(); story.append(H(s[4:], H3, 2)); i += 1; continue
        if s.startswith("## "):
            flush_all(); story.append(H(s[3:], H2, 1)); i += 1; continue
        if s.startswith("# "):
            flush_all()
            story.append(PageBreak())
            story.append(H(s[2:], H1, 0))
            if i + 1 < len(lines) and lines[i + 1].strip().startswith("%SUB "):
                story.append(Paragraph(inline(lines[i + 1].strip()[5:]), H1SUB)); i += 1
            else:
                story.append(Spacer(1, 8))
            story.append(HRule(thick=1.1, color=ACC, pad=1)); story.append(Spacer(1, 6))
            i += 1; continue
        if s.startswith("|"):
            flush_para(); flush_bul(); flush_num()
            cells = [c.strip() for c in s.strip("|").split("|")]
            table.append(cells); i += 1; continue
        if re.match(r"^[-*] ", s):
            flush_para(); flush_num(); flush_tbl()
            bullets.append(s[2:]); i += 1
            while i < len(lines) and lines[i].startswith("  ") and lines[i].strip() and not re.match(r"^\s*[-*] ", lines[i]):
                bullets[-1] += " " + lines[i].strip(); i += 1
            continue
        m = re.match(r"^(\d+)\.\s+(.*)$", s)
        if m:
            flush_para(); flush_bul(); flush_tbl()
            nums.append((m.group(1), m.group(2))); i += 1
            while i < len(lines) and lines[i].startswith("   ") and lines[i].strip() and not re.match(r"^\s*(\d+\.|[-*])\s", lines[i]):
                nums[-1] = (nums[-1][0], nums[-1][1] + " " + lines[i].strip()); i += 1
            continue
        flush_bul(); flush_num(); flush_tbl()
        para.append(s); i += 1
    flush_all()
    return story


def titlepage(title, subtitle, meta_lines, blurb):
    st = [Spacer(1, 0.55 * inch),
          HRule(thick=2.2, color=ACC, pad=2), Spacer(1, 14),
          Paragraph(inline(title), TITLE), Spacer(1, 9),
          Paragraph(inline(subtitle), STITLE), Spacer(1, 14),
          HRule(thick=0.7, color=RULE, pad=2), Spacer(1, 12)]
    st.append(boxed([Paragraph(inline(blurb), BOXBODY)]))
    st.append(Spacer(1, 16))
    for m in meta_lines:
        st.append(Paragraph(inline(m), META))
    return st


def build(outfile, title, subtitle, meta_lines, blurb, body_md, figdir=".", toc=True):
    doc = Doc(outfile, title)
    story = titlepage(title, subtitle, meta_lines, blurb)
    if toc:
        story.append(PageBreak())
        story.append(Paragraph("Contents", H2))
        t = TableOfContents()
        t.levelStyles = [
            ParagraphStyle("t0", fontName="Helvetica-Bold", fontSize=9.6, leading=13.6,
                           spaceBefore=6, textColor=ACC),
            ParagraphStyle("t1", fontName="Helvetica", fontSize=8.8, leading=12.0, leftIndent=13,
                           textColor=INK),
            ParagraphStyle("t2", fontName="Helvetica", fontSize=8.1, leading=10.8, leftIndent=28,
                           textColor=colors.HexColor("#46505e")),
        ]
        story.append(t)
    story.extend(parse(body_md, figdir))
    doc.multiBuild(story)
    return outfile
