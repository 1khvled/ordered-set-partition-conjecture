#!/usr/bin/env python3
"""Build a readable PDF from paper.md without a TeX installation."""
from __future__ import annotations

import html
import re
from pathlib import Path

from reportlab.lib import colors
from reportlab.lib.enums import TA_CENTER
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
from reportlab.lib.units import mm
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.platypus import Paragraph, Preformatted, SimpleDocTemplate, Spacer

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "ordered-big-blocks-paper.pdf"

# Arial is present on Windows and covers the Greek letters used in the paper.
font_dir = Path("C:/Windows/Fonts")
for name, filename in [
    ("PaperSans", "arial.ttf"),
    ("PaperSans-Bold", "arialbd.ttf"),
    ("PaperSans-Italic", "ariali.ttf"),
    ("PaperSans-BoldItalic", "arialbi.ttf"),
]:
    path = font_dir / filename
    if path.exists():
        pdfmetrics.registerFont(TTFont(name, str(path)))
pdfmetrics.registerFontFamily(
    "PaperSans",
    normal="PaperSans",
    bold="PaperSans-Bold",
    italic="PaperSans-Italic",
    boldItalic="PaperSans-BoldItalic",
)

SUBS = str.maketrans({"₀": "0", "₁": "1", "₂": "2", "₃": "3", "₄": "4",
                      "₅": "5", "₆": "6", "₇": "7", "₈": "8", "₉": "9",
                      "₊": "+", "₋": "-", "₌": "=", "₍": "(", "₎": ")",
                      "ᵣ": "_r", "ₘ": "_m"})
SUPS = str.maketrans({"⁰": "^0", "¹": "^1", "²": "^2", "³": "^3", "⁴": "^4",
                      "⁵": "^5", "⁶": "^6", "⁷": "^7", "⁸": "^8", "⁹": "^9",
                      "ⁿ": "^n", "ʳ": "^r", "ᶻ": "^z", "ᵐ": "^m", "ᵢ": "_i",
                      "ⱼ": "_j", "ₖ": "_k", "ₙ": "_n", "ᵩ": "q"})


def normalize(text: str) -> str:
    return (text.translate(SUBS).translate(SUPS)
            .replace("⋆", "*").replace("⋯", "...").replace("∣", "|")
            .replace("≥", ">=").replace("≤", "<=")
            .replace("⌊", "floor(").replace("⌋", ")"))


def inline(text: str) -> str:
    text = normalize(html.escape(text))
    text = re.sub(r"`([^`]+)`", r'<font name="Courier">\1</font>', text)
    text = re.sub(r"\*\*([^*]+)\*\*", r"<b>\1</b>", text)
    text = re.sub(r"(?<!\*)\*([^*]+)\*(?!\*)", r"<i>\1</i>", text)
    return text


def page_number(canvas, doc):
    canvas.saveState()
    canvas.setFont("PaperSans", 8)
    canvas.setFillColor(colors.HexColor("#666666"))
    canvas.drawRightString(A4[0] - 18 * mm, 12 * mm, f"{doc.page}")
    canvas.restoreState()


def main():
    text = (ROOT / "paper.md").read_text(encoding="utf-8")
    styles = getSampleStyleSheet()
    body = ParagraphStyle("Body", parent=styles["BodyText"], fontName="PaperSans", fontSize=10.2,
                          leading=14.2, spaceAfter=7, alignment=0)
    title = ParagraphStyle("Title", parent=styles["Title"], fontName="PaperSans-Bold", fontSize=21,
                           leading=25, alignment=TA_CENTER, textColor=colors.HexColor("#172033"),
                           spaceAfter=8)
    subtitle = ParagraphStyle("Subtitle", parent=title, fontSize=13, leading=17, spaceAfter=14)
    h1 = ParagraphStyle("H1", parent=styles["Heading1"], fontName="PaperSans-Bold", fontSize=15,
                        leading=19, textColor=colors.HexColor("#172033"), spaceBefore=13, spaceAfter=7)
    h2 = ParagraphStyle("H2", parent=styles["Heading2"], fontName="PaperSans-Bold", fontSize=12,
                        leading=16, textColor=colors.HexColor("#243b64"), spaceBefore=10, spaceAfter=6)
    formula = ParagraphStyle("Formula", parent=body, leftIndent=9, rightIndent=9, spaceBefore=3,
                             spaceAfter=8, backColor=colors.HexColor("#f3f5f8"), borderPadding=7)
    code = ParagraphStyle("Code", parent=styles["Code"], fontName="Courier", fontSize=8.5,
                          leading=11.5, leftIndent=8, backColor=colors.HexColor("#f5f5f5"),
                          borderPadding=7, spaceAfter=9)
    story = []
    para: list[str] = []
    quote: list[str] = []
    code_lines: list[str] | None = None

    def flush_para():
        if para:
            story.append(Paragraph(inline(" ".join(x.strip() for x in para)), body))
            para.clear()

    def flush_quote():
        if quote:
            content = "<br/>".join(inline(x.strip().lstrip(">").strip()) for x in quote)
            story.append(Paragraph(content, formula))
            quote.clear()

    for raw in text.splitlines():
        line = raw.rstrip()
        if code_lines is not None:
            if line.strip().startswith("```"):
                story.append(Preformatted("\n".join(code_lines), code))
                code_lines = None
            else:
                code_lines.append(line)
            continue
        if line.strip().startswith("```"):
            flush_para(); flush_quote(); code_lines = []
            continue
        if line.startswith(">"):
            flush_para()
            quote.append(line)
            continue
        flush_quote()
        if not line.strip():
            flush_para()
            continue
        if line.startswith("# "):
            flush_para(); story.append(Paragraph(inline(line[2:]), title))
        elif line.startswith("## "):
            flush_para(); story.append(Paragraph(inline(line[3:]), h1))
        elif line.startswith("### "):
            flush_para(); story.append(Paragraph(inline(line[4:]), h2))
        else:
            para.append(line)
    flush_para(); flush_quote()
    if code_lines is not None:
        story.append(Preformatted("\n".join(code_lines), code))

    doc = SimpleDocTemplate(str(OUT), pagesize=A4, rightMargin=18 * mm, leftMargin=18 * mm,
                            topMargin=17 * mm, bottomMargin=19 * mm,
                            title="Ordered Set Partitions Counted by Big Blocks",
                            author="Khaled Abdelli")
    doc.build(story, onFirstPage=page_number, onLaterPages=page_number)
    print(OUT)


if __name__ == "__main__":
    main()
