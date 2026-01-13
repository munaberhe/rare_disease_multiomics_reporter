from __future__ import annotations

import os
from pathlib import Path
from typing import Optional

import pandas as pd

try:
    from openai import OpenAI
except ImportError:
    OpenAI = None  # type: ignore[assignment]


def _load_top_variants(variants_path: Path, n: int = 5) -> pd.DataFrame:
    df = pd.read_csv(variants_path, sep="\t")
    return df.head(n)


def _load_expression_results(expression_path: Path, n: int = 10) -> Optional[pd.DataFrame]:
    if not expression_path.exists():
        return None
    df = pd.read_csv(expression_path, sep="\t")
    return df.head(n)


def _format_variant_summary(df: pd.DataFrame) -> str:
    lines = []
    for _, row in df.iterrows():
        lines.append(
            f"- {row['gene']} {row['chrom']}:{row['pos']} {row['ref']}>{row['alt']} "
            f"({row['consequence']}, af={row['af']}, score={row['score']})"
        )
    return "\n".join(lines)


def _format_expression_summary(df: Optional[pd.DataFrame]) -> str:
    if df is None or df.empty:
        return "No expression results available."

    # Prefer 'log2FC' if present, otherwise fall back to DESeq2's 'log2FoldChange'
    log2fc_col = None
    if "log2FC" in df.columns:
        log2fc_col = "log2FC"
    elif "log2FoldChange" in df.columns:
        log2fc_col = "log2FoldChange"

    lines = []
    for _, row in df.iterrows():
        if log2fc_col is not None:
            lines.append(
                f"- {row['gene']}: baseMean={row['baseMean']:.1f}, {log2fc_col}={row[log2fc_col]:.2f}"
            )
        else:
            lines.append(f"- {row['gene']}: baseMean={row['baseMean']:.1f}")
    return "\n".join(lines)


def _try_llm_report(
    phenotypes: str,
    variants_summary: str,
    expression_summary: str,
) -> Optional[str]:
    api_key = os.getenv("OPENAI_API_KEY")
    if not api_key or OpenAI is None:
        return None

    client = OpenAI(api_key=api_key)

    system_msg = (
        "You are an expert rare disease genomic analyst. "
        "Given variant summaries, expression summaries, and patient phenotypes, "
        "you write concise, structured research-style reports. "
        "Do NOT give medical advice; frame everything as preliminary, research-only analysis."
    )

    user_msg = f"""
Phenotype description:
{phenotypes}

Top candidate variants:
{variants_summary}

Expression / pathway summary (toy results):
{expression_summary}

Write a concise markdown report with sections:
1. Findings
2. Supporting Evidence
3. Limitations (research only)

Keep it to 2–4 short paragraphs overall.
"""

    response = client.chat.completions.create(
        model="gpt-4o-mini",
        messages=[
            {"role": "system", "content": system_msg},
            {"role": "user", "content": user_msg},
        ],
        temperature=0.1,
    )

    content = response.choices[0].message.content
    if not content:
        return None
    return content.strip()


def generate_report(
    phenotypes: str,
    variants_path: Path,
    expression_path: Path,
    out_path: Path,
) -> None:
    """
    Generate a markdown report combining phenotypes, variants, and expression.

    If OPENAI_API_KEY is set and the openai client is available, use the LLM
    to write the report. Otherwise, fall back to a simple template.
    """
    variants_df = _load_top_variants(variants_path, n=5)
    expr_df = _load_expression_results(expression_path, n=10)

    variants_summary = _format_variant_summary(variants_df)
    expression_summary = _format_expression_summary(expr_df)

    # Try using the LLM if configured
    report_md = _try_llm_report(phenotypes, variants_summary, expression_summary)

    if report_md is None:
        # Fallback template if no LLM
        report_md = f"""# Rare Disease Multi-Omics Report (Template)

## Phenotype Summary

{phenotypes}

## Top Candidate Variants (toy scores)

{variants_summary}

## Expression Summary (DESeq2-based, toy)

{expression_summary}

## Notes

- This report was generated without an LLM (no OPENAI_API_KEY configured), or the LLM call failed.
- All results are synthetic and for demonstration only.
"""

    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(report_md, encoding="utf-8")

