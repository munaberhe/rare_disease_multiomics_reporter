from __future__ import annotations

import json
import os
from pathlib import Path
from typing import Optional

import pandas as pd


def _format_variant_summary(variants_df: pd.DataFrame) -> str:
    """
    Turn a scored variants table into a human-readable summary.
    """
    if variants_df is None or variants_df.empty:
        return "No candidate variants were prioritized."

    lines = []
    top = variants_df.head(5)

    lines.append("Top candidate variants (sorted by score):")
    for _, row in top.iterrows():
        lines.append(
            f"- {row['gene']} {row['chrom']}:{row['pos']} "
            f"{row['ref']}>{row['alt']} "
            f"({row['consequence']}, af={row['af']}, score={row['score']})"
        )
    return "\n".join(lines)


def _format_expression_summary(expr_df: Optional[pd.DataFrame]) -> str:
    """
    Summarise gene expression / DE results if available.
    """
    if expr_df is None or expr_df.empty:
        return "No differential expression results were available."

    # Take a few top up/down genes if columns look like DESeq2 output
    cols = set(expr_df.columns)
    if {"log2FoldChange", "padj"}.issubset(cols):
        up = expr_df.sort_values("log2FoldChange", ascending=False).head(3)
        down = expr_df.sort_values("log2FoldChange", ascending=True).head(3)

        lines = ["Differential expression summary (DESeq2-style):", "", "Top up-regulated genes:"]
        for _, row in up.iterrows():
            lines.append(
                f"- {row.get('gene', row.name)}: log2FC={row['log2FoldChange']:.2f}, padj={row['padj']:.2e}"
            )

        lines.append("")
        lines.append("Top down-regulated genes:")
        for _, row in down.iterrows():
            lines.append(
                f"- {row.get('gene', row.name)}: log2FC={row['log2FoldChange']:.2f}, padj={row['padj']:.2e}"
            )

        return "\n".join(lines)

    # Fallback: just show a few gene stats
    lines = ["Expression summary (toy):"]
    sample_genes = expr_df.head(3)
    for _, row in sample_genes.iterrows():
        lines.append(str(row.to_dict()))
    return "\n".join(lines)


def _build_llm_prompt(
    phenotypes: str,
    variant_summary: str,
    expression_summary: str,
) -> str:
    """
    Construct a structured prompt for the LLM to generate a clinical-style multi-omics report.
    """
    return f"""
You are an expert clinical genomicist and rare disease specialist.

You will be given:
- Patient phenotypes
- A list of prioritized variants with simple scores
- A brief summary of differential expression results

Write a concise, structured report that:
1. Summarizes the phenotype in clinical language.
2. Integrates the variant findings and highlights the most plausible causal gene(s).
3. Comments on whether the expression data supports or contradicts the variant findings.
4. Clearly states limitations (toy data, simplified scoring) and that this is not a clinical report.

Use markdown headings (## Phenotype, ## Genomic findings, ## Expression findings, ## Interpretation / Limitations).

Patient phenotypes:
{phenotypes}

Variant summary:
{variant_summary}

Expression summary:
{expression_summary}
""".strip()


def _call_llm(prompt: str) -> str:
    """
    Call an LLM (OpenAI) if OPENAI_API_KEY is set.
    Returns the model's text, or raises an exception if something goes wrong.
    """
    api_key = os.getenv("OPENAI_API_KEY")
    if not api_key:
        raise RuntimeError("OPENAI_API_KEY not set; cannot call LLM.")

    try:
        # For OpenAI's Python SDK >= 1.x
        from openai import OpenAI

        client = OpenAI(api_key=api_key)

        response = client.chat.completions.create(
            model="gpt-4o-mini",
            messages=[
                {"role": "system", "content": "You are an expert clinical genomics report writer."},
                {"role": "user", "content": prompt},
            ],
            temperature=0.3,
        )
        return response.choices[0].message.content.strip()
    except ImportError:
        # Fallback for older SDKs (if user is using openai==0.x)
        import openai  # type: ignore

        openai.api_key = api_key
        completion = openai.ChatCompletion.create(
            model="gpt-4o-mini",
            messages=[
                {"role": "system", "content": "You are an expert clinical genomics report writer."},
                {"role": "user", "content": prompt},
            ],
            temperature=0.3,
        )
        return completion.choices[0].message["content"].strip()


def generate_report(
    phased_variants_path: Path,
    expr_results_path: Optional[Path],
    phenotypes: str,
    out_path: Path,
) -> None:
    """
    High-level function to generate a multi-omics report.

    - Reads the scored variant table.
    - Optionally reads expression / DE results.
    - Calls an LLM if OPENAI_API_KEY is available, otherwise uses a deterministic template.
    - Writes a markdown report to out_path.
    """
    # Load variants
    variants_df = pd.read_csv(phased_variants_path, sep="\t")

    # Try to load expression results (if provided and exists)
    expr_df: Optional[pd.DataFrame]
    if expr_results_path is not None and expr_results_path.exists():
        expr_df = pd.read_csv(expr_results_path, sep="\t")
    else:
        expr_df = None

    variant_summary = _format_variant_summary(variants_df)
    expression_summary = _format_expression_summary(expr_df)

    # Try LLM; if not available or fails, fall back
    use_llm = os.getenv("OPENAI_API_KEY") is not None
    if use_llm:
        prompt = _build_llm_prompt(
            phenotypes=phenotypes,
            variant_summary=variant_summary,
            expression_summary=expression_summary,
        )
        try:
            report_body = _call_llm(prompt)
            header = "# Rare Disease Multi-Omics Report (LLM-generated)\n"
            notes = (
                "\n\n> Note: This report was generated using a large language model. "
                "The underlying data are synthetic and this output is for demonstration only."
            )
            full_report = header + "\n\n" + report_body + notes
        except Exception as e:
            # Log the error in a very simple way and fall back
            fallback_header = "# Rare Disease Multi-Omics Report (LLM fallback)\n"
            fallback_body = (
                f"LLM call failed with error: {e}\n\n"
                "Below is a deterministic summary of the available data.\n\n"
                "## Phenotypes\n"
                f"{phenotypes}\n\n"
                "## Variant summary\n"
                f"{variant_summary}\n\n"
                "## Expression summary\n"
                f"{expression_summary}\n"
            )
            full_report = fallback_header + "\n\n" + fallback_body
    else:
        # Deterministic, non-LLM fallback
        full_report = (
            "# Rare Disease Multi-Omics Report (Template)\n\n"
            "## Phenotypes\n"
            f"{phenotypes}\n\n"
            "## Variant summary\n"
            f"{variant_summary}\n\n"
            "## Expression summary\n"
            f"{expression_summary}\n\n"
            "## Notes\n"
            "- This report was generated **without** an LLM (no OPENAI_API_KEY configured).\n"
            "- All results are synthetic and for demonstration only.\n"
        )

    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(full_report)
    print(f"Report written to: {out_path}")

