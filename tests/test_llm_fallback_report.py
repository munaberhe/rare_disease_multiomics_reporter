import os
from pathlib import Path

from src.rdmr.llm_report import generate_report


def test_llm_fallback_report(tmp_path, monkeypatch):
    """
    Ensure that when OPENAI_API_KEY is not set, generate_report
    produces a deterministic template-style report.
    """
    # Ensure no API key
    monkeypatch.delenv("OPENAI_API_KEY", raising=False)

    # Create tiny fake variant table
    variants_path = tmp_path / "variants.tsv"
    variants_path.write_text(
        "chrom\tpos\tref\talt\tgene\tconsequence\taf\tscore\n"
        "1\t123456\tA\tG\tBRCA1\tmissense_variant\t0.0001\t5.0\n"
    )

    # No expression file
    expr_path = tmp_path / "expr.tsv"  # doesn't exist, that's fine

    out_path = tmp_path / "report.md"

    generate_report(
        phased_variants_path=variants_path,
        expr_results_path=expr_path,
        phenotypes="short stature, developmental delay",
        out_path=out_path,
    )

    assert out_path.exists()
    text = out_path.read_text()

    assert "# Rare Disease Multi-Omics Report" in text
    assert "## Phenotypes" in text
    assert "Variant summary" in text
    assert "Expression summary" in text
    assert "This report was generated **without** an LLM" in text

