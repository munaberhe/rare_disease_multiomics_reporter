import argparse
import subprocess
from pathlib import Path

from .variant_scoring import load_variants_table, score_variants, save_scored_variants
from .llm_report import generate_report  # we'll create this next


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Rare disease multi-omics reporter (toy pipeline).",
    )
    parser.add_argument(
        "--vcf",
        type=Path,
        required=True,
        help="Path to input variant table (TSV for now, VCF-like).",
    )
    parser.add_argument(
        "--counts",
        type=Path,
        required=True,
        help="Path to RNA-seq counts matrix (TSV).",
    )
    parser.add_argument(
        "--phenotypes",
        type=str,
        required=True,
        help="Free-text phenotype description (e.g. HPO-like).",
    )
    parser.add_argument(
        "--out",
        type=Path,
        required=True,
        help="Output report path (Markdown).",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    print("=== Rare Disease Multi-Omics Reporter (MVP) ===")
    print(f"VCF/variants path : {args.vcf}")
    print(f"Counts path       : {args.counts}")
    print(f"Phenotypes        : {args.phenotypes}")
    print(f"Report outpath    : {args.out}")

    # 1) Load and score variants (Python side)
    print("\n[1/3] Loading and scoring variants...")
    variants = load_variants_table(args.vcf)
    scored = score_variants(variants)
    scored_out = Path("results/variants/scored_variants.tsv")
    save_scored_variants(scored, scored_out)
    print(f"Saved scored variants to {scored_out}")

    # 2) Run R expression analysis placeholder
    print("\n[2/3] Running R expression analysis (placeholder)...")
    expr_out_dir = Path("results/expression")
    expr_out_dir.mkdir(parents=True, exist_ok=True)
    cmd = [
        "Rscript",
        "R/01_deseq2_analysis.R",
        str(args.counts),
        str(expr_out_dir),
    ]
    try:
        result = subprocess.run(
            cmd,
            check=True,
            capture_output=True,
            text=True,
        )
        # Show R's stdout so you can see what happened
        print(result.stdout)
        if result.stderr.strip():
            print("[R stderr]")
            print(result.stderr)
    except subprocess.CalledProcessError as e:
        print("Error running R script:")
        print(e.stdout)
        print(e.stderr)
        raise

    expr_results_path = expr_out_dir / "expression_results.tsv"

    # 3) LLM / template report generation
    print("\n[3/3] Generating report (LLM or fallback template)...")
    generate_report(
        phenotypes=args.phenotypes,
        variants_path=scored_out,
        expression_path=expr_results_path,
        out_path=args.out,
    )
    print(f"Report written to {args.out}")
    print("Pipeline complete.")


if __name__ == "__main__":
    main()

