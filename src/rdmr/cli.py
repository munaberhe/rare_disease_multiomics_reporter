import argparse
import subprocess
from pathlib import Path

from .variant_scoring import load_variants_table, score_variants, save_scored_variants
from .llm_report import generate_report


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


def run_r_script(cmd: list[str]) -> None:
    """Helper to run an Rscript command and show its output."""
    try:
        result = subprocess.run(
            cmd,
            check=True,
            capture_output=True,
            text=True,
        )
        if result.stdout.strip():
            print(result.stdout)
        if result.stderr.strip():
            print("[R stderr]")
            print(result.stderr)
    except subprocess.CalledProcessError as e:
        print("Error running R script:")
        print(e.stdout)
        print(e.stderr)
        raise


def main() -> None:
    args = parse_args()
    print("=== Rare Disease Multi-Omics Reporter (MVP) ===")
    print(f"VCF/variants path : {args.vcf}")
    print(f"Counts path       : {args.counts}")
    print(f"Phenotypes        : {args.phenotypes}")
    print(f"Report outpath    : {args.out}")

    # 1) Load and score variants (Python side)
    print("\n[1/4] Loading and scoring variants...")
    variants = load_variants_table(args.vcf)
    scored = score_variants(variants)
    scored_out = Path("results/variants/scored_variants.tsv")
    save_scored_variants(scored, scored_out)
    print(f"Saved scored variants to {scored_out}")

    # 2) Run DESeq2 in R
    print("\n[2/4] Running DESeq2 (R)...")
    expr_out_dir = Path("results/expression")
    expr_out_dir.mkdir(parents=True, exist_ok=True)
    deseq_cmd = [
        "Rscript",
        "R/01_deseq2_analysis.R",
        str(args.counts),
        str(expr_out_dir),
    ]
    run_r_script(deseq_cmd)
    deseq_results_path = expr_out_dir / "deseq2_results.tsv"

    # 3) Run GO enrichment (clusterProfiler)
    print("\n[3/4] Running GO enrichment (R, clusterProfiler)...")
    enrich_cmd = [
        "Rscript",
        "R/02_pathway_analysis.R",
        str(deseq_results_path),
        str(expr_out_dir),
    ]
    run_r_script(enrich_cmd)

    # 4) LLM / template report generation
    print("\n[4/4] Generating report (LLM or fallback template)...")
    generate_report(
        phenotypes=args.phenotypes,
        variants_path=scored_out,
        expression_path=deseq_results_path,
        out_path=args.out,
    )
    print(f"Report written to {args.out}")
    print("Pipeline complete.")


if __name__ == "__main__":
    main()

