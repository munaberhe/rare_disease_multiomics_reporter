from pathlib import Path
from typing import Literal

import pandas as pd


Consequence = Literal[
    "stop_gained",
    "frameshift_variant",
    "splice_acceptor_variant",
    "splice_donor_variant",
    "missense_variant",
    "inframe_deletion",
    "inframe_insertion",
    "synonymous_variant",
    "intron_variant",
    "other",
]


def load_variants_table(path: str | Path) -> pd.DataFrame:
    """
    Load a simple TSV table of variants.

    Expected columns:
      - chrom
      - pos
      - ref
      - alt
      - gene
      - consequence
      - af  (allele frequency, float between 0 and 1)

    For this MVP we assume the file is well-formed.
    """
    path = Path(path)
    df = pd.read_csv(path, sep="\t")
    required_cols = {"chrom", "pos", "ref", "alt", "gene", "consequence", "af"}
    missing = required_cols - set(df.columns)
    if missing:
        raise ValueError(f"Missing required columns in variants table: {missing}")
    return df


def rule_based_score_row(row: pd.Series) -> float:
    """
    Very simple rule-based pathogenicity score.

    Higher = more suspicious. This is NOT clinical, just illustrative.

    Heuristics:
      - High-impact consequence (stop_gained, frameshift, splice) => +3
      - Missense / inframe => +2
      - Synonymous / intron / other => +0
      - Rarer variants get more points:
          af < 0.001  => +3
          af < 0.01   => +2
          af < 0.05   => +1
    """
    consequence = str(row["consequence"]).lower()
    af = float(row["af"])

    score = 0.0

    # Consequence impact
    high_impact = {
        "stop_gained",
        "frameshift_variant",
        "splice_acceptor_variant",
        "splice_donor_variant",
    }
    medium_impact = {
        "missense_variant",
        "inframe_deletion",
        "inframe_insertion",
    }

    if consequence in high_impact:
        score += 3.0
    elif consequence in medium_impact:
        score += 2.0
    else:
        score += 0.0  # synonymous/intron/other

    # Allele frequency
    if af < 0.001:
        score += 3.0
    elif af < 0.01:
        score += 2.0
    elif af < 0.05:
        score += 1.0

    return score


def score_variants(df: pd.DataFrame) -> pd.DataFrame:
    """
    Add a 'score' column and return variants sorted from most to least suspicious.
    """
    df = df.copy()
    df["score"] = df.apply(rule_based_score_row, axis=1)
    df = df.sort_values("score", ascending=False).reset_index(drop=True)
    return df


def save_scored_variants(df: pd.DataFrame, out_path: str | Path) -> None:
    """
    Save scored variants as TSV.
    """
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_path, sep="\t", index=False)

