# Rare Disease Multi-Omics Reporter

## Overview

This project is a toy **rare disease multi-omics reporting pipeline** that combines:

- Variant prioritisation (Python)
- RNA-seq differential expression and pathway enrichment (R / Bioconductor)
- An optional LLM-generated clinical-style summary (Python, when `OPENAI_API_KEY` is set)

The goal is to mimic how a rare disease genomics team might integrate variants, gene expression, and phenotypes into a structured report. All data are synthetic and this is **not** a clinical tool. It is intended as a portfolio project to demonstrate multi-language workflow design (Python + R) and LLM integration for reporting.

---

## Dataset

All input data are small, synthetic files committed in the repository so that the pipeline is fully self-contained.

### Variant table

- Path: `data/example_variants.tsv`
- Columns (toy example):
  - `chrom`, `pos`, `ref`, `alt`
  - `gene` – gene symbol
  - `consequence` – e.g. `missense_variant`, `stop_gained`, `frameshift_variant`
  - `af` – simulated allele frequency

This table is used to demonstrate simple, rule-based variant scoring for rare disease triage.

### RNA-seq counts

- Path: `data/counts.tsv`
- Structure:
  - One row per gene
  - Columns (example):

        gene    control1  control2  treated1  treated2
        GENE1   10        12        20        22
        GENE2   50        55        30        28
        GENE3   5         6         100       95
        GENE4   100       98        80        85
        GENE5   200       210       300       290

This small count matrix is sufficient to run DESeq2 and clusterProfiler and illustrate expression and pathway analysis.

---

## Methods

The pipeline has three main components: variant scoring (Python), expression + pathway analysis (R), and multi-omics report generation (Python with optional LLM).

### 1. Variant scoring (Python)

Key files:
- `src/rdmr/variant_scoring.py`
- orchestrated by `run_report.py`

Method summary:

1. Load `data/example_variants.tsv` into a pandas DataFrame.
2. Apply a simple, transparent scoring scheme:
   - High-impact consequences (e.g. `stop_gained`, `frameshift_variant`) receive more points.
   - Moderate-impact (`missense_variant`) receive intermediate points.
   - Synonymous or non-coding variants receive low or zero impact score.
   - Rare variants (e.g. `af < 0.001`) receive additional points.
   - Low-frequency variants (e.g. `af < 0.01`) receive some points.
   - Common variants contribute little or nothing.
3. Sum these contributions into a total variant score.
4. Sort the variants by score (descending).
5. Write the scored table to `results/variants/scored_variants.tsv`.

This approximates a basic variant triage step used in rare disease genomics pipelines.

### 2. Differential expression with DESeq2 (R)

Key file:
- `R/01_deseq2_analysis.R`

Method summary:

1. Read `data/counts.tsv` into R.
2. Construct a `DESeqDataSet` with a simple design (e.g. condition: control vs treated).
3. Estimate size factors and dispersions, using a strategy appropriate for very small toy datasets.
4. Fit the negative binomial model and perform Wald tests via DESeq2.
5. Produce a results table with columns such as:
   - `baseMean`
   - `log2FoldChange`
   - `lfcSE`
   - `pvalue`
   - `padj`
6. Save outputs:
   - Differential expression results: `results/expression/deseq2_results.tsv`
   - Volcano plot: `results/expression/volcano_deseq2.png`

This implements a standard DESeq2 analysis workflow on simulated RNA-seq counts.

### 3. GO / pathway enrichment with clusterProfiler (R)

Key file:
- `R/02_pathway_analysis.R`

Method summary:

1. Read `results/expression/deseq2_results.tsv`.
2. Select significantly differentially expressed genes, for example:
   - adjusted p-value (`padj`) < 0.05
   - absolute log2 fold change (`|log2FoldChange|`) > 1
3. Map gene identifiers to Entrez IDs using:
   - `clusterProfiler::bitr(...)`
   - `org.Hs.eg.db`
4. Perform GO Biological Process enrichment:
   - `clusterProfiler::enrichGO(..., ont = "BP")`
5. Save outputs:
   - GO enrichment table (e.g.): `results/expression/go_enrichment_airway_BP.csv`
   - Bar plot of enriched GO terms: `results/expression/go_bp_barplot_airway.png`

This step models the typical downstream interpretation of a DESeq2 result using GO enrichment.

### 4. Multi-omics report generation with optional LLM (Python)

Key files:
- `src/rdmr/llm_report.py`
- `run_report.py`

Method summary:

1. Load:
   - Scored variants from `results/variants/scored_variants.tsv`.
   - Optional DE results from `results/expression/deseq2_results.tsv` (if generated).
2. Build:
   - A variant summary: top variants with gene, position, consequence, allele frequency and score.
   - An expression summary: for example, top up- and down-regulated genes based on DESeq2 columns.
3. If `OPENAI_API_KEY` is set:
   - Construct a structured prompt that includes:
     - Patient phenotypes (passed on the command line).
     - Top candidate variants and scores.
     - Expression summary.
   - Call an OpenAI model (e.g. `gpt-4o-mini`) to generate a markdown report with headings such as:
     - Phenotype
     - Genomic findings
     - Expression findings
     - Interpretation / Limitations
   - Save an LLM-generated report with a heading such as:
     - `# Rare Disease Multi-Omics Report (LLM-generated)`
4. If `OPENAI_API_KEY` is not set or the LLM call fails:
   - Generate a deterministic, template-based report with sections:
     - Phenotypes
     - Variant summary
     - Expression summary
     - Notes
   - Explicitly state that no LLM was used and that the data are synthetic and non-clinical.
5. Write the final report to a path such as `results/reports/patient_001_report.md`.

This demonstrates how LLM-based narrative reporting can sit on top of traditional bioinformatics outputs with a robust fallback.

---

## Setup

### Python environment

From the project root:

    python -m venv .venv
    source .venv/bin/activate          # Windows PowerShell: .venv\Scripts\Activate.ps1

Install required Python packages (adapt if you have a `requirements.txt`):

    pip install pandas matplotlib openai

If you do have a `requirements.txt` file, you can instead run:

    pip install -r requirements.txt

### R / Bioconductor environment

In an R session:

    install.packages("BiocManager")
    BiocManager::install(c("DESeq2", "clusterProfiler", "org.Hs.eg.db"))
    install.packages(c("ggplot2", "dplyr", "readr"))

Make sure R can find these packages when running the scripts from the command line.

---

## How to Run

### 1. Run DESeq2 and GO enrichment (R)

From the project root:

    Rscript R/01_deseq2_analysis.R
    Rscript R/02_pathway_analysis.R

This should produce:

- `results/expression/deseq2_results.tsv`
- `results/expression/volcano_deseq2.png`
- `results/expression/go_enrichment_airway_BP.csv` (name may vary slightly)
- `results/expression/go_bp_barplot_airway.png`

### 2. Score variants and generate the multi-omics report (Python)

Activate your Python virtual environment, then run:

    python run_report.py \
      --vcf data/example_variants.tsv \
      --counts data/counts.tsv \
      --phenotypes "short stature, developmental delay" \
      --out results/reports/patient_001_report.md

This will:

1. Run variant scoring and write `results/variants/scored_variants.tsv`.
2. Use DESeq2 (and GO) results from `results/expression/` if they are present.
3. Generate a markdown report integrating phenotypes, variants and expression.
4. Use an LLM when `OPENAI_API_KEY` is set, or a deterministic template otherwise.

---

## Results

Running the full pipeline produces a set of outputs that resemble a small rare disease multi-omics analysis.

### Variant-level outputs

- `results/variants/scored_variants.tsv`
  - Columns such as:
    - `chrom`, `pos`, `ref`, `alt`
    - `gene`, `consequence`
    - `af`, `score`
  - Variants are sorted by score (highest first), highlighting the most promising candidates under the toy rules.

### Expression and pathway outputs

- `results/expression/deseq2_results.tsv`
  - DESeq2-style results with log2 fold changes, standard errors, p-values, and adjusted p-values for each gene.
- `results/expression/volcano_deseq2.png`
  - A volcano plot summarising differential expression.
- `results/expression/go_enrichment_airway_BP.csv`
  - GO Biological Process terms, adjusted p-values, and associated genes for enriched pathways (where enough genes map).
- `results/expression/go_bp_barplot_airway.png`
  - A bar plot of representative enriched GO terms.

### Report outputs

- `results/reports/patient_001_report.md`
  - In LLM mode (API key set and call successful):
    - A structured markdown report generated by the model, with clinical-style narrative.
  - In fallback mode (no API key or LLM error):
    - A deterministic template report summarising phenotypes, variants, expression and limitations.

These outputs illustrate the full journey from small synthetic tables to combined numerical/statistical results and a human-readable multi-omics report.

---

## Discussion

This project demonstrates how to combine Python, R and modern language models into a coherent rare disease multi-omics workflow:

- Python is used for:
  - Transparent, rules-based variant scoring.
  - Orchestration of inputs and outputs.
  - Markdown report generation.
  - Optional LLM integration with a clearly defined fallback.

- R / Bioconductor is used for:
  - Differential expression analysis via DESeq2.
  - GO enrichment analysis via clusterProfiler and org.Hs.eg.db.
  - Visual summaries such as volcano plots and GO bar plots.

- LLM integration:
  - Shows how an LLM can provide clinician-style summaries on top of structured genomics outputs.
  - Keeps deterministic, auditable numerical analysis separate from generative narrative text.
  - Emphasises the need for clear limitations and disclaimers in any clinical-adjacent setting.

From a portfolio perspective, this repository highlights:

- Experience with both Python and R in bioinformatics.
- Familiarity with standard tools like DESeq2 and clusterProfiler.
- Understanding of variant scoring, differential expression, and pathway analysis in rare disease genomics.
- The ability to design a small, testable pipeline that could be extended or productionised (e.g. with containers, CI, and real clinical data).

Potential future extensions include:

- Replacing toy TSVs with real VCFs and RNA-seq count matrices.
- Adding richer variant annotations (ClinVar, in silico predictive scores, constraint metrics).
- Integrating additional omics layers (CNVs, methylation, proteomics).
- Containerising the pipeline and adding continuous integration tests on a tiny dataset.
- Expanding unit tests to validate R-side outputs and cross-check consistency between R and Python layers.

Even in toy form, this project resembles workflows used in rare disease and translational genomics, making it a strong demonstration of multi-omics analysis, reporting and LLM integration for bioinformatics and pharma roles.

