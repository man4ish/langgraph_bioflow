# LangGraph BioFlow

**AI-Orchestrated Bioinformatics Workflow Supervisor (Nextflow Edition)**

---

## Overview

**LangGraph BioFlow** is an experimental orchestration framework that integrates **LangGraph** with **Nextflow** to build an **intelligent, adaptive bioinformatics workflow manager**.

It introduces **AI-assisted decision-making** into traditional bioinformatics pipelines — dynamically deciding whether to proceed, rerun, or halt based on QC outcomes.

The goal is to evolve into a **smart pipeline supervisor** capable of:

* Monitoring and controlling multi-stage computational pipelines
* Evaluating data quality in real time
* Making dynamic rerun or stop decisions based on QC metrics
* Generating concise, interpretable summaries for researchers

---

## Current Capabilities

✅ **Sequential execution**: FastQC → Alignment → Deduplication → Variant Calling
✅ **LangGraph-based control flow** with retry limits and state memory
✅ **QC-driven branching** — automatically re-runs QC if metrics fail thresholds
✅ **Nextflow integration** for modular process execution
✅ **Scalable design** — can extend to RNA-seq, GWAS, or WGS workflows

---

## Architecture

```
Raw Reads
   ↓
FASTQC → QC Decision (LangGraph)
   ├── PASS → Exome Pipeline (SAMTOOLS_SORT → PICARD_MARKDUP → GATK_HAPLOTYPECALLER)
   └── FAIL → Re-run QC (max 3 attempts)
```

Each stage is represented as a **LangGraph node**, which calls a corresponding **Nextflow workflow**.
The **QC Decision Node** analyzes the results of `qc_pipeline.nf` and decides the next action.

---

## Example Workflow Execution

```
main.py
│
├── Runs qc_pipeline.nf (FastQC)
│      ├── PASS → proceeds to exome_pipeline.nf
│      └── FAIL → retries QC (up to 3 times)
│
└── exome_pipeline.nf → SAMTOOLS_SORT → PICARD_MARKDUP → GATK_HAPLOTYPECALLER
```

Example output:

```
=== Running QC (attempt 1) ===
QC Decision: FAIL
QC failed on attempt 1. Retrying with adjusted parameters...

=== Running QC (attempt 2) ===
QC Decision: PASS
Proceeding to exome pipeline.

=== Running Exome Pipeline ===
VCF generated: results/gatk/sample.vcf

Final Workflow State:
decision='COMPLETE' step='exome' attempt=2
```

---

## Project Goals

* Extend to end-to-end workflows (Alignment → QC → Imputation → GWAS → Reporting)
* Replace rule-based QC thresholds with **LLM-driven evaluation**
* Integrate **LangGraph memory and persistence** for auditability
* Add **Nextflow monitoring hooks** for runtime introspection
* Develop a **FastAPI/Streamlit dashboard** for live graph visualization

---

## Directory Layout

```
langgraph_bioflow/
├── main.py                     # LangGraph controller orchestrating Nextflow workflows
├── workflows/
│   ├── qc_pipeline.nf          # Runs FastQC only
│   └── exome_pipeline.nf       # Runs sort → markdup → variant calling
├── modules/
│   ├── fastqc/main.nf
│   ├── samtools/main.nf
│   ├── picard/main.nf
│   └── gatk/main.nf
├── utils/
│   ├── logger.py
│   └── helpers.py
├── config/
│   └── settings.yaml
├── requirements.txt
└── README.md
```

---

## Quick Start

### 1. Clone and Install Dependencies

```bash
git clone https://github.com/yourusername/langgraph_bioflow.git
cd langgraph_bioflow
pip install -r requirements.txt
```

### 2. Install Nextflow

```bash
curl -s https://get.nextflow.io | bash
sudo mv nextflow /usr/local/bin/
nextflow -version
```

### 3. Run the Orchestrated Workflow

```bash
python main.py
```

---

## Future Enhancements

* [ ] Parse **real FastQC results** instead of random PASS/FAIL simulation
* [ ] Implement **LLM-driven QC reasoning** (LangGraph + GPT)
* [ ] Introduce **Nextflow event hooks** for progress tracking
* [ ] Store LangGraph state in **SQLite** for persistent runs
* [ ] Add **FastAPI dashboard** for interactive graph visualization
* [ ] Extend for **RNA-seq, metagenomics, and multi-omics pipelines**
