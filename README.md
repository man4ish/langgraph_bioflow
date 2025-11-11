# LangGraph BioFlow

**AI-Orchestrated Bioinformatics Workflow Supervisor (Nextflow Edition)**

---

## Overview

**LangGraph BioFlow** is an experimental framework that integrates **LangGraph** with **Nextflow** to create an intelligent, adaptive workflow supervisor for bioinformatics pipelines.

The system leverages **AI-driven decision nodes** to interpret QC metrics, control pipeline reruns, and manage downstream analyses such as alignment, QC, imputation, and GWAS.

The goal is to evolve into a **smart AI supervisor** that can:

* Monitor and control computational pipelines
* Evaluate data quality automatically
* Make dynamic rerun or stop decisions
* Generate interpretable summaries for researchers

---

## Current Features

* Node-based orchestration using **LangGraph**
* Workflow execution via **Nextflow**
* QC-based decision logic with automatic retries
* Configurable recursion/attempt limits to prevent runaway loops
* Modular structure for replacing mock nodes with real bioinformatics stages

---

## Architecture

```
Input → Alignment → QC → Decision → (Re-run if FAIL) → End
```

Each stage is represented by a **LangGraph node**, which can trigger a **Nextflow process** or another pipeline component.
The **Decision Node** interprets QC results from `qc_summary.txt` and adaptively decides whether to re-run or proceed.

---

## Example Workflow

A minimal demonstration of an adaptive QC workflow:

```
main.py → runs Nextflow pipeline (qc_pipeline.nf)
             ↓
          QC summary → PASS → END
                       FAIL → re-run (up to 3 attempts)
```

---

## Project Goals

* Extend to full bioinformatics pipeline (alignment → QC → GWAS → reporting)
* Replace rule-based QC decisions with **LLM evaluation** of QC metrics
* Integrate **LangGraph memory/state tracking** for pipeline monitoring
* Add **Nextflow event hooks** for real-time workflow introspection
* Build a **FastAPI dashboard** for live visualization of LangGraph state

---

## Project Structure

```
langgraph_bioflow/
├── main.py                   # LangGraph controller integrating Nextflow
├── workflows/
│   └── qc_pipeline.nf        # Example Nextflow pipeline
├── requirements.txt
├── config/
│   └── settings.yaml         # (Optional future config)
├── utils/
│   ├── logger.py
│   └── helpers.py
└── README.md
```

---

## Quick Start

### 1. Clone and Install

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

### 3. Run the Demo

```bash
python main.py
```

Expected output (for mock workflow):

```
Running Nextflow workflow (attempt 1)...
✅ QC summary generated: work/.../qc_summary.txt
QC Decision: FAIL
QC failed on attempt 1. Re-running with adjusted parameters (next: 2).
...
QC Decision: PASS
QC passed. Proceeding to END.
Final State: {'attempt': 2, 'decision': 'PASS'}
```

---

## Future Enhancements

* [ ] Integrate **real QC metric parsers** (FastQC, MultiQC)
* [ ] Add **LLM reasoning** for QC interpretation
* [ ] Introduce **Nextflow monitoring hooks** via API
* [ ] Use **SQLite** for persistent LangGraph state
* [ ] Support **parallel task orchestration**
* [ ] Build **FastAPI or Streamlit dashboard** for live pipeline visualization

