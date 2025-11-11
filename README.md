
# LangGraph BioFlow  
**AI-Orchestrated Bioinformatics Workflow Supervisor**

---

## Overview
**LangGraph BioFlow** is an experimental project that integrates **LangGraph** with **WDL (Workflow Description Language)** and **Docker** to build an intelligent, adaptive workflow supervisor for bioinformatics pipelines.

The system uses **LLM-driven decision nodes** to interpret QC metrics, control downstream analysis steps (alignment → QC → imputation → GWAS), and adaptively rerun or skip tasks based on real-time results.

It aims to evolve into a **smart AI supervisor** capable of:
- Monitoring computational pipelines
- Evaluating data quality automatically
- Making dynamic rerun decisions
- Generating interpretive analysis summaries

---

## Current Features
- Modular node-based workflow using **LangGraph**
- Execution of **WDL pipelines** via Cromwell + Docker
- Configurable QC thresholds (mean Q30, duplication rate, alignment rate)
- Mock nodes for each analysis stage:
  - Input registration  
  - Alignment (via WDL)  
  - QC metrics evaluation  
  - Decision logic (rerun / continue)  
  - Imputation, GWAS, and reporting  

---

## Architecture
```
Input → Alignment → QC → Decision → Imputation → GWAS → Report
```

Each stage is represented by a **LangGraph node**, which can run a **WDL task** inside Docker.  
The **Decision Node** interprets QC metrics (either rule-based or via an LLM) and conditionally controls downstream flow.

---

## Project Goals
- Integrate **Cromwell** for executing WDL pipelines inside Docker.
- Introduce **persistent state** and workflow monitoring.
- Replace rule-based QC decisions with **LLM evaluation** of metrics.
- Expand to full bioinformatics pipelines: alignment → QC → variant calling → GWAS → report generation.
- Enable human-in-the-loop approval for uncertain cases.

---

## Project Structure
```
langgraph_bioflow/
├── main.py
├── requirements.txt
├── README.md
├── config/
│   ├── workflow_config.yaml
│   └── model_config.yaml
├── nodes/
│   ├── input_manager.py
│   ├── aligner_node.py
│   ├── qc_node.py
│   ├── decision_node.py
│   ├── imputation_node.py
│   ├── gwas_node.py
│   ├── report_node.py
│   └── error_handler.py
├── utils/
│   ├── wdl_runner.py
│   ├── logger.py
│   └── shell_runner.py
└── data/
    ├── input/
    ├── output/
    └── temp/
```

---

## Quick Start

### 1. Clone and Install
```bash
git clone https://github.com/yourusername/langgraph_bioflow.git
cd langgraph_bioflow
pip install -r requirements.txt
```

### 2. Set Up Cromwell
Download and place Cromwell:
```bash
wget https://github.com/broadinstitute/cromwell/releases/download/85/cromwell-85.jar -O /usr/local/bin/cromwell.jar
chmod +x /usr/local/bin/cromwell.jar
```

Verify:
```bash
java -jar /usr/local/bin/cromwell.jar --version
```

### 3. Test Docker
Ensure Docker works:
```bash
docker run hello-world
```

### 4. Run the Demo
```bash
python main.py
```

Expected output (for mock workflow):
```
[InputManager] Starting new dataset workflow...
[AlignerNode] Launching WDL alignment workflow...
[QCNode] Collecting quality metrics...
[DecisionNode] Decision: continue
[ImputationNode] Running imputation WDL...
[GWASNode] Performing association analysis...
[ReportNode] Generating summary report...
Final Result: {'status': 'complete'}
```

---

## Future Enhancements
- [ ] Integrate **FastQC parser** for real metrics  
- [ ] Add **LLM-driven decision node** for adaptive QC reasoning  
- [ ] Store state and results in **SQLite or MongoDB**  
- [ ] Add **LangSmith tracing** for visualization  
- [ ] Create **FastAPI dashboard** for human review and workflow monitoring  


---

## Maintainer
**Manish Kumar**  
Sr. Bioinformatics Software Engineer & Data Scientist  
````
