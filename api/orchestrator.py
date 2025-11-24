import subprocess
import uuid
import time
import json
from pathlib import Path
from api.status_manager import update_status

RUNS_DIR = Path("logs/runs")
RUNS_DIR.mkdir(parents=True, exist_ok=True)

class Orchestrator:

    @staticmethod
    def run_full_pipeline(sample_id: str, pipeline: str = "exome"):

        run_id = str(uuid.uuid4())
        update_status(run_id, "STARTED", stage="initializing")

        try:
            # 1. QC Stage
            update_status(run_id, "RUNNING", stage="fastqc")
            subprocess.run(["nextflow", "run", "workflows/qc_pipeline.nf"], check=True)

            # 2. Nextflow Pipeline
            update_status(run_id, "RUNNING", stage=f"pipeline_{pipeline}")
            subprocess.run(["nextflow", "run", f"workflows/{pipeline}_pipeline.nf"], check=True)

            # 3. Pathway Enrichment
            update_status(run_id, "RUNNING", stage="pathway_enrichment")
            subprocess.run(["python", "src/analysis/pathway_enrichment_gene_list.py"], check=True)

            # 4. Gene → Pathway Visualization
            update_status(run_id, "RUNNING", stage="network_visualization")
            subprocess.run(["python", "src/analysis/gene_pathway_network_visualization_rag.py"], check=True)

            # 5. RAG Query
            update_status(run_id, "RUNNING", stage="rag")
            subprocess.run(["python", "-m", "src.rag.rag_pipeline_langchain", "--structured"], check=True)

            # 6. GNN Prediction
            update_status(run_id, "RUNNING", stage="gnn_predictor")
            subprocess.run(["python", "src/gnn_predictor/models/gnn_predictor.py"], check=True)

            update_status(run_id, "COMPLETED", stage="done")

        except subprocess.CalledProcessError as e:
            update_status(run_id, "FAILED", stage="error", error=str(e))

        return run_id
