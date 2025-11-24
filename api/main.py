# langgraph_bioflow/api/main.py
from fastapi import FastAPI
from pydantic import BaseModel
from typing import List
import random

app = FastAPI(title="LangGraph BioFlow API")

# ----------------------------
# Data Models
# ----------------------------
class StageStatus(BaseModel):
    name: str
    status: str  # "pending", "running", "completed", "failed"

class WorkflowStatus(BaseModel):
    stages: List[StageStatus]

# ----------------------------
# Endpoints
# ----------------------------

@app.get("/")
def root():
    return {"message": "LangGraph BioFlow API is running"}

# Workflow status
@app.get("/workflow/status")
def workflow_status():
    # Dummy data for testing
    stages = [
        {"name": "QC", "status": "completed"},
        {"name": "Exome Pipeline", "status": "running"},
        {"name": "Pathway Enrichment", "status": "pending"},
        {"name": "RAG", "status": "pending"},
        {"name": "GNN", "status": "pending"},
    ]
    return {"stages": stages}

# Pathway enrichment results
@app.get("/pathway/enrichment")
def pathway_enrichment():
    data = [
        {"Gene": "TP53", "Pathway": "Apoptosis", "Score": 0.95},
        {"Gene": "APP", "Pathway": "Neurodegeneration", "Score": 0.88},
        {"Gene": "SNCA", "Pathway": "Parkinson", "Score": 0.76},
        {"Gene": "LRRK2", "Pathway": "Cell Cycle", "Score": 0.67},
    ]
    return data

# RAG network edges
@app.get("/rag/network")
def rag_network():
    edges = [
        ("TP53", "Apoptosis"),
        ("APP", "Neurodegeneration"),
        ("SNCA", "Parkinson"),
        ("LRRK2", "Cell Cycle"),
    ]
    return {"edges": edges}

# GNN predictions
@app.get("/gnn/predictions")
def gnn_predictions():
    predictions = [
        {"Drug": "DrugA", "Score": 0.95},
        {"Drug": "DrugB", "Score": 0.88},
        {"Drug": "DrugC", "Score": 0.81},
    ]
    return predictions

# Workflow logs
@app.get("/workflow/logs")
def workflow_logs():
    logs = [
        "[2025-11-23 10:00] QC started",
        "[2025-11-23 10:05] QC passed",
        "[2025-11-23 10:10] Exome pipeline started",
        "[2025-11-23 10:30] Pathway enrichment completed",
        "[2025-11-23 10:45] RAG query executed",
        "[2025-11-23 11:00] GNN predictions generated",
    ]
    return {"logs": logs}
