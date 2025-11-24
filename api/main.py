from fastapi import FastAPI
from api.orchestrator import Orchestrator
from api.status_manager import read_status

app = FastAPI()

@app.get("/")
def root():
    return {"message": "LangGraph BioFlow API is running"}

@app.post("/run_pipeline")
def run_pipeline(sample_id: str, pipeline: str = "exome"):
    run_id = Orchestrator.run_full_pipeline(sample_id, pipeline)
    return {"run_id": run_id, "message": "Pipeline started"}

@app.get("/status/{run_id}")
def get_status(run_id: str):
    return read_status(run_id)
