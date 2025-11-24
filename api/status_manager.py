import json
from pathlib import Path

RUNS_DIR = Path("logs/runs")
RUNS_DIR.mkdir(exist_ok=True, parents=True)

def update_status(run_id, status, stage=None, error=None):
    status_file = RUNS_DIR / f"{run_id}.json"
    data = {
        "run_id": run_id,
        "status": status,
        "stage": stage,
        "error": error
    }
    with open(status_file, "w") as f:
        json.dump(data, f, indent=2)

def read_status(run_id):
    status_file = RUNS_DIR / f"{run_id}.json"
    if not status_file.exists():
        return {"error": "Run ID not found"}
    return json.loads(status_file.read_text())
