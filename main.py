from langgraph.graph import StateGraph, END
from pydantic import BaseModel
import random
import os
import subprocess

# --- Define State Schema ---
class WorkflowState(BaseModel):
    attempt: int = 1
    decision: str = "UNKNOWN"

# --- Define Nodes ---
def nextflow_run(state: WorkflowState):
    attempt = state.attempt
    print(f"Running Nextflow workflow (attempt {attempt})...\n")
    subprocess.run(["nextflow", "run", "workflows/qc_pipeline.nf"], check=True)
    # simulate reading QC summary
    qc_summary = os.path.join("work", "latest", "qc_summary.txt")
    # Simulate: randomly FAIL or PASS
    decision = "PASS" if random.random() > 0.7 else "FAIL"
    print(f"QC Decision: {decision}")
    return {"decision": decision, "attempt": attempt}

def qc_decision(state: WorkflowState):
    decision = state.decision
    attempt = state.attempt
    if decision == "PASS":
        print("QC passed. Proceeding to END.")
        return {"decision": "PASS", "attempt": attempt}
    else:
        print(f"QC failed on attempt {attempt}. Re-running with adjusted parameters (next: {attempt + 1}).")
        return {"decision": "FAIL", "attempt": attempt}

def qc_condition(state: WorkflowState):
    if state.decision == "PASS":
        return "done"
    elif state.attempt > 3:  # stop runaway loops
        print("Reached max attempts (3). Stopping workflow.")
        return "done"
    else:
        state.attempt += 1
        return "rerun"

# --- Build Graph ---
graph = StateGraph(WorkflowState)
graph.add_node("nextflow_run", nextflow_run)
graph.add_node("qc_decision", qc_decision)

graph.add_edge("__start__", "nextflow_run")
graph.add_edge("nextflow_run", "qc_decision")

graph.add_conditional_edges(
    "qc_decision",
    qc_condition,
    {"rerun": "nextflow_run", "done": END},
)

app = graph.compile()

# --- Run ---
if __name__ == "__main__":
    state = WorkflowState()
    final_state = app.invoke(state)
    print("Final State:", final_state)
