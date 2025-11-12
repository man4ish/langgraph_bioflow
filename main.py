from langgraph.graph import StateGraph, END
from pydantic import BaseModel
import subprocess
import os
import random

# --- Define Workflow State ---
class WorkflowState(BaseModel):
    step: str = "qc"
    decision: str = "UNKNOWN"
    attempt: int = 1


# --- Define Nodes ---
def run_qc(state: WorkflowState):
    """Runs the FastQC-only pipeline"""
    print(f"\n=== Running QC (attempt {state.attempt}) ===\n")
    subprocess.run(["nextflow", "run", "workflows/qc_pipeline.nf"], check=True)

    # Simulate QC decision (replace this with parsing QC report)
    decision = "PASS" if random.random() > 0.6 else "FAIL"
    print(f"QC Decision: {decision}")
    return {"decision": decision, "step": "qc", "attempt": state.attempt}


def run_exome_pipeline(state: WorkflowState):
    """Runs the main exome pipeline"""
    print("\n=== Running Exome Pipeline ===\n")
    subprocess.run(["nextflow", "run", "workflows/exome_pipeline.nf"], check=True)
    return {"decision": "COMPLETE", "step": "exome", "attempt": state.attempt}


def qc_condition(state: WorkflowState):
    """Decides what to do after QC"""
    if state.decision == "PASS":
        return "continue"
    elif state.attempt >= 3:
        print("Reached maximum attempts (3). Stopping workflow.")
        return "done"
    else:
        print(f"QC failed on attempt {state.attempt}. Retrying with adjusted parameters...")
        state.attempt += 1
        return "retry"


# --- Build LangGraph ---
graph = StateGraph(WorkflowState)
graph.add_node("run_qc", run_qc)
graph.add_node("run_exome_pipeline", run_exome_pipeline)

graph.add_edge("__start__", "run_qc")

graph.add_conditional_edges(
    "run_qc",
    qc_condition,
    {
        "retry": "run_qc",
        "continue": "run_exome_pipeline",
        "done": END,
    },
)

graph.add_edge("run_exome_pipeline", END)
app = graph.compile()

# --- Run Workflow ---
if __name__ == "__main__":
    state = WorkflowState()
    final_state = app.invoke(state)
    print("\n=== Final Workflow State ===")
    print(final_state)
