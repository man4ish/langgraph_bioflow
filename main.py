from typing import TypedDict
import yaml
from langgraph.graph import StateGraph, END
from nodes.nextflow_runner import nextflow_runner

# --- Define state schema ---
class BioFlowState(TypedDict, total=False):
    config: dict
    input: str
    qc_pass: bool
    wdl_output: str
    run_status: str


def load_config(path: str = "config/workflow_config.yaml") -> dict:
    with open(path, "r") as f:
        return yaml.safe_load(f)


def start_node(state: BioFlowState):
    config = state["config"]
    print("Starting BioFlow test...")

    first_input = config["input_files"][0]
    state["input"] = first_input
    print(f"Loaded input file: {first_input}")
    return state


def qc_node(state: BioFlowState):
    thresholds = state["config"]["qc_thresholds"]
    print(f"Running QC on {state['input']}")
    print(f"QC thresholds: {thresholds}")

    # Dummy pass for now
    state["qc_pass"] = True
    return state


# --- Build graph ---
graph = StateGraph(BioFlowState)
graph.add_node("start", start_node)
graph.add_node("qc", qc_node)
graph.add_node("wdl_runner", wdl_runner)

graph.add_edge("start", "qc")
graph.add_edge("qc", "wdl_runner")
graph.add_edge("wdl_runner", END)
graph.set_entry_point("start")


# --- Compile & run ---
if __name__ == "__main__":
    config = load_config()
    app = graph.compile()
    app.invoke({"config": config})
