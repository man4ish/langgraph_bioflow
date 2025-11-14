from langgraph.graph import StateGraph, END
from pydantic import BaseModel
import subprocess
import random

# ------------------------------
# Workflow State
# ------------------------------
class WorkflowState(BaseModel):
    step: str = "qc"
    decision: str = "UNKNOWN"
    attempt: int = 1
    pipeline_choice: str = None  # 'wgs', 'exome', 'rnaseq'
    gene_list: list = None
    rag_output: dict = None
    gnn_output: dict = None

# ------------------------------
# Nodes
# ------------------------------
def run_qc(state: WorkflowState):
    """Run QC (FastQC, retry if fail)"""
    print(f"\n=== Running QC (attempt {state.attempt}) ===\n")
    subprocess.run(["nextflow", "run", "workflows/qc_pipeline.nf"], check=True)
    decision = "PASS" if random.random() > 0.6 else "FAIL"
    print(f"QC Decision: {decision}")
    return {"decision": decision, "step": "qc", "attempt": state.attempt}

def choose_pipeline(state: WorkflowState):
    """Choose sequencing pipeline: WGS / Exome / RNA-Seq"""
    if not state.pipeline_choice:
        # For example, can use user input or config
        state.pipeline_choice = random.choice(["wgs", "exome", "rnaseq"])
    print(f"Chosen pipeline: {state.pipeline_choice}")
    subprocess.run(["nextflow", "run", f"workflows/{state.pipeline_choice}_pipeline.nf"], check=True)
    return {"step": "pipeline", "pipeline_choice": state.pipeline_choice}

def run_pathway_enrichment(state: WorkflowState):
    """Generate gene list and run pathway enrichment"""
    # Simulate gene list output
    state.gene_list = ["TP53", "APP", "SNCA", "LRRK2"]
    subprocess.run(["python", "pathway_enrichment_gene_list.py", "--genes"] + state.gene_list, check=True)
    return {"step": "pathway", "gene_list": state.gene_list}

def run_rag(state: WorkflowState):
    """Run RAG pipeline on gene list"""
    query = f"genes linked to disease for {state.pipeline_choice}"
    subprocess.run(["python", "src/rag_pipeline.py", "--query", query], check=True)
    state.rag_output = {"neo4j_nodes": 100, "neo4j_edges": 150}  # Example output
    return {"step": "rag", "rag_output": state.rag_output}

def run_gnn(state: WorkflowState):
    """Run GNN predictor using RAG Neo4j output"""
    subprocess.run(["python", "gnn_predictor.py", "--input", "rag_output_placeholder.json"], check=True)
    state.gnn_output = {"predictions": ["DrugA", "DrugB"], "scores": [0.95, 0.88]}
    return {"step": "gnn", "gnn_output": state.gnn_output}

def visualize_results(state: WorkflowState):
    """Visualize GNN predictions and scores"""
    subprocess.run(["python", "visualization.py", "--gnn_output", "gnn_output_placeholder.json"], check=True)
    return {"step": "visualization"}

# ------------------------------
# Conditional Logic
# ------------------------------
def qc_condition(state: WorkflowState):
    if state.decision == "PASS":
        return "continue"
    elif state.attempt >= 3:
        print("QC failed max attempts. Stopping workflow.")
        return "done"
    else:
        state.attempt += 1
        print(f"Retrying QC (attempt {state.attempt})...")
        return "retry"

# ------------------------------
# Build LangGraph
# ------------------------------
graph = StateGraph(WorkflowState)
graph.add_node("run_qc", run_qc)
graph.add_node("choose_pipeline", choose_pipeline)
graph.add_node("run_pathway_enrichment", run_pathway_enrichment)
graph.add_node("run_rag", run_rag)
graph.add_node("run_gnn", run_gnn)
graph.add_node("visualize_results", visualize_results)

# QC branch
graph.add_edge("__start__", "run_qc")
graph.add_conditional_edges(
    "run_qc",
    qc_condition,
    {
        "retry": "run_qc",
        "continue": "choose_pipeline",
        "done": END,
    },
)

# Linear downstream flow
graph.add_edge("choose_pipeline", "run_pathway_enrichment")
graph.add_edge("run_pathway_enrichment", "run_rag")
graph.add_edge("run_rag", "run_gnn")
graph.add_edge("run_gnn", "visualize_results")
graph.add_edge("visualize_results", END)

app = graph.compile()

# ------------------------------
# Run Workflow
# ------------------------------
if __name__ == "__main__":
    state = WorkflowState()
    final_state = app.invoke(state)
    print("\n=== Final Workflow State ===")
    print(final_state)
