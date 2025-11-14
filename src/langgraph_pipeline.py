# langgraph_pipeline.py
"""
Full AI-augmented Exome Workflow with LangGraph
Nodes: QC → Exome → SnpEff → Pathway Enrichment → RAG → Neo4j
Author: Manish Kumar
Date: 2025-11-14
"""

from langgraph.graph import StateGraph, END
from pydantic import BaseModel
import subprocess
import random

# Import updated modules
from pathway_enrichment_gene_list import GeneListEnrichment
from rag_runner import RAGRunner

# -------------------------------
# Workflow State
# -------------------------------
class WorkflowState(BaseModel):
    step: str = "qc"
    decision: str = "UNKNOWN"
    attempt: int = 1
    gene_list: list = []
    top_genes: list = []
    pathway_queries: list = []
    rag_output_files: list = []

# -------------------------------
# Nodes
# -------------------------------
def run_qc(state: WorkflowState):
    """Runs FastQC and decides pass/fail"""
    print(f"\n=== Running QC (attempt {state.attempt}) ===\n")
    subprocess.run(["nextflow", "run", "workflows/qc_pipeline.nf"], check=True)

    # Replace with real QC parsing
    decision = "PASS" if random.random() > 0.3 else "FAIL"
    print(f"QC Decision: {decision}")
    state.decision = decision
    return state

def run_exome_pipeline(state: WorkflowState):
    """Run exome pipeline and extract gene list from SnpEff"""
    print("\n=== Running Exome Pipeline ===\n")
    subprocess.run(["nextflow", "run", "workflows/exome_pipeline.nf"], check=True)

    # Replace with parsing actual SnpEff output
    state.gene_list = ["TP53", "BRCA1", "BRCA2", "LRRK2", "SNCA", "APP"]
    print(f"Genes extracted from SnpEff: {state.gene_list}")
    return state

def run_pathway_enrichment(state: WorkflowState):
    """Perform pathway enrichment and generate RAG queries"""
    print("\n=== Running Pathway Enrichment ===\n")
    enrichment = GeneListEnrichment(state.gene_list)
    enrichment.run_enrichment(database="Reactome_2022")

    top_pathways = enrichment.top_pathways(5)
    queries = []
    for pathway in top_pathways['Term']:
        genes = enrichment.get_pathway_genes(pathway)
        if not genes:
            continue
        query = f"Genes in the {pathway} pathway related to Parkinson's disease: {', '.join(genes)}"
        queries.append(query)

    enrichment.plot_top_pathways(10, save_path="results/enrichment/top_pathways.png")
    state.pathway_queries = queries
    print(f"Generated {len(queries)} RAG queries")
    return state

def run_rag_pipeline(state: WorkflowState):
    """Run RAG pipeline for all pathway queries"""
    print("\n=== Running RAG Pipeline ===\n")
    if not state.pathway_queries:
        raise ValueError("No pathway queries provided")
    
    rag = RAGRunner()
    state.rag_output_files = rag.run_queries(state.pathway_queries)
    return state

def load_into_neo4j(state: WorkflowState):
    """Load RAG outputs into Neo4j"""
    print("\n=== Loading into Neo4j ===\n")
    for file in state.rag_output_files:
        subprocess.run(["python3", "src/rag_to_neo4j.py", "--input", file], check=True)
    print("Neo4j ingestion completed")
    return state

# -------------------------------
# Conditions
# -------------------------------
def qc_condition(state: WorkflowState):
    """Decide whether to retry QC or continue"""
    if state.decision == "PASS":
        return "continue"
    elif state.attempt >= 3:
        print("Reached max QC attempts. Stopping workflow.")
        return "done"
    else:
        print(f"QC failed on attempt {state.attempt}. Retrying...")
        state.attempt += 1
        return "retry"

# -------------------------------
# Build LangGraph
# -------------------------------
graph = StateGraph(WorkflowState)
graph.add_node("run_qc", run_qc)
graph.add_node("run_exome_pipeline", run_exome_pipeline)
graph.add_node("run_pathway_enrichment", run_pathway_enrichment)
graph.add_node("run_rag_pipeline", run_rag_pipeline)
graph.add_node("load_into_neo4j", load_into_neo4j)

# QC branching
graph.add_edge("__start__", "run_qc")
graph.add_conditional_edges(
    "run_qc",
    qc_condition,
    {
        "retry": "run_qc",
        "continue": "run_exome_pipeline",
        "done": END,
    }
)

# Linear flow for remaining nodes
graph.add_edge("run_exome_pipeline", "run_pathway_enrichment")
graph.add_edge("run_pathway_enrichment", "run_rag_pipeline")
graph.add_edge("run_rag_pipeline", "load_into_neo4j")
graph.add_edge("load_into_neo4j", END)

# Compile workflow
app = graph.compile()

# -------------------------------
# Run Workflow
# -------------------------------
if __name__ == "__main__":
    state = WorkflowState()
    final_state = app.invoke(state)
    print("\n=== Final Workflow State ===")
    print(final_state)
