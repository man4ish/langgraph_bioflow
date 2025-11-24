import streamlit as st
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
import requests
from time import sleep

# ----------------------------
# Config
# ----------------------------
DEMO_MODE = True  # Set False in production
API_BASE_URL = "http://127.0.0.1:8000"  # FastAPI backend

# ----------------------------
# Apply Light Theme Styling
# ----------------------------
st.set_page_config(page_title="LangGraph BioFlow", layout="wide")

st.markdown(
    """
    <style>
    body { background-color: #f5f5f5; color: #000000; }
    .stProgress > div > div > div > div { background-color: #4CAF50; }
    .stDataFrame th { background-color: #4CAF50; color: white; }
    </style>
    """,
    unsafe_allow_html=True
)

# ----------------------------
# Sidebar
# ----------------------------
st.sidebar.title("LangGraph BioFlow")
workflow_stage = st.sidebar.radio(
    "Start from stage:",
    ["Pipeline Choice", "Pathway Enrichment", "RAG", "GNN"]
)
st.sidebar.markdown("---")
st.sidebar.header("Pipeline Options")
pipeline_choice = st.sidebar.selectbox(
    "Choose pipeline (WGS / Exome / RNA-Seq):",
    ["WGS", "Exome", "RNA-Seq"]
)
st.sidebar.button("Run Workflow")

# ----------------------------
# Tabs
# ----------------------------
tab1, tab2, tab3, tab4, tab5 = st.tabs(
    ["Workflow Status", "Pathway Enrichment", "RAG Network", "GNN Predictions", "Logs"]
)

# ----------------------------
# Tab 1: Workflow Status
# ----------------------------
with tab1:
    st.header("Workflow Progress")
    status_placeholder = st.empty()
    progress_bar = st.progress(0)

    if DEMO_MODE:
        for i in range(1, 101):
            sleep(0.03)
            progress_bar.progress(i)
            status_placeholder.markdown(f"**Current Stage:** {workflow_stage} | Step {i}% completed")
    else:
        try:
            resp = requests.get(f"{API_BASE_URL}/workflow/status")
            status_data = resp.json()
            for i, stage in enumerate(status_data.get("stages", []), start=1):
                progress_bar.progress(int((i / len(status_data['stages'])) * 100))
                status_placeholder.markdown(f"**Stage:** {stage['name']} | Status: {stage['status']}")
        except Exception as e:
            st.error(f"Error fetching workflow status: {e}")

# ----------------------------
# Tab 2: Pathway Enrichment
# ----------------------------
with tab2:
    st.header("Pathway Enrichment Results")
    if DEMO_MODE:
        enrichment_data = pd.DataFrame({
            "Gene": ["TP53", "APP", "SNCA", "LRRK2"],
            "Pathway": ["Apoptosis", "Neurodegeneration", "Parkinson", "Cell Cycle"],
            "Score": [0.95, 0.88, 0.76, 0.67]
        })
    else:
        try:
            resp = requests.get(f"{API_BASE_URL}/pathway/enrichment")
            enrichment_data = pd.DataFrame(resp.json())
        except Exception as e:
            st.error(f"Error fetching enrichment data: {e}")
            enrichment_data = pd.DataFrame()

    st.dataframe(enrichment_data.style.set_properties(**{
        'background-color': '#ffffff',
        'color': '#000000',
        'border-color': '#4CAF50'
    }), use_container_width=True)
    st.download_button("Download Enrichment JSON",
                       data=enrichment_data.to_json(orient="records"),
                       file_name="enrichment_results.json")

# ----------------------------
# Tab 3: RAG Network
# ----------------------------
with tab3:
    st.header("RAG Knowledge Graph")
    G = nx.Graph()
    if DEMO_MODE:
        G.add_edges_from([("Gene1", "PathwayA"), ("Gene2", "PathwayB"), ("Gene3", "PathwayC")])
    else:
        try:
            resp = requests.get(f"{API_BASE_URL}/rag/network")
            edges = resp.json().get("edges", [])
            G.add_edges_from(edges)
        except Exception as e:
            st.error(f"Error fetching RAG network: {e}")

    fig, ax = plt.subplots(figsize=(7, 5))
    nx.draw(G, with_labels=True, node_color='lightgreen', node_size=1400, font_size=12, ax=ax)
    st.pyplot(fig)

# ----------------------------
# Tab 4: GNN Predictions
# ----------------------------
with tab4:
    st.header("GNN Drug Predictions")
    if DEMO_MODE:
        predictions = pd.DataFrame({
            "Drug": ["DrugA", "DrugB", "DrugC"],
            "Score": [0.95, 0.88, 0.81]
        })
    else:
        try:
            resp = requests.get(f"{API_BASE_URL}/gnn/predictions")
            predictions = pd.DataFrame(resp.json())
        except Exception as e:
            st.error(f"Error fetching GNN predictions: {e}")
            predictions = pd.DataFrame()

    st.dataframe(predictions.style.set_properties(**{
        'background-color': '#ffffff',
        'color': '#000000',
        'border-color': '#4CAF50'
    }), use_container_width=True)
    st.download_button("Download GNN Predictions",
                       data=predictions.to_json(orient="records"),
                       file_name="gnn_predictions.json")

# ----------------------------
# Tab 5: Logs
# ----------------------------
with tab5:
    st.header("Workflow Logs")
    if DEMO_MODE:
        log_text = """
        [2025-11-23 10:00] QC started
        [2025-11-23 10:05] QC passed
        [2025-11-23 10:10] Pipeline: WGS started
        [2025-11-23 10:30] Pathway enrichment completed
        [2025-11-23 10:45] RAG query executed
        [2025-11-23 11:00] GNN predictions generated
        """
    else:
        try:
            resp = requests.get(f"{API_BASE_URL}/workflow/logs")
            log_entries = resp.json().get("logs", [])
            log_text = "\n".join(log_entries)
        except Exception as e:
            st.error(f"Error fetching logs: {e}")
            log_text = ""
    st.text_area("Live Logs", value=log_text, height=200)
