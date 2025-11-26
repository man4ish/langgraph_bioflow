import streamlit as st
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
import requests
from time import sleep
from PIL import Image

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
# Functional Header Bar + Soft Background Styling
# ----------------------------

# Add soft pastel background to header row
st.markdown("""
<div style="
    display:flex;
    justify-content:space-between;
    align-items:center;  /* vertical align center */
    padding:2px 0px;
">
    <div style='font-size:26px; font-weight:bold; color:#2E4053;'>
        Workflow Orchestration Dashboard
    </div>
    <div style='width:250px;'>
        <input type="text" placeholder="Search..."
               style="width:100%; padding:8px; border-radius:6px; border:1px solid #ccc;" />
    </div>
</div>
""", unsafe_allow_html=True)


# ----------------------------
# Sidebar
# ----------------------------
# Load your logo
logo = Image.open("../assets/Langgraph_bioflow.png")
st.sidebar.image(logo, width=150)

# Sidebar title and description
st.sidebar.markdown("## LangGraph BioFlow")
st.sidebar.markdown("AI-Orchestrated Bioinformatics Workflow Supervisor")
st.sidebar.markdown("---")

# Workflow stage selection
workflow_stage = st.sidebar.radio(
    "Start from stage:",
    ["Pipeline Choice", "Pathway Enrichment", "RAG", "GNN"]
)

# Pipeline options
st.sidebar.header("Pipeline Options")
pipeline_choice = st.sidebar.selectbox(
    "Choose pipeline (WGS / Exome / RNA-Seq):",
    ["WGS", "Exome", "RNA-Seq"]
)

# Run workflow button
st.sidebar.button("Run Workflow")

# ----------------------------
# Tabs
# ----------------------------
tab1, tab2, tab3, tab4, tab5 = st.tabs(
    ["Workflow Status", "Pathway Enrichment", "RAG Network", "GNN Predictions", "Logs"]
)
tab_height = 50  # fixed height in px, adjust as needed

# ----------------------------
# Tab 1: Workflow Status
# ----------------------------
with tab1:
    st.markdown(f"<div style='height:{tab_height}px; overflow-y:auto; padding-right:10px;'>", unsafe_allow_html=True)

    st.markdown("""
    <div style="
        background-color:#C8E6C9;  
        padding:5px 10px; 
        border-radius:5px; 
        display:inline-block;
        font-size:18px;
        font-weight:bold;
        color:#2E4053;
        margin-bottom:5px;">
        Workflow Progress
    </div>
    """, unsafe_allow_html=True)

    st.markdown("Monitor the pipeline execution and QC progress")
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

    st.markdown("</div>", unsafe_allow_html=True)

# ----------------------------
# Tab 2: Pathway Enrichment
# ----------------------------
with tab2:
    st.markdown(f"<div style='height:{tab_height}px; overflow-y:auto; padding-right:10px;'>", unsafe_allow_html=True)

    st.markdown("""
    <div style="
        background-color:#A3CFE5; 
        padding:5px 10px; 
        border-radius:5px; 
        display:inline-block;
        font-size:18px;
        font-weight:bold;
        color:#2E4053;
        margin-bottom:5px;">
        Pathway Enrichment Results
    </div>
    """, unsafe_allow_html=True)

    st.markdown("Explore enriched pathways for your gene sets")
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

    st.markdown("</div>", unsafe_allow_html=True)

# ----------------------------
# Tab 3: RAG Network
# ----------------------------
with tab3:
    st.markdown(f"<div style='height:{tab_height}px; overflow-y:auto; padding-right:10px;'>", unsafe_allow_html=True)

    st.markdown("""
    <div style="
        background-color:#FFE0B2;  
        padding:5px 10px; 
        border-radius:5px; 
        display:inline-block;
        font-size:18px;
        font-weight:bold;
        color:#2E4053;
        margin-bottom:5px;">
        RAG Knowledge Graph
    </div>
    """, unsafe_allow_html=True)

    st.markdown("Visualize gene-pathway-drug relationships using the RAG network")
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

    st.markdown("</div>", unsafe_allow_html=True)

# ----------------------------
# Tab 4: GNN Predictions
# ----------------------------
with tab4:
    st.markdown(f"<div style='height:{tab_height}px; overflow-y:auto; padding-right:10px;'>", unsafe_allow_html=True)

    st.markdown("""
    <div style="
        background-color:#D1C4E9;  
        padding:5px 10px; 
        border-radius:5px; 
        display:inline-block;
        font-size:18px;
        font-weight:bold;
        color:#2E4053;
        margin-bottom:5px;">
        GNN Drug Predictions
    </div>
    """, unsafe_allow_html=True)

    st.markdown("View predicted drugs and their scores from the GNN model")
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

    st.markdown("</div>", unsafe_allow_html=True)

# ----------------------------
# Tab 5: Logs
# ----------------------------
with tab5:
    st.markdown(f"<div style='height:{tab_height}px; overflow-y:auto; padding-right:10px;'>", unsafe_allow_html=True)

    st.markdown("""
    <div style="
        background-color:#CFD8DC;  
        padding:5px 10px; 
        border-radius:5px; 
        display:inline-block;
        font-size:18px;
        font-weight:bold;
        color:#2E4053;
        margin-bottom:5px;">
        Workflow Logs
    </div>
    """, unsafe_allow_html=True)

    st.markdown("View detailed logs of the workflow execution and status")
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

    st.markdown("</div>", unsafe_allow_html=True)


st.markdown("---")  # Horizontal separator

st.markdown("""
<div style="
    min-height:40vh;  /* minimum height for content */
    display:flex;
    flex-direction:column;
    justify-content:flex-end;">
    <div style="text-align:center; color:#555555; font-size:12px;">
        Developed by Manish Kumar | © 2025
    </div>
</div>
""", unsafe_allow_html=True)
