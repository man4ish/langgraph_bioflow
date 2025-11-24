import streamlit as st
import requests
import time

API_URL = "http://localhost:8000"

st.title("LangGraph BioFlow – Real-Time Dashboard")

run_id = st.text_input("Enter Run ID:")

if st.button("Check Status"):
    if run_id:
        while True:
            status = requests.get(f"{API_URL}/status/{run_id}").json()

            st.write(status)

            if status["status"] in ("COMPLETED", "FAILED"):
                st.success(f"Pipeline {status['status']}")
                break

            time.sleep(2)
            st.experimental_rerun()
