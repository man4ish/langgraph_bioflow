import subprocess
import os

def run_wdl(wdl_path, inputs_path=None, cromwell_path=None, docker_image=None, output_dir="./data/output"):
    """
    Runs a WDL file using Cromwell in Docker mode.
    """
    if not os.path.exists(cromwell_path):
        raise FileNotFoundError(f"Cromwell not found at {cromwell_path}")

    cmd = [
        "java", "-jar", cromwell_path, "run", wdl_path,
        "--inputs", inputs_path if inputs_path else "{}"
    ]

    print(f"[WDLRunner] Executing: {' '.join(cmd)}")

    try:
        subprocess.run(cmd, check=True)
        print("[WDLRunner] Workflow completed successfully.")
        return {"status": "success", "output_dir": output_dir}
    except subprocess.CalledProcessError as e:
        print(f"[WDLRunner] Workflow failed: {e}")
        return {"status": "failed"}
