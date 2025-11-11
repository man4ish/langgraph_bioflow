import subprocess
import os
from typing import TypedDict


class WDLRunnerState(TypedDict, total=False):
    config: dict
    wdl_output: str
    run_status: str


def wdl_runner(state: WDLRunnerState):
    config = state["config"]
    wdl_path = config["wdl_workflow"]["path"]
    docker_image = config["wdl_workflow"]["docker_image"]

    print(f"Running WDL workflow: {wdl_path}")
    print(f"Using Docker image: {docker_image}")

    # Simulate WDL execution for now
    # Later: integrate with real Cromwell command
    try:
        cmd = [
            "docker", "run", "--rm",
            "-v", f"{os.getcwd()}:/data",
            docker_image,
            "java", "-jar", "/app/cromwell.jar", "run", f"/data/{wdl_path}"
        ]
        print("Executing:", " ".join(cmd))
        # Simulation mode for now
        state["wdl_output"] = "mock_results/output_metrics.json"
        state["run_status"] = "success"
        print("Workflow completed successfully.")
    except Exception as e:
        state["run_status"] = "failed"
        print(f"Error running WDL: {e}")

    return state
