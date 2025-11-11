import subprocess
import os
from typing import TypedDict


class NextflowRunnerState(TypedDict, total=False):
    config: dict
    nf_output: str
    run_status: str


def nextflow_runner(state: NextflowRunnerState):
    config = state["config"]
    nf_script = config["nextflow_workflow"]["script"]
    docker_image = config["nextflow_workflow"]["docker_image"]

    print(f"Running Nextflow pipeline: {nf_script}")
    print(f"Using Docker image: {docker_image}")

    try:
        cmd = [
            "nextflow", "run", nf_script,
            "-with-docker", docker_image,
            "-work-dir", config["nextflow_workflow"]["workdir"]
        ]
        print("Executing:", " ".join(cmd))

        subprocess.run(cmd, check=True)
        state["nf_output"] = "results/qc_summary.txt"
        state["run_status"] = "success"
        print("Nextflow workflow completed successfully.")

    except subprocess.CalledProcessError as e:
        state["run_status"] = "failed"
        print(f"Error running Nextflow: {e}")

    return state
