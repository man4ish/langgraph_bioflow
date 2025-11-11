from utils.wdl_runner import run_wdl
import yaml

def aligner_node(state):
    print("[AlignerNode] Launching WDL alignment workflow...")
    with open("config/workflow_config.yaml") as f:
        config = yaml.safe_load(f)

    cromwell = config["paths"]["cromwell_jar"]
    wdl_file = config["paths"]["wdl_files"]["alignment"]
    docker_img = config["paths"]["docker_image"]
    output_dir = config["paths"]["output_dir"]

    result = run_wdl(wdl_file, cromwell_path=cromwell, docker_image=docker_img, output_dir=output_dir)
    return {"alignment_done": True, "alignment_status": result["status"]}
