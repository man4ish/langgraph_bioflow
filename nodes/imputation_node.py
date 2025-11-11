from utils.wdl_runner import run_wdl
import yaml

def imputation_node(state):
    print("[ImputationNode] Running imputation WDL...")
    with open("config/workflow_config.yaml") as f:
        config = yaml.safe_load(f)

    cromwell = config["paths"]["cromwell_jar"]
    wdl_file = config["paths"]["wdl_files"]["imputation"]
    docker_img = config["paths"]["docker_image"]
    output_dir = config["paths"]["output_dir"]

    result = run_wdl(wdl_file, cromwell_path=cromwell, docker_image=docker_img, output_dir=output_dir)
    return {"imputation_done": True, "imputation_status": result["status"]}
