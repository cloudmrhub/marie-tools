#!/bin/bash

# Check if correct number of arguments provided
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <MCR_ROOT> <RUNNER_NAME> <INPUT_FILE>"
    echo "Example: $0 /opt/matlabruntime/R2023a MARIE_runner data/inputs/inp_Duke_StadiumTriangular.json"
    exit 1
fi

MCR_ROOT=$1
RUNNER_NAME=$2
INPUT_FILE=$3

# Determine which runner script to use based on the runner name
case "$RUNNER_NAME" in
    MARIE_runner)
        RUNNER_SCRIPT="./run_MARIE_runner.sh"
        ;;
    MRGF_runner)
        RUNNER_SCRIPT="./run_MRGF_runner.sh"
        ;;
    *)
        echo "Error: Unknown runner '$RUNNER_NAME'"
        echo "Supported runners: MARIE_runner, MRGF_runner"
        exit 1
        ;;
esac

# Check if the runner script exists
if [ ! -f "$RUNNER_SCRIPT" ]; then
    echo "Error: Runner script '$RUNNER_SCRIPT' not found"
    exit 1
fi

# Execute the appropriate runner
echo "Running $RUNNER_NAME with input: $INPUT_FILE"
exec "$RUNNER_SCRIPT" "$MCR_ROOT" "$INPUT_FILE"

