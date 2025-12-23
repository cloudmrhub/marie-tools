#!/bin/bash

# Check if correct number of arguments provided
if [ "$#" -lt 3 ] || [ "$#" -gt 4 ]; then
    echo "Usage: $0 <MCR_ROOT> <RUNNER_NAME> <INPUT_FILE> [S3_BUCKET]"
    echo "Example: $0 /opt/matlabruntime/R2023a MARIE_runner data/inputs/inp_Duke_StadiumTriangular.json"
    echo "Example with S3: $0 /opt/matlabruntime/R2023a MARIE_runner data/inputs/inp_Duke_StadiumTriangular.json s3-bucket-name"
    exit 1
fi

MCR_ROOT=$1
RUNNER_NAME=$2
INPUT_FILE=$3
S3_BUCKET=$4

SOLUTIONS_DIR="data/solutions"

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
"$RUNNER_SCRIPT" "$MCR_ROOT" "$INPUT_FILE"
RUNNER_EXIT_CODE=$?

# Check if runner executed successfully
if [ $RUNNER_EXIT_CODE -ne 0 ]; then
    echo "Error: Runner failed with exit code $RUNNER_EXIT_CODE"
    exit $RUNNER_EXIT_CODE
fi

echo "Runner completed successfully"

# Upload solutions to S3 if bucket is provided
if [ -n "$S3_BUCKET" ]; then
    # Upload files to S3
    S3_PATH="s3://${S3_BUCKET}/${RUNNER_NAME}_solutions/"
    echo ""
    echo "======================================"
    echo "Uploading solutions to S3"
    echo "======================================"
    echo "Source directory: $SOLUTIONS_DIR"
    echo "Files to upload:"
    ls -la "$SOLUTIONS_DIR"
    echo "S3 destination:   $S3_PATH"
    echo "======================================"
    echo ""

    if aws s3 sync "$SOLUTIONS_DIR" "$S3_PATH"; then
        echo ""
        echo "Solutions uploaded successfully to: $S3_PATH"

        # List uploaded files
        echo ""
        echo "Uploaded files:"
        aws s3 ls "$S3_PATH/" --recursive
    else
        echo "Warning: S3 upload failed"
        exit 1
    fi
else
    echo "No S3 bucket specified, skipping upload"
fi

