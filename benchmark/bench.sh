#!/bin/bash

# --- Configuration ---
RUNS=10
DATA_FILE="$(pwd -P)/../data/reads.fasta"

# Commands to benchmark
# We silence standard output (> /dev/null) to keep the terminal clean
CMD_MINIA="pixi run minia -in $DATA_FILE"
CMD_RAMILASS="pixi run ramilass $DATA_FILE bench --debug --fuse"

# Check for GNU time
TIME_CMD="/usr/bin/time"
if [ ! -f "$TIME_CMD" ]; then
    echo "Error: GNU time not found at $TIME_CMD."
    echo "Please install it (e.g., 'sudo apt install time' or 'brew install gnu-time')."
    exit 1
fi

# --- Benchmarking Function ---
run_benchmark() {
    local tool_name="$1"
    local cmd_string="$2"
    local tmp_file="${tool_name}_metrics.txt"

    echo "---------------------------------------------------"
    echo "Starting benchmark for: $tool_name ($RUNS runs)"
    echo "${cmd_string}"
    echo "---------------------------------------------------"

    # Clear previous metrics file
    > "$tmp_file"

    total_start_time=$(date +%s%N)

    for ((i=1; i<=RUNS; i++)); do
        # Print progress bar
        echo -ne "Run $i/$RUNS... \r"

        # Run command wrapped in /usr/bin/time
        # -f "%e %M" outputs: Elapsed_Seconds Max_RSS_KB
        # -o appends stats to our tmp file
        # stdout is sent to /dev/null to hide tool output
        $TIME_CMD -f "%e %M" -o "$tmp_file" --append $cmd_string > /dev/null 2>&1
    done

    total_end_time=$(date +%s%N)
    echo -e "\nDone."

    # --- Calculate Statistics using awk ---
    # We read the tmp_file which contains lines like: "0.54 24000"
    awk -v name="$tool_name" '
    BEGIN {
        sum_time=0; max_mem=0; sum_mem=0;
    }
    {
        # $1 is time (seconds), $2 is memory (KB)
        sum_time += $1;
        sum_mem += $2;
        if ($2 > max_mem) max_mem = $2;
    }
    END {
        avg_time = sum_time / NR;
        avg_mem_mb = (sum_mem / NR) / 1024;
        max_mem_mb = max_mem / 1024;

        print "\n📊 Results for " name ":";
        printf "  Total Wall Time (%d runs): %.2f seconds\n", NR, sum_time;
        printf "  Average Time per run:     %.4f seconds\n", avg_time;
        printf "  Average Peak Memory:      %.2f MB\n", avg_mem_mb;
        printf "  Max Peak Memory:          %.2f MB\n", max_mem_mb;
    }' "$tmp_file"

    # Clean up
    rm "$tmp_file"
}

# --- Execution ---

echo "Benchmark initialized. Running $RUNS iterations per tool."

# 1. Run Minia
run_benchmark "Minia" "$CMD_MINIA"

echo ""

# 2. Run Ramilass (Your Assembly)
run_benchmark "Ramilass" "$CMD_RAMILASS"

echo "---------------------------------------------------"
echo "Benchmark Complete."