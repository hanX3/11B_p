#!/usr/bin/env bash

NRUNS=${1:-8}
MACRO=${2:-macros/run.mac}
DATA_DIR=../data
LOG_DIR=${DATA_DIR}/logs_parallel_$(date +%Y%m%d_%Hh%Mm%Ss)

mkdir -p "$DATA_DIR"
mkdir -p "$LOG_DIR"

echo "Run number : $NRUNS"
echo "Macro      : $MACRO"
echo "Data dir   : $DATA_DIR"
echo "Log dir    : $LOG_DIR"
echo

for i in $(seq 1 "$NRUNS"); do
    echo "Starting job $i / $NRUNS at $(date '+%F %T')"

    ./HB "$MACRO" > "${LOG_DIR}/run_${i}.log" 2>&1 &

    sleep 3
done

echo
echo "All jobs submitted. Waiting for them to finish..."
wait

echo
echo "All jobs finished. Merging reaction ROOT files..."

hadd -f "${DATA_DIR}/reaction_all.root" "${DATA_DIR}"/reaction_*.root
hadd -f "${DATA_DIR}/event_all_${TAG}.root"    "${DATA_DIR}"/event_*.root
hadd -f "${DATA_DIR}/track_all_${TAG}.root"    "${DATA_DIR}"/track_*.root

echo "Done."
echo "Merged file: ${DATA_DIR}/reaction_merged.root"
