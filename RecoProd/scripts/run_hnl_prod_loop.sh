#!/bin/bash

ABORT_ON_ERRORS=1
for s in $(seq 100 119); do
    echo "Running for seed $s"
    ./run_hnl_prod.sh $s 5000
    RESULT=$?
    if [ $RESULT -ne 0 ]; then
        echo "ERROR while running for seed $s. Aborting the loop."
        if [ $ABORT_ON_ERRORS -eq 1 ]; then
            echo "Aborting the loop."
            exit $RESULT
        fi
    fi
done
