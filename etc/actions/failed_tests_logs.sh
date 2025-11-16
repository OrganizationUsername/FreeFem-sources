#!/usr/bin/env bash

# This script checks if tests have failed. If yes, it prints the content of the
# failed tests log, and returns 1.
#
# It allows keeping the '-i' option for 'make check', print the failed
# tests, and only then, return exit 1 to interrupt the CI pipeline.

FAILED_TESTS=$(grep ":test-result: FAIL" examples/*/*.trs)

if [ -z "$FAILED_TESTS" ]; then
    exit 0;
fi

NB_FAILED_TESTS=$(echo "$FAILED_TESTS" | wc -l)
echo "$NB_FAILED_TESTS tests failed:"

echo "$FAILED_TESTS" | while IFS= read -r line; do
    path="${line%::test-result: FAIL}"
    printf  "\n******** ${path%.*}.log ********\n\n"
    tail -n 100 ${path%.*}.log
done

exit 1
