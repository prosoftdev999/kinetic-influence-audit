#!/bin/sh
set -u
mkdir -p /logs/verifier
printf '0' > /logs/verifier/reward.txt
if pytest -q /tests/test_verify.py --ctrf=/logs/verifier/ctrf.json; then
  printf '1' > /logs/verifier/reward.txt
fi
exit 0
