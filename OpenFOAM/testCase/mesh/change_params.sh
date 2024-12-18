#!/bin/bash -e

find  -name decomposeParDict | xargs sed -i 's/numberOfSubdomains [0-9]*;/numberOfSubdomains 8;/g'
find  -name run_solver.sh | xargs sed -i 's/processors=[0-9]*/processors=8/g'