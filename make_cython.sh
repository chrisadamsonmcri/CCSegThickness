#!/bin/bash

FILES="Streamline2DCythonSetup.py"

for i in $FILES
do
	python3 $i build_ext --inplace
done
