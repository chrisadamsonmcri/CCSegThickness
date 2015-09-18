#!/bin/bash

FILES="Streamline2DCythonSetup.py"

for i in $FILES
do
	python $i build_ext --inplace
done
