#!/bin/bash

mex -O interp2q_linear_fast_c.c
mex -O otsu2_c.c
mex -O -largeArrayDims setup_linear_equations_2d_spacing.c
mex -O numeric_gradient2d_c.c
mex -O mystream2c_vector_bounding_box.c
mex -O cummax.c
