#!/bin/bash

#OUTFILE=../python-dist.zip

#rm -f $OUTFILE
#zip -9r $OUTFILE CCSegManedit CCSegPipe.py CCSegProcess CCSegStatResultDisplay CCSegStatTest CCSegUtils.py LaplaceThicknessMethod.py LKTracker.py Streamline2DCython.c Streamline2DCython.pyx Streamline2DCythonSetup.py Streamline2DCython.so ART

TODAY=`date +%Y%m%d`

OUTFILE=/home/addo/AARNETownCloud/cc_seg/releases/${TODAY}.tar.xz

#rm -fr matlab
#mkdir -p matlab

#cp ../cc_seg_stattest_display_results.m ../cc_seg_tube_axes_display.m ../cc_seg_tube_load_cdata.m ../cc_seg_tube_p_display.m ../tube_surface.m ../../exportfig.m ../ideal_callosum_interp_points_curved.mat matlab

rm -fr doc
mkdir -p doc

cp ~/AARNETownCloud/cc_seg/doc/user_guide.pdf doc

rm -f ${OUTFILE}
tar cp CCSegThicknessToCSV CCSegManedit CCSegPipe*.py CCSegProcess CCSegStatResultDisplay CCSegStatTest CCSegUtils.py CCSegSubjectSelection.py LaplaceThicknessMethod.py LKTracker.py Streamline2DCython.c Streamline2DCython.pyx Streamline2DCythonSetup.py Streamline2DCython.so Otsu.py ART data doc FLIRT.py make_cython.sh matlab | 7z a -txz -m0=lzma2 -mx=9 -si ${OUTFILE}

#rm -fr matlab
rm -fr doc
