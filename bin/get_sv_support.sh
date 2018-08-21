#!/bin/bash
# =========================================================
#
# Copyright (C) 2016-2018, Nuno A. Fonseca (nuno dot fonseca at gmail
# dot com)
#
# This is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# if not, see <http://www.gnu.org/licenses/>.
#
#
# =========================================================
# usage: get_sv_support.sh sv_dir fusion_aliquot_dir distance
SV_DIR=$1
FUS_DIR=$2
DISTANCE=$3
set -e
## dependencies
# R
which R &> /dev/null || (echo "Missing dependency: R needs to be installed" && exit 1)
which sv_genefusions_overlap.R &> /dev/null || (echo "Missing dependency: sv_genefusions_overlap.R needs to be installed and in the PATH" && exit 1)

set +eux
# 
# 
if [ "$DISTANCE-" == "-" ]; then
    echo "ERROR: Missing arguments"
    echo "Usage: get_sv_support.sh sv_dir fusion_aliquot_dir distance"
    echo "sv_dir: directory containing bedpe files for each aliquot with the structural variants (filename format <aliquot_id>.sv.bedpe)."
    echo "fusion_aliquot_dir: directory with the files containing the putative gene fusions for each aliquot  (filename format <aliquot_id>.sum.tsv)."
    echo "distance: maximum allowed distance (in bp) between an overlapping SV and gene fusion"
    echo "Output files will be placed in the folder dist_<DISTANCE>"
    exit 1
fi
if [ ! -e $SV_DIR ]; then
    echo "ERROR: sv_dir $SV_DIR not found!"
    exit 1
fi
if [ ! -e $FUS_DIR ]; then
    echo "ERROR: fusion_aliquot_dir $FUS_DIR not found!"
    exit 1
fi

if [ "$DISTANCE" == "-" ]; then
    echo "ERROR: usage: get_sv_support.sh sv_dir fusion_aliquot_dir  distance"
    exit 1
fi

ofolder=dist_$DISTANCE
rm -rf $ofolder
mkdir -p $ofolder
ids=$(ls --color=never $SV_DIR/*.sv.bedpe|cut -f 4 -d/|sed "s/.sv.bedpe//g")
#epty bedpe
head -n 1 $(ls --color=never $SV_DIR/*.sv.bedpe|head -n 1) > empty_sv_bedpe
set -e
for id in $ids; do
    sv_file=$SV_DIR/$id.sv.bedpe
    if [ ! -e $sv_file ]; then
	sv_file=empty_sv_bedpe
    fi
    if [ -e $FUS_DIR/$id.sum.tsv.gz ]; then
	echo cmd: 	sv_genefusions_overlap.R $sv_file $FUS_DIR/$id.sum.tsv.gz $DISTANCE n $ofolder/$id

	sv_genefusions_overlap.R $sv_file $FUS_DIR/$id.sum.tsv.gz $DISTANCE n $ofolder/$id
    else
	echo Not found $FUS_DIR/$id.sum.tsv.gz
    fi
done
#merge_sv_fusion_matches_files.R  $ofolder $PIPE ${PIPE}_${ofolder}.tsv
#gzip -f ${PIPE}_${ofolder}.tsv
#echo Created ${PIPE}_${ofolder}.tsv.gz
#rm -rf $ofolder
exit 0
