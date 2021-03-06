#
## sv2gf  [![Dockerhub](https://img.shields.io/docker/automated/jrottenberg/ffmpeg.svg)](https://hub.docker.com/r/nunofonseca/sv2gf/tags/) [![License](http://img.shields.io/badge/license-GPL%203-brightgreen.svg?style=flat)](http://www.gnu.org/licenses/gpl-3.0.html)


Match gene fusions with structural variants.

# Requirements

Linux operating system with R (3.4 or above) installed.

# Installation

Simply copy the scripts in the bin folder to a folder in the PATH.

# How to run

## Input

Structural variant (SV) calls and gene fusion calls for a given sample. Examples are provided in the example folder.

SVs file (bedpe format) should contain the following columns: "chrom1","start1","end1","chrom2","start2","end2","sv_id","pe_support","strand1","strand2","svclass","svmethod".

The gene fusions file is expected to contain the following columns: "FusionGene","KnownGene1","KnownGene2","GeneId1","GeneId2","Strand","Chromosome1","Breakpoint1","Chromosome2","Breakpoint2","FrameShift","FusionJunctionSequence","SplicePattern","Number of Supporting Reads".

## Running

usage: sv_genefusions_overlap.R <sv.bedpe> <fusions.tsv> <distance in bp>  n <out_file>

sv.bedpe: file with the SVs
fusions.tsv: gene fusion file
distance: maximum allowed distance between an overlapping SV and gene fusion
out
out_file: output file with the gene fusions and closer SV (within the given distance) and SV support type. For an example check the file example/out

