#!/bin/bash

# 1 = output bin3c tsv

sed -E "s/(.*) (contig:)(.*) (.*) (.*)\t(.*)/\3\t\6/" $1