#!/bin/bash
#we will search for isoforms of this contig
contig_to_search=$1;
reference=$2;
#including only putative isoforms
part=`echo $contig_to_search | sed -n "s/^\(.*\_c[0-9]\+\)\_seq[0-9]\+/\1/p"`; egrep "$part\_" $reference | grep -v $contig_to_search;