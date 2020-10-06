#!/bin/bash
bedtools intersect -wo -a idr.optimal_peak.narrowPeak.gz -b ../bQTL/junD.txt.bed | bedtools sort -i - | cut -f1,2,3,14,15 | sort |uniq > bQTL_intersections/junD.bQTL.bed
bedtools intersect -wo -a idr.optimal_peak.narrowPeak.gz -b ../bQTL/nfkB.txt.bed | bedtools sort -i - | cut -f1,2,3,14,15 | sort |uniq > bQTL_intersections/nfkB.bQTL.bed
bedtools intersect -wo -a idr.optimal_peak.narrowPeak.gz -b ../bQTL/oct1.txt.bed | bedtools sort -i - |  cut -f1,2,3,14,15 | sort| uniq  > bQTL_intersections/oct1.bQTL.bed
bedtools intersect -wo -a idr.optimal_peak.narrowPeak.gz -b ../bQTL/pu1.txt.bed | bedtools sort -i - | cut -f1,2,3,14,15 | sort |uniq  > bQTL_intersections/pu1.bQTL.bed
bedtools intersect -wo -a idr.optimal_peak.narrowPeak.gz -b ../bQTL/stat1.txt.bed | bedtools sort -i - | cut -f1,2,3,14,15 | sort| uniq > bQTL_intersections/stat1.bQTL.bed
