#!/bin/bash
for f in /oak/stanford/groups/akundaje/projects/atlas/histone_chip/caper_out/f58a0e98-3f70-4821-b380-fac5ddb37d9c/call-align/K562.merged.ENCSR000APK.gz.bam.bpnet.plus.bw /oak/stanford/groups/akundaje/projects/atlas/histone_chip/caper_out/f58a0e98-3f70-4821-b380-fac5ddb37d9c/call-align/K562.merged.ENCSR000APK.gz.bam.bpnet.minus.bw /oak/stanford/groups/akundaje/projects/atlas/histone_chip/caper_out/f58a0e98-3f70-4821-b380-fac5ddb37d9c/call-align_ctl/K562.merged.control.ENCSR000AKY.gz.bam.bpnet.unstranded.bw /oak/stanford/groups/akundaje/projects/atlas/histone_chip/caper_out/f58a0e98-3f70-4821-b380-fac5ddb37d9c/call-align_ctl/K562.merged.control.ENCSR000AKY.gz.bam.bpnet.plus.bw /oak/stanford/groups/akundaje/projects/atlas/histone_chip/caper_out/f58a0e98-3f70-4821-b380-fac5ddb37d9c/call-align_ctl/K562.merged.control.ENCSR000AKY.gz.bam.bpnet.minus.bw
do
cp $f .
done
