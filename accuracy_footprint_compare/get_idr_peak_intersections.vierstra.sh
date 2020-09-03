bedtools intersect -wo -a /srv/scratch/annashch/chrombpnet/accuracy_footprint_compare/vierstra/K562.DS16924.interval.all.fps.0.0001.bed -b /srv/scratch/annashch/chrombpnet/k562_dnase/bpnet/interpret/K562.dnase.idr.optimal_peak.narrowPeak | cut -f1,2,3,6,7,8,15 | uniq > vierstra_idr_intersections/K562.vierstra.idr.intersections.bed 


#bedtools intersect -wo -a /srv/scratch/annashch/chrombpnet/accuracy_footprint_compare/vierstra/HEPG2.DS24838.interval.all.fps.0.01.bed -b /srv/scratch/annashch/chrombpnet/hepg2_dnase/interpret/HEPG2.dnase.idr.optimal_peak.narrowPeak | cut -f1,2,3,6,7,8,15 | uniq > vierstra_idr_intersections/HEPG2.vierstra.idr.intersections.bed 

#/srv/scratch/annashch/chrombpnet/accuracy_footprint_compare/vierstra/GM12865.DS12436.interval.all.fps.0.01.bed.gz
#/srv/scratch/annashch/chrombpnet/accuracy_footprint_compare/vierstra/IMR90.DS13219.interval.all.fps.0.01.bed.gz
