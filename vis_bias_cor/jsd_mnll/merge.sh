python merge.py --prefix baseline bpnet tobias bpnet.tobias \
       --files hepg2.dnase.baseline.txt hepg2.dnase.bpnet.txt hepg2.dnase.tobias.txt hepg2.dnase.bpnet.tobias.txt  \
       --outf hepg2.dnase

python merge.py --prefix baseline bpnet tobias bpnet.tobias \
       --files k562.dnase.baseline.txt k562.dnase.bpnet.txt k562.dnase.tobias.txt k562.dnase.bpnet.tobias.txt  \
       --outf k562.dnase

python merge.py --prefix baseline tobias \
       --files hepg2.atac.baseline.txt hepg2.atac.tobias.txt \
       --outf hepg2.atac

python merge.py --prefix baseline tobias \
       --files gm12878.atac.baseline.txt gm12878.atac.tobias.txt \
       --outf gm12878.atac
