wget -O GSE165554_RAW.tar "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE165554&format=file"

tar -xvf GSE165554_RAW.tar

mkdir outs
mkdir outs/g1
mkdir outs/g2
mkdir seurat

mv GSM5039270_G1_CellRanger_outs_barcodes.tsv.gz outs/g1
mv GSM5039270_G1_CellRanger_outs_features.tsv.gz outs/g1
mv GSM5039270_G1_CellRanger_outs_matrix.mtx.gz outs/g1

mv GSM5039270_G2_CellRanger_outs_barcodes.tsv.gz outs/g2
mv GSM5039270_G2_CellRanger_outs_features.tsv.gz outs/g2
mv GSM5039270_G2_CellRanger_outs_matrix.mtx.gz outs/g2

mv GSM5039270_scSeq.rds.gz seurat
