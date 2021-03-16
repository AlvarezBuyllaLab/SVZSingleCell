wget -O GSE165554_RAW.tar "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE165554&format=file"

mkdir outs
mkdir outs/g1
mkdir outs/g2
mkdir seurat

tar -xvf GSE165554_RAW.tar GSM5039270_G1_CellRanger_outs_barcodes.tsv.gz GSM5039270_G1_CellRanger_outs_features.tsv.gz GSM5039270_G1_CellRanger_outs_matrix.mtx.gz -C /outs/G1
tar -xvf GSE165554_RAW.tar GSM5039270_G2_CellRanger_outs_barcodes.tsv.gz GSM5039270_G2_CellRanger_outs_features.tsv.gz GSM5039270_G2_CellRanger_outs_matrix.mtx.gz -C /outs/G2
tar -xvf GSE165554_RAW.tar GSM5039270_scSeq.rds.gz -C /seurat
