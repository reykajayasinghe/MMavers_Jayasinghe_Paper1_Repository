#!/bin/bash

#https://pyscenic.readthedocs.io/en/latest/installation.html#command-line-interface

###############################################################################
# NOTE: The expression data needs to have genes in columns and samples in rows
###############################################################################



sample=$1
run=$2

## Set path variables
## Set path variables
transfacPath='/data/allTFs_hg38.txt'
transDB1='/data/hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather'
transDB2='/data/hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather'
annotTABLE='/data/motifs-v9-nr.hgnc-m0.001-o0.0.tbl'

#outdir=${sample}_results
#inputdata=${sample}.rawdata.tsv
inputdata=${sample}.loom
exprout=${sample}_${run}_expr_mat.adjacencies.tsv
regulonsout=${sample}_${run}_regulons.csv
aucout=${sample}_${run}_auc_mtx.loom

## Docker commands
# START
docker run --rm -t -v /diskmnt/Datasets/mmy_scratch/Dipersio/MyDrug/Analysis2/Rerun_October2023/Downstream_Analysis/scenic/data:/data aertslab/pyscenic:0.12.1 pyscenic grn \
        --num_workers 8 \
        -o /data/${exprout} \
	/data/${inputdata} \
	${transfacPath} \
	--cell_id_attribute CellID\
	--gene_attribute Gene
 
docker run --rm -t -v /diskmnt/Datasets/mmy_scratch/Dipersio/MyDrug/Analysis2/Rerun_October2023/Downstream_Analysis/scenic/data:/data aertslab/pyscenic:0.12.1 pyscenic ctx \
        /data/${exprout} \
      	${transDB1} \
	${transDB2} \
        --annotations_fname ${annotTABLE} \
        --expression_mtx_fname /data/${inputdata} \
        --mode "dask_multiprocessing" \
        --output /data/${regulonsout} \
        --num_workers 8 \
        --cell_id_attribute CellID\
        --gene_attribute Gene

docker run --rm -t -v /diskmnt/Datasets/mmy_scratch/Dipersio/MyDrug/Analysis2/Rerun_October2023/Downstream_Analysis/scenic/data:/data aertslab/pyscenic:0.12.1 pyscenic aucell \
        /data/${inputdata} \
        /data/${regulonsout} \
        -o /data/${aucout} \
        --num_workers 8


# END
