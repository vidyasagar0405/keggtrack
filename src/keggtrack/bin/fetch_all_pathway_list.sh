#!/usr/bin/env bash

ORG="hsa"
OUTDIR="kgml/${ORG}"

mkdir -p ${OUTDIR}

wget --output-document ${OUTDIR}/pathway_list_raw.tsv https://rest.kegg.jp/link/${ORG}/pathway/
cut -f1 ${OUTDIR}/pathway_list_raw.tsv | sort -u > ${OUTDIR}/pathway_list.tsv

for ID in $(cat ${OUTDIR}/pathway_list.tsv);
do
    wget --output-document ${OUTDIR}/"${ID}".kgml http://rest.kegg.jp/get/"${ID}"/kgml
    sleep 0.3333333333333333
    echo -e "Fetched ${ID} pathway"
done
