featureCounts -a Annotations/atRTD3_TS_21Feb22_transfix.gtf \
      -o counts_genes.txt \
      -T 8 \
      -s 0 \
      -t exon \
      -g gene_id \
      merged_alignments.bam
