#Generate sam file(当比对模式是singleend的时候，test.fq表示唯一的fastq)
bwa mem hg19.fasta test.fq -R '@RG\tID:test\tLB:test\tPL:ILLUMINA\tSM:test' > file_out2.sam
#Reorder Sam file
java -jar picard.jar ReorderSam I=file_out2.sam O=file_out3_reorder.sam REFERENCE=hg19.fasta
#SAM trans to BAM
samtools view -bS file_out3_reorder.sam -o file_out4_sam2bam.bam
#Sort sam file
java -jar picard.jar SortSam INPUT=file_out4_sam2bam.bam OUTPUT=file_out5_bam_sort.bam SORT_ORDER='coordinate'
#Build index
samtools index file_out5_bam_sort.bam
#Call SNP indel
java -jar GenomeAnalysisTK.jar -T 'UnifiedGenotyper' -glm 'BOTH' -R hg19.fasta -I file_out5_bam_sort.bam -D dbsnp_138.hg19.vcf -o file_out6.vcf -stand_call_conf 10 -stand_emit_conf 30 -minIndelCnt 5 -minIndelFrac 0.05
#Get depth summary result
java -jar GenomeAnalysisTK.jar -R hg19.fasta -T 'DepthOfCoverage' -ct 0 -ct 1 -ct 5 -ct 10 -ct 15 -ct 20 -ct 25 -ct 30 -I file_out5_bam_sort.bam -o test.depth

