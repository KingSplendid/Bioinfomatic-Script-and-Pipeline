#Generate sam file(当比对模式是singleend的时候，test.fq表示唯一的fastq)
soft/bwa-0.7.15/bwa mem ref/hg19.fasta test.fq -R '@RG\tID:test\tLB:test\tPL:ILLUMINA\tSM:test' > file_out2.sam
#Reorder Sam file
java -jar /home/lijunhui/software/bwa_test/soft/picard.jar ReorderSam I=file_out2.sam O=file_out3_reorder.sam REFERENCE=ref/hg19.fasta
#SAM trans to BAM
/home/lijunhui/software/bwa_test/soft/samtools-1.3.1/samtools view -bS file_out3_reorder.sam -o file_out4_sam2bam.bam
#Sort sam file
java -jar /home/lijunhui/software/bwa_test/soft/picard.jar SortSam INPUT=file_out4_sam2bam.bam OUTPUT=file_out5_bam_sort.bam SORT_ORDER='coordinate'
#Build index
/home/lijunhui/software/bwa_test/soft/samtools-1.3.1/samtools index file_out5_bam_sort.bam
#Call SNP indel
java -jar /home/lijunhui/software/bwa_test/soft/GenomeAnalysisTK.jar -T 'UnifiedGenotyper' -glm 'BOTH' -R ref/hg19.fasta -I file_out5_bam_sort.bam -D /home/lijunhui/software/bwa_test/soft/dbsnp_138.hg19.vcf -o file_out6.vcf -stand_call_conf 10 -stand_emit_conf 30 -minIndelCnt 5 -minIndelFrac 0.05
#Get depth summary result
java -jar /home/lijunhui/software/bwa_test/soft/GenomeAnalysisTK.jar -R ref/hg19.fasta -T 'DepthOfCoverage' -ct 0 -ct 1 -ct 5 -ct 10 -ct 15 -ct 20 -ct 25 -ct 30 -I file_out5_bam_sort.bam -o test.depth

