"""
Snakefile to perform trimming, alignment and collect statistics on fastQ files
generated from bisulfite-treated cfDNA libraries.
"""


################################################################################
# Config file and setting parameters
################################################################################


samples = ['test']
AUTOSOMALCHROMO="chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22"
BBDUK='resources/bbmap/bbduk.sh'
BISMARK='resources/bismark/Bismark-0.24.0/bismark'
RMDUPS='resources/bismark/Bismark-0.24.0/deduplicate_bismark'
METHREF='resources/hg19'
METHEXT='resources/bismark/Bismark-0.24.0/bismark_methylation_extractor'
ADAPTERS='resources/bbmap/resources/adapters.fa'

rule all:
    input:
        expand('sample_output/tissues_of_origin/{sample}.tsv', sample = samples)
        
    

################################################################################
################################################################################
#
#
# SAMPLE PROCESSING STEPS
#
#
################################################################################
################################################################################

################################################################################
# Processing raw fastQ
################################################################################
rule trim:
    input:
        r1 = 'data/{sample}_R1.fastq.gz',
        r2 = 'data/{sample}_R2.fastq.gz'
    output:
        r1p = 'sample_output/trim/{sample}_R1_trim.fastq',
        r2p = 'sample_output/trim/{sample}_R2_trim.fastq'
    threads: 2
    log: 'sample_output/logs/trim/{sample}.trim.log'
    shell:
        """
        {BBDUK} in1={input.r1} \
                in2={input.r2} \
                out1={output.r1p} \
                out2={output.r2p} \
                -Xmx1g -threads={threads} \
                ref={ADAPTERS} \
                tbo tpe maq=10 entropy=0.25 &>{log}
        """

################################################################################
# Alignment to hg19
################################################################################
rule alignment:
	input:
		r1p = 'sample_output/trim/{sample}_R1_trim.fastq',
		r2p = 'sample_output/trim/{sample}_R2_trim.fastq'
	output:
		bam = 'sample_output/aligned/raw_aligned/{sample}.bam'
	log: 'sample_output/logs/alignment/{sample}.alignment.log'
	threads: 2
	params:
		outdir = 'sample_output/aligned/raw_aligned/'
	shell:
		"""
		{BISMARK} --genome {METHREF} \
					--parallel {threads} \
					--quiet \
					-o {params.outdir} \
					-1 {input.r1p} \
					-2 {input.r2p}
		mv {params.outdir}{wildcards.sample}_R1_trim_bismark_bt2_pe.bam {output.bam}
		mv {params.outdir}{wildcards.sample}_R1_trim_bismark_bt2_PE_report.txt {log}
		"""

rule filter_bam_get_statistics:
	input:
		bam ='sample_output/aligned/raw_aligned/{sample}.bam'
	output:
		sorted_bam = temp('sample_output/aligned/raw_aligned/{sample}.sorted.bam'),
		bismark_dup = temp('sample_output/aligned/raw_aligned/{sample}.deduplicated.bam'),
		mapped_all_chr='sample_output/aligned/all_chr/{sample}_mapped_all_chr.bam',
		mapped_all_chr_bai = 'sample_output/aligned/all_chr/{sample}_mapped_all_chr.bam.bai',
		mapped_autosomal='sample_output/aligned/autosomal/{sample}_mapped_autosomal.bam',
		mapped_autosomal_bai = 'sample_output/aligned/autosomal/{sample}_mapped_autosomal.bam.bai',
		mapped_chrM = 'sample_output/aligned/chrM/{sample}_mapped_chrM.bam',
		mapped_chrM_bai = 'sample_output/aligned/chrM/{sample}_mapped_chrM.bam.bai',
		dup_stats = 'sample_output/deduplication/{sample}.txt',
		name_sorted = 'sample_output/aligned/autosomal/{sample}_mapped_autosomal_namesorted.bam'
	threads: 1
	log: 'sample_output/logs/deduplication/{sample}.log'
	params:
		mapQ='10',
		outdir = 'sample_output/aligned/raw_aligned/'
	shell:
		"""
		samtools view -f3 -F 256,512 {input.bam} -h -o - | samtools sort -n -@ {threads} - -o {output.sorted_bam}
		{RMDUPS} -p -o {wildcards.sample} --output_dir {params.outdir} --bam {output.sorted_bam} &> {log}
		num_reads_analyzed=$(grep 'Total number of alignments analysed' {log} | cut -f2)
		num_reads_kept=$(grep 'Total count of deduplicated leftover sequences' {log} | cut -f7 -d' ')
		echo -e "total_reads\treads_kept" > {output.dup_stats}
		echo -e "$num_reads_analyzed\t$num_reads_kept" >> {output.dup_stats}
		samtools sort -@ {threads} {output.bismark_dup} -o {output.mapped_all_chr}
		samtools index {output.mapped_all_chr}
		samtools view -@ {threads} -b -q {params.mapQ} {output.mapped_all_chr} {AUTOSOMALCHROMO} -o {output.mapped_autosomal}
		samtools index {output.mapped_autosomal}
		samtools view -@ {threads} -b -q {params.mapQ} {output.mapped_all_chr} chrM -o {output.mapped_chrM}
		samtools index {output.mapped_chrM}
		samtools sort -@ {threads} -n {output.mapped_autosomal} -o {output.name_sorted}
		"""

################################################################################
# Methylation extraction
# note: requires methylation references to have been prepared.
################################################################################
rule methylation_extraction:
	input:
		mapped_autosomal= 'sample_output/aligned/autosomal/{sample}_mapped_autosomal.bam',
		name_sorted = 'sample_output/aligned/autosomal/{sample}_mapped_autosomal_namesorted.bam'
	output:
		CpG_bg='sample_output/methylation_extraction/{sample}.bedGraph',
		CpG_bg_gz = temp('sample_output/methylation_extraction/{sample}_mapped_autosomal_namesorted.bedGraph.gz'),
		CpG_bismark = 'sample_output/methylation_extraction/{sample}.bismark.cov.gz',
		mbias = 'sample_output/methylation_extraction/mbias/{sample}.M-bias.txt',
		CHGOB = 'sample_output/methylation_extraction/CHG/{sample}_CHG_OB.txt.gz',
		CHGOT = 'sample_output/methylation_extraction/CHG/{sample}_CHG_OT.txt.gz',
		CHHOB = 'sample_output/methylation_extraction/CHH/{sample}_CHH_OB.txt.gz',
		CHHOT = 'sample_output/methylation_extraction/CHH/{sample}_CHH_OT.txt.gz',
		CpGOB = 'sample_output/methylation_extraction/CpG/{sample}_CpG_OB.txt.gz',
		CpGOT = 'sample_output/methylation_extraction/CpG/{sample}_CpG_OT.txt.gz',
		log = 'sample_output/logs/methylation_extraction/{sample}.methylation_extraction.log'
	threads: 8
	params:
		outdir = 'sample_output/methylation_extraction/'
	shell:
		"""
		extra_params="--gzip --ignore 10 --ignore_r2 5"
		echo $extra_params
		{METHEXT} --parallel {threads} \
					-p \
					--bedGraph \
					--genome {METHREF} \
					-o {params.outdir} \
					$extra_params \
					{input.name_sorted}
		gunzip -c {params.outdir}{wildcards.sample}_mapped_autosomal_namesorted.bedGraph.gz | sort-bed - > {output.CpG_bg}
		mv {params.outdir}{wildcards.sample}_mapped_autosomal_namesorted.bismark.cov.gz {output.CpG_bismark}
		mv {params.outdir}{wildcards.sample}_mapped_autosomal_namesorted_splitting_report.txt {output.log}
		mv {params.outdir}{wildcards.sample}_mapped_autosomal_namesorted.M-bias.txt {output.mbias}
		mv {params.outdir}CHG_OT_{wildcards.sample}_mapped_autosomal_namesorted.txt.gz {output.CHGOT}
		mv {params.outdir}CHG_OB_{wildcards.sample}_mapped_autosomal_namesorted.txt.gz {output.CHGOB}
		mv {params.outdir}CHH_OT_{wildcards.sample}_mapped_autosomal_namesorted.txt.gz {output.CHHOT}
		mv {params.outdir}CHH_OB_{wildcards.sample}_mapped_autosomal_namesorted.txt.gz {output.CHHOB}
		mv {params.outdir}CpG_OT_{wildcards.sample}_mapped_autosomal_namesorted.txt.gz {output.CpGOT}
		mv {params.outdir}CpG_OB_{wildcards.sample}_mapped_autosomal_namesorted.txt.gz {output.CpGOB}
		"""

rule binned_methylation:
	input:
		CpG_bg='sample_output/methylation_extraction/{sample}.bedGraph',
		golden_bed = '/gpfs/commons/groups/landau_lab/alcheng/helping_rougvielab/resources/methylmatrix/regions'
	output:
		CpG_bg_tmp=temp('sample_output/aligned/autosomal/{sample}_mapped_autosomal_CpG.bedGraph.tmp'),
		intersected_golden_sample_tmp = temp('sample_output/{sample}.intersect.golden.tmp'),
		binned_CpG_DMR_golden = 'sample_output/binned_samples/golden_markers/{sample}'
	shell:
		"""
		tail -n +2 {input.CpG_bg} | bedtools sort -i - > {output.CpG_bg_tmp}
		bedtools intersect -wo -a {input.golden_bed} -b {output.CpG_bg_tmp} -sorted |
			awk '$6-$5==1 {{print $0}}' | awk 'NF{{NF-=1}};1' > {output.intersected_golden_sample_tmp}
		Rscript Bin/aggregate_over_regions.R {output.intersected_golden_sample_tmp} {output.binned_CpG_DMR_golden}
		"""

rule tissue_of_origin:
	input:
		binned_CpG_DMR_golden = 'sample_output/binned_samples/golden_markers/{sample}',
		reference_methylomes ='/gpfs/commons/groups/landau_lab/alcheng/helping_rougvielab/resources/methylmatrix/MethylMatrix_binned',
		lookup = '/gpfs/commons/groups/landau_lab/alcheng/helping_rougvielab/resources/methylmatrix/lookup_table.txt'
	output:
		mp = 'sample_output/tissues_of_origin/{sample}.tsv'
	shell:
		"""
		Rscript Bin/tissues_of_origin.R \
			{input.binned_CpG_DMR_golden} \
			{input.reference_methylomes} \
			{output.mp} \
			{input.lookup} \
			{wildcards.sample}
		"""