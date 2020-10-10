ref=/public/agis/huangsanwen_group/fengshuangshuang/huyong/Ref_genome/V4/S_lycopersicum_chromosomes.4.00.fa

reads_path=/public/agis/huangsanwen_group/fengshuangshuang/huyong/20200527_tomato_haploid_dmp
reads_TS_1_1=${reads_path}/TS-1/100505_I331_FC61K8DAAXX_L2_LYChxvRAADIAAPE_1.fq.gz
reads_TS_1_2=${reads_path}/TS-1/100505_I331_FC61K8DAAXX_L2_LYChxvRAADIAAPE_2.fq.gz
reads_TS_43_1=${reads_path}/TS-43/TS43_L7_I103.R1.clean.fastq.gz
reads_TS_43_2=${reads_path}/TS-43/TS43_L7_I103.R2.clean.fastq.gz
reads=(${reads_TS_1_1} ${reads_TS_1_2} ${reads_TS_43_1} ${reads_TS_43_2})
created_path=/public/agis/huangsanwen_group/fengshuangshuang/huyong/tem_file
SM=(I331 0 I103)
ID=(TS-1 0 TS-43)

i=0
while((${i}<${#reads[*]}))
do
        bwa mem -t 5 -R "@RG\tID:${ID[i]}\tSM:${SM[i]}\tPL:illumina" ${ref} ${reads[i]} ${reads[i+1]} > ${created_path}/${reads[i]:83:4}.aligned_reads.sam
        samtools sort ${created_path}/${reads[i]:83:4}.aligned_reads.sam -o ${created_path}/${reads[i]:83:4}.sort.bam
        samtools rmdup ${created_path}/${reads[i]:83:4}.sort.bam ${created_path}/${reads[i]:83:4}.sort.rmduped.bam
        samtools index ${created_path}/${reads[i]:83:4}.sort.rmduped.bam
        /public/agis/huangsanwen_group/fengshuangshuang/huyong/software/gatk-4.1.7.0/gatk --java-options "-Xms11g" HaplotypeCaller -R ${ref} -I ${created_path}/${reads[i]:83:4}.sort.rmduped.bam -O ${created_path}/${reads[i]:83:4}.g.vcf.gz -ERC GVCF
        let 'i=i+2'
done

i=0
/public/agis/huangsanwen_group/fengshuangshuang/huyong/software/gatk-4.1.7.0/gatk --java-options "-Xmx50g" CombineGVCFs -R ${ref} --variant ${created_path}/${reads[i]:83:4}.g.vcf.gz --variant ${created_path}/${reads[i+2]:83:4}.g.vcf.gz --variant ${created_path}/1255.g.vcf.gz --variant ${created_path}/1452.g.vcf.gz --variant ${created_path}/1746.g.vcf.gz --variant ${created_path}/198_.g.vcf.gz -O ${created_path}/cohort.g.vcf.gz
/public/agis/huangsanwen_group/fengshuangshuang/huyong/software/gatk-4.1.7.0/gatk --java-options "-Xmx50g" GenotypeGVCFs -R ${ref} -V ${created_path}/cohort.g.vcf.gz -O ${created_path}/cohort.output.vcf.gz
/public/agis/huangsanwen_group/fengshuangshuang/huyong/software/gatk-4.1.7.0/gatk --java-options "-Xmx50g" SelectVariants -R ${ref} -V ${created_path}/cohort.output.vcf.gz --select-type-to-include SNP -O ${created_path}/output.vcf.gz
/public/agis/huangsanwen_group/fengshuangshuang/huyong/software/gatk-4.1.7.0/gatk --java-options VariantsToTable -R ${ref} -O ${created_path}/output.vcf.gz -F CHROM -F POS -F REF -F ALT -GF GT -O ${created_path}/SNP.output.table
~