ref=/public/agis/huangsanwen_group/fengshuangshuang/huyong/Ref_genome/V4/S_lycopersicum_chromosomes.4.00.fa
bwa index ${ref}
/public/agis/huangsanwen_group/fengshuangshuang/huyong/software/samtools-1.10/samtools faidx ${ref}
java -jar /public/agis/huangsanwen_group/suntianshu/software/picard-tools-1.124/picard.jar CreateSequenceDictionary REFERENCE=${ref} OUTPUT=/public/agis/huangsanwen_group/fengshuangshuang/huyong/Ref_genome/V4/S_lycopersicum_chromosomes.4.00.dict
echo 'index done' >> /public/agis/huangsanwen_group/fengshuangshuang/huyong/sourced/Output_done
##Build index files

reads_path=/public/agis/huangsanwen_group/fengshuangshuang/huyong/20200527_tomato_haploid_dmp
reads_1255_1=${reads_path}/1255_FDSW202387373-1r_1.clean.fq.gz
reads_1255_2=${reads_path}/1255_FDSW202387373-1r_2.clean.fq.gz
reads_1452_1=${reads_path}/1452_FDSW202387374-1r_1.clean.fq.gz
reads_1452_2=${reads_path}/1452_FDSW202387374-1r_2.clean.fq.gz
reads_1746_1=${reads_path}/1746_FDSW202387375-1r_1.clean.fq.gz
reads_1746_2=${reads_path}/1746_FDSW202387375-1r_2.clean.fq.gz
reads_198_1=${reads_path}/198_FDSW202387372-1r_1.clean.fq.gz
reads_198_2=${reads_path}/198_FDSW202387372-1r_2.clean.fq.gz
reads=(${reads_1255_1} ${reads_1255_2} ${reads_1452_1} ${reads_1452_2} ${reads_1746_1} ${reads_1746_2} ${reads_198_1} ${reads_198_2})
created_path=/public/agis/huangsanwen_group/fengshuangshuang/huyong/tem_file
#Definition of path strings

SM=(FDSW202387373 0 FDSW202387374 0 FDSW202387375 0 FDSW202387372)
ID=(1255 0 1452 0 1746 0 198)
i=0
while((${i}<${#reads[*]}))
do
    bwa mem -t 11 -R "@RG\tID:${ID[i]}\tSM:${SM[i]}\tPL:illumina" ${ref} ${reads[i]} ${reads[i+1]} > ${created_path}/${reads[i]:83:4}.aligned_reads.sam

    /public/agis/huangsanwen_group/fengshuangshuang/huyong/software/samtools-1.10/samtools sort ${created_path}/${reads[i]:83:4}.aligned_reads.sam -o ${created_path}/${reads[i]:83:4}.sort.bam
    /public/agis/huangsanwen_group/fengshuangshuang/huyong/software/samtools-1.10/samtools rmdup ${created_path}/${reads[i]:83:4}.sort.bam ${created_path}/${reads[i]:83:4}.sort.rmduped.bam
    /public/agis/huangsanwen_group/fengshuangshuang/huyong/software/samtools-1.10/samtools index ${created_path}/${reads[i]:83:4}.sort.rmduped.bam

    /public/agis/huangsanwen_group/fengshuangshuang/huyong/software/gatk-4.1.7.0/gatk --java-options "-Xms11g" HaplotypeCaller -R ${ref} -I ${created_path}/${reads[i]:83:4}.sort.rmduped.bam -O ${created_path}/${reads[i]:83:4}.g.vcf.gz -ERC GVCF
    let "i=i+2"
done
#Maping,sort,rmove duplicate,build index after remove duplicate,call mutution

i=0
/public/agis/huangsanwen_group/fengshuangshuang/huyong/software/gatk-4.1.7.0/gatk --java-options "-Xmx50g" CombineGVCFs -R ${ref} --variant ${created_path}/${reads[i]:83:4}.g.vcf.gz --variant ${created_path}/${reads[i+2]:83:4}.g.vcf.gz --variant ${created_path}/${reads[i+4]:83:4}.g.vcf.gz --variant ${created_path}/${reads[i+6]:83:4}.g.vcf.gz -O ${created_path}/cohort.g.vcf.gz
#Runing time extremely long,due to 198_.g.vcf.gz file damaged possibly.

/public/agis/huangsanwen_group/fengshuangshuang/huyong/software/gatk-4.1.7.0/gatk --java-options "-Xmx50g" GenotypeGVCFs -R ${ref} -V ${created_path}/cohort.g.vcf.gz -O ${created_path}/output.vcf.gz
/public/agis/huangsanwen_group/fengshuangshuang/huyong/software/gatk-4.1.7.0/gatk --java-options "-Xmx50g" SelectVariants -R ${ref} -V ${created_path}/output.vcf.gz --select-type-to-include SNP -O ${created_path}/SNP.output.vcf.gz
/public/agis/huangsanwen_group/fengshuangshuang/huyong/software/gatk-4.1.7.0/gatk --java-options "-Xmx50g" VariantsToTable -R ${ref} -V ${created_path}/SNP.output.vcf.gz  -F CHROM -F POS -F REF -F ALT -GF GT -O ${created_path}/SNP.output.table
#Combine gvcf,convert to vcf,export SNP as a table file
