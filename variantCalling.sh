echo "Mapping using Bowtie2"

# Bowtie2 Index reference
bowtie2-build ${ref} ${refidx}

# Bowtie2 alignment 
for file_1 in ${read}/*_1.fq.gz; do 
	file_2="${file_1%_1.fq.gz}_2.fq.gz" 
	base_name=$(basename "$file_1" _1.fq.gz)    
      	OUTPUT="${base_name}.sam"   
      	bowtie2 -x ${refidx} -1 "$file_1" -2 "$file_2" -S "${aligned_read}/${OUTPUT}" --local -p 4 
done

echo "Convert sam to bam and sort"

for file in ${aligned_read}/*.sam; do 
    base_name=$(basename "$file" .sam) 
    samtools sort "$file" -o "${bam}/${base_name}_sorted.bam" 
done

echo "Add or replace read groups"

i=1 
for file in ${bam}/*_sorted.bam; do 
    base_name=$(basename "$file" _sorted.bam) 
    java -jar /path/to/tools/picard/picard.jar AddOrReplaceReadGroups \ 
    INPUT="$file" \ 
    OUTPUT="${bam}/${base_name}_sorted_rg.bam" \ 
    RGID=${i} \ 
    RGLB=${i} \ 
    RGPL=Illumina \ 
    RGPU=1 \ 
    RGSM="${base_name}" 
   ((i+=1)) 
done

for file in ${bam}/*_sorted_rg.bam; do 
    samtools index ${file} 
done

echo "Mark duplicates"

for file in ${bam}/*_sorted_rg.bam; do
	  base_name=$(basename "$file" .bam)
	  java -jar /path/to/tools/picard/picard.jar MarkDuplicates \
	  INPUT="${file}" \
	  OUTPUT="${bam}/${base_name}_dup.bam" \
	  METRICS_FILE="${bam}/${base_name}.metric" \
	  CREATE_INDEX=true
done

echo "Call variants - gatk HaplotypeCaller"

for file in ${bam}/*_sorted_rg_dup.bam; do 
    base_name=$(basename "${file}" _sorted_rg_dup.bam)  
    gatk HaplotypeCaller -R ${ref} -I ${file} -O "${results}/${base_name}.gvcf.gz" -ERC GVCF
done

echo "Step VI: Combine .gvcf files"

samples=$(find . -name "*.gvcf.gz" | sed 's/^/--variant /') 
gatk CombineGVCFs $(echo $samples) -O "${results}/combined.gvcf.gz" -R ${ref}

echo "Call variants - gatk GenotypeGVCF"

for file in ${results}/*.gvcf.gz; do
	base_name=$(basename "${file}" .gvcf.gz)
	gatk GenotypeGVCFs -R ${ref} -V ${file} -O "${results}/${base_name}.vcf.gz"
done

echo "Varaint filtration"

# Extract only SNPs 
gatk SelectVariants -R ${ref} -V ${results}/combined.vcf.gz --select-type SNP -O ${results}/raw_snps.vcf.gz 

# gatk Variantfiltration 
gatk VariantFiltration \ 
    -R ${ref} \ 
    -V ${results}/raw_snps.vcf.gz \ 
    -O ${results}/filtered_snps.vcf \ 
    --filter-name "QD_filter" -filter "QD < 2.0" \ 
    --filter-name "QUAL_filter" -filter "QUAL < 30.0" \ 
    --filter-name "FS_filter" -filter "FS > 60.0" \ 
    --filter-name "MQ_filter" -filter "MQ < 30.0" \ 
    --filter-name "SOR_filter" -filter "SOR > 3.0" \ 
    --filter-name "MQRankSum_filter" -filter "MQRankSum < -12.5" \ 
    --filter-name "ReadPosRankSum_filter" -filter "ReadPosRankSum < -8.0" \ 
    -genotype-filter-expression "DP < 5" \ 
    -genotype-filter-name "LowDP_filter" \ 
    -genotype-filter-expression "GQ < 10" \ 
    -genotype-filter-name "GQ_filter"

# Select Variants that PASS filters
gatk SelectVariants --exclude-filtered -V ${results}/filtered_snps.vcf -O ${results}/analysis_pass.vcf
cat ${results}/analysis_pass.vcf|grep -v -E "DP_filter|GQ_filter" > ${results}/Final-snps-filtered.vcf 
