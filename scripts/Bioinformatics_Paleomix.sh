#Processing of sequences using Paleomics workflow

##1.Pre-mapping process
###Merging reads for all lanes for R1 and R2 sequences set
while read sample; do
    cat "$sample"_*_R2_001.fastq.gz > "$sample"_R2.fastq.gz
 done < ID2.txt

 while read sample; do
    cat "$sample"_*_R1_001.fastq.gz > "$sample"_R1.fastq.gz
 done < ID1.txt

###Removal of adapters and collapse for paired end reads
for infile in $(pwd)/*.fastq.gz
do
bname=$(basename $infile)
echo $bname
bname2=$(echo $bname | cut -f1-7 -d-| cut -f1-2 -d_) #to be adjusted to file names

AdapterRemoval --file1 "$bname2"_R1.fastq.gz --file2  "$bname2"_R2.fastq.gz --mm 3 --minlength 30 --basename $bname2 --trimns --trimqualities --minquality 30 --adapter1 AGATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNNNATCTCGTATGCCGTCTTCTGCTTG --adapter2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTNNNNNNNNGTGTAGATCTCGGTGGTCGCCGTATCATT --collapse;

###Merging outputs 
cat $bname.collapsed $bname.pair1.truncated $bname.pair2.truncated > $bname.merged.fq

done

###Removing poly tails
for file in $(pwd)/*_merged.fq
do
bname=$(basename $file)
echo $bname
bname2=$(echo $bname | cut -f1-7 -d- | cut -f1-2 -d_) #to be adjusted to file names
echo $bname2

echo Step 1. Removing poly A tails
fastq-grep -v "AAAAA$" $file > kmer_$bname
echo Step 2. Removing reverse complemented A tails
fastq-grep -v "^TTTTT" kmer_$bname > kmer2_$bname 
echo Step 3. Removing poly G tails
fastq-grep -v "GGGGG$" kmer2_$bname > kmer3_$bname
echo Step 4. Removing reverse complemented C tails
fastq-grep -v "^CCCCC" kmer3_$bname > kmer4_$bname 

echo Step 5. Removing rememnants adapter sequence 1 = AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC
fastq-grep -v "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC" kmer4_$bname > adap1_kmer2_$bname
echo Step 6. Removing remnants adapter sequence 2 = ATCTCGTATGCCGTCTTCTGCTTG
fastq-grep -v "ATCTCGTATGCCGTCTTCTGCTTG" adap1_kmer2_$bname > ${bname2}.pp.rmdup.fq
echo Step 10. Calculating read length distribution and outputting file
awk '{if(NR%4==2) print length($1)}' ${bname2}.pp.rmdup.fq | sort -n | uniq -c > ${bname2}.pp.rmdup.fq.read_length.txt

done

##2.Performing competitive mapping
for file in /path/to/your/files/*.rmdup.fq
do
bname=$(basename "$file" .rmdup.fq)
db=Database_name_all_mitogenomes

bowtie2 -x $db -U $file --very-sensitive --threads 20 --no-unal | samtools view -bq 30 -o ${bname}.bam

samtools sort -O BAM -o sort_${bname}.bam ${bname}.bam
done

##3.Comparing mapping performances of Bowtie2 and BWA for separate mapping
###a.Mapping with Bowtie2
for file in /path/to/your/files/*.rmdup.fq
do
bname=$(basename "$file" .rmdup.fq)

db=Database_name_species_mitogenomes
bowtie2 -x $db -U $file --very-sensitive --threads 20 --no-unal | samtools view -bq 30 -o ${bname}.bam
samtools sort -O BAM -o sort_${bname}.bam ${bname}.bam

###b.Mapping with BWA
echo $bname
bwa aln -l 1024 -n 0.001 -t 10 file_species_mitogenome.fasta $file | bwa samse file_species_mitogenome.fasta - $file | samtools view -F 4 -q 30 -@ 10 -uS - | samtools sort -@ 10 -o /path/sort_bwa_${bname}.bam
done

##4. Post-processing and statistics
###Removal of duplicated with Picard
for file in sort*
do
output_file="no_dupl_${file}"
metrics_file="markdupl_${file}.txt"

java -jar /path/to/picard.jar MarkDuplicates \
   I="$file" \
   O="$output_file" \
   M="$metrics_file" \
   REMOVE_DUPLICATES= true
done

###Evaluating mapping statistics
for file in "$output_file"
do
    bname=$(basename "$file" .bam)

    echo $file
    bamcov -m -w0 $file | paste > bamcov_histo_${bname}.txt
    bamcov $file | paste > bamcov_table_${bname}.txt
done

##5.Generating consensus sequences for phylogenetic analyses
for file in "$output_file"
do
angsd -doFasta 2 -doCounts 1 -i $file -out $file.consensus.fa
done

gunzip *.consensus.fq.gz

for file in *.fa.fa
 do
    filename="${file%.fa.fa}"

###Renaming the headers within each consensus sequence
    seqtk seq -a "$file" | awk -v filename="$filename" '/^>/{print ">" filename "." ++i; next} 1' > "${filename}_renamed.fasta"
done
