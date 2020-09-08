#! /bin/bash

#*****************************************************************************
#   popHyRAD: Pipeline to perform population genetic analyses from HyRAD data
#   Authors: J. Gauthier
#
#  This program is free software: you can redistribute it and/or modify
#  it under the terms of the GNU Affero General Public License as
#  published by the Free Software Foundation, either version 3 of the
#  License, or (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU Affero General Public License for more details.
#
#  You should have received a copy of the GNU Affero General Public License
#  along with this program.  If not, see <http://www.gnu.org/licenses/>.
#*****************************************************************************




function help {
	echo "popHyRAD: Pipeline to perform population genetic analyses from HyRAD data"
	echo "Usage: ./popHyRAD.sh -r loci_from_ipyrad -m [OPTIONS]"
}

while getopts ":r:" opt; do
	case $opt in
		r)
			echo "use reference=$OPTARG" >&2
			reference=$OPTARG
			;;
		m)
			echo "use mapping tool=$OPTARG" >&2
			mapping_tool=$OPTARG
			;;
	esac
done



cat "$reference" | tr '\n' '\t' | tr '|' '\n' | grep "/" | sed -e 's/^\t//g' | awk '{print NR"\t"$0}' > temp_file
while read a ; do name=`echo "$a" | awk '{print $1}'` ; echo "$a" | tr '\t' '\n' | awk '{print ">locus_"'$name'"_"$0}' | tr ' ' '\n' | grep "." | grep "_" -B 1 >> temp_all_seq_locus.fasta ; done < temp_file
sed -e 's/-//g' temp_all_seq_locus.fasta > "$reference"_allseq_locus.fasta
rm temp_*


grep ">" "$reference"_allseq_locus.fasta | awk -F "_" '{print $1"_"$2"_"}' | sed -e 's/>//g' > list_locus
fastalength "$reference"_allseq_locus.fasta > "$reference"_allseq_locus.length
while read a ; do grep "$a" "$reference"_allseq_locus.length | sort -nrk 2 | head -n 1 | awk '{print $1}' >> "$reference"_allseq_locus.longest ; done < list_locus
./scripts/fastaselect.pl "$reference"_allseq_locus.longest "$reference"_allseq_locus.fasta > "$reference"_locus_uniq.fasta


bwa index "$reference"_locus_uniq.fasta


for i in `ls *_R1.fastq`
		do
		sample=`echo $i | sed -e 's/_R1.fastq//g'`
		if [ "$mapping_tool" == "mem" ]
		then
			bwa mem "$reference"_locus_uniq.fasta "$sample"_R1.fastq "$sample"_R2.fastq > "$nameR2"_on_ref.sam
		elif [ "$mapping_tool" == "aln" ]
		then
			bwa aln "$reference"_locus_uniq.fasta "$sample"_R1.fastq > "$i"_aln_R1.sai
			bwa aln "$reference"_locus_uniq.fasta "$sample"_R2.fastq > "$i"_aln_R2.sai
			bwa sampe "$reference"_locus_uniq.fasta "$i"_aln_R1.sai "$i"_aln_R2.sai "$sample"_R1.fastq "$sample"_R1.fastq> "$i"_on_ref.sam
			rm *.sai
		fi
		name=`echo "$i"_on_ref.sam | sed -e 's/\.sam//g'`
		samtools sort "$i"_on_ref.sam -o temp_1_sorted.bam
		samtools view -bF 4 temp_1_sorted.bam > temp_1_sorted_keep.bam
		java -jar -mx128G /home/jeremy/local/envpicard/share/picard-2.20.2-0/picard.jar AddOrReplaceReadGroups I=temp_1_sorted_keep.bam O=temp_1_sorted_keep_rg.bam ID=["$sample"] RGLB=[id] PL=[pl] PU=[pu] SM=["$sample"]
		GenomeAnalysisTK -T RealignerTargetCreator -I temp_1_sorted_keep_rg.bam -R "$reference"_locus_uniq.fasta -o temp.intervals
		GenomeAnalysisTK -T IndelRealigner -I temp_1_sorted_keep_rg.bam -R "$i".fasta -targetIntervals temp.intervals -o temp_1_sorted_keep_rg_realign.bam
		java -jar -mx128G /home/jeremy/local/envpicard/share/picard-2.20.2-0/picard.jar MarkDuplicates I=temp_1_sorted_keep_rg_realign.bam O=temp_1_sorted_keep_pcrdup.bam M=marks REMOVE_DUPLICATES=true
		mapDamage -i temp_1_sorted_keep_pcrdup.bam -r "$reference"_locus_uniq.fasta --rescale --merge-reference-sequences
		cp ./results_temp_1_sorted_keep_pcrdup/temp_1_sorted_keep_pcrdup.rescaled.bam ./"$sample"_sorted_keep_pcrdup.rescaled.bam
		rm temp_*
		done

freebayes -f "$reference"_locus_uniq.fasta *.rescaled.bam > snpset_freebayes.vcf



