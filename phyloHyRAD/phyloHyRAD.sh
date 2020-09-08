#! /bin/bash

#*****************************************************************************
#   phyloHyRAD: Pipeline to perform phylogenetic analyses from HyRAD data
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
	echo "phyloHyRAD: Pipeline to perform phylogenetic analyses from HyRAD data"
	echo "Usage: ./phyloHyRAD.sh -r loci_from_ipyrad -m [OPTIONS]"
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


nameREF=`echo "$reference"_allseq_locus.fasta | sed -e 's/.fasta//g'`
bwa index "$reference"_allseq_locus.fasta
for j in `ls *_R1.fastq`
	do
	nameR2=`echo $j | sed -e 's/_R1.fastq//g'`
	if [ "$mapping_tool" == "mem" ]
	then
		bwa mem "$reference"_allseq_locus.fasta "$j" "$nameR2"_R2.fastq > "$nameR2"_on_ref.sam
	elif [ "$mapping_tool" == "aln" ]
	then
		bwa aln "$reference"_allseq_locus.fasta "$sample"_R1.fastq > "$i"_aln_R1.sai
		bwa aln "$reference"_allseq_locus.fasta "$sample"_R2.fastq > "$i"_aln_R2.sai
		bwa sampe "$reference"_allseq_locus.fasta "$i"_aln_R1.sai "$i"_aln_R2.sai "$sample"_R1.fastq "$sample"_R1.fastq> "$i"_on_ref.sam
	fi
	grep -v "^@" "$nameR2"_on_ref.sam | awk '{print $3}' | sort | uniq -c > "$nameR2"_on_ref.sam_nb_reads
	rm *.sai
	done

grep ">" "$reference"_allseq_locus.fasta | awk -F "_" '{print $1"_"$2"_"}' | sed -e 's/>//g' > list_locus
for i in `ls *.sam_nb_reads` ; do while read a ; do grep "$a" "$i" | sort -nrk 1 | head -n 1 | awk '{print $2}' >> "$i"_best_probe ; done < list_locus ; done

rm *.sam

for i in `ls *nb_reads_best_probe`
		do
		sample=`echo $i | sed -e 's/_on_ref.sam_nb_reads_best_probe//g'`
		./scripts/fastaselect.pl "$reference"_allseq_locus.fasta "$i" > "$i".fasta
		bwa index "$i".fasta
		if [ "$mapping_tool" == "mem" ]
		then
			bwa mem "$i".fasta "$sample"_R1.fastq "$sample"_R2.fastq > "$nameR2"_on_specific_ref.sam
		elif [ "$mapping_tool" == "aln" ]
		then
			bwa aln "$i".fasta "$sample"_R1.fastq > "$i"_aln_R1.sai
			bwa aln "$i".fasta "$sample"_R2.fastq > "$i"_aln_R2.sai
			bwa sampe "$i".fasta "$i"_aln_R1.sai "$i"_aln_R2.sai "$sample"_R1.fastq "$sample"_R1.fastq> "$i"_on_specific_ref.sam
			rm *.sai
		fi
		name=`echo "$i"_on_specific_ref.sam | sed -e 's/\.sam//g'`
		samtools sort "$i"_on_specific_ref.sam -o temp_1_sorted.bam
		samtools view -bF 4 temp_1_sorted.bam > temp_1_sorted_keep.bam
		java -jar -mx128G /home/jeremy/local/envpicard/share/picard-2.20.2-0/picard.jar AddOrReplaceReadGroups I=temp_1_sorted_keep.bam O=temp_1_sorted_keep_rg.bam ID=["$sample"] RGLB=[id] PL=[pl] PU=[pu] SM=["$sample"]
		GenomeAnalysisTK -T RealignerTargetCreator -I temp_1_sorted_keep_rg.bam -R "$i".fasta -o temp.intervals
		GenomeAnalysisTK -T IndelRealigner -I temp_1_sorted_keep_rg.bam -R "$i".fasta -targetIntervals temp.intervals -o temp_1_sorted_keep_rg_realign.bam
		java -jar -mx128G /home/jeremy/local/envpicard/share/picard-2.20.2-0/picard.jar MarkDuplicates I=temp_1_sorted_keep_rg_realign.bam O=temp_1_sorted_keep_pcrdup.bam M=marks REMOVE_DUPLICATES=true
		mapDamage -i temp_1_sorted_keep_pcrdup.bam -r "$i".fasta --rescale --merge-reference-sequences
		cp ./results_temp_1_sorted_keep_pcrdup/temp_1_sorted_keep_pcrdup.rescaled.bam ./"$sample"_sorted_keep_pcrdup.rescaled.bam
		samtools mpileup -A -uf "$i".fasta "$sample"_sorted_keep_pcrdup.rescaled.bam | bcftools call -c | vcfutils.pl vcf2fq > temp.fastq
		seqtk seq -q30 -n N -a temp.fastq > "$name"_consensus.fasta
		rm temp_*
		done


grep ">" "$reference"_allseq_locus.fasta | awk -F "_" '{print $1"_"$2"_"}' | sort | uniq > temp_list_loci

while read a
	do
	locus=`echo $a | sed -e 's/>//g'`
	for i in `ls *consensus.fasta`
		do
		sample=`echo $i | sed -e 's/_consensus.fasta//g'`
		grep -A 1 "$a" "$i" | sed 's/>.*/>'"$sample"'/' >> "$locus".fasta
		done
	mafft --auto "$locus".fasta > "$locus"_align.fasta
	rm "$locus".fasta
	done < temp_list_loci

find . -size 0 -delete

AMAS.py concat -i *_align.fasta -f fasta -d dna
mkdir all_loci_separated
mv *_align.fasta ./all_loci_separated


