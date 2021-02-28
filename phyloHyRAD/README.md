# phyloHyRAD is a pipeline to perform phylogenetic analyses from HyRAD data. 


## Dependances

* Python and BioPython
* [samtools](http://www.htslib.org/)
* [IpyRAD](https://ipyrad.readthedocs.io)
* [bwa](http://bio-bwa.sourceforge.net/)
* [picard](https://broadinstitute.github.io/picard/)
* gatk3 
* [mapDamage2.0](https://ginolhac.github.io/mapDamage)
* [mafft](https://mafft.cbrc.jp/alignment/software/)
* [amas](https://github.com/marekborowiec/AMAS)


## Usage


Detailled description:

Step 1 : Reference loci assembly

This step is performed on probes using [IpyRAD](https://ipyrad.readthedocs.io). Probes are produced using a classical ddRAD approach. [IpyRAD](https://ipyrad.readthedocs.io) is run on the probes in a classical way (see [IpyRAD](https://ipyrad.readthedocs.io) manual). As recommended, several "% of identity" values can be tested to identify the best clustering parameters according to the divergence between the species used as probes.

Step 2 : phyloHyRAD pipeline

Command:

```
./phyloHyRAD.sh -r ipyrad_output.loci -m 

```

All the following steps are included in the phyloHyRAD pipeline: 

1.Reference loci catalog formatting. 

Once the IpyRAD analyses on the probes are successful, one of the output, finishing by ".loci", contains the sequences used as reference.
This file needs to be indicated in the phyloHyRAD command after the option -r.
The pipeline will carry out the necessary modifications.

2.Mapping

The second step is the mapping of the reads (R1 and R2) from each historical sample on the reference catalog.
R1 and R2 fastq files need to be located in the current directory with a classical name: *sample\_R1.fastq* and *sample\_R2.fastq*
These reads need to be previously cleaned to remove adaptors and bases with low quality and synchronized. 
Two mapping method are implemented in the pipeline: BWA-MEM algorithm (option -m mem) and BWA-backtrack aln/sampe (option -m aln).

3.Mapping cleaning

Various classical mapping filtering are implemented including: the indel realignment and the removal of PCR duplicates using [picard](https://broadinstitute.github.io/picard/) suite and the deamination correction using [mapDamage2.0](https://ginolhac.github.io/mapDamage).

4.Consensus

For each sample and each loci a consensus is generated from the mapping file.
This consensus is performed using [samtools mpileup / bcftools / vcfutils.pl](https://samtools.github.io/bcftools/howtos/consensus-sequence.html) and, by default, consider samples as diploid, each heterozygote position is called using the IUPAC code.

5.Loci combination and alignment

The previous consensuses are combined and aligned using [mafft](https://mafft.cbrc.jp/alignment/software/).
At the end all alignments are concatenated to a global alignment. 
Before this final step, alignments are stored in a folder named *all\_loci\_separated* and
can be filtered and manually checked at the user's convenience.

Two final files are generated:

* *concatenated.out*: the global alignment in fasta format
* *partitions.txt*: the partition file

These files can directly be used in all phylogeny tools.





