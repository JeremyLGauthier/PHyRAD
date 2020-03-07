# PHyRAD: loci reconstruction from HyRAD sequencing for phylogeny


PHyRAD is a pipeline to perform phylogenetic analyses from HyRAD data. 

3 steps:

* probes reconstruction using [IpyRAD](https://ipyrad.readthedocs.io/) 
* sample loci reconstruction, combinaition and alignment
* phylogeny

**Reference:**   
Gauthier J, Pajlokiv M,  Neuenschwander S, Kaila L, Schmid S, rlando L, Alvarez N (2020). **Museomics identifies genetic erosion in two butterfly species across the 20th century in Finland.** under review in _Molecular Ecology Resources_.

## Dependances

* Python and BioPython
* [samtools](http://www.htslib.org/)
* [IpyRAD](https://ipyrad.readthedocs.io)
* [bwa](http://bio-bwa.sourceforge.net/)
* [picard](https://broadinstitute.github.io/picard/)
* gatk3 
* [freebayes](https://github.com/ekg/freebayes)
* mafft 
* [amas](https://github.com/marekborowiec/AMAS)


## Usage
Step1

IpyRAD loci assembly on probes

Step2


```
./PhyRAD.sh -r ipyrad_output.loci 

```
