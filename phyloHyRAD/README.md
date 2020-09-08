# phyloHyRAD is a pipeline to perform phylogenetic analyses from HyRAD data. 

3 steps:

* probes reconstruction using [IpyRAD](https://ipyrad.readthedocs.io/) 
* sample loci reconstruction, combination and alignment
* phylogeny


## Dependances

* Python and BioPython
* [samtools](http://www.htslib.org/)
* [IpyRAD](https://ipyrad.readthedocs.io)
* [bwa](http://bio-bwa.sourceforge.net/)
* [picard](https://broadinstitute.github.io/picard/)
* gatk3 
* [mapDamage2](https://ginolhac.github.io/mapDamage)
* mafft 
* [amas](https://github.com/marekborowiec/AMAS)


## Usage
Step1

IpyRAD loci assembly on probes

Step2


```
./phyloHyRAD.sh -r ipyrad_output.loci -m 

```
