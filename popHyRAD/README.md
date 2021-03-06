# popHyRAD: small variant genotyping for HyRAD sequencing


popHyRAD is a pipeline to perform population genetic analyses from HyRAD data. 

3 steps:

* probes reconstruction using [IpyRAD](https://ipyrad.readthedocs.io/) 
* sample reads mapping and cleaning
* SNP calling using [freebayes](https://github.com/ekg/freebayes)

**Reference:**   
Gauthier J, Pajlokiv M, Neuenschwander S, Kaila L, Schmid S, Orlando L, Alvarez N (2020). **Museomics identifies genetic erosion in two butterfly species across the 20th century in Finland.**  _Molecular Ecology Resources_.
(https://doi.org/10.1111/1755-0998.13167)


## Dependances

* Python and BioPython
* [samtools](http://www.htslib.org/)
* [IpyRAD](https://ipyrad.readthedocs.io)
* [bwa](http://bio-bwa.sourceforge.net/)
* [picard](https://broadinstitute.github.io/picard/)
* gatk3
* [mapDamage2](https://ginolhac.github.io/mapDamage)
* [freebayes](https://github.com/ekg/freebayes)


## Usage
Step1

IpyRAD loci assembly on probes

Step2


```
./popHyRAD.sh -r ipyrad_output.loci -m 

```
