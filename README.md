# Dove
Variant analysis tool kit written in python3 for downstream vcf analysis.

## Online Variant Anotation
Annotates given vcf file using Ensembl VEP API.

## Multi Genome Analysis
Filters index sample's variants with parents and siblings based on given inheritance pattern.

## Variant Filter
Filters annotation file by given parameters.

## Quick Install
* Linux&Mac  
> sudo pip3 install --index-url https://test.pypi.org/simple/ dove
* Windows
> pip install --index-url https://test.pypi.org/simple/ dove

## Example Uses

* OVA examples  
>  dove OVA -i input.vcf -o output.csv

* MGA examples  
>  dove MGA -i index.csv -m mother.csv -f father.csv -a trio -p AR -o output.csv
>  dove MGA -d pedigree.csv -a siblings -o output.csv

* VF generic filter examples
>  dove VF -i index.csv -z 1/1 -f 0.01 -r 20 -b protein_coding -o output.csv

* VF Specific filter examples
>  dove VF -t filterfile.txt

* dove examplefile
>  dove examplefile --help
>  dove examplefile -e both

* dove updateomim
>  dove updateomim --help
>  dove updateomim -k apikey
