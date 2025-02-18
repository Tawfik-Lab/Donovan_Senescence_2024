# Data sources

Below are the datasets used in this analysis, including publication references and direct access links.

- **GSE154659**  
  - Paper DOI: [https://doi.org/10.1016/j.neuron.2020.07.026](https://doi.org/10.1016/j.neuron.2020.07.026)
  - Source: [GEO Accession GSE154659](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE154659)
  - Species: *Mus musculus*

- **GSE155622**  
  - Paper DOI: [https://doi.org/10.1038/s41422-021-00479-9](https://doi.org/10.1038/s41422-021-00479-9)
  - Source: [GEO Accession GSE155622](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE155622)
  - Species: *Mus musculus*

- **GSE249746**  
  - Paper DOI: [https://doi.org/10.1101/2023.03.17.533207](https://doi.org/10.1101/2023.03.17.533207)
  - Source: [GEO Accession GSE249746](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE249746)
  - Species: *Homo sapiens*

- **North 2019**  
  - Paper DOI: [https://doi.org/10.1093/brain/awz063](https://doi.org/10.1093/brain/awz063)
  - Sources:
      - [Supplementary Data](https://academic.oup.com/brain/article/142/5/1215/5382382#supplementary-data)
      - [UT Dallas Pain Neuro Lab](https://apps.utdallas.edu/bbs/painneurosciencelab/sensoryomics/hdrgclinical/)
  - Species: *Homo sapiens*
  - Note: If `wget` fails, please download this data manually using the provided links.


# Build instructions

This section provides instructions for building and cleaning the project.


## Build All Targets
```sh
cd build
make
```

## Build individual targets

```sh
cd build
make downloads
make datasets
make figures
```

## Clean up
```sh
make clean
```
