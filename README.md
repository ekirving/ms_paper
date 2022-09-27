# Polygenic selection analysis for Multiple Sclerosis (MS)
This repository contains code for the ancestry stratified polygenic selection analyses from 
[Genetic Risk for Multiple Sclerosis Originated in Pastoralist Steppe Populations](https://).

![Figure 2](./figure/Figure_5.png?raw=true)

If you reuse any of this code then please cite the preprint:
> Barrie, W.&ast;, Yang, Y.&ast;, Attfield, K.E.&ast;, Irving-Pease, E.&ast;, Scorrano, G.&ast;, Jensen, L.T.&ast;, 
> Armen, A.P., Dimopoulos, E.A., Stern, A., Refoyo-Martinez, A., Ramsøe, A., Gaunitz, C., Demeter, F., 
> Jørkov, M.L.S., Møller, S.B., Springborg, B., Klassen, L., Hyldgård, I.M., Wickmann, N., Vinner, L., 
> Korneliussen, T.S., Allentoft, M.E., Sikora, M., Kristiansen, K., Rodriguez, S., Nielsen, R., Iversen, A.K.N., 
> Lawson, D.J., Fugger, L., Willerslev, E., 2022. Genetic risk for Multiple Sclerosis originated in Pastoralist Steppe 
> populations. *bioRxiv* 2022.09.23.509097. https://doi.org/10.1101/2022.09.23.509097

## Installation
Download the code: 
```bash
git clone git@github.com:ekirving/ms_paper.git && cd ms_paper/
```

The easiest way to install all the dependencies is with the [conda package manager](https://docs.conda.io/en/latest/).

```bash
conda env create --name ms --file environment.yaml
```

Then activate the environment:
```bash
conda activate ms
```

## Running the code

To reproduce all the analyses, simply run:

```bash
snakemake all
```

## Author

Evan K. Irving-Pease, [GLOBE Institute](https://globe.ku.dk/), University of Copenhagen 

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details
