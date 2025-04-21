# RLadies-Morelia: Reprohack 2025

-   Artículo: [Landscape of mobile genetic elements and their antibiotic resistance cargo in prokaryotic genomes](https://academic.oup.com/nar/article/50/6/3155/6552054)
-   Github del articulo: https://git.embl.de/khedkar/promge
-   Fecha: lunes 28 de abril, 2025
-   Integrantes: Evelia Coss, Marisol Navarro, Diana Barceló y Johana Castelán.

## Material

Libro de [Quarto](https://iscb-rsg-mexico.github.io/Reprohack2025_RLadiesMorelia_RSGMexico/docs/index.html) con informacion. Dataset del articulo en [Git](https://git.embl.de/khedkar/promge)

## Figuras propuestas

En el Github se encuentra el codigo accesible de las siguientes figuras:

| Figura              | Script original  | Script modificado (en español)  |  Input  |
| :---------------- | :------: | ----: |  ----: |
| Python Hat        |   True   | 23.99 |
| SQL Hat           |   True   | 23.99 |
| Codecademy Tee    |  False   | 19.99 |
| Codecademy Hoodie |  False   | 42.99 |



+-------------+--------------------------------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------------------------------+-----------------------------------------------------------------------------------------------+
| Figura      | Script original                                                                            | Script modificado (en español)                                                                                                             | Input                                                                                         |
+=============+============================================================================================+============================================================================================================================================+===============================================================================================+
| Figura 2    | [`R/Figure2.R`](https://git.embl.de/khedkar/promge/-/blob/main/R/Figure2.R?ref_type=heads) | [`Figura2_modificado.R`](https://github.com/ISCB-RSG-MEXICO/Reprohack2025_RLadiesMorelia_RSGMexico/blob/main/scripts/Figura2_modificado.R) | ``` r                                                                                         |
|             |                                                                                            |                                                                                                                                            | mge_pg <- read_tsv("data/processed_data/mge_bins_per_genome_final.txt.gz", col_names = T)     |
|             |                                                                                            |                                                                                                                                            | mge_solitary <- read_tsv("data/processed_data/solitary_mge_bins_final.txt.gz", col_names = F) |
|             |                                                                                            |                                                                                                                                            | rec_class <- read_tsv("data/raw_data/recombinase.list.gz", col_names = F)                     |
|             |                                                                                            |                                                                                                                                            | glist <- read_tsv("data/raw_data/genome_status_supplementary_tableS2.txt.gz", col_names = T)  |
|             |                                                                                            |                                                                                                                                            | ```                                                                                           |
+-------------+--------------------------------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------------------------------+-----------------------------------------------------------------------------------------------+
| Figura 3    | [`R/Figure3.R`](https://git.embl.de/khedkar/promge/-/blob/main/R/Figure3.R?ref_type=heads) | [`Figura3_modificado.R`](https://github.com/ISCB-RSG-MEXICO/Reprohack2025_RLadiesMorelia_RSGMexico/blob/main/scripts/Figura3_modificado.R) | ``` r                                                                                         |
|             |                                                                                            |                                                                                                                                            | tax <- read_tsv("data/raw_data/species_with_atleast_2genomes.list.gz", col_names=F)           |
|             |                                                                                            |                                                                                                                                            | db <-read_tsv("data/processed_data/mge_bins_per_genome_final.txt.gz", col_names = T)          |
|             |                                                                                            |                                                                                                                                            | gs <-read_tsv("data/raw_data/genome_size.txt.gz", col_names = T)                              |
|             |                                                                                            |                                                                                                                                            | class_tree <- read.tree("data/raw_data/progenomes2_class_tree.nwk")                           |
|             |                                                                                            |                                                                                                                                            | glist <- read_tsv("data/raw_data/genome_status_supplementary_tableS2.txt.gz", col_names = T)  |
|             |                                                                                            |                                                                                                                                            | ```                                                                                           |
+-------------+--------------------------------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------------------------------+-----------------------------------------------------------------------------------------------+
| Figura 4    | [`R/Figure4.R`](https://git.embl.de/khedkar/promge/-/blob/main/R/Figure4.R?ref_type=heads) |                                                                                                                                            |                                                                                               |
+-------------+--------------------------------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------------------------------+-----------------------------------------------------------------------------------------------+
| Figura 5    | [`R/Figure5.R`](https://git.embl.de/khedkar/promge/-/blob/main/R/Figure5.R?ref_type=heads) |                                                                                                                                            |                                                                                               |
+-------------+--------------------------------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------------------------------+-----------------------------------------------------------------------------------------------+

------------------------------------------------------------------------

Curso hecho con amor 💜.
