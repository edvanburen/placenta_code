# "Single-Cell RNA-Sequencing Reveals Placental Response under Environmental Stress" by Van Buren et al

## Contact Emails

Eric Van Buren: [evb\@hsph.harvard.edu], Hae-Ryung Park: [hae-ryung_park\@urmc.rochester.edu]

## Code Repository

This repository contains the code needed to reproduce analysis and results from the manuscript "Single-Cell RNA-Sequencing Reveals Placental Response under Environmental Stress" by Van Buren, Azzara, Rangel-Moreno, de la Luz Garcia-Hernandez, Murphy, Cohen, Lin, and Park. The corresponding data are available on Zenodo at https://zenodo.org/records/10258020. 

# read_data_to_submit.R

This file contains the R code needed to create the final data object used in all analyses and figures, including the preprocessing and data integration steps used. The final object is provided in the <code>final_Seurat_object.RData</code> object on Zenodo, so users not need to run this code themselves.

# manuscript_plots_to_submit.R

This file contains the code needed to produce all figures and results which rely on single-cell RNA-sequencing data used in the manuscript. Before running, users should download the <code>final_Seurat_object.RData</code> object from Zenodo and install required R packages.

