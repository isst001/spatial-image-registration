# Spatial Image Registration Workflow for Seurat / Xenium Data

This repository provides an R-based workflow for semi-automated image registration of spatial transcriptomics data onto a reference contour image.

The pipeline is designed for **Seurat-based spatial objects** (e.g., Xenium / MERSCOPE style data) and includes manual region selection, manual anchor selection, coordinate rotation, thin-plate spline (TPS) warping, export of registered coordinates, and visualization of mapped cell types.

---

# Features

* standardizes cell IDs for downstream processing
* extracts molecule coordinates for selected genes
* interactively selects a region of interest
* rotates spatial coordinates before registration
* collects homologous anchors manually from the spatial plot
* reads contour anchors from a CSV file
* fits thin-plate spline models for image registration
* warps molecule and centroid coordinates into contour space
* exports registered coordinates and QC metrics
* visualizes selected cell types on the reference contour

---

# Repository Structure

```
spatial-image-registration/
├── R/
│   ├── image_registration_public_release.R
│   ├── run_example.R
│   └── utils_io.R
├── inst/
│   └── extdata/
│       ├── contour_img_example.png
│       └── contour_anchors_example.csv
├── data/
│   ├── input/
│   ├── anchors/
│   └── outputs/
├── docs/
│   ├── workflow_overview.md
│   └── zenodo_release_checklist.md
├── README.md
├── LICENSE
├── .gitignore
└── CITATION.cff
```

---

# Requirements

## R version

Tested with:

```
R ≥ 4.3
```

## R packages

Required packages:

* Seurat
* dplyr
* ggplot2
* imager
* Morpho
* sp
* purrr
* stringr
* fields
* tibble

Install missing packages:

```r
install.packages(c(
"dplyr",
"ggplot2",
"imager",
"Morpho",
"sp",
"purrr",
"stringr",
"fields",
"tibble"
))
```

Install **Seurat** separately following the official instructions.

---

# Input Requirements

## 1️⃣ Seurat Object

The Seurat object should contain:

* spatial image entry

  ```
  object@images[[slide_id]]
  ```

* molecule coordinates

  ```
  object@images[[slide_id]]$molecules
  ```

* centroid coordinates

  ```
  object@images[[slide_id]]$centroids
  ```

* metadata column defining cell type or cluster identity

Default column used in the script:

```
RNA_snn_res.14
```

This can be changed in the configuration:

```
cell_type_column
```

---

## 2️⃣ Contour Image

A PNG image used as the reference anatomical contour.

Example location:

```
inst/extdata/contour_img_example.png
```

---

## 3️⃣ Contour Anchors

CSV file containing anchor points on the contour image.

Format:

```
x,y
1710,887
1884,257
...
```

These anchors must correspond to the **same anatomical landmarks** that will be clicked in the spatial dataset.

---

# Usage

1️⃣ Place Seurat object:

```
data/input/seurat_object.rds
```

2️⃣ Place contour files:

```
inst/extdata/contour_img_example.png
inst/extdata/contour_anchors_example.csv
```

3️⃣ Edit configuration in:

```
R/image_registration_public_release.R
```

4️⃣ Run:

```r
source("R/image_registration_public_release.R")

result <- run_image_registration(config)
```

---

# Example Run Script

```
R/run_example.R
```

Run with:

```r
source("R/run_example.R")
```

---

# Output Files

Generated in:

```
data/outputs/
```

Outputs include:

* rotated spatial preview image
* warped overlay image
* anchor QC plot
* warped molecule coordinate table
* warped centroid coordinate table
* centroid metadata table
* anchor residual error table
* cell type spatial plots

Spatial anchors clicked interactively are saved to:

```
data/anchors/
```

---

# Reproducibility

Recommended publication workflow:

1. finalize analysis code
2. upload repository to GitHub
3. create a tagged release (e.g. `v1.0.0`)
4. archive release with **Zenodo**
5. cite the DOI in your manuscript

---

# Citation

If you use this code, please cite the software release described in **CITATION.cff**.

---

# License

Released under the **MIT License**.

See the LICENSE file for details.
