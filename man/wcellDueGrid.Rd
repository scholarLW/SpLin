\name{wcellDueGrid}
\alias{wcellDueGrid}
\alias{wcellduegrid}

\title{
	Window segmentation function.
}

\description{
	In accordance with the predefined grid density, the signal matrix derived from the spatial transcriptomic chip is partitioned into grids and subsequently divided into multiple smaller windows.
}

\usage{
	wcellDueGrid(rds, Assay = 'Spatial', grid_density = 30)
}

\arguments{
  \item{rds}{       
     The Seurat object file processed using the R package SpLin.
  }
  \item{Assay}{       
    Assay used for calculating the expression levels of cells within windows during window segmentation.
  }
  \item{grid_density}{       
    Grid density for splitting the signal matrix must be a positive integer between 10 and 1000.
  } 
}
\author{
	Wen Luo(design, implementation, testing)
}
\examples{
	rds = 'HE/SpLin_EXP_output.rds'
	grid_density = 30
	library(SpLin)
	wcellDueGrid(rds, grid_density = grid_density)
}