\name{markPlot}
\alias{markPlot}

\title{
	Linearized visualization of gene expression levels.
}

\description{
	Used for visualizing the distribution of gene expression levels in the straightened coordinate system.
}

\usage{
	markPlot(rds, markerGenes = c('Reg3b','Dcn','Hmgcs2'), SpatialFeaturePlotstriptextsize = 20, SpatialFeaturePlotpointsize = 0.1, SpatialFeaturePlotpointalpha = 0.5, pdf.width = 16, pdf.height = 4, png.width = 50, png.height = 8, png.dpi = 300, sv = 'V3')
}

\arguments{
  \item{rds}{       
     The Seurat object file processed using the R package SpLin.
  }
  \item{markerGenes}{       
     A list of genes of interest, with the input format being a vector.
  }
  \item{SpatialFeaturePlotstriptextsize}{       
     Color bar text size.
  }
  \item{SpatialFeaturePlotpointsize}{       
     The size of the scatter points.
  }
  \item{SpatialFeaturePlotpointalpha}{       
     The transparency of the scatter points, ranging from 0 to 1.
  }
  \item{pdf.width}{       
     The width of the PDF figure.
  }
  \item{pdf.height}{       
     The height of the PDF figure.
  } 
  \item{png.width}{       
     The width of the PNG figure.
  }
  \item{png.height}{       
     The height of the PNG figure.
  }  
  \item{png.dpi}{       
     The nominal resolution in ppi which will be recorded in the bitmap file, if a positive integer.  Also used for 'units' other than the default, and to convert points to pixels.
  }
  \item{sv}{       
	The Seurat object of class, V5 assay or V3 assay. Generally, versions of Seurat prior to V5 (excluding V5) are considered V3 assay. The parameter value can be either V3 or V5, and it defaults to V3 here.
  }      
}
\author{
	Wen Luo(design, implementation, testing)
}
\examples{
	rds = 'HE/SpLin_EXP_output.rds'
	library(SpLin)
	markPlot(rds, markerGenes = c('Igkc', 'Muc2', 'Lypd8', 'Igha', 'Actb', 'Zg16', 'Saa1', 'S100a6'), SpatialFeaturePlotpointsize = 1.5)
}