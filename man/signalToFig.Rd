\name{signalToFig}
\alias{signalToFig}
\alias{signaltofig}
\alias{signal}
\alias{sTF}

\title{
	signal to figure.
}
\description{
	We present a function to convert signal matrices into visual figures for analysis.
}
\usage{
	signalToFig(rds, colN = NULL, pointcellsize = NULL, dpi = 300) 	
}
\arguments{
  \item{rds}{       
     A Seurat object file containing a spatial signal matrix with spatial coordinates, derived from 10X Genomics or Stereo-seq sequencing.
  }
  \item{colN}{       
     A class name used to determine group colors, which must come from a column in the metadata of the Seurat object file, and defaults to NULL.
  }
   \item{pointcellsize}{       
     The size of spatial scatter points.
  }
  \item{dpi}{       
     The nominal resolution in ppi which will be recorded in the bitmap file, if a positive integer.  Also used for 'units' other than the default, and to convert points to pixels.	
  }
}
\author{
	Wen Luo(design, implementation, testing)
}
\examples{
	rds = 'modified_data.rds'
	colN = NULL
	library(SpLin)
	signalToFig(rds, colN = colN)
}
