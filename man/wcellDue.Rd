\name{wcellDue}
\alias{wcellDue}
\alias{wcelldue}

\title{
	Window segmentation function.
}

\description{
	Given a step size, segment the data into windows according to the specified method based on the straightened data.
}

\usage{
	wcellDue(rds, Assay = 'Spatial', method = NULL, windowSV = 250, windowSH = 10)
}

\arguments{
  \item{rds}{       
     The Seurat object file processed using the R package SpLin.
  }
  \item{Assay}{       
    Assay used for calculating the expression levels of cells within windows during window segmentation.
  }
  \item{method}{       
    The method for window segmentation, which includes two options: "Vertical" (segmentation perpendicular to the x-axis) and "Horizontal" (segmentation perpendicular to the y-axis).
  }
  \item{windowSV}{       
    The window size (perpendicular to the xnew direction) should be within the range of 100 to 1000. For 10X Genomics data, a value of 250 is recommended; for Stereo-seq data, a value of 500 is recommended.
  }
  \item{windowSH}{       
    The window size (perpendicular to the ynew direction) should be within the range of 10 to 100.
  } 
}
\author{
	Wen Luo(design, implementation, testing)
}
\examples{
	rds = 'HE/SpLin_EXP_output.rds'
	windowSV = 250
	library(SpLin)
	wcellDue(rds, windowSV = windowSV)
}