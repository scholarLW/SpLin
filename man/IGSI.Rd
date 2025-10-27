\name{IGSI}
\alias{IGSI}
\alias{igsi}

\title{
	Hotspot region and hotspot center analysis function.
}

\description{
	The hotspot region and hotspot center analysis function is primarily used to identify areas or central zones that are enriched with certain functional patterns due to biological processes or external influences. It aims to provide data support for interpreting biological processes.
}

\usage{
	IGSI(rds, GMF, topFreq = 0.05, nbin = 10, alpha = 0.001, SD = 3, topSegFreq = 0.4, TopCellFreq = 0, eps = 200, grid_density = 50, kneighbor = 10, corr = 0.5, piontsize = 0.5, sv = 'V3')
}

\arguments{
  \item{rds}{       
	The Seurat object file processed using moduleGene from the R package SpLin.
  }
  \item{GMF}{       
	The functional gene module result file processed using moduleGene from the R package SpLin.
  }  
  \item{topFreq}{       
	Select cells with scores in the top XX% for downstream analysis, with the range being between 0 and 1.
  }
  \item{nbin}{       
	Number of bins of aggregate expression levels for all analyzed features.
  }
  \item{alpha}{       
	significance levels for the test to accept change-points.
  }
  \item{SD}{       
	The number of SDs between means to keep a split if undo.splits="sdundo".
  } 
  \item{topSegFreq}{       
	Select segments with scores ranking in the top XX% for downstream analysis, with the range being between 0 and 1.
  } 
  \item{TopCellFreq}{       
	The proportion threshold for TopCell falling into a region.
  } 
  \item{eps}{       
	Size (radius) of the epsilon neighborhood. Can be omitted if 'x' is a frNN object.
  }
  \item{grid_density}{       
	Grid density for krige function, which is a positive integer.
  }   
  \item{kneighbor}{       
	The number of neighboring points of a point in Module Co-Activity Spots analysis. The numerical range is between 1 and 20.
  }  
  \item{corr}{       
	The similarity threshold between a point and its neighboring points in Module Co-Activity Spots analysis. The numerical range is between 0.5 and 1.
  }      
  \item{piontsize}{       
	The size of the scatter points.
  }  
  \item{sv}{       
	The Seurat object of class, V5 assay or V3 assay. Generally, versions of Seurat prior to V5 (excluding V5) are considered V3 assay. The parameter value can be either V3 or V5, and it defaults to V3 here.
  }      
}
\author{
	Wen Luo(design, implementation, testing)
}
\examples{
	rds = 'HE/IGSI/winCell/SpLin_EXP_output_withWindow.RDS'
	GMF = 'HE/IGSI/winCell/Modules/Gene/moduleLabel.txt'
	library(SpLin)
	IGSI(rds, GMF, topSegFreq = 0.4, piontsize = 2.5)
}