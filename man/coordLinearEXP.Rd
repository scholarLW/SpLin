\name{coordLinearEXP}
\alias{coordLinearEXP}
\alias{coordlinearexp}
\alias{coordlinear}

\title{
	Linearize and straighten the spatial coordinates of the signal matrix.
}

\description{
	We construct an ordered bin set (typically quadrilaterals, but may also include polygons) using the ordered “Nei” and “Wai” signal points, and then perform spatial linear transformations and other processes to ultimately achieve the linearization and straightening of the spatial coordinates of the signal matrix.
}

\usage{
	coordLinearEXP(rds, json, nlim = 0, drop = TRUE, triangleMerge = TRUE, triangle.probs = 0.2)
}

\arguments{
  \item{rds}{       
     The Seurat object file processed using the R package SpLin.
  }
  \item{json}{       
     This is a Labelme JSON result file for images containing the labels "Nei" and "Wai", or additionally including the "Sign" label after adjustment.
  }
  \item{nlim}{       
     The maximum number of cells allowed within the range formed by three adjacent points in the “Wai” label, used for filtering the points in the “Wai” label.
  }
  \item{drop}{       
     Whether to output the image of the Bin.
  }
  \item{triangleMerge}{       
     For triangular-shaped Bins, whether to merge them with adjacent Bins.
  }
  \item{triangle.probs}{       
     The quantile probability of Bins determined to be triangular in shape.
  }
}
\author{
	Wen Luo(design, implementation, testing)
}
\examples{
	rds = 'SignalMatrix/Update_Manual_modified_data.rds'
	json = 'SignalMatrix/Adjusted_Output_Manual.json'
	nlim = 0
	drop = TRUE  
	triangleMerge = TRUE 
	triangle.probs = 0.2 
	library(SpLin)
	coordLinearEXP(rds, json, nlim = nlim, drop = drop, triangleMerge = triangleMerge, triangle.probs = triangle.probs)
}