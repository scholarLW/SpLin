\name{getPolygonPionts}
\alias{getPolygonPionts}
\alias{getpolygonpionts}

\title{
	Get the signal points on the edges of the polygon.
}

\description{
	Map the pixel coordinates back to the high-resolution stained image and obtain the color pixel coordinates of the stained film corresponding to the regional location.
}

\usage{
	getPolygonPionts(rds, imageFile, json, idjson, Signal = FALSE)
}

\arguments{
  \item{rds}{       
     A Seurat object file containing a spatial signal matrix with spatial coordinates, derived from 10X Genomics or Stereo-seq sequencing.
  }
  \item{imageFile}{       
     The resized image, or the image generated from the signal matrix.
  }
  \item{json}{       
     This is a Labelme JSON result file for images containing the labels "Nei" and "Wai", or additionally including the "Sign" label.
  }
  \item{idjson}{       
     The dimension information file for the image, typically named "Dimensions.txt".
  }
  \item{Signal}{       
     Whether it is an image formed by the signal matrix.
  }
}
\author{
	Wen Luo(design, implementation, testing)
}
\examples{
	rds = 'modified_data.rds'
	image = 'SignalMatrix/signalMatrix.png'
	json = 'SignalMatrix/signalMatrix.json' 
	idjson = 'SignalMatrix/image_dimensions.json'
	library(SpLin)
	getPolygonPionts(rds, image, json, idjson, Signal = TRUE)
}