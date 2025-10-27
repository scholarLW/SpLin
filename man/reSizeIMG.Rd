\name{reSizeIMG}
\alias{reSizeIMG}
\alias{resizeimg}
\alias{reimg}

\title{
	Resize an image and save its dimensions.
}

\description{
	Resize an image and save its dimensions.
}

\usage{	
	reSizeIMG(ori, resi = NULL, dpi = 300, it = 'Signal')
}

\arguments{
  \item{ori}{       
     The original image file.
  }
  \item{resi}{       
     The resized image file, which defaults to NULL if it is set as 'Signal'.
  }
  \item{dpi}{       
     The nominal resolution in ppi which will be recorded in the bitmap file, if a positive integer.  Also used for 'units' other than the default, and to convert points to pixels.	
  }  
  \item{it}{       
     Type of the original image, can be 'ssDNA', 'HE', or 'Signal'.
  }
}
\author{
	Wen Luo(design, implementation, testing)
}
\examples{
	inputfile = 'SignalMatrix/signalMatrix.png'
	library(SpLin)
	reSizeIMG(inputfile)
}