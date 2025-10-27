\name{rdsMerge}
\alias{rdsMerge}
\alias{rdsmerge}
\alias{mergerds}

\title{
    Merge two Seurat object files.
}

\description{
    The script can perform lasso operations on Seurat object files using StereoMap and merge them with Seurat object files that have already undergone cell clustering and annotation before the lasso operation, resulting in a Seurat object file with clustering annotations after lasso.
    Note: Since the integration analysis process may change cell names, to ensure a successful merge, please first check if there are any differences in cell names between the two files.
}

\usage{
    rdsMerge(lassoRDS, celltypeRDS)
}

\arguments{
    \item{lassoRDS}{
        After StereoMap, obtain the expression matrix of the lasso-selected region using the SAW software to form a Seurat object file with binXX or cellbin (with the suffix .RDS, such as DSS12.lasso.bin50.RDS).
    }
    \item{celltypeRDS}{
        If cell clustering and annotation have been performed on the sample before lasso selection, you can provide the annotated Seurat object file (with the suffix .RDS, such as DSS12.celltype.rds).
    }
}

\author{
    Wen Luo (design, implementation, testing)
}

\examples{
    cell = './test/DSS12.celltype.rds'
    lasso = './test/DSS12.lasso.bin50.RDS'
    library(SpLin)
    rdsMerge(lasso, cell)
}