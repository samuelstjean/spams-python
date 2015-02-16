\name{spams.im2col_sliding}
\alias{spams.im2col_sliding}
\title{
spams.im2col_sliding
}
\description{
  A replacement for matlab function im2col for mode %% 'sliding'.
}
\usage{
B <- spams.im2col_sliding(A,m,n,RGB = FALSE);
}
\arguments{
  \item{A}{%
    double M x N matrix
    \preformatted{%
      contains the image pixels. For a color image,
      - if rgb pixels are consecutive in rows (N = image-ncols * 3),
        RGB must be set to TRUE.
 %%     - if colors are stored in 3 "planes" (M = 3 * image-nrows),
	RGB must be FALSE
      }
    }
    \item{m}{%
      number of rows of a sample
    }
    \item{n}{%
      number of columns of a sample
    }
    \item{RGB}{%
      TRUE for a colour image with RGB pixels consecutive.
    }
}
\value{
  \item{B}{%
%%    result matrix of size [m*n,(M -m + 1) * (N - n + 1)]
  }
}
\author{
Jean-Paul CHIEZE, 2011
}
