\name{matrixpls.boot}
\alias{matrixpls.boot}
\title{Bootstrapping of matrixpls function}
\usage{
  matrixpls.boot(data, R = 500, ...,
    parallel = c("no", "multicore", "snow"),
    ncpus = getOption("boot.ncpus", 1L),
    stopOnError = FALSE)
}
\arguments{
  \item{data}{Matrix or data frame containing the raw
  data.}

  \item{R}{Number of bootstrap samples to draw.}

  \item{...}{All other parameters are passed through to
  \code{matrixpls}.}

  \item{stopOnError}{A logical indicating whether
  boostrapping should be continued when error occurs in a
  replication.}

  \item{parallel}{ The type of parallel operation to be
  used (if any).  If missing, the default is taken from the
  option \code{"boot.parallel"} (and if that is not set,
  \code{"no"}).  }

  \item{ncpus}{ integer: number of processes to be used in
  parallel operation: typically one would chose this to the
  number of available CPUs.  }
}
\value{
  An object of class \code{\link[boot]{boot}}.
}
\description{
  \code{matrixpls.boot} is a convenience method that
  implements bootstrapping of \code{matrixpls} with
  \code{\link[boot]{boot}} method of the \code{boot}
  package.
}
\seealso{
  \code{\link[boot]{boot}}
}

