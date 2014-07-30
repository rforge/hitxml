\name{HITXML-package}
\alias{HITXML-package}
\alias{HITXML}
\docType{package}

\title{
Package to create XML plans for HIT.
}

\description{
This package contains code that enable to create simple XML-based irradiation plans. 
}

\details{
\tabular{ll}{
Package: \tab HITXML\cr
Type: \tab Package\cr
Version: \tab 0.6\cr
Date: \tab 2012-04-17\cr
License: \tab GPL3\cr
}

The planning process is executed as a command line script distributed with this package and should be started using \code{Rscript} rather then using the individual R functions of the package in a custom environment (as this is highly likely to fail due to the immatureness of the code...). The use is explained in the 'example' section below. The package contains proprietary data from HIT (libC) AND MUST NOT BE DISTRIBUTED BEYOND DKFZ E0409 or HIT despite the GPL license! The authors refuse any responsibility in this case.

In case you encounter an error message like 'bad interpreter /usr/bin/Rscript^M' the Unix line-ends (default) have been accidentally changed to Windows style. You can fix this by 'dos2unix' for all '.R'-files. 

If run successful, output will be a request and a plan file with a path ready to be copied onto the control machine as well as a pdf showing the field shape etc.

For the latter, we recommend a pdf viewer capable of reloading updating files (i.e. evince, not Acrobat reader).

Suggestion for prompt under Windows $P$_$G$S

The fluence to dose conversion uses tabulated stopping power data for liquid water from ICRU49/73 reports.

HIT_XML uses the txtplot function from Claus E. Andersen, DTU Risoe.

!!IMPORTANT!!: If your plan should not work at HIT, esp. when successfully applied before or with only one ion species double check with accelerator team if MD% checksums for the control system have been updated after any kind of system maintainance. This is the most likely cause for HITXML plans to fail.
}

\author{
Steffen Greilich, Sarah Stefanowicz

Maintainer: Steffen Greilich <s.greilich@dkfz.de>
}

\references{
Nix
}

\keyword{ HIT plan generation }
\seealso{
XML
}
\examples{
\dontrun{
## To start the planning tool, use Rscript from the command line
# As the location of the script can vary with your installation
# you can use the HITXML function which return the location
# and invokes the script directly
Rscript -e "HITXML::HITXML()"
# Unfortunately this does not work under Windows where it will only
# return the file and its location which you then have to manually 
# insert in the command, e.g.
Rscript C:\Users\greilich\Documents\R\win-library\2.14\HITXML\exec\HIT_XML.R
}
}
