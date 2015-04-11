\name{HITXML-package}
\alias{HITXML-package}
\alias{HITXML}
\docType{package}

\title{
Package to create XML plans for HIT.
}

\description{Code to create simple XML-based irradiation plans and forward-calculate plans or beam records. The functionality is consider to be used withs sscript distributed with this package (in /exec). The package is tailored to use with the scanning ion beam at the Heidelberg Ion-beam Therapy center (HIT) but might be useful for generic application as well. }

\details{
Command line scripts are distributed with this package and can be started using \code{Rscript} or R as an interpreter (Linux, Max). The use is explained in the 'example' section below. The package needs data on accelerator setting that are proprietary data to HIT (libC) and are therefore NOT distributed with the package. Users of DKFZ / HIT contact s.greilich@dkfz.de for more information.

In case you encounter an error message like 'bad interpreter /usr/bin/Rscript^M' the Unix line-ends (default) have been accidentally changed to Windows style. You can fix this by 'dos2unix' for all '.R'-files. 

If run successful, output will be a request and a plan file with a path ready to be copied onto the control machine as well as a pdf showing the field shape etc.

For the latter, we recommend a pdf viewer capable of reloading updating files (i.e. evince, Sumatra, ... *not* Acrobat reader).

Suggestion for prompt under Windows $P$_$G$S

The fluence to dose conversion uses tabulated stopping power data for liquid water from ICRU49/73 reports from the R-package libamtrack.

HIT_XML uses the txtplot function from Claus E. Andersen, DTU Risoe.

!!IMPORTANT!!: If your plan should not work at HIT, esp. when successfully applied before or with only one ion species double check with accelerator team if MD checksums for the control system have been updated after any kind of system maintainance. This is the most likely cause for HITXML plans to fail.}

\author{
Steffen Greilich, Sarah Stefanowicz and others

Maintainer: Steffen Greilich <s.greilich@dkfz.de>
}

\keyword{ Plan generation }
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

# There are two scripts available at the moment:
HIT_XML.R       - generates xml-based planes
HIT_XML_PBR.R   - forward-caculate fields from plans or beam records (PBR files)
}
}
