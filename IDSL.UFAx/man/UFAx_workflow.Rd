\name{UFAx_workflow}
\alias{UFAx_workflow}
\title{
UFAx Workflow
}
\description{
This function runs the exhaustive chemical enumeration part of the IDSL.UFAx pipeline.
}
\usage{
UFAx_workflow(spreadsheet)
}
\arguments{
  \item{spreadsheet}{
IDSL.UFAx parameter spreadsheet
}
}
\value{
The MolecularFormulaAnnotationTable is saved in the assigned folder in the parameter spreadsheet.
}
\note{
You should load the IDSL.UFA package to run the IDSL.UFAx functions.
}
\seealso{
See this link to illusterate an example \url{https://github.com/idslme/IDSL.UFAx/wiki/Example-for-IDSL.UFAx}
}