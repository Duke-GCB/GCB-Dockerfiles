FROM continuumio/miniconda
MAINTAINER Alejandro Barrera <alejandro.barrera@duke.edu>

RUN conda install -c bioconda deeptools

# Default command to execute at startup of the container
ENV DEEPTOOLS_SUBCOMMANDS="multiBamSummary multiBigwigSummary correctGCBias bamCoverage bamCompare bigwigCompare computeMatrix plotCorrelation plotPCA plotFingerprint bamPEFragmentSize computeGCBias plotCoverage plotHeatmap plotProfile"
CMD echo "deepTools subcommands:\n\t$DEEPTOOLS_SUBCOMMANDS"
