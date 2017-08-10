bwa-samtools
============

This Dockerfile contains bwa with samtools.

The tool is bwa, but since it produces SAM output, we include samtools and a wrapper script in the Docker image to produce BAM output.
