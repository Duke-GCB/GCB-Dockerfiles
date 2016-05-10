#!/usr/bin/env python2.7

# Written by Thomas Konneker, adapted by Ian McDowell and Alejandro Barrera
# This script is intended to crawl through bed files to filter
# out reads that are too concentrated within a specified base-pair
# window.  This is presumably because they are PCR artifacts.
# See Boyle et. al 2011, Genome Research
# Currently this program does not 'window' in the traditional sense
# of tiling out the genome and then evaluating those
# windows. Instead starts a window at the next place there is an entry
# that follows the previous window.

# Modified 7-1-2014 to include requirement for 2 reads minimum for removal of putative artifact.
# Modified 4-25-2016 to reformat code including option to specify a BED file instead of reading stdin.

import sys
from optparse import OptionParser
from collections import Counter


class BedEntry:
    def __init__(self, chromosome, start, end, name,
                 score, strand, fullBedLine):
        self.chromosome = str(chromosome)
        self.startposition = int(start)
        self.endposition = int(end)
        self.name = name
        self.score = score
        self.strand = strand
        self.fullBedLine = fullBedLine


class BedReader:
    """BED file reader object. For parsing and yielding bed entries from a file."""

    def __init__(self, inFile):
        self.inFile = inFile

    def spitLines(self):
        for line in self.inFile:
            fullBedLine = line.strip('\n\r')
            wholeline = fullBedLine.split('\t')
            bedLine = BedEntry(wholeline[0], wholeline[1], wholeline[2],
                               wholeline[3], wholeline[4], wholeline[5], fullBedLine)
            yield bedLine


def windowChooser(position_list, cutoff):
    """Evaluate whether a list of bed is too concentrated at a single position, given a concentration threshold.

    Returns: boolean.
    """
    position_list, cutoff = position_list, cutoff
    if position_list:
        count = Counter(position_list)
        peak = float(count.most_common()[0][1])
        if peak / len(position_list) >= cutoff and peak >= 2:
            return 0
        else:
            return 1
    else:
        return 1


def windowizer(bedIn, windowsize, cutoff):
    """Compute bedlines that fall under the cutoff criteria for concentration of reads to a single base pair within a window."""

    windowend = 0
    current_chr = "chr1"
    full_bed_list = []
    position_list = []
    for line in bedIn.spitLines():
        if line.startposition > windowend:
            windowstart = line.startposition
            windowend = windowstart + windowsize
            if windowChooser(position_list, cutoff):
                for item in full_bed_list:
                    print item
            position_list = [line.startposition]
            full_bed_list = [line.fullBedLine]
        elif line.startposition <= windowend:
            if line.chromosome == current_chr:
                position_list.append(line.startposition)
                full_bed_list.append(line.fullBedLine)
            else:
                current_chr = line.chromosome
                windowstart = line.startposition
                windowend = windowstart + windowsize
                if windowChooser(position_list, cutoff):
                    for item in full_bed_list:
                        print item
                position_list = [line.startposition]
                full_bed_list = [line.fullBedLine]
        else:
            if windowChooser(position_list, cutoff):
                for item in full_bed_list:
                    print item


def main():
    parser = OptionParser()
    parser.add_option('-w', '--windowsize', type='int',
                      dest='windowsize', default=31,
                      help='threshold for trimming reads that have \
              the same mapping start position')
    parser.add_option('-c', '--cutoff', type="float",
                      dest='cutoff', default=0.70,
                      help="threshold for concentration at a single base \
              within a window to cutoff")
    parser.add_option('-i', '--infile', type="str", dest='infile',
                      help="Input BED file to be scanned and filtered")

    (options, args) = parser.parse_args()

    windowsize = options.windowsize
    cutoff = options.cutoff
    infile_name = options.infile

    # parse command line
    if infile_name:
        infile = open(infile_name)
    else:
        infile = sys.stdin
    entry = BedReader(infile)
    windowizer(entry, windowsize, cutoff)

    if infile_name is not sys.stdin:
        infile.close()


if __name__ == '__main__':
    main()