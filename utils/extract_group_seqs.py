#!/usr/bin/env python
# -*- coding: utf-8 -*-

def parse_args():
    '''Parse arguments given to script'''

    import argparse

    parser = argparse.ArgumentParser(description="Given two genomes, get both \
        common and specific regions to both genomes following a kmer-based algorithm")
    parser.add_argument("--max-seeds", metavar="maximum seeds per kmer", 
                        dest="max_seeds", default=20)
    parser.add_argument("-m, --mapping", dest="mfile", metavar="mapping file")
    parser.add_argument("-f, --fasta", dest="fasta", metavar="fasta file")
    parser.add_argument("-o, --output", dest="outfile", metavar="output file")
    parser.add_argument("-g, --group", dest="group", metavar="group to extract")

    args = parser.parse_args()

    return args

def LoadMapping(mfile, group):
    '''This function reads a mapping file and returns a dictionary where
    provided groups are keys and a respective sequence with its coordinates
    are values inside another dict Each group has an associated list of tuples,
    where each tuple contains the identifier of the sequence and another 
    tuple of coordinates where that sequence has aligned with this group

    {
    Group1 : { seq1: (0,100), seq2 : (200,300) }
    }

    In this case it is returned a dict with sequences asociated to Group1.
    These sequences are seq1 from 0 to 100 and seq2 from 200 to 300

    groups are provided as a set
    
    '''

    mapping = {}
    with open(mfile) as infile:
        for line in infile:
            line = line.rstrip()
            g, seqs  = line.split("\t")
            if g == group:
                seqs = seqs.split(";")
                for seq in seqs:
                    gi, start, end = seq.split(":")
                    mapping[gi] = (int(start), int(end))

    return mapping

def FastaParser(handle):
    """THIS FUNCTION IS A COPY PASTE FROM BIOPYTHON PARSER
    It is included so this software does not need a biopython installation

    Generator function to iterate over Fasta records (as string tuples).
    For each record a tuple of two strings is returned, the FASTA title
    line (without the leading '>' character), and the sequence (with any
    whitespace removed). The title line is not divided up into an
    identifier (the first word) and comment or description.
    >>> with open("Fasta/dups.fasta") as handle:
    ...     for values in SimpleFastaParser(handle):
    ...         print(values)
    ...
    ('alpha', 'ACGTA')
    ('beta', 'CGTC')
    ('gamma', 'CCGCC')
    ('alpha (again - this is a duplicate entry to test the indexing code)', 'ACGTA')
    ('delta', 'CGCGC')
    """
    # Skip any text before the first record (e.g. blank lines, comments)
    while True:
        line = handle.readline()
        if line == "":
            return  # Premature end of file, or just empty?
        if line[0] == ">":
            break

    while True:
        if line[0] != ">":
            raise ValueError(
                "Records in Fasta files should start with '>' character")
        title = line[0:].rstrip()
        lines = []
        line = handle.readline()
        while True:
            if not line:
                break
            if line[0] == ">":
                break
            lines.append(line.rstrip())
            line = handle.readline()

        # Remove trailing whitespace, and any internal spaces
        # (and any embedded \r which are possible in mangled files
        # when not opened in universal read lines mode)
        yield title, "".join(lines).replace(" ", "").replace("\r", "")

        if not line:
            return  # StopIteration

    assert False, "Should not reach this line"

def GetGI(record_string):
    '''This function return the gi number from a fasta header with format
    gi|xxxxxxx|. If fasta header is not in this format it raises an error and
    exits program with an error message'''
    import sys
    
    if record_string[1:4] == "gi|":
        header = record_string.split("|")
        gi = header[1]
        return gi
    else:
        sys.exit("Error: there is one or more headers that don't fit the:  \
    >gi|xxxxxxxxx|.... format\
    Error in record: {}".format(record_string))

def ExtractSeqs(parser, mapping):
    '''GENERATOR: This function takes a fasta parser of sequences to be extracted
    according to groups they belong to and a mapping dictionary with
    groups to be extracted containing coordinates to be extracted for
    each sequence and yields extracted sequence by coordinates'''

    # For each sequence of the original fasta file
    for header, seq in parser:
        # Get its gi
        gi = GetGI(header)
        # Check if that gi is part of the group being extracted
        if gi in mapping:
            # In that case extract that sequence
            print mapping
            print mapping[gi]
            start_coord, end_coord = mapping[gi]
            extracted = seq[start_coord : end_coord]
            yield (header, extracted)

def WriteSeq(handle, header, seq, ruler=100):
    '''This function takes a fasta header and sequence and writes
    sequence in fasta format with ruler as a value for formating
    lines'''
    
    handle.write(header + "\n")
    rul = 0
    for i in xrange(len(seq)):
        if rul < ruler:
            handle.write(seq[i])
            rul += 1
        else:
            handle.write("\n" + seq[i])
            rul = 0
    handle.write("\n")

def main():

    args = parse_args()

    mapping = LoadMapping(args.mfile, args.group)
    with open(args.fasta) as handle:
        parser = FastaParser(handle)
        with open(args.outfile, "w") as handle:
            for header, seq in ExtractSeqs(parser, mapping):
                print header
                WriteSeq(handle, header, seq)

if __name__ == "__main__":
    main()