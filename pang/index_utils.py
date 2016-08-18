#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Module for indexing nucleotide sequence for pangmer

def BuildIndex(k):
    '''This function initializes an index with all possible k-mers of
    length k as keys and empty lists as values. It is initialized with key
    <<start_offset : 0 >>. This value points which coordinate start counting
    from. First sequence indexed allways will start from 0, but next sequences
    added to index should start where last sequence ended.

    For example. Let's say there are sequecenes A and B with lengths 1000  and
    500. Sequence A will be indexed from 0 to 1000, and B from 1000 to 1700 so
    each sequence always have its own coordinates

    '''

    from itertools import product
    from array import array
 
    print "Building index of {}-mers...".format(k)
    index = {
              "start_offset" : 0
            }
    for kmer in product("ATCG", repeat = k):
        kmer = "".join(kmer)
        # Initialize with an empty value since not all possible kmers are going
        # to be present in indexed sequences due to DNA nature
        index[kmer] = None

    print "Done!"
    return index

def IndexSequence(sequence, k, index):
    '''This function takes an index of k-mer keys mapping to start coordinates
    in a sequence and updates it with new coordinates for sequence provided
    
    # start_offset is used when a new sequence is going to be added to the
    # index, so each sequence has its unique coordinates beginning always 
    # from the last coordinate of the previous sequence

    # At he same time reverse complement is indexed too for each k-mer, with
    # reverse coordinates too, so when a k-mer is taken in scanned sequence
    # it can seed and extend an alignment in the forward and reverse strand
    # at the same time

    '''
    from array import array
    #from seq_utils import ReverseComplement

    complement = {
                   "A":"T",
                   "a":"t",
                   "T":"A",
                   "t":"a",
                   "G":"C",
                   "g":"c",
                   "C":"G",
                   "c":"g",
                 }

    sequence_length = len(sequence)
    # Iterate over each k-mer of sequence, so the last kmer to be taken starts
    # len(sequence) - (k - 1)
    sequence_end = sequence_length - (k - 1)
    # Get the start offset to be taken into account for indexing this sequence
    start_offset = index["start_offset"]
    for i in xrange(0, sequence_end):
        kmer = sequence[i : i + k]
        # If that kmr has no ambiguous nucleotides
        if kmer in index:
            # Get also the reverse complement kmer
            rc_kmer = ReverseComplement(kmer, complement)
            # If that kmer is empty
            # AQUI HAY QUE IMPLEMENTAR UN DEFAULTDICT
            if not index[kmer]:
                index[kmer] = array("I", [i + start_offset])
            else:
                # If it is already present
                # Add to index the position where that k-mer starts plus the offset
                index[kmer].append(i + start_offset)

            # Now do the same for the reverse complement
            if not index[rc_kmer]:
                index[rc_kmer] = array("I", [i + start_offset])
            else:
                index[rc_kmer].append(i + start_offset)

    # Update the start offset with the length of this sequence
    index["start_offset"] += sequence_length
    
    return index