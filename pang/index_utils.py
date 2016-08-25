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
            # If that kmer is empty
            # AQUI HAY QUE IMPLEMENTAR UN DEFAULTDICT
            if not index[kmer]:
                index[kmer] = array("I", [i + start_offset])
            else:
                # If it is already present
                # Add to index the position where that k-mer starts plus the offset
                index[kmer].append(i + start_offset)

    # Now index the reverse strand 
    # (So we can detect inversions)
    sequence = sequence[::-1]
    for i in xrange(0, sequence_end):
        kmer = sequence[i : i + k]
        # If that kmr has no ambiguous nucleotides
        if kmer in index:
            # If that kmer is empty
            # AQUI HAY QUE IMPLEMENTAR UN DEFAULTDICT
            if not index[kmer]:
                index[kmer] = array("I", [i + start_offset])
            else:
                # If it is already present
                # Add to index the position where that k-mer starts plus the offset
                index[kmer].append(i + start_offset)

    # Now index the reverse_complement strand 
    # (So the complemetary strand is aligned too)
    sequence = "".join([complement[nt] if nt in complement else "N" for nt in sequence])
    for i in xrange(0, sequence_end):
        kmer = sequence[i : i + k]
        # If that kmr has no ambiguous nucleotides
        if kmer in index:
            # If that kmer is empty
            # AQUI HAY QUE IMPLEMENTAR UN DEFAULTDICT
            if not index[kmer]:
                index[kmer] = array("I", [i + start_offset])
            else:
                # If it is already present
                # Add to index the position where that k-mer starts plus the offset
                index[kmer].append(i + start_offset)
                
    # Now index the "Forward complement"
    # (So we can detect inversions in the complemetary strand)
    sequence = sequence[::-1]
    for i in xrange(0, sequence_end):
        kmer = sequence[i : i + k]
        # If that kmr has no ambiguous nucleotides
        if kmer in index:
            # If that kmer is empty
            # AQUI HAY QUE IMPLEMENTAR UN DEFAULTDICT
            if not index[kmer]:
                index[kmer] = array("I", [i + start_offset])
            else:
                # If it is already present
                # Add to index the position where that k-mer starts plus the offset
                index[kmer].append(i + start_offset)

    # Update the start offset with the length of this sequence
    index["start_offset"] += (sequence_length)
    
    return index

def ReindexRecord(header, k, index, index_map, new_seq):
    ''' After a record has been aligned, this function takes sequences from that
    record that do not produce core alignments, reindex them, and update the
    index keeping updated the index_map too, in order to now which are the
    coordinates in <<index>> corresponding to new record sequences
    '''
    # Next calc which coordinate this record sequences will be indexed from
    start_record = index["start_offset"]
    # New seqs will store new seqs to be written to CORE GENOME
    index = IndexSequence(new_seq, k, index)
    # End record will coincide with new start_offset value
    end_record = index["start_offset"] - 1
    # Update index_map records list with new record
    index_map[0].append(header)
    # Update in same position of coords list new coords
    index_map[1].append( (start_record, end_record) )
    # Return index to get the updated version
    return index