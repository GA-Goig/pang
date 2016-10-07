#!/usr/bin/env python
# -*- coding: utf-8 -*-


def parse_args():
    '''Parse arguments given to script'''

    import argparse

    parser = argparse.ArgumentParser(description="Given two genomes, get both "\
        "common and specific regions to both genomes following a kmer-based algorithm")
    
    parser.add_argument("-d", dest="genome_dir", metavar="Genome dir", required=True, 
        help="Directory containing fasta or multifasta files to be clustered and a cinfo file")
    
    parser.add_argument("-r", dest="recursive", action="store_true", 
        help="Use for check all directories within directory provided by << -d >>")
    
    parser.add_argument("-g", metavar="jumping Gap between k-mers. DEFAULT: 6", dest="G", 
        default=6, help="Gap between k-mers. From inf to -k+1. (e.g G = 0 if no"\
        " gaps are desired. G = -K+1 if an overlap of 1 base between k-mers is desired")
    
    parser.add_argument("-k", metavar="k-mer length. DEFAULT: 12", dest="k", default=12, 
        help="k-mer length. Avoid use k > 12 in desktop computers")
    
    parser.add_argument("-f", metavar="min alignment length. DEFAULT: 48", dest="F", default=48,
        help="Distance checked by k-mer jumping from a seed k-mer to be considered as" \
        " a significant alignment")
    
    parser.add_argument("-j", metavar="max distance to combine fragments. DEFAULT: 25", 
        dest="J", default=25, help="Maximum distance to join significant alignments "\
        "from k-mer jumping")
    
    parser.add_argument("--max-seeds", metavar="maximum seeds per kmer. DEFAULT: 20", 
                        dest="max_seeds", default=20, help="Maximum seeds of a k-mer"\
                        " in order to seed alignments. This limit allows for a dramatic "\
                        "speedup in repetitive regions")
    
    parser.add_argument("-l, --length", dest="L", metavar="min cluster length. DEFAULT: 500", 
        default=500, help="Minimum length of a sequence in order to be considered as a new Cluster")
    
    parser.add_argument("-n", dest="N", metavar="max percentage of N per cluster. DEFAULT: 0.8",
     default=0.8, help="Maximum percentage of ambiguous nucleotides in clusters. If a new"\
     " cluster exceeds this maximum, it is discarded")

    args = parser.parse_args()

    return args

def Align(k_start, seed_coordinate, index, sequence, k, G):
    '''After a seeding match, get kmers from sequence moving k+G bases between
    each one, and "moving" k+G bases from the seed_coordinate too. 
    For each n-kmer, if it has a coordinate in index that coincides with 
    seed_coordinate + (k+G)*n, they are considered to be contiguous in both
    sequences, the scanned one and the indexed one. This function looks for all
    contiguous k-mers and returns the length of the extension produced'''
    from pang.seq_utils import GappedKmerGenerator
    from pang.binary_search import binary_search as bs
    gapped_kmer_gen = GappedKmerGenerator(sequence, k_start, k, G)
    current_index_coord = seed_coordinate
    kmer = gapped_kmer_gen.next()
    jump = k + G
    while kmer in index:
        kmer_index_coords = index[kmer] # Take coords for next k+G gapped kmer
        current_index_coord += jump # Move k+G bases of index
        # When both coordinates don't coincide they are not contiguous
        # binary search of current_index_coord in kmer_index_coords
        if bs(kmer_index_coords, current_index_coord):
            # If both coordinates coincide, get and check next gapped kmer
            kmer = gapped_kmer_gen.next()
        else:
            current_index_coord -= jump
            alignment_length = current_index_coord + k - seed_coordinate
            return alignment_length
    else:
        # I kmer not in index could be both, that is a kmer with ambiguous
        # nucleotides, or that kmer == 0 since the end of sequence has been
        # reached. In both cases, return the alignment length produced up to
        # this point.
        alignment_length = current_index_coord + k - seed_coordinate
        return alignment_length

def ExtendSeeds(k_start, seed_coordinates, index, sequence, k, G, F):
    '''This function takes all coordinates where a kmer of the scanned sequence
    is found in the indexed sequence and tries to make alignments starting from
    each seeding_coordinate. It returns a tuple of coordinates for the shortest
    aligned region of the scanned sequence and a list containing tuples of 
    coordinates for each aligment produced in the indexed sequence

    Parameter F allows to set a cutoff of a minimum extension for an aligment
    to be considered

    We keep the shortest alignment in order to continue scanning from the
    shortest fragment aligned to avoid pass some other k-mers that could seed
    another interesting alignments
    '''

    indexed_alignments = [] # To keep all aligments produced in indexed sequence
    shortest_alignment = int # int always > any value
    alignment = False # To check if any aligment > F has been produced
    alignments = [(-1,-1)] # To store coords for those regions already aligned
    for seed_coordinate in seed_coordinates:
        if CheckSeed(seed_coordinate, alignments):
            alignment_length = Align(k_start, seed_coordinate, index, sequence, k, G)
            if alignment_length >= F:
                if alignment_length < shortest_alignment:
                    shortest_alignment = alignment_length # Keep shortest alignment
                indexed_start = seed_coordinate
                indexed_end = seed_coordinate + alignment_length - 1
                indexed_alignments.append( (indexed_start, indexed_end) )
                alignments.append((indexed_start, indexed_end))
                alignment = True

    if alignment:
        scanned_start = k_start
        scanned_end = k_start + shortest_alignment - 1
        scanned_alignment = (scanned_start, scanned_end)
        return (scanned_alignment, indexed_alignments)  
    else:
        return 0

def CheckSeed(seed, alignments):
    '''This function takes a seed coord and a list of tuples of coordinates
    with those alignments that have been already produced. It eliminates
    seeds that fall in regions already aligned:

    |----(Seed 1)----------(Seed 2)----------------------|
         |---------------------------------------------->| Alignment from Seed 1
         
         Then eliminate Seed 2 TO AVOID THIS:
                           |---------------------------->| Alignment from Seed 2
                                                           (shortest_alignment)
    
    This function returns False if Seed falls in a region already aligned
    or True otherwise
    '''

    for alignment_coords in alignments:
        start = alignment_coords[0]
        end = alignment_coords[1]
        if seed >= start and seed <= end:
            return False

    return True


def SeedAndExtend(sequence, index, k, G, F, max_seeds, k_start=0):

    '''This function takes k-mers from a sequence and looks if they can seed 
    alignments with the "indexed sequence". If so, it tries to extend those
    alignments. If alignments of length > F are produced, then coordinates of
    that alignments are added to the "core coordinates" and continues scanning
    from the end of the shortest_alignment produced in the scanned sequence

    It returns two lists for scanned and indexed sequence that contains tuples
    of pair of coordinates (start, end) that correspond to core regions in both
    sequences'''
    
    from pang.seq_utils import KmerGenerator, SkipAmbiguous

    alignment_coordinates = []
    non_ambiguous={"A", "T", "G", "C"}
    min_non_ambiguous = k # Parameter for SkipAmbiguous
    kmer_gen = KmerGenerator(sequence, k_start, k) # Get kmers overlaping by 1

    kmer = kmer_gen.next() # Get first kmer
    while True:
        if kmer in index:
            seed_coordinates = index[kmer] # Get coordinates of that kmer in index
            if seed_coordinates and len(seed_coordinates) <= max_seeds: 
                # Try to obtain aligments for each seed
                alignments = ExtendSeeds(k_start, seed_coordinates, index, sequence,
                                         k, G, F)

                if alignments:
                    alignment_coordinates.append(alignments)
                    scanned_alignment = alignments[0]
                    # k_start now is last coordinate of the shortest alignment
                    k_start = scanned_alignment[1] + 1
                    # Start again the generator with updated k_start
                    kmer_gen = KmerGenerator(sequence, k_start, k)
                    kmer = kmer_gen.next() # And get new kmer
                else: # If not alignment produced
                    kmer = kmer_gen.next()
                    k_start += 1
            else:
                kmer = kmer_gen.next() # If not seed coordinates
                k_start += 1
        else:
            # If kmer not in index, check if it is != 0, in that case
            # kmer is a string not in index, therefore is a kmer containing
            # ambiguous DNA bases. In that case, pass over the ambiguous sequence
            # until at least a number of contiguous nucleotides defined by 
            # min_non_ambiguous are found and continue checking from that point
            if kmer != 0:
                # Check position after contiguous non_ambiguous nucleotides
                k_start = SkipAmbiguous(sequence, k_start, non_ambiguous, 
                                        min_non_ambiguous)
                # Start again the generator with updated k_start
                kmer_gen = KmerGenerator(sequence, k_start, k)
                kmer = kmer_gen.next() # And get new kmer
            # If kmer is == 0 then the generator has finished and the end of
            # the sequence has been reached. In that case, return
            else:
                return alignment_coordinates

def NewMappingGroup(mapping_file, gi):
    '''This function updates the mapping_file after new sequences have ben added
    to the core genome file'''
    
    with open(mapping_file, "a") as handle:
        handle.write(str(CURRENT_CLUSTER) + "\t" + gi + "\n")

def UpdateMappingGroups(mapping_file, groups, gi):
    '''This function updates the mapping file by adding current gi to each group
    it has produced an alignment with'''
    import shutil

    mapping_file_tmp = mapping_file + ".tmp"
    with open(mapping_file) as handle:
        # Open a temporary file for update the current mapping_file info
        with open(mapping_file_tmp, "w") as tmp:
            for line in handle:
                group, gis = line.split("\t")
                # Write an exact line if that group is not implied
                if group not in groups:
                    tmp.write(line)
                else:
                    # If not, add ;gi to the end of the line
                    tmp.write(line.rstrip() + ";" + gi + "\n")
    # And mv the overwrite the older version with the updated one
    shutil.move(mapping_file_tmp, mapping_file)

def AlignRecords(fasta, index, index_map, k, G, F, J, L, N, pangenome, mapping, max_seeds):
    '''This function iterates over each record in the fasta file and performs
    the SeedAndExtend function over each one, returning coordinates of "core"
    alignments for each record and indexed sequence in a dict where record
    title is the key'''
    import sys
    from pang.parse_utils import FastaParser, GetGI, GetRecordGroups
    from pang.coordinates_utils import SortCoordinates, JoinFragments
    from pang.coordinates_utils import  MapAlignments, GetNewCoreSeqs
    from pang.index_utils import ReindexRecord

    # USE global CORE_TITLE, CURRENT_CLUSTER
    global CORE_TITLE, CURRENT_CLUSTER

    title = "FILE_EMPTY"
    with open(fasta) as handle:
        for title, seq in FastaParser(handle):
            # Get the gi from current record
            gi = GetGI(title)
            if seq:
                # First produce alignments with indexed sequences and retrieve
                # coordinates of these alignments
                alignment_coordinates = SeedAndExtend(seq, index, k, G, F, max_seeds)
                if alignment_coordinates:
                    # Sort alignment coordinates
                    sorted_coordinates = SortCoordinates(alignment_coordinates)
                    # And Join fragments that are as close as parameter J
                    joined_coords = JoinFragments(sorted_coordinates, J)
                    # Get new sequences that do not produce core alignments
                    # This new sequences will be added to a new GROUP so first
                    # update the global variable with the new group
                    new_seqs = GetNewCoreSeqs(joined_coords, seq, L, N)
                    for new_seq, new_seq_start, new_seq_end in new_seqs:
                        CURRENT_CLUSTER += 1
                        cluster = "Cluster_" + str(CURRENT_CLUSTER)
                        # Format header to CORE_TITLE format
                        header = CORE_TITLE + cluster
                        # Add new sequences to index and update index_map
                        index = ReindexRecord(header, k, index, index_map, new_seq)
                        # Update pangenome dictionary with new_core_seq
                        pangenome[header] = new_seq
                        # Update mapping_dict with new groups
                        gi_coords = (0, len(new_seq), gi, "strand", 
                                     new_seq_start, new_seq_end)
                        mapping[cluster] = [gi_coords]


                else:
                    # If no alignment is produced, scanned_core_coords is evaluated
                    # as False since it contains an empty list of coordinates
                    # in that case all sequence is new
                    CURRENT_CLUSTER += 1
                    cluster = "Cluster_" + str(CURRENT_CLUSTER)
                    header = CORE_TITLE + cluster
                    index = ReindexRecord(header, k, index, index_map, seq)
                    # Update pangenome with a new_core_seq
                    pangenome[header] = seq
                    # Update mapping dict for that new_seq
                    gi_coords = (0, len(seq)-1, gi, "strand", 0, len(seq)-1)
                    mapping[cluster] = [gi_coords]
                    # If all sequence is new, there is no need to look which records
                    # each alignment maps to, so continue with following iteration of
                    # for loop
                    continue

            else: # If there is one empty seq
                sys.exit("Error: one or more scanned records are empty")

            # For sequences that were already aligned in the core genome
            # check which records they have been aligned with
            # Get records that have produced alignments

            mapped_alignments = MapAlignments(joined_coords, index_map)
            # If any core alignment has been produced
            if mapped_alignments:
                for scan_coords, clusters_mapped in mapped_alignments:
                    scan_start = scan_coords[0]
                    scan_end = scan_coords[1]
                    for cluster_mapped in clusters_mapped:
                        cluster = cluster_mapped[0]
                        c_start = cluster_mapped[1]
                        c_end = cluster_mapped[2]
                        mapping[cluster].append(
                            (c_start, c_end, gi, "strand", scan_start, scan_end )
                            )

 
    if title == "FILE_EMPTY": # If no title, seq returned by parser title
                              # remains <<FILE_EMPTY>>
        sys.exit("Error: scanned fasta file seems to be empty")

    return pangenome, mapping

def InitCore(core_info, genome_dir_path):
    '''This function initializes all files necessary to start build the core
    genome within a given directory. It gets the bigger fasta file as de first 
    indexed-reference file and creates a ".core" file containing that genome
    as the first record of the core genome belonging to Group 1.
    It also creates the .mapping file which will be initialize with Group 1 
    in first field and gi from that first reference in second field.

    ir returns .core and .mapping filenames with absolute path

    '''
    import os
    from pang.parse_utils import GetGI, FastaParser
    
    genus = core_info[0]
    species = core_info[1]
    core_title = core_info[2]

    # Create the genus_species.core string to create/stat the file
    core_file = genus + "_" + species + ".core"
    mapping_file = genus + "_" + species + ".mapping"
    # Add absolute path to each file
    core_file = os.path.normpath(os.path.join(genome_dir_path, core_file))
    mapping_file = os.path.normpath(os.path.join(genome_dir_path, mapping_file))
  
    return core_file, mapping_file

def BuildCore(genome_dir_path, k, G, F, J, L, N, index, max_seeds):
    '''This function takes a genome dir -that sould be a directory containing
    genome files in fasta format for a given specie- to build the set of 
    sequences that will form the core genome

    It creates:
    -a .mapping file that as a tab delimited text file containing a group number
    in the first field and all gi's pretaining to that group in the second field

    -a .core file, which is actually a fasta formated file containing all
    sequences that form the core genome

    Then it calls all necessary functions in order to align genomes and build
    the actual core

    It uses a core.info file that should have been created previously by the script
    that splits genomes from refseq .fna files in respective species folders.

    core.info has a row with common header for every sequence that will form the
    core, and another row with the number of different core groups
    '''

    from os.path import normpath, realpath
    from os.path import join as pjoin
    from pang.pang_IO import LoadInfo, WritePangenome, WriteMapping, ListFasta
    import sys
    import glob

    genome_dir_path = realpath(genome_dir_path)
    # Load the CORE_TITLE and CURRENT_CLUSTER, from core.cinfo file
    info = normpath(pjoin(genome_dir_path, "core.info"))
    if glob.glob(info):
        core_info = LoadInfo(info)
        # Making it global variables
        global CURRENT_CLUSTER, CORE_TITLE
        CORE_TITLE = core_info[2]
        CURRENT_CLUSTER = core_info[3]

        # Get genomes list from genome_dir_path as a deque object
        genomes_list = ListFasta(genome_dir_path)
     
        # Initialize necessary files in the directory as .core or .mapping
        core_file, mapping_file = InitCore(core_info, genome_dir_path)
        
        # Build first iteration index and index_map
        # index, index_map = IndexRecords(core_file, k)

        # Build a dict with core headers has keys and core seqs as values
        pangenome = {}
        # Build a dict with groups as keys and gi's as values
        mapping = {}

        # Initialize index map
        index_map = [ [] , [] ]

        # Perform alignments of records for each of the remaining genomes in list
        for genome in genomes_list:
            print "Calculating pangenome for {}...".format(genome)
            # Get absolute path to each genome
            genome = normpath(pjoin(genome_dir_path, genome))
            # And update pangenome and mapping dicts
            pangenome, mapping = \
            AlignRecords(genome, index, index_map, k, G, F, J, L, N, pangenome, mapping,
                         max_seeds)
        # Once all records have been aligned, write final core and mapping from
        # pangenome and mapping dicts
        WritePangenome(pangenome, core_file)
        WriteMapping(mapping, mapping_file )

def ProcessGenomesDir(genomes_dir, k, G, F, J, L, N, max_seeds):
    '''This function takes a genomes_dir where refseq records are stored
    in separated folders by species and calls BuildCore in each one'''
    import os
    from array import array
    from time import time
    from pang.index_utils import BuildIndex
    
    genomes_dir = os.path.realpath(genomes_dir)
    start_run_time = time()
    # first create the empty index of k-mers
    index = BuildIndex(k)

    processed_dirs = 0
    total_dirs = len(os.listdir(genomes_dir))
    for species_dir in os.listdir(genomes_dir):
        time_start = time()
        species_dir = os.path.normpath(os.path.join(genomes_dir, species_dir))
        BuildCore(species_dir, k, G, F, J, L, index, max_seeds)
        processed_dirs += 1
        
        # When finished, set index again to empty
        # arrays in order to be used by next pangenome
        index = index.fromkeys(index, array("I"))
        
        runtime = time() - start_run_time
        time_end = time() - time_start
        print "Pangenome calculated in {} seconds".format(time_end)
        print "Processed {} of {} species in {} seconds".format(processed_dirs, total_dirs, runtime)
        
        index["start_offset"] = 0

def ProcessDir(genome_dir, k, G, F, J, L, N, max_seeds):
    import os
    from time import time
    from array import array
    from pang.index_utils import BuildIndex

    genome_dir = os.path.realpath(genome_dir)
    start_time = time()
    index = BuildIndex(k)
    BuildCore(genome_dir, k, G, F, J, L, N, index, max_seeds)
    index = index.fromkeys(index, array("I"))
    time_end = time() - start_time
    print "Pangenome calculated in {} seconds".format(time_end)

def main():
    args = parse_args()
    genome_dir = args.genome_dir
    # scanned = args.scanned
    G = int(args.G)
    k = int(args.k)
    F = int(args.F)
    J = int(args.J)
    L = int(args.L)
    N = float(args.N)
    max_seeds = int(args.max_seeds)
    recursive = args.recursive
    
    if recursive:
        ProcessGenomesDir(genome_dir, k, G, F, J, L, N, max_seeds)
    elif not recursive:
        ProcessDir(genome_dir, k, G, F, J, L, N, max_seeds)
    else:
        assert False

if __name__ == "__main__":
    main()


