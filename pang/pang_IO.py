def WritePangenome(pgnome_dict, pgnome_file):
     # Open with append, to add new sequences
     with open(pgnome_file, "w") as outfile:
        for record in pgnome_dict:
            new_seq = pgnome_dict[record]
            outfile.write(record + "\n")
            # The ruler is to write sequences
            # formated with a maximum of 100
            # columns
            ruler = 0
            for i in xrange(len(new_seq)):
                if ruler < 101:
                    outfile.write(new_seq[i])
                    ruler += 1
                else:
                    outfile.write("\n" + new_seq[i])
                    ruler = 0
            outfile.write("\n")

def WriteMapping(mapping_dict, mapping_file):
    with open(mapping_file, "w") as outfile:
        for cluster in mapping_dict:
            alignments = mapping_dict[cluster]
            for alignment in alignments:
                c_start, c_end, acc, strand, seq_start, seq_end = alignment
                score = len(alignments)
                string = "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
                    cluster, c_start, c_end, acc, score, strand, seq_start, seq_end
                    )
                outfile.write(string)

def WriteNewCoreSeqs(new_core_seq, core_file, header):
    '''Writes new sequences added to core genome in the core genome file
  
    This function takes some GLOBAL variables:

    current_group is an int value of the new group to be generated to map te gi
    identifier where new sequences come from

    CORE_TITLE is a string of the common part of the header of all sequences
    comprising the core genome
    
    For example, building a Mycobacterium tuberculosis core genome

    >taxid|1773| Mycobacterium_tuberculosis Group:1

    Group 1 will correspond to some gi in an additional file called mapping_file
    So in this example, if Group 1 is mapping to gi's X, Y, Z, you can tell this
    sequence is a Mycobacterium tuberculosis sequence and, more specifically, 
    common to sequences with gene identifiers X, Y, Z
    
    So the common part, that is stored <<CORE_TITLE>> would be:

    >taxid|1773| Mycobacterium tuberculosis Group:

    '''
    
    with open(core_file, "a") as handle:
        WriteSeq(handle, header, new_core_seq)

def LoadInfo(core_info):
    '''This function loads global variables from core.cinfo file.
    It loads the common header for all sequences that will form the core genome
    and the number of groups that '''

    with open(core_info) as infile:
        GENUS = infile.readline().rstrip()
        GENUS = GENUS.split("=")[1]

        SPECIES = infile.readline().rstrip()
        SPECIES = SPECIES.split("=")[1]
        # Read first line that corresponds to CORE_TITLE
        CORE_TITLE = infile.readline().rstrip()
        # And take the second field that is the core title itself
        CORE_TITLE = CORE_TITLE.split("=")[1]
        # Read following line that corresponds to GROUPS
        CURRENT_GROUP = infile.readline().rstrip()
        # Second field corresponds to number of groups
        CURRENT_GROUP = CURRENT_GROUP.split("=")[1]

    
    return (GENUS, SPECIES, CORE_TITLE, int(CURRENT_GROUP))

def ListFasta(genome_dir_path, compressed=False):
    '''This function gets all files within genome_dir and keeps only those with
    fasta, fna or fa extension, returning a list with them sorted descending 
    size in bytes

    This way it is possible to take the biggest genome file in order to use it
    as first reference-indexed genome to begin genome alignments.

    It is intended to work with directories that hold several genome files, from
    2 or 3 to thousands of them, plus some necessary files like .mapping or 
    .cinfo
    '''
    import re
    import os
    from os.path import normpath, join
    
    # List all files within genome_dir
    files = os.listdir(genome_dir_path)
    # Compile regexp to match only fasta|fna|fa extension files
    if compressed:
        extension = re.compile("\.gz$")
    else:
        extension = re.compile("\.(fasta|fna|fa)$")
        # Keep only files with that extension
    fasta_files = filter(extension.search, files)
    
    fasta_files = [ normpath(join(genome_dir_path, fasta)) for fasta in fasta_files ]
    # Sort fasta_files in-place by ascending size so the last element will be 
    # the larger genome
    fasta_files.sort(key=os.path.getsize)

    return fasta_files

def WriteSeq(handle, header, seq, ruler=101):
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