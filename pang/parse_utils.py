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

def GetGiList(record_string):
    '''This function returns the gi list from a fasta header with a format
    used by the database built by FastCore. Given following header:

    >gi|xxxxxxxxx|ref|zzzzzzzzz|@Genus species@comprising:gi1;gi2;gi3

    it returns the list [gi1, gi2, gi3]

    or raises an error if header do not fit the format
    '''
    import sys
    try:
        header = record_string.split("@")
        gi_list = header[2]
        gi_list = gi_list.split(":")[1]
        gi_list = gi_list.split(";")
        return gi_list
    except:
        sys.exit("One or more records do not fit FastCore header format:\
    >gi|xxxxxxxxx|ref|zzzzzzzzz|@Genus species@comprising:gi1;gi2;gi3\
    Error in record: {}".format(record_string))

def GetRecordGroups(records):
    '''This function returns a list of groups extracted from a list of records
    headers that fit the following format:

        >taxid|1773| Mycobacterium_tuberculosis Group:1
        >taxid|1773| Mycobacterium_tuberculosis Group:56

    so it would return [1,56]
    '''
    groups = []
    for record in records:
        record 
        # First divide three parts of the header
        parts = record.split("@")
        # Take the third one, which is <<Group:X>>
        group = parts[2]
        # Then take only the number X
        group = group.split(":")[1]
        # And add it to final group list
        groups.append(group)
    
    return groups