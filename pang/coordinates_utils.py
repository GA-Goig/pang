def SortCoordinates(coordinates):
    '''This function simply sort a list of pairs os coordinates in tuples
    by ascending value of the start coordinate (first number of each tuple)

    It is a naive function just to avoid making same calculation more than one
    time in each program run 
    '''

    return sorted(coordinates, key=lambda x:x[0])

def GetNewCoreSeqs(scanned_sorted_coords, seq):
    '''GENERATOR: Given a list of sorted/fragment_joined coordinates of core 
    alignments, and the sequence that produced those alginments, Yield sequences
    that DO NOT BELONG TO THE CORE, (opposite seqs to coordinates provided), so 
    they could be reindexed and added as new core sequences
    
    Coordinates of regions aligned marked with ( )  

              FIRST                                      REMAINING
    Input =  |------(-------)--------(---)-------(-----)--------|
    Output = |-----|         |------|     |-----|       |-------| 

    '''
    
    seq_length = len(seq)

    # Check first if there are new seqs in the beginning  of the seq
    first_coords = scanned_sorted_coords[0] # Take first coordinates
    first_start = first_coords[0]
    # Check if first aligned coord is at the beginning of scanned seq
    if first_start == 0:
        pass
    else: # If first aligned coord is not at the beginning of seq
        # take first part of sequence as new_seq (as in example above)
        new_seq = seq[0 : first_start]
        yield (new_seq, (0, first_start))     

    for i in xrange(len(scanned_sorted_coords) - 1):
        tuple_A = scanned_sorted_coords[i] # First pair of coordinates (previous)
        tuple_B = scanned_sorted_coords[i+1] # Next pair of coordinates
        previous_end = tuple_A[1]
        next_start = tuple_B[0]
        new_seq = seq[previous_end : next_start]
        yield (new_seq, (previous_end, next_start))

    # Take the remaining sequence after last pair of coordinate if it didn't produce
    # any alignment
    # Check if the end of sequence is within a region aligned
    last_coords = scanned_sorted_coords[-1]
    end_coord = last_coords[1]
    if end_coord < seq_length: # Then the end of seq is not within an alignment
        new_seq = seq[end_coord : seq_length]
        yield (new_seq, (end_coord, seq_length))

def MapCoordinates(index_map, start_coord, end_coord):
    '''This function takes the index_map and a pair of start end coordinates
    and returns all records that coordinates belong to

    Index = |--------|    |-------|    |-------|    |-------|
            0   R1   10  20  R2   30  40   R3  50  50   R4  60
  
                |--------------------------|
                5                          45

    Given start_coord = 5, end_coord = 45, result will be [R1, R2, R3]
    '''
    records = []

    # Iterate over each pair of coordinates of index_map
    map_records = index_map[0]
    map_coords = index_map[1]
    for i in range(len(map_coords)):
        map_start = map_coords[i][0] # Start index coordinate for record i
        map_end = map_coords[i][1] # End index coordinate for record i
        if start_coord >= map_start: 
            if end_coord <= map_end:
                # Alignment happens just within a record
                records.append(map_records[i])
                return records
        else:
            # Start_coord was lower than current map_start, but greater than previuos 
            records.append(map_records[i - 1])
            # Check if also end_coord is lower than map_end
            if end_coord <= map_end:
                # In that case this is the last record to be appended,
                records.append(map_records[i])
                return records
            else:
                # If end_coord is not yet lower than map_end, still checking
                continue

def MapRecords(aligned_index, index_map):
    '''This function checks which regions of the indexed sequence(s) produced
    core alignments and returns a list with all records implicated in those
    alignments
    
    aligned_index is a list of tuples containing coordinates in the form 
    (start, end) so each pair is checked to see which records of the indexed
    sequence(s) comprise and added to final list
    '''
    records_aligned = set()
    for coordinates in aligned_index:
        start = coordinates[0]
        end = coordinates[1]
        records_matched = MapCoordinates(index_map, start, end)
        for record in records_matched:
            records_aligned.add(record)

    return records_aligned
  
def JoinFragments(sorted_coords, J):
    '''This function takes a list of sorted tuples of coordinates (start, end) 
    and join fragments that are as close as defined by user with parameter J'''


    # Join fragments that are no more distant that G value
    i = 0
    while i < len(sorted_coords) - 1:
        # Take contiguous tuples of coordinates
        tuple_A = sorted_coords[i]
        tuple_B = sorted_coords[i+1]

        # Take coordinates of each tuple
        start_coord_A = tuple_A[0]
        end_coord_A = tuple_A[1]
        start_coord_B = tuple_B[0]
        end_coord_B = tuple_B[1]
        distance = start_coord_B - end_coord_A # Then calculate distance         
        if distance <= J: # If distance is lower than value G
            new_start = start_coord_A
            new_end = end_coord_B
            # Update second tuple with new coordinates and delete first tuple            
            sorted_coords[i] = (new_start, new_end)   
            del sorted_coords[i+1]  
        else:
            i += 1

    return sorted_coords