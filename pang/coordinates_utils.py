def SortCoordinates(coordinates):
    '''This function simply sort a list of pairs os coordinates in tuples
    by ascending value of the start coordinate (first number of each tuple)

    It is a naive function just to avoid making same calculation more than one
    time in each program run 
    '''

    return sorted(coordinates, key=lambda x:x[0])

def GetNewCoreSeqs(joined_coords, seq, L):
    '''GENERATOR: Given a list of sorted/fragment_joined coordinates of core 
    alignments, and the sequence that produced those alginments, Yield sequences
    that DO NOT BELONG TO THE CORE, (opposite seqs to coordinates provided), so 
    they could be reindexed and added as new core sequences
    
    Coordinates of regions aligned marked with ( )  

              FIRST                                      REMAINING
    Input =  |------(-------)--------(---)-------(-----)--------|
    Output = |-----|         |------|     |-----|       |-----q--| 
    
    Yield sequence that have at least length L
    '''
    
    seq_length = len(seq)
    # joined_coords are tuples pairs of scanned_coordinates with an
    # asociated list of indexed coordinates they aligned with, get
    # only scanned coordinates in order to get new sequences
    scanned_coords = [coords[0] for coords in joined_coords]
    # Check first if there are new seqs in the beginning  of the seq
    first_coords = scanned_coords[0] # Take first coordinates
    first_start = first_coords[0]
    # Check if first aligned coord is at the beginning of scanned seq
    if first_start == 0:
        pass
    else: # If first aligned coord is not at the beginning of seq
        # take first part of sequence as new_seq (as in example above)
        new_seq = seq[0 : first_start]
        if len(new_seq) >= L:
            yield (new_seq, (0, first_start))     

    for i in xrange(len(scanned_coords) - 1):
        tuple_A = scanned_coords[i] # First pair of coordinates (previous)
        tuple_B = scanned_coords[i+1] # Next pair of coordinates
        previous_end = tuple_A[1]
        next_start = tuple_B[0]
        new_seq = seq[previous_end : next_start]
        if len(new_seq) >= L:
            yield (new_seq, (previous_end, next_start))

    # Take the remaining sequence after last pair of coordinate if it didn't produce
    # any alignment
    # Check if the end of sequence is within a region aligned
    last_coords = scanned_coords[-1]
    end_coord = last_coords[1]
    if end_coord < seq_length: # Then the end of seq is not within an alignment
        new_seq = seq[end_coord : seq_length]
        if len(new_seq) >= L:
            yield (new_seq, (end_coord, seq_length))

def MapCoordinates(index_map, coords):
    '''This function takes the index_map and a pair of start end coordinates
    and returns all records that coordinates belong to

    Index = |--------|    |-------|    |-------|    |-------|
            0   R1   10  20  R2   30  40   R3  50  50   R4  60
  
                |--------------------------|
                5                          45

    Given start_coord = 5, end_coord = 45, result will be [R1, R2, R3]
    '''
    print "index_map = {}".format(index_map)
    print "coords = {}".format(coords)

    records = []
    start_coord = coords[0]
    end_coord = coords[1]

    # Iterate over each pair of coordinates of index_map
    map_records = index_map[0]
    map_coords = index_map[1]
    for i in range(len(map_coords)):
        map_start = map_coords[i][0] # Start index coordinate for record i
        map_end = map_coords[i][1] # End index coordinate for record i
        if start_coord < map_end:
            record = map_records[i]
            if end_coord <= map_end:
                records.append(record)
                return records  # Last record to be included
                                # no need to still checking
            else:
                records.append(record)
    # coordinates at are not actually in the index are provided, and this
    # should not occur either

    assert False

def MapAlignments(joined_coords, index_map):
    '''This function checks which regions of the index have been aligned
    by each region of the scanned_sequence. It takes a joined_coords
    list containing tuples of: A - coordinates of a region from the 
    scanned sequence and B - index coordinates list that region aligns with 

    For example: [ (0,100), [(400,500), (600,700)]] 

    In this case region from 0 to 100 in scanned sequence produced two
    alignments, from 400 to 500 and from 600 to 700 in the index, so
    this could be a duplication. Say those coordinates from the index
    belong to two different groups, G1 and G2. In that case
    this function should return ( (0,100), ["G1", "G2"]) So we now that
    region 0-100 from current sequence is in group1 and group2.

    '''

    records_map = []
    for scan_coords, index_coords in joined_coords:
        records_matched = []
        for coords in index_coords:
            records = (MapCoordinates(index_map, coords))
            for record in records:
                if record not in records_matched:
                    records_matched.append(record)

        records_map.append( (scan_coords, records_matched) )

    return records_map
  
def JoinFragments(sorted_coords, J):
    '''This function takes a list of sorted tuples of coordinates (start, end) 
    and join fragments that are as close as defined by user with parameter J'''
    from copy import copy

    scoords = copy(sorted_coords)
    # Join fragments that are no more distant that G value
    i = 0
    while i < len(scoords) - 1:
        # Take contiguous t uples of coordinates
        tuple_A = scoords[i][0]
        tuple_B = scoords[i+1][0]
        index_coords_A = scoords[i][1]
        index_coords_B = scoords[i+1][1]

        # Take coordinates of each tuple
        start_coord_A = tuple_A[0]
        end_coord_A = tuple_A[1]
        start_coord_B = tuple_B[0]
        end_coord_B = tuple_B[1]
        distance = start_coord_B - end_coord_A # Then calculate distance         
        if distance <= J: # If distance is lower than value J
            new_start = start_coord_A
            new_end = end_coord_B
            new_index_coords = index_coords_A + index_coords_B
            # Update second tuple with new coordinates and delete first tuple            
            scoords[i] = ( (new_start, new_end), new_index_coords)   
            del scoords[i+1]  
        else:
            i += 1

    return scoords
