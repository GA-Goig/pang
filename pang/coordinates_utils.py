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
    that DO NOT BELONG TO THE CORE, so they could be reindexed and added as new
    core sequences
    '''

    for i in xrange(len(scanned_sorted_coords) - 1):
        tuple_A = scanned_sorted_coords[i] # First pair of coordinates (previous)
        tuple_B = scanned_sorted_coords[i+1] # Next pair of coordinates
        previous_end = tuple_A[1]
        next_start = tuple_B[0]
        yield seq[previous_end : next_start]

def MapCoordinates(index_map, start_coord, end_coord):
    '''This function takes the index_map and a pair of start end coordinates
    and returns all records that coordinates belong to

    Index = |--------|    |-------|    |-------|    |-------|
            0   R1   10  20  R2   30  40   R3  50  50   R4  60
  
                |--------------------------|
                5                          45

    Given start_coord = 5, end_coord = 45, result sould be [R1, R2, R3]
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
    core alignments and return a list with all records implicated in those
    alignments
    
    aligned_index is a list of tuples containing coordinates in the form 
    (start, end) so each pair is checked to see which records of the indexed
    sequence(s) comprise and added to final list
    '''
    records_aligned = []
    records_matched = False
    for coordinates in aligned_index:
        start = coordinates[0]
        end = coordinates[1]
        records_matched = MapCoordinates(index_map, start, end)
        records_aligned.append(records_matched)

    return records_matched
  
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