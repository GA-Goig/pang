from bisect import bisect_left

def binary_search(sorted_array, x, lo=0):  
    if sorted_array:
        hi = len(sorted_array) # hi defaults to len(a)   
        pos = bisect_left(sorted_array, x, lo, hi)  # find insertion position
        if pos != hi and sorted_array[pos] == x:
            return True
        else:
            return False
    else:
        return False