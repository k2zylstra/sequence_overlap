#! /bin/python3
import Test
from bitarray import bitarray
import sys

# bit string of values deliminated by 9 consecutive 0's
# h is a header table containing any and all locations of naturally occuring
#   deliminators
def sort(S, h, Ssize):
    # create bitvector of ones. Doesnt matter that I am doubling the memory.
    #   will fix this in the C implementation if I try and do it that way
    delim_size = 9
    cp_bitvec = bitarray(Ssize - 9)
    cp_bitvec.setall(False)
    cp_bitvec.append('111111111') # this is the delimiter portion of the array

    shift_amount = 0
    tmp = 0

    val1 = bitarray(22)
    val1.setall(False)
    val2 = bitarray(22)
    val2.setall(False)

    is_val1 = True

    # need to change this so that I am always looking at the second value except for the first time
    for i in range(Ssize):
        x = S & cp_bitvec
        if x == cp_bitvec:
            # delimiter start or end
            shift_amount += 0
            cp_bitvec << delim_size
            if is_val1:
                is_val1 = False
                #check if val1 is less than key value and if val2 is greater than key value
            else:
                is_val1 = True
        if val1:
            val1 |= (bitarray(False) * 21).append(True) << shift_amount
        else:
            val2 |= (bitarray(False) * 21).append(True) << shift_amount
                
        cp_bitvec << 1
        shift_amount += 1


def linsearch(B, a):
    pass

def si_counter(Bs, Be, a):
    first = linsearch(Bs, a[1])
    last = linsearch(Be, a[0])
    count = first - last
    return count


def si_counter(Bs, Be, a):
    first = linsearch(Bs, a[1])
    last = linsearch(Be, a[0])
    count = first - last
    return count

# A is a 2d array where each object is a tuple start and end value
#
# B is a 2D array that is stored so the the first item is a list
#   of the starts and the second item is a list of the ends
def perinterval_counter(A, B):
    B[0] = sorted(B[0])
    B[1] = sorted(B[1])
    C = []
    for i in range(len(A)):
        C.append(si_counter(B[0], B[1], A[i]))
    return C

def prefixsum(C):
    R = [C[0]]
    for i in range(1, len(C)):
        R.append(R[i-1] + C[i])
    return R

def enumerator(A, B):
    B[0] = sorted(B[0])
    B[1] = sorted(B[1])
    C = perinterval_counter(A, B)
    R = prefixsum(C)
    E = []
    start = 0
    for i in range(len(A)):
        end = R[i]
        f = binary_search(B[0], A[i][1])
        while end - start > 0:
            if A[i][0] <= B[1][f] and A[i][1] >= B[0][f]:
                E.append([A[i], [B[0][f], B[1][f]]])
                start += 1
            f -= 1
    
    return E

def main():
    #Test.test_bsearch(binary_search)
    e = Test.test_enumerator(enumerator)
    print(e)

if __name__ == "__main__":
    main()