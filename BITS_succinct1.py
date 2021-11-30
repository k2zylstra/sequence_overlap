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
    cp_bitvec = bitarray(Ssize)
    for i in range(Ssize):
        

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