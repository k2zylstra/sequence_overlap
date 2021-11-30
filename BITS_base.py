#! /bin/python3
import Test


#only difference between these two is the first if statement that 
#   is arguably not necessary; update it is worse
def binary_search(B, a):
    j = len(B) // 2
    i = 0
    while j >= 1:
        if B[i] > a:
            i = i - j
        elif B[i] < a:
            i = i + j
        else:
            break 
        j = j // 2
    return i

def si_counter(Bs, Be, a):
    first = binary_search(Bs, a[1])
    last = binary_search(Be, a[0])
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