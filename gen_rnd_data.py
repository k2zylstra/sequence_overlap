import math
import random
import time
import os
import string
import numpy

def generate_sequence(size):
    ran = round(random.random())
    sample = "ACGT"
    s = ''.join((random.choice(sample)) for x in range(size))
    return s

def create_AandB(size_chrom, num_tocreate):
    A = []
    B = [[],[]]
    for i in range(num_tocreate):
        A.append([int(random.random() * size_chrom), None])
        B[0].append(int(random.random() * size_chrom))
    for i in range(num_tocreate):
        A[i][1] = int(random.random() * size_chrom)
        B[1].append(int(random.random() * size_chrom))
    
    return (A, B)

def test_enumerator(enumerator):
    A, B = create_AandB(5000, 300)
    ret = enumerator(A, B)
    return ret

def test_bsearch(bsearch):
    rs = [99, 1000, 10000, 100000, 200000]
    repeat = 999
    depth_list = []
    for k in range(len(rs)):
        r = rs[k]
        tl = math.log(r, 2)
        depth_list.append([0,0,tl, 0])
        for j in range(repeat):
            Bs = []
            for i in range(r):
                Bs.append(int(random.random()* r))
            Bs.sort()

            Be = []
            for i in range(len(Bs)):
                Be.append(Bs[i] + int(random.random() * r))
            Be.sort()

            a_s = int(random.random() * r)
            ae = int(random.random() * r + a_s)

            start = time.perf_counter_ns()
            indexe, depthe = bsearch(Bs, ae) 
            indexs, depths = bsearch(Be, a_s)
            end = time.perf_counter_ns()
            
            depth_list[k][1] += depthe
            depth_list[k][0] += depths
            depth_list[k][3] += end - start
        print("Length:\t", r, "\n")
        print("tree depth:\t\t", depth_list[k][2])
        print("average end depths:\t", depth_list[k][1] / repeat)
        print("average start depths:\t", depth_list[k][0] / repeat)
        print("average time:\t\t", depth_list[k][3] / repeat)
        
        print("==========")
        print("==========\n")

def main():
    s = generate_sequence(1000)


if __name__ == "__main__":
    main()