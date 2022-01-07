import itertools

def generate_strings(s_size, alphabet):
    l = map(''.join, itertools.product(alphabet, repeat=s_size))
    return l    

def get_all_probs(base_int_size):
    probs = []
    x = generate_strings(base_int_size, "01")
    x = list(x)
    S = len(x)

    for m in range(base_int_size):
        c = 0
        s =  ''.join('0' for i in range(m))
        for i in range(S):
            if s in x[i]:
                c += 1
        probs.append(c/S)
    
    return probs

def get_prob(base_int_size, delim_size):
    probs = []
    x = generate_strings(base_int_size, "01")
    x = list(x)
    S = len(x)

    c = 0
    s = ''.join('0' for i in range(delim_size))
    for i in range(S):
        if s in x[i]:
            c += 1
    print(c)
    return c/S

def main():
    #get_all_probs(22)
    print(get_prob(9,2))

if __name__ == "__main__":
    main()