import itertools

def generate_strings(s_size, alphabet):
    l = map(''.join, itertools.product(alphabet, repeat=s_size))
    return l    

x = generate_strings(22, "01")
x = list(x)
S = len(x)

for m in range(22):
    c = 0
    s = ''.join('0' for i in range(m))
    for i in range(S):
        if s in x[i]:
            c += 1

    print(c/S)