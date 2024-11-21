# Authors : Sankalp Mittal (395001), Ilia Badanin (350775), Vasco Frazão (396229)


import math
import random
import numpy as np

def isPrime(n):
    if n <= 1:
        return False
    if n <= 3:
        return True
    if n % 2 == 0 or n % 3 == 0:
        return False
    i = 5
    while i * i <= n:
        if n % i == 0 or n % (i + 2) == 0:
            return False
        i += 6
    return True


# fine unless we have prime over p_{20'000}
def next_prime(n):
    while True:
        if isPrime(n):
            return n
        n += 1


def npsolver(weights, n, m, b, t, c):
    max_wt = max(t, c)
    min_siz = max(max_wt * n + 1, n * n)
    prime = next_prime(min_siz + 1)
    
    set_of_vals = np.arange(prime)

    weights = np.array(weights)
    
    num_runs = math.ceil(math.log(n))
    for _ in range(num_runs):
        alpha = np.random.randint(1,prime+1, (n,n))
        gammas = np.random.choice(set_of_vals[1:], size=n * max_wt + 1, replace=False)
        r_vals =  np.linalg.det((gammas[:, None, None].astype(float))**weights[None, ...]*(weights != -1)).astype(int)
    
        P_matrix = gammas[:,None] ** np.arange(len(gammas))[None,:]
        
        c = np.linalg.solve(P_matrix, r_vals)
    
        if c[b] % prime != 0:
            print("yes")
            return 0
    
    print("no")
    return 0

def main():
    n, m, b, t, c = map(int, input().split())
    weights = [[-1] * n for _ in range(n)]
    for _ in range(m):
        u, v, w = map(int, input().split())
        weights[u][v] = w
    return npsolver(weights, n, m, b, t, c)


if __name__ == "__main__":
    main()
