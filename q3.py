# Authors : Sankalp Mittal (395001), Tamar Alphaidze (393635)

import math
import sympy
import numpy as np
import random

def isPrime(n):
    if n <= 1:
        return False
    if n <= 3:
        return True
    if n%2 == 0 or n%3 == 0:
        return False
    i = 5
    while i*i <= n:
        if n%i == 0 or n%(i+2) == 0:
            return False
        i += 6
    return True

def next_prime(n):
    while True:
        if isPrime(n):
            return n
        n += 1

def get_rand_val(set_of_vals):
    return random.choice(set_of_vals)

def get_rand_vector(set_of_vals, m):
    return random.sample(set_of_vals, m)

def mod_matrix_inv(matrix, prime):
    n = len(matrix)
    inv_matrix = np.zeros((n,n))
    # TODO : Write the function since the smypy one is not working
    det_K = np.linalg.det(matrix).round()
    det_K_mod = int(((det_K%prime)+prime)%prime)
    # print(det_K)
    det_inv = pow(det_K_mod,-1, prime)
    # print(det_inv)
    adj_matrix = np.linalg.inv(matrix)
    for i in range(n):
        for j in range(n):
            adj_matrix[i][j] = (((adj_matrix[i][j]*det_K).round() % prime) + prime)%prime
    for i in range(n):
        for j in range(n):
            inv_matrix[i][j] = (det_inv*adj_matrix[i][j])%prime
    return inv_matrix

def main():
    n, m, b, t, c = map(int, input().split())
    weights = [[-1] * n for _ in range(n)]
    for _ in range(m):
        u, v, w = map(int, input().split())
        weights[u][v] = w

    max_wt = max(t, c)
    min_siz = max(max_wt * n + 1, n * n)
    prime = next_prime(min_siz+1)
    # print(prime)
    set_of_vals = list(range(prime))

    num_runs = math.ceil(math.log(n))
    isTrue = False
    for _ in range(num_runs):
        alpha = [[get_rand_val(set_of_vals[1:]) for _ in range(n)] for _ in range(n)]
        gammas = get_rand_vector(set_of_vals[1:], n * max_wt + 1) # Sliced to prevent 0

        # #Testing
        # alpha = [[0,5,0],[9,7,1],[2,7,1]]
        # gammas = [2,1,6,8,4,10,9,5,7,3]
        # #Test end

        # print(alpha)
        # print(gammas)

        y = sympy.symbols('y')
        H_matrix = [[0]*n for _ in range(n)]
        for i in range(n):
            for j in range(n):
                if weights[i][j] != -1:
                    H_matrix[i][j] = alpha[i][j] * (y**weights[i][j])
        # print(H_matrix)
        H_matrix = sympy.Matrix(H_matrix)
        H_det = H_matrix.det()
        # print(H_det)

        r_vals = []
        for gam in gammas:
            r_vals.append(H_det.subs(y, gam)%prime)
        s_vector = np.array(r_vals)
        P_matrix = np.zeros((len(gammas), len(gammas)))
        for i in range(len(gammas)):
            for j in range(len(gammas)):
                P_matrix[j][i] = pow(gammas[j], i,prime)
        # print(P_matrix)


        P_inv = mod_matrix_inv(P_matrix, prime)
        # print(P_inv)

        # should_be_identity = P_inv @ P_matrix
        # for i in range(len(should_be_identity)):
        #     for j in range(len(should_be_identity[i])):
        #         should_be_identity[i][j] = should_be_identity[i][j]%prime
        # print(should_be_identity)

        c = P_inv.dot(s_vector)
        for i in range(len(c)):
            c[i] = c[i]%prime
        # print(c)
        if c[b] != 0:
            # print(c[b])
            print("yes")
            isTrue = True
            break
    if not isTrue:
        print("no")
    return 0

if __name__ == "__main__":
    main()