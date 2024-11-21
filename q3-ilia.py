# Authors : Sankalp Mittal (395001), Ilia Badanin (350775), Vasco Frazão (396229)


import math
import random


def determinant_mod(matrix, prime):
    """
    Gaussian elimination with modular arithmetic to compute the determinant.
    All calculations are performed modulo `prime`.
    """
    # Get the size of the matrix
    n = len(matrix)

    # Create a copy of the matrix to avoid mutating the input
    mat = [[element % prime for element in row] for row in matrix]

    # Initialize determinant as 1
    det = 1

    for i in range(n):
        # Find pivot for column i
        pivot = i
        for j in range(i + 1, n):
            if mat[j][i] > mat[pivot][i]:
                pivot = j

        # If pivot is zero, determinant is zero
        if mat[pivot][i] == 0:
            return 0

        # Swap rows if needed
        if pivot != i:
            mat[i], mat[pivot] = mat[pivot], mat[i]
            det = (-det) % prime  # Swapping rows flips the sign of the determinant

        # Multiply determinant by the pivot element
        det = (det * mat[i][i]) % prime

        # Normalize row i (modular division by mat[i][i])
        inv = pow(mat[i][i], -1, prime)  # Modular inverse of the pivot
        for j in range(i, n):
            mat[i][j] = (mat[i][j] * inv) % prime

        # Eliminate column i for rows below
        for j in range(i + 1, n):
            factor = mat[j][i]
            for k in range(i, n):
                mat[j][k] = (mat[j][k] - factor * mat[i][k]) % prime

    return det


def linsolve(A, b, prime):
    """
    Solves the system of linear equations A * x = b using Gaussian elimination
    modulo a prime number.

    Parameters:
    A (list of lists): Coefficient matrix (n x n).
    b (list): Right-hand side vector (n).
    prime (int): A prime number for modulo operations.

    Returns:
    list: Solution vector x (n).
    """
    n = len(A)

    # Forward elimination: Reduce A to upper triangular form
    for i in range(n):
        # Find the pivot row
        pivot = i
        for j in range(i + 1, n):
            if abs(A[j][i]) > abs(A[pivot][i]):
                pivot = j

        # Swap rows in A and b if pivot changes
        if pivot != i:
            A[i], A[pivot] = A[pivot], A[i]
            b[i], b[pivot] = b[pivot], b[i]

        # Ensure the pivot is non-zero modulo prime
        if A[i][i] % prime == 0:
            raise ValueError("Matrix is singular or not invertible modulo prime.")

        # Eliminate entries below the pivot
        for j in range(i + 1, n):
            # Compute the factor to zero out A[j][i]
            factor = (A[j][i] * pow(A[i][i], -1, prime)) % prime  # Modular division
            for k in range(i, n):
                A[j][k] = (A[j][k] - factor * A[i][k]) % prime
            b[j] = (b[j] - factor * b[i]) % prime

    # Back substitution: Solve for x in Ux = c
    x = [0] * n
    for i in range(n - 1, -1, -1):
        x[i] = b[i]
        for j in range(i + 1, n):
            x[i] = (x[i] - A[i][j] * x[j]) % prime
        x[i] = (x[i] * pow(A[i][i], -1, prime)) % prime  # Modular division

    return x


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


def hdet(gamma, alpha, weights, n, prime):
    H_matrix = [
        [
            (alpha[i][j] * pow(gamma, weights[i][j], prime)) % prime
            if weights[i][j] != -1 else 0
            for j in range(n)
        ]
        for i in range(n)
    ]
    return determinant_mod(H_matrix, prime)


def get_P(r_vals, gammas, prime):
    # return [[pow(g, i, prime) for i in range(len(gammas))] for g in gammas]
    return [[pow(g % prime, i, prime) for i in range(len(gammas))] for g in gammas]



def solver(weights, n, m, b, t, c):
    max_wt = max(t, c)
    min_siz = max(max_wt * n + 1, n * n)
    prime = next_prime(min_siz + 1)

    set_of_vals = list(range(prime))

    num_runs = math.ceil(math.log(n))
    for _ in range(num_runs):
        alpha = [[random.choice(set_of_vals[1:]) for _ in range(n)] for _ in range(n)]
        gammas = random.sample(set_of_vals[1:], n * max_wt + 1)
        # print(gammas)

        r_vals = [hdet(gamma, alpha, weights, n, prime) for gamma in gammas]

        P_matrix = get_P(r_vals, gammas, prime)

        c = linsolve(P_matrix, r_vals,prime)    

        # print(c)

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
        weights[u][v] = 0 if w == c else 1
    return solver(weights, n, m, b, t, c)


if __name__ == "__main__":
    main()
