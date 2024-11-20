# Authors : Sankalp Mittal (395001), Ilia Badanin (350775), Vasco Frazão (396229)


import math
import random


def determinant(matrix):
    """
    37.3ms vs 1.26ms (np.linalg.det) for 100x100
    Gaussian elimination O(n^3), almost the same as LU-decomposition by time
    Returns correct determinant
    might have an overflow, but it happens only for 1000x1000 matrices
    """
    # Get the size of the matrix
    n = len(matrix)

    # Create a copy of the matrix to avoid mutating the input
    mat = [row[:] for row in matrix]

    # Initialize determinant as 1
    det = 1

    for i in range(n):
        # Find pivot for column i
        pivot = i
        for j in range(i + 1, n):
            if abs(mat[j][i]) > abs(mat[pivot][i]):
                pivot = j

        # If pivot is zero, determinant is zero
        if mat[pivot][i] == 0:
            return 0

        # Swap rows if needed
        if pivot != i:
            mat[i], mat[pivot] = mat[pivot], mat[i]
            det *= -1  # Swapping rows flips the sign of the determinant

        # Multiply determinant by the pivot element
        det *= mat[i][i]

        # Normalize row i
        for j in range(i + 1, n):
            mat[i][j] /= mat[i][i]

        # Eliminate column i for rows below
        for j in range(i + 1, n):
            for k in range(i + 1, n):
                mat[j][k] -= mat[j][i] * mat[i][k]

    return det


def solve(A, b):
    """
    Solves the system of linear equations A * x = b using Gaussian elimination.

    Parameters:
    A (list of lists): Coefficient matrix (n x n).
    b (list): Right-hand side vector (n).

    Returns:
    list: Solution vector x (n).
    """
    n = len(A)

    # Forward elimination: Reduce to upper triangular form
    for i in range(n):
        # Find the pivot
        pivot = i
        for j in range(i + 1, n):
            if abs(A[j][i]) > abs(A[pivot][i]):
                pivot = j

        # Swap rows in A and b
        if pivot != i:
            A[i], A[pivot] = A[pivot], A[i]
            b[i], b[pivot] = b[pivot], b[i]

        # Check for singular matrix
        if A[i][i] == 0:
            raise ValueError("Matrix is singular or nearly singular.")

        # Eliminate entries below the pivot
        for j in range(i + 1, n):
            factor = A[j][i] / A[i][i]
            for k in range(i, n):
                A[j][k] -= factor * A[i][k]
            b[j] -= factor * b[i]

    # Back substitution: Solve for x in Ux = c
    x = [0] * n
    for i in range(n - 1, -1, -1):
        x[i] = b[i]
        for j in range(i + 1, n):
            x[i] -= A[i][j] * x[j]
        x[i] /= A[i][i]

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


def hdet(gamma, alpha, weights, n):
    H_matrix = [
        [
            alpha[i][j] * (gamma ** weights[i][j]) * (weights[i][j] != -1)
            for j in range(n)
        ]
        for i in range(n)
    ]
    return determinant(H_matrix)


def get_P(r_vals, gammas, prime):
    return [[pow(g, i, prime) for i in range(len(gammas))] for g in gammas]


def solver(weights, n, m, b, t, c):
    max_wt = max(t, c)
    min_siz = max(max_wt * n + 1, n * n)
    prime = next_prime(min_siz + 1)

    set_of_vals = list(range(prime))

    num_runs = math.ceil(math.log(n))
    for _ in range(num_runs):
        alpha = [[random.choice(set_of_vals[1:]) for _ in range(n)] for _ in range(n)]
        gammas = random.sample(set_of_vals[1:], n * max_wt + 1)

        r_vals = [hdet(gamma, alpha, weights, n) % prime for gamma in gammas]

        P_matrix = get_P(r_vals, gammas, prime)

        c = solve(P_matrix, r_vals)

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
    return solver(weights, n, m, b, t, c)


if __name__ == "__main__":
    main()
