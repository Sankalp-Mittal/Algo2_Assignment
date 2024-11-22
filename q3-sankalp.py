# Authors : Sankalp Mittal (395001), Ilia Badanin (350775), Vasco Frazão (396229)


import math
import random


def determinant(matrix, prime):
    """
    Gaussian elimination with improved handling of zero pivots and numerical stability.
    Returns the correct determinant.
    """
    # Get the size of the matrix
    n = len(matrix)

    # Create a copy of the matrix to avoid mutating the input
    mat = [row[:] for row in matrix]
    for i in range(n):
        for j in range(n):
            mat[i][j] = mat[i][j] % prime
            
    # Initialize determinant as 1
    det = 1

    for i in range(n):
        # Find pivot for column i
        pivot = i
        for j in range(i + 1, n):
            if mat[j][i] > mat[pivot][i]:
                pivot = j

        # If pivot is zero, determinant is zero
        if mat[pivot][i] == 0:  # Check for very small pivots
            return 0

        # Swap rows if needed
        if pivot != i:
            mat[i], mat[pivot] = mat[pivot], mat[i]
            det = -det % prime  # Swapping rows flips the sign of the determinant

        # Multiply determinant by the pivot element
        mat[i][i] = mat[i][i] % prime # Field
        det = det*mat[i][i] % prime

        # Normalize row i (avoid division by zero)
        if mat[i][i]:  # Ensure it's not a very small value
            inv = pow(mat[i][i], -1, prime)  # Modular inverse of the pivot
            for j in range(i + 1, n):
                mat[i][j] = (mat[i][j] * inv)%prime

        # Eliminate column i for rows below
        for j in range(i + 1, n):
            if mat[j][i]:  # Skip rows where the element is too small
                for k in range(i + 1, n):
                    mat[j][k] = (mat[j][k] - (mat[j][i] * mat[i][k])) % prime

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
            if A[j][i] > A[pivot][i]:
                pivot = j

        # Swap rows in A and b if pivot changes
        if pivot != i:
            A[i], A[pivot] = A[pivot], A[i]
            b[i], b[pivot] = b[pivot], b[i]

        # Ensure the pivot is non-zero
        if A[i][i] % prime == 0:
            raise ValueError("Matrix is singular or not invertible modulo prime.")

        # Eliminate entries below the pivot
        for j in range(i + 1, n):
            # Compute the factor to zero out A[j][i]
            factor = (A[j][i] * pow(A[i][i], -1, prime)) % prime
            for k in range(i, n):
                A[j][k] = (A[j][k] - factor * A[i][k]) % prime
            b[j] = (b[j] - factor * b[i]) % prime

    # Back substitution: Solve for x in Ux = c
    x = [0] * n
    for i in range(n - 1, -1, -1):
        x[i] = b[i]
        for j in range(i + 1, n):
            x[i] = (x[i] - A[i][j] * x[j]) % prime
        x[i] = (x[i] * pow(A[i][i], -1, prime)) % prime

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
            (alpha[i][j] * pow(gamma, weights[i][j], prime) * (weights[i][j] != -1)) % prime
            for j in range(n)
        ]
        for i in range(n)
    ]
    return determinant(H_matrix, prime)


def get_P(r_vals, gammas, prime):
    return [[pow(g, i, prime) for i in range(len(gammas))] for g in gammas]


def solver(weights, n, m, b, t, c):
    max_wt = 1
    min_siz = max(max_wt * n + 1, n * n)
    prime = next_prime(min_siz + 1)

    set_of_vals = list(range(prime))

    num_runs = math.ceil(math.log(n)) + 1
    for _ in range(num_runs):
        alpha = [[random.choice(set_of_vals[1:]) for _ in range(n)] for _ in range(n)]
        gammas = random.sample(set_of_vals[1:], n * max_wt + 1)

        r_vals = [hdet(gamma, alpha, weights, n, prime) for gamma in gammas]

        P_matrix = get_P(r_vals, gammas, prime)

        x = linsolve(P_matrix, r_vals,prime)

        for power in range(len(x)):
            if x[power]:
                if b == power*c + (n-power)*t:
                    print("yes")
                    return 0
        
    print("no")
    return 0


def main():
    n, m, b, t, c = map(int, input().split())
    weights = [[-1] * n for _ in range(n)]
    for _ in range(m):
        u, v, w = map(int, input().split())
        weights[u][v] = 1 if w == c else 0
    return solver(weights, n, m, b, t, c)


if __name__ == "__main__":
    main()
