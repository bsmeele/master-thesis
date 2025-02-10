from sympy import symbols, Matrix, zeros

def main():
    debug = True

    # Matrix dimensions
    M, N = 3, 3

    # Setting up symbols
    Rswl1, Rswl2, Rsbl1, Rsbl2, Rwl, Rbl = symbols('Rswl1 Rswl2 Rsbl1 Rsbl2 Rwl Rbl')
    R = [[symbols(f'R({i+1}_{j+1})') for j in range(N)] for i in range(M)]
    Vappwl1 = [symbols(f'Vappwl1({i+1})') for i in range(M)]
    Vappwl2 = [symbols(f'Vappwl2({i+1})') for i in range(M)]
    Vappbl1 = [symbols(f'Vappbl1({j+1})') for j in range(N)]
    Vappbl2 = [symbols(f'Vappbl2({j+1})') for j in range(N)]
    Vwl = [symbols(f'Vwl({i+1}_{j+1})') for i in range(M) for j in range(N)]
    Vbl = [symbols(f'Vbl({i+1}_{j+1})') for i in range(M) for j in range(N)]

    # V = Matrix.vstack(Matrix(Vwl), Matrix(Vbl))
    
    # Rswl1, Rswl2, Rsbl1, Rsbl2, Rwl, Rbl = 3, float("inf"), float("inf"), 5, 3, 2
    # R = [[10, 15, 20], [25, 30, 35], [40, 45, 50]]
    # R = [[i+1 + j+1 for j in range(N)] for i in range(M)]
    # Vappwl1 = [5, 7, 9]
    # Vappwl1 = [i+1 for i in range(M)]
    # Vappwl2 = [0 for i in range(M)]
    # Vappbl1 = [0 for j in range(N)]
    # Vappbl2 = [0 for j in range(N)]

    A = zeros(M*N, M*N)
    for i in range(M):
        for j in range(N):
            if j == 0:
                A[i*N + j, i*N + j] = 1/Rswl1 + 1/R[i][j] + 1/Rwl
                A[i*N + j, i*N + j+1] = -1/Rwl
            elif j == N-1:
                A[i*N + j, i*N + j] = 1/Rswl2 + 1/R[i][j] + 1/Rwl
                A[i*N + j, i*N + j-1] = -1/Rwl
            else:
                A[i*N + j, i*N+j] = 1/R[i][j] + 2/Rwl
                A[i*N + j, i*N + j-1] = -1/Rwl
                A[i*N + j, i*N + j+1] = -1/Rwl

    if debug:
        print("A:")
        for row in A.tolist():
            print(row)
        print()

    B = zeros(M*N, M*N)
    for i in range(M):
        for j in range(N):
            B[i*N + j, i*N + j] = -1/R[i][j]

    if debug:
        print("B:")
        for row in B.tolist():
            print(row)
        print()

    C = -B.copy()

    if debug:
        print("C:")
        for row in C.tolist():
            print(row)
        print()

    D = zeros(M*N, M*N)
    for i in range(M):
        for j in range(N):
            if i == 0:
                D[i*N + j, i*N + j] = -1/Rsbl1 + -1/R[i][j] + -1/Rbl
                D[i*N + j, (i+1)*N + j] = 1/Rbl
            elif i == M-1:
                D[i*N + j, i*N + j] = -1/Rsbl2 + -1/R[i][j] + -1/Rbl
                D[i*N + j, (i-1)*N + j] = 1/Rbl
            else:
                D[i*N + j, i*N + j] = -1/R[i][j] + -2/Rbl
                D[i*N + j, (i-1)*N + j] = 1/Rbl
                D[i*N + j, (i+1)*N + j] = 1/Rbl
    
    if debug:
        print("D:")
        for row in D.tolist():
            print(row)
        print()

    G_AB = Matrix.hstack(A, B)
    G_CD = Matrix.hstack(C, D)
    G_ABCD = Matrix.vstack(G_AB, G_CD)

    if debug:
        print("G_ABCD:")
        for row in G_ABCD.tolist():
            print(row)
        print()

    E = zeros(2*M*N, 1)
    for i in range(M):
        for j in range(N):
            if j == 0:
                E[i*N + j] = Vappwl1[i]/Rswl1
            elif j == N-1:
                E[i*N + j] = Vappwl2[i]/Rswl2
    for i in range(M):
        for j in range(N):
            if i == 0:
                E[i*N + j + M*N] = -Vappbl1[j]/Rsbl1
            elif i == M-1:
                E[i*N + j + M*N] = -Vappbl2[j]/Rsbl2
    
    if debug:
        print("E:")
        for row in E.tolist():
            print(row)
        print()

    # res = G_ABCD * V
    # for i in range(2*N*M):
    #     print(res[i], ' = ', E[i])
    # print()

    V_solved = G_ABCD.solve(E)

    if debug:
        print("V:")
        for row in V_solved:
            print(row)
        print()

    Iout = []
    for j in range(N):
        Ioutj = 0
        for i in range(M):
            Ioutj += (V_solved[i*N + j] - V_solved[i*N + j + M*N])/R[i][j]
        Iout.append(Ioutj)
    
    for i in Iout:
        print(i)
    print()


if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        print(e)
    
    input("Press Enter to exit...")

