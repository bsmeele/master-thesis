from sympy import symbols, Matrix, zeros

def main():
    # Matrix dimensions
    M, N = 3, 3

    # Setting up symbols
    Rsense, Rcol, Rrow = symbols('Rsense Rcol Rrow')
    G = [[symbols(f'g{i+1}{j+1}') for j in range(N)] for i in range(M)]
    V = [[symbols(f'V{i+1}{j+1}') for j in range(N)] for i in range(M)]
    Va = [[symbols(f'Va{i+1}{j+1}') for j in range(N)] for i in range(M)]
    Vin = [symbols(f'Vin{i+1}') for i in range(M)]

    # M, N = 3, 3
    # Rsense, Rcol, Rrow = 5, 2, 3
    # G = [[1/10, 1/15, 1/20], [1/25, 1/30, 1/35], [1/40, 1/45, 1/50]]
    # # G = [[0, 0, 0], [0, 0, 0], [1/10, 0, 0]]
    # V = [[((i+1)*(j+1)) for j in range(N)] for i in range(M)]
    # Va = [[((i+1)*(j+1))  for j in range(N)] for i in range(M)]
    # Vin = [5, 7, 9]

    # M, N = 2, 2
    # Rsense, Rcol, Rrow = 5, 2, 3
    # # G = [[1/10, 1/15, 1/20], [1/25, 1/30, 1/35], [1/40, 1/45, 1/50]]
    # G = [[1/11, 1/17], [1/29, 1/37]]
    # Vin = [7, 5]

    print("G:")
    for row in G:
        print(row)
    print()

    print("Vin:")
    for v in Vin:
        print(v)
    print()
    
    # ----- Step 1: Formulate Column Linear Systems -----
    A = []
    for j in range(N):
        Aj = []
        for row in range(M):
            row_elem = []
            for i in range(M):
                if i == row:
                    row_elem.append(-1)
                elif i > row:
                    row_elem.append((i - row) * G[i][j] * Rcol)
                else:
                    row_elem.append(0)
            Aj.append(row_elem)
        A.append(Matrix(Aj))
    
    for j in range(N):
        print("A", j, ":")
        for row in A[j].tolist():
            print(row)
        print()
    print()

    Vcol = []
    for j in range(N):
        Vcolj = []
        for i in range(M):
            Vcolj.append([V[i][j]])
        Vcol.append(Matrix(Vcolj))
    
    VAcol = []
    for j in range(N):
        VAcolj = []
        for i in range(M):
            VAcolj.append([Va[i][j]])
        VAcol.append(Matrix(VAcolj))
    
    J = Matrix([[(Rsense + i * Rcol)] for i in range(M).__reversed__()])

    print("J:")
    for j in J:
        print(j)
    print()
    
    K = []
    for j in range(N):
        K.append(Matrix([G[i][j] for i in range(M)]).T)

    print("K:")
    for j in range(N):
        print("K", j, ":")
        print(K[j])
        print()
    print()
    
    Iout = []
    for j in range(N):
        Ioutj = 0
        for i in range(M):
            Ioutj += G[i][j] * V[i][j]
        Iout.append(Ioutj)
    
    # -Aj * Vcolj = VAcolj - J * Ioutj
    for j in range(N):
        print(f"Column {j}:")
        for i in range(M):
            print((-A[j] * Vcol[j])[i], " = ", (VAcol[j] - J * Iout[j])[i])
        print()
    
    # (-Aj + J * Kj) * Vcolj = VAcolj
    for j in range(N):
        print(f"Column {j}:")
        for i in range(M):
            print(((-A[j] + J * K[j])*Vcol[j])[i], " = ", VAcol[j][i])
        print()

    for j in range(N):
        tmp = J * K[j]
        for row in tmp.tolist():
            print(row)
        print()
    
    # ----- Step 2: Merge Column Linear Systems -----
    COLmat_matrices = [(-A[j] + J * K[j]) for j in range(N)]
    COLmat = direct_sum(*COLmat_matrices)

    print("COLmat:")
    for row in COLmat.tolist():
        print(row)
    print()
    
    Gmat = direct_sum(*K)

    print("Gmat:")
    for row in Gmat.tolist():
        print(row)
    print()
    
    CVcol = Matrix.vstack(*Vcol)
    CVAcol = Matrix.vstack(*VAcol)

    # COLmat * CVcol = CVAcol
    for i in range(M*M):
        print((COLmat * CVcol)[i], " = ", CVAcol[i])
    print()

    # Gmat * CVcol = (Iout)^T
    for i in range(M):
        print((Gmat * CVcol)[i], " = ", Iout[i])
    print()

    # ----- Step 3: Formulate Row Linear Systems -----
    B = []
    for i in range(M):
        Bi = []
        for row in range(N):
            row_elem = []
            for j in range(N):
                if (row <= j):
                    row_elem.append(G[i][j] * (row + 1))
                else:
                    row_elem.append(G[i][j] * (j + 1))
            Bi.append(row_elem)
        B.append(Matrix(Bi) * Rrow)

    for i in range(M):
        print("B", i, ":")
        for row in B[i].tolist():
            print(row)
        print()
    print()

    Vrow = []
    for i in range(M):
        Vrowi = []
        for j in range(N):
            Vrowi.append([V[i][j]])
        Vrow.append(Matrix(Vrowi))

    VArow = []
    for i in range(M):
        VArowi = []
        for j in range(N):
            VArowi.append([Va[i][j]])
        VArow.append(Matrix(VArowi))

    Vrowin = [Matrix([Vin[i] for j in range(N)]) for i in range(M)]

    # Bi * Vrowi = Vrowini - VArowi
    for i in range(M):
        print(f"Row {i}:")
        for j in range(N):
            print((B[i] * Vrow[i])[j], " = ", (Vrowin[i] - VArow[i])[j])
        print()

    # ----- Step 4: Merge Row Linear Systems -----
    ROWmat = direct_sum(*B)

    print("ROWmat:")
    for row in ROWmat.tolist():
        print(row)
    print()

    CVrow = Matrix.vstack(*Vrow)
    CVArow = Matrix.vstack(*VArow)
    CVrowin = Matrix.vstack(*Vrowin)

    # CVrowin - ROWmat * CVrow = CVArow
    for j in range(N*N):
        print((CVrowin - ROWmat * CVrow)[j], " = ", CVArow[j])
    print()

    # ----- Step 5: Eliminate Internal Variables -----
    ROWmatAint = zeros(N*M, N*M)
    CVrowinA = zeros(N*M, 1)
    for j in range(N):
        for i in range(M):
            ROWmatAint[j*M + i, :] = ROWmat[i*N + j, :]
            CVrowinA[j*M + i, :] = CVrowin[i*N + j, :]
    ROWmatA = zeros(N*M, N*M)
    for j in range(N):
        for i in range(M):
            ROWmatA[:, j*M + i] = ROWmatAint[:, i*N + j]

    print("ROWmatA:")
    for row in ROWmatA.tolist():
        print(row)
    print()

    print("CVrowINA:")
    for v in CVrowinA:
        print(v)
    print()

    # CVrowINA - ROWmatA * CVcol = CVAcol
    for i in range(N*M):
        print((CVrowinA - ROWmatA * CVcol)[i], " = ", CVAcol[i])
    print()

    # tmp = COLmat + ROWmatA
    # tmp_inv = tmp.inv()

    # with open("tmp_inv.txt", "w") as f:
    #     for i in range(N*N):
    #         for j in range(M*M):
    #             print(f"({i} {j}):", file=f)
    #             print(tmp_inv[i, j], file=f)
    #             print(file=f)

    # for i in range(N*N):
    #     for j in range(M*M):
    #         print(f"({i} {j}):")
    #         print(tmp_inv[i, j])
    #         print()

    tmp = COLmat + ROWmatA
    print("COLmat + ROWmatA:")
    for row in tmp.tolist():
        print(row)
    print()

    print("(COLmat + ROWmatA)^-1:")
    for row in tmp.inv().tolist():
        print(row)
    print()

    NETmat = Gmat * (COLmat + ROWmatA).inv()

    print("NETmat:")
    for row in NETmat.tolist():
        print(row)
    print()
    
    # CVrowina - ROWmata * CVcol = CVAcol
    # for i in range(N*N):
    #     print((CVrowina - ROWmata * CVcol)[i], " = ", CVAcol[i])

    # (COLmat + ROWmatA) * CVcol = CVrowina
    # for i in range(N*N):
    #     print(((COLmat + ROWmata) * CVcol)[i], " = ", CVrowina[i])

    # (Iout)^T = NETmat * CVrowina
    # for i in range(M):
    #     print(Iout[i], " = ", (NETmat * CVrowina)[i])

    # ----- Step 6: Reduce Matrix Dimension -----
    NETmatC = zeros(N, M)
    for i in range(M):
        for j in range(N):
            NETmatC[:, i] += NETmat[:, j*M + i]

    print("NETmatC:")
    for row in NETmatC.tolist():
        print(row)
    print()

    Iout = NETmatC * Matrix(Vin)

    print("Iout:")
    for i in Iout:
        print(i.evalf())
    print()


def direct_sum(*matrices):
    # Calculate the total size of the resulting matrix
    total_rows = sum(m.rows for m in matrices)
    total_cols = sum(m.cols for m in matrices)
    
    # Create an empty matrix of the final size, filled with zeros
    result = zeros(total_rows, total_cols)
    
    # Position each matrix in its block-diagonal location
    row_offset = 0
    col_offset = 0
    for m in matrices:
        rows, cols = m.shape
        result[row_offset:row_offset + rows, col_offset:col_offset + cols] = m
        row_offset += rows
        col_offset += cols
    
    return result

if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        print(e)
    
    input("Press Enter to exit...")
