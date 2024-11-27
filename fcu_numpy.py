import numpy as np
import time

def main():
    # Matrix dimensions
    M, N = 32, 32

    print_flag = False

    total_time = 0

    for i in range(100):
        Rmin = 100.0
        Rmax = 1000.0
        Vdd = 5.0
        
        Rcol = 2.0
        Rrow = 3.0
        Rsense = 5.0

        
        if print_flag:
            print(f"Rcol: {Rcol}")
            print(f"Rrow: {Rrow}")
            print(f"Rsense: {Rsense}\n")

        R = np.random.rand(M, N)
        R = (R + 1.0) * 0.5 * (Rmax - Rmin) + Rmin
        R[0, 0] = 10.0
        R[0, 1] = 15.0
        R[0, 2] = 20.0
        R[1, 0] = 25.0
        R[1, 1] = 30.0
        R[1, 2] = 35.0
        R[2, 0] = 40.0
        R[2, 1] = 45.0
        R[2, 2] = 50.0

        if print_flag:
            print("R:")
            print(R, "\n")

        G = np.divide(1, R)

        if print_flag:
            print("G:")
            print(G, "\n")

        Vin = np.random.rand(M)
        Vin = np.where(Vin > 0.5, Vdd, 0)  # Set values > 0.5 to Vdd, others to 0
        Vin[0] = 5.0
        Vin[1] = 7.0
        Vin[2] = 9.0

        if print_flag:
            print("Vin:")
            print(Vin)

        start_time = time.time()

        # Step 1: Formulate Column Linear Systems
        A = []
        for j in range(N):
            Aj = []
            for row in range(M):
                row_elem = []
                for i in range(M):
                    if i == row:
                        row_elem.append(-1)
                    elif i > row:
                        row_elem.append((i - row) * G[i, j] * Rcol)
                    else:
                        row_elem.append(0)
                Aj.append(row_elem)
            A.append(np.array(Aj, dtype=float))

        if print_flag:
            for j in range(N):
                print("A", j, ":")
                print(A[j])
                print()
            print()

        J = np.array([[(Rsense + i * Rcol)] for i in range(M)][::-1], dtype=float)

        if print_flag:
            print("J:")
            print(J)
            print()

        K = [G[:, j].reshape(1, -1) for j in range(N)]

        if print_flag:
            print("K:")
            for j in range(N):
                print("K", j, ":")
                print(K[j])
                print()
            print()

        # Step 2: Merge Column Linear Systems
        COLmat_matrices = [(-A[j] + J @ K[j]) for j in range(N)]
        COLmat = direct_sum(*COLmat_matrices)

        if print_flag:
            print("COLmat:")
            print(COLmat)
            print()

        Gmat = direct_sum(*K)

        if print_flag:
            print("Gmat:")
            print(Gmat)
            print()

        # Step 3: Formulate Row Linear Systems
        B = []
        for i in range(M):
            Bi = []
            for row in range(N):
                row_elem = []
                for j in range(N):
                    if row <= j:
                        row_elem.append(G[i, j] * (row + 1))
                    else:
                        row_elem.append(G[i, j] * (j + 1))
                Bi.append(row_elem)
            B.append(np.array(Bi, dtype=float) * Rrow)

        if print_flag:
            for i in range(M):
                print("B", i, ":")
                print(B[i])
                print()
            print()

        Vrowin = [np.array([[Vin[i]] for j in range(N)], dtype=float) for i in range(M)]

        # Step 4: Merge Row Linear Systems
        ROWmat = direct_sum(*B)

        if print_flag:
            print("ROWmat:")
            print(ROWmat)
            print()

        CVrowin = np.vstack(Vrowin)

        # Step 5: Eliminate Internal Variables
        ROWmatAint = np.zeros((N * M, N * M), dtype=float)
        CVrowinA = np.zeros((N * M, 1), dtype=float)
        for j in range(N):
            for i in range(M):
                ROWmatAint[j * M + i, :] = ROWmat[i * N + j, :]
                CVrowinA[j * M + i, :] = CVrowin[i * N + j, :]

        ROWmatA = np.zeros((N * M, N * M), dtype=float)
        for j in range(N):
            for i in range(M):
                ROWmatA[:, j * M + i] = ROWmatAint[:, i * N + j]

        if print_flag:
            print("ROWmatA:")
            print(ROWmatA)
            print()

        if print_flag:
            print("CVrowINA:")
            print(CVrowinA)
            print()

        if print_flag:
            tmp = COLmat + ROWmatA
            print("COLmat + ROWmatA:")
            print(tmp)
            print()

        if print_flag:
            print("(COLmat + ROWmatA)^-1:")
            print(np.linalg.inv(tmp))
            print()

        NETmat = Gmat @ np.linalg.inv(COLmat + ROWmatA)

        if print_flag:
            print("NETmat:")
            print(NETmat)
            print()

        # Step 6: Reduce Matrix Dimension
        NETmatC = np.zeros((N, M), dtype=float)
        for i in range(M):
            for j in range(N):
                NETmatC[:, i] += NETmat[:, j * M + i]

        if print_flag:
            print("NETmatC:")
            print(NETmatC)
            print()

        Iout = NETmatC @ Vin

        end_time = time.time()

        print(f"Execution time: {(end_time - start_time) * 1000} (ms)")
        total_time += end_time - start_time

        if print_flag:
            print("Iout:")
            print(Iout)
            print()
    
    print(f"Average time: {total_time*1000/100} (ms)")

def direct_sum(*matrices):
    # Calculate the total size of the resulting matrix
    total_rows = sum(m.shape[0] for m in matrices)
    total_cols = sum(m.shape[1] for m in matrices)
    
    # Create an empty matrix of the final size, filled with zeros
    result = np.zeros((total_rows, total_cols), dtype=float)
    
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
