import matplotlib.pyplot as plt

def main():
    file_path = 'out.txt'
    file_path2 = 'out2.txt'

    t = []  # Time
    V = []  # Votage
    I = []  # Current
    N = []  # Oxygen vacancies (state variable)
    T = []  # Temperature
    R_schottky = []
    R_disc = []
    R_plug = []
    R_series = []
    V_schottky = []
    V_discplugserial = []
    phibn = []
    V_solve_bottom = []
    V_solve_top = []
    R = []
    V2 = []
    R0 = []
    R1 = []
    R2 = []
    R3 = []
    R4 = []
    R5 = []
    R6 = []
    R7= []
    R8 = []
    R9 = []
    R10 = []
    R11 = []
    R12 = []
    R13 = []
    R14 = []
    R15 = []
    R16 = []
    R17 = []
    R18 = []
    R19 = []
    R20 = []

    with open(file_path, 'r') as file:
        for line in file:
            numbers = line.split()
            if len(numbers) >= 15:
                try:
                    t.append(float(numbers[0]))
                    V.append(float(numbers[1]))
                    I.append(abs(float(numbers[2])))
                    N.append(float(numbers[3]))
                    T.append(float(numbers[4]))
                    R_schottky.append(float(numbers[5]))
                    R_disc.append(float(numbers[6]))
                    R_plug.append(float(numbers[7]))
                    R_series.append(float(numbers[8]))
                    V_schottky.append(float(numbers[9]))
                    V_discplugserial.append(float(numbers[10]))
                    phibn.append(float(numbers[11]))
                    V_solve_bottom.append(float(numbers[12]))
                    V_solve_top.append(float(numbers[13]))
                    R.append(float(numbers[14]))

                except ValueError:
                    print(f"Skipping invalid line: {line.strip()}")

    with open(file_path2, 'r') as file:
        for line in file:
            numbers = line.split()
            if len(numbers) >= 22:
                try:
                    V2.append(float(numbers[0]))
                    R0.append(float(numbers[1]))
                    R1.append(float(numbers[2]))
                    R2.append(float(numbers[3]))
                    R3.append(float(numbers[4]))
                    R4.append(float(numbers[5]))
                    R5.append(float(numbers[6]))
                    R6.append(float(numbers[7]))
                    R7.append(float(numbers[8]))
                    R8.append(float(numbers[9]))
                    R9.append(float(numbers[10]))
                    R10.append(float(numbers[11]))
                    R11.append(float(numbers[12]))
                    R12.append(float(numbers[13]))
                    R13.append(float(numbers[14]))
                    R14.append(float(numbers[15]))
                    R15.append(float(numbers[16]))
                    R16.append(float(numbers[17]))
                    R17.append(float(numbers[18]))
                    R18.append(float(numbers[19]))
                    R19.append(float(numbers[20]))
                    R20.append(float(numbers[21]))

                except ValueError:
                    print(f"Skipping invalid line: {line.strip()}")

    plt.figure(figsize=(10, 7))

    # plt.subplot(2, 2, 1)
    # plt.plot(t, V, '-', markersize=2)
    # plt.title('Voltage waveform')
    # plt.xlabel('Time (s)')
    # plt.ylabel('Voltage (V)')
    # plt.grid(True)

    plt.subplot(2, 2, 1)
    plt.plot(V, I, '-', markersize=2)
    plt.yscale('log')
    plt.xlim(min(V), max(V))
    plt.ylim(1e-8, max(I)*2)
    plt.title('I-V Characteristic')
    plt.xlabel('Voltage (V)')
    plt.ylabel('Current (A)')
    plt.grid(True)

    plt.subplot(2, 2, 2)
    plt.plot(t, N, '-', markersize=2)
    plt.yscale('log')
    plt.xlim(0, 6)
    plt.title('Oxygen vacancy concentration of the disc')
    plt.xlabel('Time (s)')
    plt.ylabel('Disc Concentration')
    plt.legend()
    plt.grid(True)

    plt.subplot(2, 2, 3)
    plt.plot(t, T, '-', markersize=2)
    plt.xlim(0, 6)
    # plt.ylim(0, 2000)
    plt.title('Temperature')
    plt.xlabel('Time (s)')
    plt.ylabel('Temperature (K)')
    plt.legend()
    plt.grid(True)

    # plt.subplot(2, 2, 4)
    # plt.plot(t, R_schottky, '-', label="Schottky", color='blue')
    # plt.plot(t, R_disc, '-', label="Disc", color='orange')
    # plt.plot(t, R_plug, '-', label="Plug", color='yellow')
    # plt.plot(t, R_series, '-', label="Series", color='purple')
    # plt.plot(t, R, '-', label="Total", color='green')
    # plt.yscale('log')
    # plt.xlim(0, 6)
    # plt.ylim(0, 1e6)
    # plt.title('Resistance components')
    # plt.xlabel('Time (s)')
    # plt.ylabel('Resistance (Ohm)')
    # plt.legend()
    # plt.grid(True)

    plt.subplot(2, 2, 4)
    plt.plot(V2, R0, '-', label="N=0.008", color='blue')
    plt.plot(V2, R5, '-', label="N=5", color='orange')
    plt.plot(V2, R10, '-', label="N=10", color='yellow')
    plt.plot(V2, R15, '-', label="N=15", color='purple')
    plt.plot(V2, R20, '-', label="N=20", color='green')
    # plt.yscale('log')
    plt.ylim(1500, 2000)
    plt.title('Resistance per oxygen vacancy concentration')
    plt.xlabel('Votlage (V)')
    plt.ylabel('Resistance (Ohm)')
    plt.legend()
    plt.grid(True)

    # plt.subplot(2, 2, 3)
    # plt.plot(t, V, '-', label='Applied', color='blue')
    # plt.plot(t, V_schottky, '-', label='Schottky', color='orange')
    # plt.plot(t, V_discplugserial, '-', label='DiscPlugSerial', color='purple')
    # plt.xlim(0, 6)
    # plt.title('Voltage')
    # plt.xlabel('Time (s)')
    # plt.ylabel('Voltage (V)')
    # plt.legend()
    # plt.grid(True)

    # plt.subplot(2, 2, 4)
    # plt.plot(t, phibn, '-', label='phibn', color='blue')
    # plt.plot(t, V_schottky, '-', label='Vschottky', color='orange')
    # plt.xlim(0, 6)
    # plt.ylim(-0.2, 0.2)
    # plt.title('phibn')
    # plt.xlabel('Time (s)')
    # plt.ylabel('phibn')
    # plt.legend()
    # plt.grid(True)

    # plt.subplot(2, 2, 4)
    # plt.plot(t, V_solve_bottom, '-', label='Solve bottom', color='blue')
    # plt.plot(t, V_solve_top, '-', label='Solve top', color='orange')
    # plt.title('V solve')
    # plt.xlabel('Time (s)')
    # plt.ylabel('Voltage (V)')
    # plt.legend()
    # plt.grid(True)


    # Adjust layout and display the plots
    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        print(e)
    
    input("Press Enter to exit...")
