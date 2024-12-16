import matplotlib.pyplot as plt

def main():
    file_path = 'out.txt'

    t = []  # Time
    V = []  # Votage
    I = []  # Current
    N = []  # Oxygen vacancies (state variable)
    T = []  # Temperature
    R_schottky = []
    R_disc = []
    R_plug = []
    R_series = []

    with open(file_path, 'r') as file:
        for line in file:
            numbers = line.split()
            if len(numbers) >= 9:
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
    plt.title('Temperature')
    plt.xlabel('Time (s)')
    plt.ylabel('Temperature (K)')
    plt.legend()
    plt.grid(True)

    plt.subplot(2, 2, 4)
    plt.plot(t, R_schottky, '-', label="Schottky", color='blue')
    plt.plot(t, R_disc, '-', label="Disc", color='orange')
    plt.plot(t, R_plug, '-', label="Plug", color='yellow')
    plt.plot(t, R_series, '-', label="Series", color='purple')
    plt.yscale('log')
    plt.xlim(0, 6)
    plt.ylim(0, 1e6)
    plt.title('Resistance components')
    plt.xlabel('Time (s)')
    plt.ylabel('Resistance (Ohm)')
    plt.legend()
    plt.grid(True)

    # Adjust layout and display the plots
    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        print(e)
    
    input("Press Enter to exit...")
