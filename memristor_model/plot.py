import matplotlib.pyplot as plt

def main():
    file_path = 'out.txt'

    t = []
    V = []
    I = []
    N = []

    with open(file_path, 'r') as file:
        for line in file:
            numbers = line.split()
            if len(numbers) >= 4:
                try:
                    t.append(float(numbers[0]))
                    V.append(float(numbers[1]))
                    I.append(abs(float(numbers[2])))
                    N.append(float(numbers[3]))
                except ValueError:
                    print(f"Skipping invalid line: {line.strip()}")

    # Create the plots
    plt.figure(figsize=(7, 7))

    # Plot 1: First and second numbers
    plt.subplot(2, 2, 1)
    plt.plot(t, V, '-', markersize=2)
    plt.title('Voltage waveform')
    plt.xlabel('Time (s)')
    plt.ylabel('Voltage (V)')
    plt.grid(True)

    # Plot 2: Second and third numbers
    plt.subplot(2, 2, 2)
    plt.plot(V, I, '-', markersize=2)
    plt.yscale('log')
    plt.xlim(min(V), max(V))
    plt.ylim(1e-8, max(I)*2)
    plt.title('I-V Characteristic')
    plt.xlabel('Voltage (V)')
    plt.ylabel('Current (A)')
    plt.grid(True)

    # Plot 3: First and last numbers
    plt.subplot(2, 2, 3)
    plt.plot(t, N, '-', markersize=2)
    plt.yscale('log')
    plt.title('Oxygen vacancy concentration of the disc')
    plt.xlabel('Time (s)')
    plt.ylabel('Disc Concentration')
    plt.legend()
    plt.grid(True)

    # Leave the bottom-right plot empty
    plt.subplot(2, 2, 4)
    plt.axis('off')  # Turn off the axis for the empty plot

    # Adjust layout and display the plots
    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        print(e)
    
    input("Press Enter to exit...")
