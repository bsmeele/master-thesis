import matplotlib.pyplot as plt
import csv

def main():
    # file_path = 'out_files/out.txt'
    file_path = 'out.txt'
    # file_path = 'out_16x16_bottom_left.txt'
    # file_path = 'out_16x16_top_right.txt'
    # file_path = 'out_16x16_top_right_access_transistors.txt'
    file_path4 = 'out_files/all-signals-no-var.csv'

    t = []  # Time
    V = []  # Votage
    I = []  # Current
    N = []  # Oxygen vacancies (state variable)
    T = []  # Temperature
    V_schottky = []
    V_discplugserial = []
    R_schottky = []
    R_disc = []
    R_plug = []
    R_series = []
    R_total = []

    c_t_T = []
    c_T = []
    c_t_be = []
    c_be = []
    c_t_ion = []
    c_ion = []
    c_t_resistorR0 = []
    c_resistorR0 = []
    c_t_schottkytunnel = []
    c_schottkytunnel = []
    c_t_AE = []
    c_AE = []
    c_t_idt0 = []
    c_idt0 = []
    c_t_A = []
    c_A = []
    c_t_E_ion = []
    c_E_ion = []
    c_t_Ischottkytunnel = []
    c_Ischottkytunnel = []
    c_t_Nchange = []
    c_Nchange = []
    c_t_Nreal = []
    c_Nreal = []
    c_t_Rdisc = []
    c_Rdisc = []
    c_t_Rline = []
    c_Rline = []
    c_t_Rplug = []
    c_Rplug = []
    c_t_Rtheff = []
    c_Rtheff = []
    c_t_Treal = []
    c_Treal = []
    c_t_W0 = []
    c_W0 = []
    c_t_W00 = []
    c_W00 = []
    c_t_cvo = []
    c_cvo = []
    c_t_dWamax = []
    c_dWamax = []
    c_t_dWamin = []
    c_dWamin = []
    c_t_eps_eff = []
    c_eps_eff = []
    c_t_epsphib_eff = []
    c_epsphib_eff = []
    c_t_epsprime = []
    c_epsprime = []
    c_t_gamma = []
    c_gamma = []
    c_t_phiBn = []
    c_phiBn = []
    c_t_psi = []
    c_psi = []
    c_t_IAE = []
    c_IAE = []
    c_t_N_gnd_flow = []
    c_N_gnd_flow = []
    c_t_OE = []
    c_OE = []
    c_t_T_gnd_flow = []
    c_T_gnd_flow = []
    c_t_be_OE_flow = []
    c_be_OE_flow = []
    c_t_resistorR0_be_flow = []
    c_resistorR0_be_flow = []
    c_t_schottkytunnel_resistorR0_flow = []
    c_schottkytunnel_resistorR0_flow = []
    c_t_p = []
    c_p = []
    c_t_N = []
    c_N = []
    c_t_ion_gnd_flow = []
    c_ion_gnd_flow = []

    with open(file_path, 'r') as file:
        for line in file:
            numbers = line.split()
            if len(numbers) >= 12:
                try:
                    t.append(float(numbers[0]))
                    V.append(float(numbers[1]))
                    I.append(float(numbers[2]))
                    N.append(float(numbers[3]))
                    T.append(float(numbers[4]))
                    V_schottky.append(float(numbers[5]))
                    V_discplugserial.append(float(numbers[6]))
                    R_schottky.append(float(numbers[7]))
                    R_disc.append(float(numbers[8]))
                    R_plug.append(float(numbers[9]))
                    R_series.append(float(numbers[10]))
                    R_total.append(float(numbers[11]))

                except ValueError:
                    print(f"Skipping invalid line: {line.strip()}")

    with open(file_path4, mode='r', newline='') as file:
        reader = csv.reader(file)
        for row in reader:
            if len(row) >= 76:
                try:
                    c_t_T.append(float(row[0]))
                    c_T.append(float(row[1]))
                    c_t_be.append(float(row[2]))
                    c_be.append(float(row[3]))
                    c_t_ion.append(float(row[4]))
                    c_ion.append(float(row[5]))
                    c_t_resistorR0.append(float(row[6]))
                    c_resistorR0.append(float(row[7]))
                    c_t_schottkytunnel.append(float(row[8]))
                    c_schottkytunnel.append(float(row[9]))
                    c_t_AE.append(float(row[10]))
                    c_AE.append(float(row[11]))
                    c_t_idt0.append(float(row[12]))
                    c_idt0.append(float(row[13]))
                    c_t_A.append(float(row[14]))
                    c_A.append(float(row[15]))
                    c_t_E_ion.append(float(row[16]))
                    c_E_ion.append(float(row[17]))
                    c_t_Ischottkytunnel.append(float(row[18]))
                    c_Ischottkytunnel.append(float(row[19]))
                    c_t_Nchange.append(float(row[20]))
                    c_Nchange.append(float(row[21]))
                    c_t_Nreal.append(float(row[22]))
                    c_Nreal.append(float(row[23]))
                    c_t_Rdisc.append(float(row[24]))
                    c_Rdisc.append(float(row[25]))
                    c_t_Rline.append(float(row[26]))
                    c_Rline.append(float(row[27]))
                    c_t_Rplug.append(float(row[28]))
                    c_Rplug.append(float(row[29]))
                    c_t_Rtheff.append(float(row[30]))
                    c_Rtheff.append(float(row[31]))
                    c_t_Treal.append(float(row[32]))
                    c_Treal.append(float(row[33]))
                    c_t_W0.append(float(row[34]))
                    c_W0.append(float(row[35]))
                    c_t_W00.append(float(row[36]))
                    c_W00.append(float(row[37]))
                    c_t_cvo.append(float(row[38]))
                    c_cvo.append(float(row[39]))
                    c_t_dWamax.append(float(row[40]))
                    c_dWamax.append(float(row[41]))
                    c_t_dWamin.append(float(row[42]))
                    c_dWamin.append(float(row[43]))
                    c_t_eps_eff.append(float(row[44]))
                    c_eps_eff.append(float(row[45]))
                    c_t_epsphib_eff.append(float(row[46]))
                    c_epsphib_eff.append(float(row[47]))
                    c_t_epsprime.append(float(row[48]))
                    c_epsprime.append(float(row[49]))
                    c_t_gamma.append(float(row[50]))
                    c_gamma.append(float(row[51]))
                    c_t_phiBn.append(float(row[52]))
                    c_phiBn.append(float(row[53]))
                    c_t_psi.append(float(row[54]))
                    c_psi.append(float(row[55]))
                    c_t_IAE.append(float(row[56]))
                    c_IAE.append(float(row[57]))
                    c_t_N_gnd_flow.append(float(row[58]))
                    c_N_gnd_flow.append(float(row[59]))
                    c_t_OE.append(float(row[60]))
                    c_OE.append(float(row[61]))
                    c_t_T_gnd_flow.append(float(row[62]))
                    c_T_gnd_flow.append(float(row[63]))
                    c_t_be_OE_flow.append(float(row[64]))
                    c_be_OE_flow.append(float(row[65]))
                    c_t_resistorR0_be_flow.append(float(row[66]))
                    c_resistorR0_be_flow.append(float(row[67]))
                    c_t_schottkytunnel_resistorR0_flow.append(float(row[68]))
                    c_schottkytunnel_resistorR0_flow.append(float(row[69]))
                    c_t_p.append(float(row[70]))
                    c_p.append(float(row[71]))
                    c_t_N.append(float(row[72]))
                    c_N.append(float(row[73]))
                    c_t_ion_gnd_flow.append(float(row[74]))
                    c_ion_gnd_flow.append(float(row[75]))
                except ValueError:
                    print(f"Skipping invalid line: {line.strip()}")
    
    equal_t = True
    for i in range(0, len(c_t_T)):
        if (c_t_T[i] != c_t_be[i]):
            equal_t = False
            break
        if (c_t_T[i] != c_t_ion[i]):
            equal_t = False
            break
        if (c_t_T[i] != c_t_resistorR0[i]):
            equal_t = False
            break
        if (c_t_T[i] != c_t_schottkytunnel[i]):
            equal_t = False
            break
        if (c_t_T[i] != c_t_AE[i]):
            equal_t = False
            break
        if (c_t_T[i] != c_t_idt0[i]):
            equal_t = False
            break
        if (c_t_T[i] != c_t_A[i]):
            equal_t = False
            break
        if (c_t_T[i] != c_t_E_ion[i]):
            equal_t = False
            break
        if (c_t_T[i] != c_t_Ischottkytunnel[i]):
            equal_t = False
            break
        if (c_t_T[i] != c_t_Nchange[i]):
            equal_t = False
            break
        if (c_t_T[i] != c_t_Nchange[i]):
            equal_t = False
            break
        if (c_t_T[i] != c_t_Nreal[i]):
            equal_t = False
            break
        if (c_t_T[i] != c_t_Rdisc[i]):
            equal_t = False
            break
        if (c_t_T[i] != c_t_Rline[i]):
            equal_t = False
            break
        if (c_t_T[i] != c_t_Rplug[i]):
            equal_t = False
            break
        if (c_t_T[i] != c_t_Rtheff[i]):
            equal_t = False
            break
        if (c_t_T[i] != c_t_Treal[i]):
            equal_t = False
            break
        if (c_t_T[i] != c_t_W0[i]):
            equal_t = False
            break
        if (c_t_T[i] != c_t_W00[i]):
            equal_t = False
            break
        if (c_t_T[i] != c_t_cvo[i]):
            equal_t = False
            break
        if (c_t_T[i] != c_t_dWamax[i]):
            equal_t = False
            break
        if (c_t_T[i] != c_t_dWamin[i]):
            equal_t = False
            break
        if (c_t_T[i] != c_t_eps_eff[i]):
            equal_t = False
            break
        if (c_t_T[i] != c_t_epsphib_eff[i]):
            equal_t = False
            break
        if (c_t_T[i] != c_t_epsprime[i]):
            equal_t = False
            break
        if (c_t_T[i] != c_t_gamma[i]):
            equal_t = False
            break
        if (c_t_T[i] != c_t_phiBn[i]):
            equal_t = False
            break
        if (c_t_T[i] != c_t_psi[i]):
            equal_t = False
            break
        if (c_t_T[i] != c_t_IAE[i]):
            equal_t = False
            break
        if (c_t_T[i] != c_t_N_gnd_flow[i]):
            equal_t = False
            break
        if (c_t_T[i] != c_t_OE[i]):
            equal_t = False
            break
        if (c_t_T[i] != c_t_T_gnd_flow[i]):
            equal_t = False
            break
        if (c_t_T[i] != c_t_be_OE_flow[i]):
            equal_t = False
            break
        if (c_t_T[i] != c_t_resistorR0_be_flow[i]):
            equal_t = False
            break
        if (c_t_T[i] != c_t_schottkytunnel_resistorR0_flow[i]):
            equal_t = False
            break
        if (c_t_T[i] != c_t_p[i]):
            equal_t = False
            break
        if (c_t_T[i] != c_t_N[i]):
            equal_t = False
            break
        if (c_t_T[i] != c_t_ion_gnd_flow[i]):
            equal_t = False
            break

    if equal_t:
        print("Ts are equal")
    else:
        print("Ts are not equal")

    I_abs = [abs(e) for e in I]
    c_IAE_abs = [abs(e) for e in c_IAE]
    c_AE = [-e for e in c_AE]
    c_T = [e * 1e3 for e in c_T]
    c_Rseries = [e + 650 for e in c_Rline]
    c_Rschottky = [(v - (rd + rp + rs) * i) / i if i != 0 else 0 for v, rd, rp, rs, i in zip(c_AE, c_Rdisc, c_Rplug, c_Rseries, c_IAE)]

    i = 0
    j = 0
    t_dif = []
    I_dif = []
    while True:
        if i >= len(V) or j >= len(c_AE)-1:
            break
        if c_t_AE[j] <= t[i] and c_t_AE[j+1] >= t[i]:
            t_dif.append(t[i])
            I_dif.append(abs((c_IAE[j] + ((c_IAE[j+1] - c_IAE[j])/(c_t_AE[j+1] - c_t_AE[j])) * (t[i] - c_t_AE[j])) - I[i]))
            i += 1
        else:
            j += 1

    I_dif_reg = [(i_dif / i) if (i != 0) else 0 for i_dif, i in zip(I_dif, I)]

    plt.figure(figsize=(10, 7))

    # plt.subplot(2, 2, 1)
    # plt.plot(t, V, '-', markersize=2)
    # plt.title('Voltage waveform')
    # plt.xlabel('Time (s)')
    # plt.ylabel('Voltage (V)')
    # plt.grid(True)

    plt.subplot(2, 2, 1)
    plt.plot(V, I_abs, '-', markersize=2, label='C++', color='blue')
    plt.plot(c_AE, c_IAE_abs, markersize=2, label='Cadence', color='orange')
    plt.yscale('log')
    # plt.xlim(min(V), max(V))
    plt.xlim(-1.5, 1.5)
    plt.ylim(1e-8, max(I)*2)
    plt.title('I-V Characteristic')
    plt.xlabel('Voltage (V)')
    plt.ylabel('Current (A)')
    plt.legend()
    plt.grid(True)

    plt.subplot(2, 2, 2)
    plt.plot(t, N, '-', markersize=2, label='C++', color='blue')
    plt.plot(c_t_Nreal, c_Nreal, '-', markersize=2, label='Cadence', color='orange')
    plt.yscale('log')
    plt.xlim(0, 6)
    plt.title('Oxygen vacancy concentration of the disc')
    plt.xlabel('Time (s)')
    plt.ylabel('Disc Concentration')
    plt.legend()
    plt.grid(True)

    plt.subplot(2, 2, 3)
    plt.plot(t, T, '-', markersize=2, label='C++', color='blue')
    plt.plot(c_t_T, c_T, '-', markersize=2, label='Cadence', color='orange')
    plt.xlim(0, 6)
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
    # plt.plot(t, R_total, '-', label="Total", color='green')
    # plt.yscale('log')
    # plt.xlim(0, 6)
    # plt.ylim(0, 1e6)
    # plt.title('Resistance components')
    # plt.xlabel('Time (s)')
    # plt.ylabel('Resistance (Ohm)')
    # plt.legend()
    # plt.grid(True)

    plt.subplot(2, 2, 4)
    plt.plot(t_dif, I_dif, '-', label='absolute', color='blue')
    plt.plot(t_dif, I_dif_reg, '-', label='relative', color='orange')
    plt.yscale('log')
    plt.xlim(0, 6)
    # plt.ylim(1e-12, 1e-3)
    plt.ylim(1e-12, 1e2)
    plt.title('Current difference between Cadence and C++ simulation')
    plt.xlabel('Time (s)')
    plt.ylabel('Current difference (A)')
    plt.legend()
    plt.grid(True)

    # plt.subplot(2, 2, 4)
    # plt.plot(t_dif, I_dif_reg, '-')
    # plt.yscale('log')
    # plt.xlim(0, 6)
    # # plt.ylim(1e-12, 1e-3)
    # plt.title('Regularized current difference between Cadence and C++ simulation')
    # plt.xlabel('Time (s)')
    # plt.ylabel('Current difference (A)')
    # plt.grid(True)

    # plt.subplot(2, 2, 4)
    # # plt.plot(t, R_disc, '-', label="c++", color='blue')
    # # plt.plot(c_t_Rdisc, c_Rdisc, '-', label="Cadence", color='orange')
    # # plt.plot(t, R_plug, '-', label="c++", color='blue')
    # # plt.plot(c_t_Rplug, c_Rplug, '-', label="Cadence", color='orange')
    # # plt.plot(t, R_series, '-', label="c++", color='blue')
    # # plt.plot(c_t_Rline, c_Rseries, '-', label="Cadence", color='orange')
    # plt.plot(t, Rtheff, '-', label="c++", color='blue')
    # plt.plot(c_t_Rtheff, c_Rtheff, '-', label="Cadence", color='orange')
    # plt.yscale('log')
    # plt.xlim(0, 6)
    # # plt.ylim(0, 1e6)
    # plt.title('Resistance components')
    # plt.xlabel('Time (s)')
    # plt.ylabel('Resistance (Ohm)')
    # plt.legend()
    # plt.grid(True)

    # plt.subplot(2, 2, 4)
    # plt.plot(t, V, '-', label='Applied', color='blue')
    # plt.plot(t, V_schottky, '-', label='Schottky', color='orange')
    # plt.plot(t, V_discplugserial, '-', label='DiscPlugSerial', color='purple')
    # plt.xlim(0, 6)
    # plt.title('Voltage')
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
