import os
import re

import numpy as np
import matplotlib.pyplot as plt

def save_bin(path, array):
    shape = np.array(array.shape, dtype=np.float32)  # Store shape as int64
    # print(shape)
    with open(path, "wb") as f:
        # shape = array.shape
        shape.tofile(f)  # Save shape first
        array.tofile(f)  # Save data


def load_bin(path, dtype=np.float64, shape_len=2):
    with open(path, "rb") as f:
        # Read shape first
        shape = np.fromfile(f, dtype=np.int64, count=shape_len).astype(int)
        # shape = np.fromfile(f, dtype=np.float32, count=shape_len).astype(int)
       
        # Read data
        data = np.fromfile(f, dtype=dtype)
        shape = tuple(shape)

        # Reshape
        array = data.reshape(shape)
       
    return array


def main():
    # base_folder = "data/bin/batch_1"
    # dir_pattern = r"row_\d+-\d+_col_\d+-\d+"
    base_folder = "RRAM_validation_data/64x64/save_c_arrays"
    dir_pattern = r"batch_\d+"

    crossbar_size = 64

    mac_err = []
    mac_err_reg = []
    mem_err = np.zeros((crossbar_size, crossbar_size))
    mem_err_reg = np.zeros((crossbar_size, crossbar_size))
    lrs_err = np.zeros((crossbar_size, crossbar_size))
    lrs_err_reg = np.zeros((crossbar_size, crossbar_size))
    hrs_err = np.zeros((crossbar_size, crossbar_size))
    hrs_err_reg = np.zeros((crossbar_size, crossbar_size))
    mem_avg = np.zeros((crossbar_size, crossbar_size))
    weight_sum = np.zeros((crossbar_size, crossbar_size))
    cnt = np.zeros(crossbar_size)
    lrs_cnt = np.zeros((crossbar_size, crossbar_size))
    hrs_cnt = np.zeros((crossbar_size, crossbar_size))

    num_files = 0

    for subdir in os.listdir(base_folder):
        subdir_path = os.path.join(base_folder, subdir)

        print(f"{num_files} files processed".ljust(30), end='\r')

        if os.path.isdir(subdir_path) and re.match(dir_pattern, subdir):
            input_data = load_bin(os.path.join(subdir_path, "input.bin"), dtype=np.int64, shape_len=2)
            weight_data = load_bin(os.path.join(subdir_path, "weight.bin"), dtype=np.int64, shape_len=2)

            mac_data = load_bin(os.path.join(subdir_path, "mac.bin"), dtype=np.float64, shape_len=2)
            mem_data = load_bin(os.path.join(subdir_path, "mem.bin"), dtype=np.float64, shape_len=3)
            mac_data_out = load_bin(os.path.join(subdir_path, "out_MAC.bin"), dtype=np.float32, shape_len=2)
            mem_data_out = load_bin(os.path.join(subdir_path, "output.bin"), dtype=np.float32, shape_len=3)

            for batch in range(len(mac_data)):
                for col in range(len(mac_data[batch])):
                    # mac_err.append(mac_data[batch][col]*1e-6 - mac_data_out[batch][col])
                    # mac_err_reg.append(abs(mac_data[batch][col]*1e-6 - mac_data_out[batch][col]) / (mac_data[batch][col]*1e-6))
                    # mac_err.append(abs(mac_data[batch][col]*1e-6 - mac_data_out[batch][col]))
                    mac_err_reg.append(abs(mac_data[batch][col] - mac_data_out[batch][col]) / (mac_data[batch][col]))
                    mac_err.append(abs(mac_data[batch][col] - mac_data_out[batch][col]))
            
            for batch in range(len(mem_data)):
                for row in range(len(mem_data[batch])):
                    if input_data[batch][row]:
                        for col in range(len(mem_data[batch][row])):
                            # mem_err[row][col] += abs(mem_data[batch][row][col]*1e-6 - mem_data_out[batch][row][col])
                            # mem_err_reg[row][col] += abs(mem_data[batch][row][col]*1e-6 - mem_data_out[batch][row][col]) / (mem_data[batch][row][col]*1e-6)
                            # mem_avg[row][col] += mem_data[batch][row][col]*1e-6
                            mem_err[row][col] += abs(mem_data[batch][row][col] - mem_data_out[batch][row][col])
                            mem_err_reg[row][col] += abs(mem_data[batch][row][col] - mem_data_out[batch][row][col]) / (mem_data[batch][row][col])
                            mem_avg[row][col] += mem_data[batch][row][col]

                            if weight_data[row][col] > 0:
                                # lrs_err[row][col] += abs(mem_data[batch][row][col]*1e-6 - mem_data_out[batch][row][col])
                                # lrs_err_reg[row][col] += abs(mem_data[batch][row][col]*1e-6 - mem_data_out[batch][row][col]) / (mem_data[batch][row][col]*1e-6)
                                lrs_err[row][col] += abs(mem_data[batch][row][col] - mem_data_out[batch][row][col])
                                lrs_err_reg[row][col] += abs(mem_data[batch][row][col] - mem_data_out[batch][row][col]) / (mem_data[batch][row][col])
                                lrs_cnt[row][col] += 1
                            else:
                                # hrs_err[row][col] += abs(mem_data[batch][row][col]*1e-6 - mem_data_out[batch][row][col])
                                # hrs_err_reg[row][col] += abs(mem_data[batch][row][col]*1e-6 - mem_data_out[batch][row][col]) / (mem_data[batch][row][col]*1e-6)
                                hrs_err[row][col] += abs(mem_data[batch][row][col] - mem_data_out[batch][row][col])
                                hrs_err_reg[row][col] += abs(mem_data[batch][row][col] - mem_data_out[batch][row][col]) / (mem_data[batch][row][col])
                                hrs_cnt[row][col] += 1
                        cnt[row] += 1
        num_files += 1

    print("Finished reading data")

    for row in range(len(mem_err)):
        for col in range(len(mem_err[0])):
            mem_avg[row][col] = mem_avg[row][col] / cnt[row]
            mem_err[row][col] = mem_err[row][col] / cnt[row]
            mem_err_reg[row][col] = mem_err_reg[row][col] / cnt[row]
            lrs_err[row][col] = lrs_err[row][col] / lrs_cnt[row][col] 
            lrs_err_reg[row][col] = lrs_err_reg[row][col] / lrs_cnt[row][col]
            hrs_err[row][col] = hrs_err[row][col] / hrs_cnt[row][col] 
            hrs_err_reg[row][col] = hrs_err_reg[row][col] / hrs_cnt[row][col]

    counts, bins = np.histogram(mac_err, bins='auto')   
    counts_reg, bins_reg = np.histogram(mac_err_reg, bins='auto')

    max_bin_index = np.argmax(counts)
    max_bin_index_reg = np.argmax(counts_reg)

    # # Get the edges of the bin with the highest count
    bin_with_highest_count = (bins[max_bin_index], bins[max_bin_index + 1])
    bin_with_highest_count_reg = (bins_reg[max_bin_index_reg], bins_reg[max_bin_index_reg + 1])

    # # Get the count of entries in the bin with the highest count
    max_count = counts[max_bin_index]
    max_count_reg = counts_reg[max_bin_index_reg]

    print(f"Bin with the highest count: {bin_with_highest_count}")
    print(f"Count of entries in this bin: {max_count}")
    
    print(f"Relative bin with the highest count: {bin_with_highest_count_reg}")
    print(f"Count of entries in this bin: {max_count_reg}")
    
    plt.figure(1)
    plt.hist(bins[:-1], bins, weights=counts)
    plt.xlabel('Error')
    plt.ylabel('Count')
    plt.title('Error histogram')

    plt.figure(2)
    plt.hist(bins_reg[:-1], bins_reg, weights=counts_reg)
    plt.xlabel('Relative error')
    plt.ylabel('Count')
    plt.title('Relative error histogram')

    # plt.figure(figsize=(10, 7))

    # plt.subplot(2, 2, 1)
    plt.figure(3)
    plt.imshow(mem_err, cmap='hot', interpolation='nearest')
    plt.colorbar(label='Error magnitude (A)')
    plt.title('Crossbar Error Map')
    plt.xlabel('Crossbar Column')
    plt.ylabel('Crossbar Row')

    # plt.subplot(2, 2, 2)
    plt.figure(4)
    plt.imshow(mem_err_reg, cmap='hot', interpolation='nearest')
    plt.colorbar(label='Error magnitude')
    plt.title('Crossbar Relative Error Map')
    plt.xlabel('Crossbar Column')
    plt.ylabel('Crossbar Row')

    # plt.subplot(2, 2, 3)
    plt.figure(5)
    plt.imshow(lrs_err, cmap='hot', interpolation='nearest')
    plt.colorbar(label='Error magnitude')
    plt.title('Crossbar LRS error map')
    plt.xlabel('Crossbar Column')
    plt.ylabel('Crossbar Row')
    
    # plt.subplot(2, 2, 3)
    plt.figure(6)
    plt.imshow(lrs_err_reg, cmap='hot', interpolation='nearest')
    plt.colorbar(label='Relative error magnitude')
    plt.title('Crossbar LRS relative error map')
    plt.xlabel('Crossbar Column')
    plt.ylabel('Crossbar Row')

    # plt.subplot(2, 2, 4)
    plt.figure(7)
    plt.imshow(hrs_err, cmap='hot', interpolation='nearest')
    plt.colorbar(label='Error magnitude')
    plt.title('Crossbar HRS Error Map')
    plt.xlabel('Crossbar Column')
    plt.ylabel('Crossbar Row')
    
    # plt.subplot(2, 2, 4)
    plt.figure(8)
    plt.imshow(hrs_err_reg, cmap='hot', interpolation='nearest')
    plt.colorbar(label='Relative error magnitude')
    plt.title('Crossbar HRS relative error map')
    plt.xlabel('Crossbar Column')
    plt.ylabel('Crossbar Row')
    
    # plt.tight_layout()
    plt.show()


def other_main():
    base_folder = "RRAM_test_data/batch_66"
    dir_pattern = r"sim_\d+"

    crossbar_size = 128
    
    lrs = np.zeros((crossbar_size, crossbar_size))
    hrs = np.zeros((crossbar_size, crossbar_size))
    lrs_cnt = np.zeros((crossbar_size, crossbar_size))
    hrs_cnt = np.zeros((crossbar_size, crossbar_size))

    for subdir in os.listdir(base_folder):
        subdir_path = os.path.join(base_folder, subdir)

        if os.path.isdir(subdir_path) and re.match(dir_pattern, subdir):
            input_data = load_bin(os.path.join(subdir_path, "input.bin"), dtype=np.int64, shape_len=1)
            weight_data = load_bin(os.path.join(subdir_path, "weight.bin"), dtype=np.int64, shape_len=2)

            mac_data_out = load_bin(os.path.join(subdir_path, "out_MAC.bin"), dtype=np.float32, shape_len=1)
            mem_data_out = load_bin(os.path.join(subdir_path, "out.bin"), dtype=np.float32, shape_len=2)
            
            for row in range(len(mem_data_out)):
                if input_data[row]:
                    for col in range(len(mem_data_out[row])):
                        if (mem_data_out[row][col] > 1 or mem_data_out[row][col] < -1):
                            print(mem_data_out[row][col])
                        if weight_data[row][col] > 0:
                            lrs[row][col] += mem_data_out[row][col]
                            lrs_cnt[row][col] += 1
                        else:
                            hrs[row][col] += mem_data_out[row][col]
                            hrs_cnt[row][col] += 1

    for row in range(len(lrs)):
        for col in range(len(lrs[0])):
            lrs[row][col] = lrs[row][col] / lrs_cnt[row][col]
            hrs[row][col] = hrs[row][col] / hrs_cnt[row][col]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

    # Heatmap for matrix1
    cax1 = ax1.imshow(lrs, cmap='hot', aspect='auto')
    ax1.set_title('LRS')
    fig.colorbar(cax1, ax=ax1)

    # Heatmap for matrix2
    cax2 = ax2.imshow(hrs, cmap='hot', aspect='auto')
    ax2.set_title('HRS')
    fig.colorbar(cax2, ax=ax2)

    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    try:
        main()
        # other_main()
    except Exception as e:
        print(e)
    
    input("Press Enter to exit...")
