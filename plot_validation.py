import os
import re

import numpy as np
import matplotlib.pyplot as plt

def save_bin(path,array):
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
        # print(shape)
        # Reshape
        array = data.reshape(shape)
       
    return array


def main():
    base_folder = "data/bin/batch_1"
    dir_pattern = r"row_\d+-\d+_col_\d+-\d+"

    mac_err = []
    mem_err = np.zeros((32, 32))
    mem_err_reg = np.zeros((32, 32))
    lrs_err = np.zeros((32, 32))
    lrs_err_reg = np.zeros((32, 32))
    hrs_err = np.zeros((32, 32))
    hrs_err_reg = np.zeros((32, 32))
    mem_avg = np.zeros((32, 32))
    weight_sum = np.zeros((32, 32))
    cnt = np.zeros(32)
    lrs_cnt = np.zeros((32, 32))
    hrs_cnt = np.zeros((32, 32))

    for subdir in os.listdir(base_folder):
        subdir_path = os.path.join(base_folder, subdir)

        if os.path.isdir(subdir_path) and re.match(dir_pattern, subdir):
            input_data = load_bin(os.path.join(subdir_path, "input.bin"), dtype=np.int64, shape_len=2)
            weight_data = load_bin(os.path.join(subdir_path, "weight.bin"), dtype=np.int64, shape_len=2)

            mac_data = load_bin(os.path.join(subdir_path, "mac.bin"), dtype=np.float64, shape_len=2)
            mem_data = load_bin(os.path.join(subdir_path, "mem.bin"), dtype=np.float64, shape_len=3)
            mac_data_out = load_bin(os.path.join(subdir_path, "out_MAC.bin"), dtype=np.float32, shape_len=2)
            mem_data_out = load_bin(os.path.join(subdir_path, "output.bin"), dtype=np.float32, shape_len=3)

            for batch in range(len(mac_data)):
                for col in range(len(mac_data[batch])):
                    mac_err.append(mac_data[batch][col]*1e-6 - mac_data_out[batch][col])
                    # mac_err.append(abs(mac_data[batch][col]*1e-6 - mac_data_out[batch][col]))
            
            for batch in range(len(mem_data)):
                for row in range(len(mem_data[batch])):
                    if input_data[batch][row]:
                        for col in range(len(mem_data[batch][row])):
                            mem_err[row][col] += abs(mem_data[batch][row][col]*1e-6 - mem_data_out[batch][row][col])
                            mem_err_reg[row][col] += abs(mem_data[batch][row][col]*1e-6 - mem_data_out[batch][row][col]) / (mem_data[batch][row][col]*1e-6)
                            mem_avg[row][col] += mem_data[batch][row][col]*1e-6

                            if weight_data[row][col] > 0:
                                lrs_err[row][col] += abs(mem_data[batch][row][col]*1e-6 - mem_data_out[batch][row][col])
                                lrs_err_reg[row][col] += abs(mem_data[batch][row][col]*1e-6 - mem_data_out[batch][row][col]) / (mem_data[batch][row][col]*1e-6)
                                lrs_cnt[row][col] += 1
                            else:
                                hrs_err[row][col] += abs(mem_data[batch][row][col]*1e-6 - mem_data_out[batch][row][col])
                                hrs_err_reg[row][col] += abs(mem_data[batch][row][col]*1e-6 - mem_data_out[batch][row][col]) / (mem_data[batch][row][col]*1e-6)
                                hrs_cnt[row][col] += 1
                        cnt[row] += 1

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
    plt.hist(bins[:-1], bins, weights=counts)
    plt.xlabel('Error (uA)')
    plt.ylabel('Count')
    plt.title('Histogram of error')
    plt.show()

    # plt.figure(figsize=(10, 7))

    # plt.subplot(2, 2, 1)
    # plt.imshow(mem_err, cmap='hot', interpolation='nearest')
    # plt.colorbar(label='Error magnitude')
    # plt.title('Crossbar Error Map')
    # plt.xlabel('Crossbar Column')
    # plt.ylabel('Crossbar Row')

    # plt.subplot(2, 2, 2)
    # plt.imshow(mem_err_reg, cmap='hot', interpolation='nearest')
    # plt.colorbar(label='Error magnitude')
    # plt.title('Crossbar Relative Error Map')
    # plt.xlabel('Crossbar Column')
    # plt.ylabel('Crossbar Row')

    # plt.subplot(2, 2, 3)
    # plt.imshow(lrs_err, cmap='hot', interpolation='nearest')
    # # plt.imshow(lrs_err_reg, cmap='hot', interpolation='nearest')
    # plt.colorbar(label='Error magnitude')
    # plt.title('Crossbar LRS Error Map')
    # plt.xlabel('Crossbar Column')
    # plt.ylabel('Crossbar Row')

    # plt.subplot(2, 2, 4)
    # plt.imshow(hrs_err, cmap='hot', interpolation='nearest')
    # # plt.imshow(hrs_err_reg, cmap='hot', interpolation='nearest')
    # plt.colorbar(label='Error magnitude')
    # plt.title('Crossbar HRS Error Map')
    # plt.xlabel('Crossbar Column')
    # plt.ylabel('Crossbar Row')
    
    # plt.tight_layout()
    # plt.show()


if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        print(e)
    
    input("Press Enter to exit...")
