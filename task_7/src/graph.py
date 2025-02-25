import matplotlib.pyplot as plt
import csv
import numpy as np

def plot_from_csv(file_path):
    with open(file_path, newline='') as csvfile:
        reader = csv.reader(csvfile)
        headers = next(reader)
        data = list(reader)
    
    n_values = [int(row[0]) for row in data]
    D_left = [float(row[1]) for row in data]
    D_mid = [float(row[2]) for row in data]
    D_trapz = [float(row[3]) for row in data]
    D_simpson = [float(row[4]) for row in data]
    D_nc = [float(row[5]) for row in data]
    D_gs = [float(row[6]) for row in data]
    
    plt.figure(figsize=(10, 5))
    
    plt.subplot(1, 2, 1)
    plt.yscale('log')
    plt.plot(n_values, D_left, marker='o', label='D_left')
    plt.plot(n_values, D_mid, marker='s', label='D_mid')
    plt.plot(n_values, D_trapz, marker='^', label='D_trapz')
    plt.plot(n_values, D_simpson, marker='d', label='D_simpson')
    plt.xlabel('n')
    plt.ylabel('D values (log scale)')
    plt.title('Graph 1: First Four D Values')
    plt.legend()
    plt.grid()
    
    plt.subplot(1, 2, 2)
    plt.yscale('log')
    plt.plot(n_values, D_nc, marker='o', label='D_nc', color='r')
    plt.plot(n_values, D_gs, marker='s', label='D_gs', color='g')
    plt.xlabel('n')
    plt.ylabel('D values (log scale)')
    plt.title('Graph 2: Last Two D Values')
    plt.legend()
    plt.grid()
    
    plt.tight_layout()
    plt.show()

def plot_second_csv(file_path):
    with open(file_path, newline='') as csvfile:
        reader = csv.reader(csvfile)
        headers = next(reader)
        data = list(reader)
    
    n_values = [int(row[0]) for row in data]
    D_nc_compound = [float(row[1]) for row in data]
    D_nc_def = [float(row[2]) for row in data]
    
    plt.figure(figsize=(8, 5))
    plt.yscale('log')
    plt.plot(n_values, D_nc_compound, marker='o', label='D_nc_compound', color='b')
    plt.plot(n_values, D_nc_def, marker='s', label='D_nc_def', color='m')
    plt.xlabel('n')
    plt.ylabel('D values (log scale)')
    plt.title('Graph 3: D_nc_compound vs D_nc_def')
    plt.legend()
    plt.grid()
    plt.show()


def main():
    file_path_1 = 'out/tonethree.csv'
    file_path_2 = 'out/pentwo.csv'
    plot_from_csv(file_path_1)
    plot_second_csv(file_path_2)

if __name__ == "__main__":
    main()
