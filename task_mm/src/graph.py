# SPDX-License-Identifier: Apache-2.0

import matplotlib.pyplot as plt
import csv
import glob
import sys

def read_csv(file_path):
    with open(file_path, newline='') as csvfile:
        reader = csv.reader(csvfile)
        headers = next(reader)
        data = list(reader)
    n_values = [int(row[0]) for row in data]
    labels = headers[1:]
    values = {label: [float(row[i + 1]) for row in data] for i, label in enumerate(labels)}
    return n_values, values

def plot_all_csv():
    csv_files = glob.glob("results/*.csv")
    if not csv_files:
        print("Error: No csv's found in 'results/'")
        sys.exit(1)

    ref_n_values, _ = read_csv(csv_files[0])
    
    data_list = []
    for file in csv_files:
        n_values, values = read_csv(file)
        if n_values != ref_n_values:
            print(f"Error: files have diffrent n-s")
            sys.exit(1)
        data_list.append({'file': file, 'n_values': n_values, 'values': values})

    categories = ['avg', '95p', '1p']
    markers = ['o', 's', '^', 'd', 'x', 'v', '<', '>', 'p']

    for scale in ['log']: # ['linear', 'log']:
        for category in categories:
            plt.figure(figsize=(10, 6))
            if scale == 'log':
                plt.yscale('log')
            marker_idx = 0
            for entry in data_list:
                file_name = entry['file']
                n_values = entry['n_values']
                values = entry['values']
                for label, val in values.items():
                    if category in label:
                        plot_label = f"{file_name} - {label}"
                        plt.plot(n_values, val, marker=markers[marker_idx % len(markers)], label=plot_label)
                        marker_idx += 1
            plt.xlabel('n')
            plt.ylabel(f'Elapsed Time ({category}) [{scale} scale]')
            plt.title(f'Elapsed Time vs n {category.upper()} [{scale} scale]')
            plt.legend()
            plt.grid()
            plt.show()

if __name__ == "__main__":
    plot_all_csv()

