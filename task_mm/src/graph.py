# SPDX-License-Identifier: Apache-2.0

import matplotlib.pyplot as plt
import csv

def plot_new_csv(file_path):
    with open(file_path, newline='') as csvfile:
        reader = csv.reader(csvfile)
        headers = next(reader)
        data = list(reader)

    n_values = [int(row[0]) for row in data]
    labels = headers[1:]
    values = {label: [float(row[i + 1]) for row in data]
              for i, label in enumerate(labels)}

    categories = ['avg', '95p', '1p']
    markers = ['o', 's', '^', 'd', 'x', 'v', '<', '>', 'p']

    for scale in ['linear', 'log']:
        for i, category in enumerate(categories):
            plt.figure(figsize=(10, 6))
            if scale == 'log':
                plt.yscale('log')
            for j, (label, val) in enumerate(values.items()):
                if category in label:
                    plt.plot(n_values, val, marker=markers[j % len(markers)], label=label)
            plt.xlabel('n')
            plt.ylabel(f'Execution Time ({category}) ({scale} scale)')
            plt.title(f'Execution Time vs n for {category.upper()} Metrics ({scale} scale)')
            plt.legend()
            plt.grid()
            plt.show()

if __name__ == "__main__":
    file_path = 'out.csv'
    plot_new_csv(file_path)
