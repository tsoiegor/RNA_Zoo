import pandas as pd
import numpy as np
import os
import argparse
import matplotlib.pyplot as plt
from tqdm import tqdm

parser = argparse.ArgumentParser(prog='Calculate_readLoss')
parser.add_argument('-b', '--barcode2count')
parser.add_argument('-o', '--out_path')
parser.add_argument('-n', '--n_points', type=int)
parser.add_argument('-l', '--max_count_limit', type=int)

args = parser.parse_args()


PATH_TO_barcode2count = args.barcode2count
OUT_PATH = args.out_path
N_POINTS = args.n_points


barcode2count = pd.read_csv(PATH_TO_barcode2count, skiprows=1, names=['barcode', 'count'])

min_count, max_count = barcode2count['count'].min(), barcode2count['count'].max()
max_count = min(max_count, args.max_count_limit)
step = (max_count - min_count)/N_POINTS
reads_total = barcode2count['count'].sum()


reads_left_df = pd.DataFrame(columns=['count_border', 'reads_left_frac'])
for count_step in tqdm(np.arange(min_count, max_count, step)):
    count_border = count_step
    reads_left = barcode2count['count'][barcode2count['count'] >= count_step].sum()
    reads_left_df.loc[len(reads_left_df)] = [count_border, reads_left/reads_total*100]

reads_left_df.to_csv(os.path.join(OUT_PATH, 'readsN_by_count_border.csv'))

plt.scatter(x = reads_left_df['count_border'], y = reads_left_df['reads_left_frac'], s=5)
plt.grid(axis='both')
plt.xlabel('count_border')
plt.ylabel('reads_left_%')
plt.savefig(os.path.join(OUT_PATH, 'readsN_by_count_border.png'), format='png', dpi=600, bbox_inches='tight')