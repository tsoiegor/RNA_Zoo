import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import argparse
import os

parser = argparse.ArgumentParser(prog='barcode_freq_dist')
parser.add_argument('-b', '--barcode')
parser.add_argument('-o', '--out_path')
args = parser.parse_args()

barcode_df = pd.read_csv(args.barcode, names=['barcodes'])
value_count_DF = pd.DataFrame(barcode_df.value_counts())

sns.displot(data=value_count_DF, x='count', kde=True)
plt.savefig(os.path.join(args.out_path, "barcode_freq_dist.png"), dpi=600, bbox_inches='tight', format='png')

