import pandas as pd
import argparse

parser = argparse.ArgumentParser(prog="gtf_stats")
parser.add_argument('annotation')

args = parser.parse_args()


annotation = pd.read_csv(args.annotation, sep='\t', header=None, comment='#')
print(annotation.iloc[:, 2].value_counts())
