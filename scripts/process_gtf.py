import pandas as pd
import numpy as np
from tqdm import tqdm
import csv
import argparse

def modifications_needed(annotation_path: str):
    annotation = pd.read_csv(annotation_path, sep='\t', names=['scaffold', 'produced_by', 'feature_type', 'start', 'end', 'dot', 'strand', 'frame', 'desc'], comment='#')
    genes= annotation[annotation['feature_type'] == 'gene']
    if genes.desc.apply(lambda x: "biotype" in x).sum().item() == len(genes):
        return False
    else:
        return True


class process_annotation():
    def __init__(self, annotation_path, add_type=True):
        self.annotation_path = annotation_path
        self.annotation = pd.read_csv(self.annotation_path, sep='\t', names=['scaffold', 'produced_by', 'feature_type', 'start', 'end', 'dot', 'strand', 'frame', 'desc'], comment='#')
        self.add_type= add_type
    def desc2fullDesc(self):
        desc2fullDescDict = {}
        unique_desc_full = self.annotation.desc[self.annotation.desc.apply(lambda x: '_id' in x)].unique()
        for desc in tqdm(self.annotation.desc):
            if "_id" in desc:
                desc2fullDescDict[desc] = desc
            else:
                for full_desc in unique_desc_full:
                    if f'"{desc}"' in full_desc:
                        desc2fullDescDict[desc] = full_desc
        return desc2fullDescDict
    def process_description(self):
        desc2fullDescDict = self.desc2fullDesc()
        annotation = self.annotation
        annotation['desc'] = annotation['desc'].apply(lambda x: desc2fullDescDict[x])
        if self.add_type: annotation['desc'] = (annotation['desc'] + ' gene_biotype "protein_coding";').str.replace('\"\"', '')
        return annotation

if __name__ == '__main__':
    parser = argparse.ArgumentParser(prog='process_annotation')
    parser.add_argument('-i', '--ingtf')
    parser.add_argument('-o', '--outgtf')
    args = parser.parse_args()
    if modifications_needed(args.ingtf):
        outgtf = process_annotation(args.ingtf).process_description()
        outgtf.to_csv(args.outgtf, index=False, sep='\t', header=False, quoting=csv.QUOTE_NONE)
    else:
        os.symlink(args.ingtf, args.outgtf)
