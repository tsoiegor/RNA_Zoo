import pandas as pd
import numpy as np
import os
import pysam
from tqdm import tqdm
import argparse
import seaborn as sns
import matplotlib.pyplot as plt
import csv
from itertools import chain
from collections import defaultdict, Counter

class MixedGenomeAlignment():
    def __init__(self, bam: str, barcode_fq_csv: str, genomes_dir: str, out_path: str, valid_barcodes: list):
        self.bam = bam
        self.barcode_fq_csv = barcode_fq_csv
        self.genomes_dir = genomes_dir
        self.out_path = out_path
        self.valid_barcodes = valid_barcodes
        
    def make_Name2seq_dict(self):
        print('********** make_Name2seq_dict **********', flush=True)
        Name2seq_dict = {}
        with open(self.barcode_fq_csv, mode='r') as infile:
            reader = csv.reader(infile)
            for rows in reader:
                read_name = rows[0]
                barcode = rows[1]
                Name2seq_dict[read_name] = barcode
        print('[make_Name2seq_dict]: DONE')
        return Name2seq_dict
        
    @staticmethod
    def get_spBarcodeDict_keys(genomes_dir: str) -> list:
        print('********** get_spBarcodeDict_keys **********', flush=True)
        genomes = os.listdir(genomes_dir)
        print('[get_spBarcodeDict_keys]: DONE', flush=True)
        return [genome.split('.')[0] for genome in genomes]
    
    @staticmethod    
    def chrName2spName(chrName: str) -> str:
        return "_".join(chrName.split('_')[0:2])
    
    
    
    def make_SpBarcodeDict(self):
        '''SpBarcodeDict = {}
        keys = self.get_spBarcodeDict_keys(self.genomes_dir)
        for key in keys:
            SpBarcodeDict[key] = []
        print(f'Initialized: {SpBarcodeDict}', flush=True)
        print('[make_SpBarcodeDict]: SpBarcodeDict initialized...', flush=True)'''
        bam = pysam.AlignmentFile(self.bam, 'rb', threads=10)
        print('[make_SpBarcodeDict]: Bam file opened...', flush=True)
        Name2seq = self.make_Name2seq_dict()
        counter = defaultdict(Counter)
        print('[make_SpBarcodeDict]: starting iterations over bam file...', flush=True)
        errors = 0
        try:
            for align in tqdm(bam.fetch(), total=bam.mapped):
                if align.mapping_quality >= 20:
                    try:
                        barcode = Name2seq[align.query_name.replace('/2', '/1')][0:20]
                        if barcode in self.valid_barcodes:
                            sp_name = self.chrName2spName(bam.get_reference_name(align.reference_id))
                            counter[sp_name][barcode] += 1
                        else:
                            continue
                    except:
                        print(f'ERROR: refID:{align.reference_id} -> ref:{bam.get_reference_name(align.reference_id)}', flush=True)
                        print(f'queryName: {align.query_name}', flush=True)
                        print(f'sp name: {sp_name}, barcode: {barcode}', flush=True)
                        errors += 1
                else:
                    continue
        except Exception as e:
            print(f"Error processing BAM file: {e}", flush=True)
            raise e
        finally:
            bam.close()
            print('[make_SpBarcodeDict]: bam file closed', flush=True)
        
        print(f'[make_SpBarcodeDict]: DONE, error_count: {errors}', flush=True)
        return counter
            
    def calculate_distribution(self):
        sp_counter = self.make_SpBarcodeDict()
        output_df = pd.DataFrame(columns=['barcode', 'M'], index=None, )
        print('[calculate_distribution]: SpBarcodeDict loaded...', flush=True)
        M_stats = []
        for barcode in tqdm(set(chain.from_iterable(sp_counter.values()))):
            counts = [c[barcode] for c in sp_counter.values() if barcode in c.keys()]
            total = sum(counts)
            if total == 0: continue
            M = max(counts) / total
            output_df.loc[len(output_df)] = [barcode, M]
            M_stats.append(M)
        output_df.to_csv(os.path.join(self.out_path, 'barcode2M_stat.csv'))
        print('[calculate_distribution]: DONE', flush=True)
        return M_stats

    def Visualize(self):
        M_Statistics_list = self.calculate_distribution()
        print('[Visualize]: M_Statistics_list loaded...', flush=True)
        M_stats_df = pd.DataFrame({'M_statistics': M_Statistics_list})
        M_stats_df.to_csv(os.path.join(self.out_path, 'M_stats_df.csv'))
        sns.displot(data=M_stats_df, kind='hist', kde=True, x = 'M_statistics')
        plt.savefig(os.path.join(self.out_path, "M_Statistics_distribution.png"), dpi=600, bbox_inches='tight', format='png')
        print('[Visualize]: DONE', flush=True)
        
if __name__ == '__main__':
    parser = argparse.ArgumentParser(prog='MixedGenomeAlignment')
    parser.add_argument('-b', '--bam')
    parser.add_argument('-f', '--barcode_fq_csv')
    parser.add_argument('-g', '--genomes_dir')
    parser.add_argument('-o', '--out_path')
    parser.add_argument('-v', '--valid_barcodes')
    
    args = parser.parse_args()
    
    print('CHECK: validating arguments...', flush=True)
    assert os.path.exists(args.genomes_dir)
    assert os.path.isfile(args.bam)
    assert os.path.exists(args.out_path)
    assert os.path.isfile(args.barcode_fq_csv)
    assert os.path.isfile(args.valid_barcodes)
    print('CHECK: arguments are valid', flush=True)
    
    print('[VALID BARCODES]: Processing provided valid barcodes into python list...', flush=True)
    SAMPLE = 1000
    valid_barcodes = pd.read_csv(args.valid_barcodes, names=['barcodes'])['barcodes'].sample(n=SAMPLE, random_state=42).to_list()
    print('[VALID BARCODES]: DONE', flush=True)

    MixedGenomeAlignment(args.bam, args.barcode_fq_csv, args.genomes_dir, args.out_path, valid_barcodes).Visualize()
