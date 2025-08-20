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
<<<<<<< HEAD
import sys
=======
>>>>>>> refs/remotes/origin/main

class MixedGenomeAlignment():
    def __init__(self, bam: str, barcode_fq_csv: str, genomes_dir: str, out_path: str, valid_barcodes: list):
        self.bam = bam
        self.barcode_fq_csv = barcode_fq_csv
        self.genomes_dir = genomes_dir
        self.out_path = out_path
<<<<<<< HEAD
        self.valid_barcodes = set(valid_barcodes)
=======
        self.valid_barcodes = valid_barcodes
>>>>>>> refs/remotes/origin/main
        
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
    
    
<<<<<<< HEAD
    def BamSpDistribution(self):
        bam = pysam.AlignmentFile(self.bam, 'rb', threads=10)
        print('[BamSpDistribution]: Bam file opened...', flush=True)
        print('[BamSpDistribution]: starting iterations over bam file...', flush=True)
        errors = 0
        counter = Counter()
        try:
            for align in tqdm(bam.fetch(), total=bam.mapped):
                if align.mapping_quality >= 20:
                    try:
                        sp_name = self.chrName2spName(bam.get_reference_name(align.reference_id))
                        counter[sp_name] += 1
                    except:
                        print(f'ERROR: refID:{align.reference_id} -> ref:{bam.get_reference_name(align.reference_id)}', flush=True)
                        print(f'queryName: {align.query_name}', flush=True)
                        errors += 1
                else:
                    continue
        except Exception as e:
            print(f"Error processing BAM file: {e}", flush=True)
            raise e
        finally:
            bam.close()
            print('[BamSpDistribution]: bam file closed', flush=True)
        print(f'[BamSpDistribution]: DONE, error_count: {errors}', flush=True)
        print(counter, flush=True)
        return counter
    
    def make_SpBarcodeDict(self):
=======
    
    def make_SpBarcodeDict(self):
        '''SpBarcodeDict = {}
        keys = self.get_spBarcodeDict_keys(self.genomes_dir)
        for key in keys:
            SpBarcodeDict[key] = []
        print(f'Initialized: {SpBarcodeDict}', flush=True)
        print('[make_SpBarcodeDict]: SpBarcodeDict initialized...', flush=True)'''
>>>>>>> refs/remotes/origin/main
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
<<<<<<< HEAD
    
    def barcode2sp(self):
        sp_counter = self.make_SpBarcodeDict()
        print('[barcode2sp]: SpBarcodeDict loaded...', flush=True)
        barcode2sp = pd.DataFrame(columns=['barcode', 'species'])
        species = [sp for sp in sp_counter.keys()]
        for barcode in tqdm(set(chain.from_iterable(sp_counter.values()))):
            counts = [c[barcode] if barcode in c.keys() else 0 for c in sp_counter.values()]
            sp = species[np.argmax(np.array(counts))]
            barcode2sp.loc[len(barcode2sp)] = [barcode, sp]
        barcode2sp.to_csv(os.path.join(self.out_path, 'barcode2sp.csv'))
        print('[barcode2sp]: DONE', flush=True)
        
=======
            
>>>>>>> refs/remotes/origin/main
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
<<<<<<< HEAD
    parser.add_argument('method', type=str)
=======
>>>>>>> refs/remotes/origin/main
    parser.add_argument('-b', '--bam', type=str)
    parser.add_argument('-f', '--barcode_fq_csv', type=str)
    parser.add_argument('-g', '--genomes_dir', type=str)
    parser.add_argument('-o', '--out_path', type=str)
    parser.add_argument('-v', '--valid_barcodes', type=str)
    parser.add_argument('-s', '--barcode_samples', type=int)
    
    args = parser.parse_args()
    
    print('CHECK: validating arguments...', flush=True)
    assert os.path.exists(args.genomes_dir)
    assert os.path.isfile(args.bam)
    assert os.path.exists(args.out_path)
    assert os.path.isfile(args.barcode_fq_csv)
    assert os.path.isfile(args.valid_barcodes)
    print('CHECK: arguments are valid', flush=True)
    
    print('[VALID BARCODES]: Processing provided valid barcodes into python list...', flush=True)
    if args.barcode_samples == None:
        valid_barcodes = pd.read_csv(args.valid_barcodes, names=['barcodes'])['barcodes'].to_list()
    else:
        SAMPLE = args.barcode_sample
        valid_barcodes = pd.read_csv(args.valid_barcodes, names=['barcodes'])['barcodes'].sample(n=SAMPLE, random_state=42).to_list()
<<<<<<< HEAD
    print('[VALID BARCODES]: DONE', flush=True)
    if args.method == 'visualize':
        MixedGenomeAlignment(args.bam, args.barcode_fq_csv, args.genomes_dir, args.out_path, valid_barcodes).Visualize()
    elif args.method == 'barcode2sp':
        MixedGenomeAlignment(args.bam, args.barcode_fq_csv, args.genomes_dir, args.out_path, valid_barcodes).barcode2sp()
    elif args.method == 'BamSpDistribution':
        MixedGenomeAlignment(args.bam, args.barcode_fq_csv, args.genomes_dir, args.out_path, valid_barcodes).BamSpDistribution()
    else:
        print('ERROR: wrong positional argument')
        sys.exit(1)
=======
    valid_barcodes = pd.read_csv(args.valid_barcodes, names=['barcodes'])['barcodes'].to_list()
    print('[VALID BARCODES]: DONE', flush=True)

    MixedGenomeAlignment(args.bam, args.barcode_fq_csv, args.genomes_dir, args.out_path, valid_barcodes).Visualize()
>>>>>>> refs/remotes/origin/main
