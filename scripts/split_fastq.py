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
from multiprocessing import Pool
import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from pandarallel import pandarallel

pandarallel.initialize(nb_workers=11)

class splitFastq():
    def __init__(self, 
                out_path : str, 
                barcode_fq_csv: str, 
                barcode2sp: str, 
                fastq: str,
                processes=10, 
                ):
        self.out_path = out_path
        self.barcode_fq_csv = barcode_fq_csv
        self.barcode2sp = barcode2sp
        self.processes = processes
        self.fastq = fastq
        
    def make_barcode2readName_dict(self) -> dict:
        print('********** make_barcode2readName_dict **********', flush=True)
        barcode2readName_dict = defaultdict(list)
        with open(self.barcode_fq_csv, mode='r') as infile:
            reader = csv.reader(infile)
            for rows in reader:
                read_name = rows[0]
                barcode = rows[1][0:20]
                barcode2readName_dict[barcode].append(read_name)
        print('[make_barcode2readName_dict]: DONE', flush=True)
        return barcode2readName_dict
    
    def make_sp2barcodes(self):
        def write_barcode_lists(row, out_path=self.out_path):
            filename = row.values[0] + '_barcodes.csv'
            barcodes = pd.DataFrame({'barcode': row.values[1]})
            out_file = os.path.join(out_path, filename)
            barcodes.to_csv(out_file, header=False, index=False)
        barcode2sp = pd.read_csv(self.barcode2sp, index_col=0)
        print(f'[make_sp2barcodes]: barcode2sp df loaded, columns: {barcode2sp.columns}', flush=True)
        sp2barcodes = barcode2sp.groupby(by='species').agg(lambda x: ','.join(x).split(',')).reset_index()
        print(f'[make_sp2barcodes]: sp2barcodes groupby df loaded, columns: {sp2barcodes.columns}', flush=True)
        sp2barcodes.apply(write_barcode_lists, axis=1)
        print(f'[make_sp2barcodes]: DONE', flush=True)
        return sp2barcodes
    
    def sp2barcodes_to_sp2readNames(self) -> pd.DataFrame:
        sp2barcodes = self.make_sp2barcodes()
        barcode2readName_dict = self.make_barcode2readName_dict()
        print('[sp2barcodes_to_sp2readNames]: sp2barcodes and barcode2readName_dict loaded', flush=True)
        sp2barcodes['barcode'] = sp2barcodes['barcode'].apply(lambda barcodes: [x for xs in [barcode2readName_dict[barcode] for barcode in barcodes] for x in xs if type(xs) == list])
        return sp2barcodes
    
    def writeFastq(self):
        def make_fastq(row):
            out_file = row.iloc[0]
            readNames_set = set(row.iloc[1])
            with open(self.fastq, "r") as infastq, open(out_file, "w") as outfastq:
                for record in tqdm(SeqIO.parse(infastq, "fastq"), total=770731614):
                    if record.id.replace('/2', '/1') in readNames_set:
                        outfastq.write(f"@{record.id}\n")
                        outfastq.write(f"{record.seq}\n")
                        outfastq.write("+\n")
                        outfastq.write(f"{''.join([chr(i+33) for i in record.letter_annotations['phred_quality']])}\n")
        sp2readNames = self.sp2barcodes_to_sp2readNames()
        sp2readNames['species'] = sp2readNames.species.apply(lambda x: os.path.join(self.out_path, x+'.fastq')).to_list()
        print('[writeFastq]: writing fastq...', flush=True)
        sp2readNames.parallel_apply(make_fastq, axis=1)
        print('DONE', flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    
    parser.add_argument('-b', '--barcode2sp', type=str)
    parser.add_argument('-f', '--fastq', type=str)
    parser.add_argument('-o', '--out_path', type=str)
    parser.add_argument('-c', '--barcode_fq_csv', type=str)
    
    args = parser.parse_args()
    print('arguments are loaded...', flush=True)
    
    splitFastq(args.out_path, args.barcode_fq_csv, args.barcode2sp, args.fastq).writeFastq()
    
    