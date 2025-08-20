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
        barcode2readName_dict = {}
        with open(self.barcode_fq_csv, mode='r') as infile:
            reader = csv.reader(infile)
            for rows in reader:
                read_name = rows[0]
                barcode = rows[1][0:20]
                barcode2readName_dict[barcode] = read_name
        print('[make_barcode2readName_dict]: DONE')
        return barcode2readName_dict

    @staticmethod
    def barcodes2readnames(barcodes: list, barcode2readName_dict: dict) -> list:
        return [barcode2readName_dict[barcode] for barcode in barcodes]
    
    def make_sp2barcodes(self):
        def write_barcode_lists(row, out_path=self.out_path):
            filename = row.values[0] + '_barcodes.csv'
            barcodes = pd.DataFrame({'barcode': row.values[1]})
            out_file = os.path.join(out_path, filename)
            barcodes.to_csv(out_file, header=False, index=False)
        barcode2sp = pd.read_csv(self.barcode2sp)
        sp2barcodes = barcode2sp.groupby(by='species').agg(lambda x: ','.join(x).split(',')).reset_index()
        sp2barcodes.apply(write_barcode_lists, axis=1)
    
    def sp2barcodes_to_sp2readNames(self) -> pd.DataFrame:
        barcode2readName_dict = self.make_barcode2readName_dict()
        sp2barcodes = self.make_sp2barcodes()
        p = Pool(self.processes)
        sp2barcodes['barcode'] = p.map(self.barcodes2readnames(), sp2barcodes['barcode'])
        return sp2barcodes
    
    def writeFastq(self):
        def make_fastq(row):
            out_file = row[0]
            readNames_set = set(row[1])
            with open(self.fastq, "r") as infastq, open(out_file, "w") as outfastq:
                for record in SeqIO.parse(infastq, "fastq"):
                    if record.id.replace('/1', '/2') in readNames_set:
                        outfastq.write(f"@{record.id}\n")
                        outfastq.write(f"{record.seq}\n")
                        outfastq.write("+\n")
                        outfastq.write(f"{record.letter_annotations['phred_quality']}\n")
        sp2readNames = self.sp2barcodes_to_sp2readNames(self)
        p = Pool(self.processes)
        p.map(lambda df: df.apply(make_fastq, axis=1), sp2readNames)
        print('DONE', flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    
    parser.add_argument('-b', '--barcode2sp', type=str)
    parser.add_argument('-f', '--fastq', type=str)
    parser.add_argument('-o', '--out_path', type=str)
    parser.add_argument('-c', '--barcode_fq_csv', type=str)
    
    args = parser.parse_args()
    
    splitFastq(args.out_path, args.barcode_fq_csv, args.barcode2sp, args.fastq).writeFastq()
    
    