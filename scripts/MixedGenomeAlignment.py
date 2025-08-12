import pandas as pd
import numpy as np
import os
from Bio import SeqIO
import gzip
import pysam
from tqdm import tqdm
import argparse
import seaborn as sns
import matplotlib.pyplot as plt

class MixedGenomeAlignment():
    def __init__(self, bam: str, barcode_fq_list: list, genomes_dir: str):
        self.bam = bam
        self.barcode_fq_list = barcode_fq_list
        self.genomes_dir = genomes_dir
        
    def make_Name2seqDF(self):
        name2seq = pd.DataFrame(columns=['read_name', 'barcode'])
        for fastq in self.barcode_fq_list:
            tmp_df = pd.DataFrame(columns=['read_name', 'barcode'])
            with gzip.open(fastq) as fq:
                for record in SeqIO.parse(fq, 'fastq'):
                    tmp_df.loc[len(tmp_df)] = [record.id, str(record.seq)]
            name2seq = pd.concat([name2seq, tmp_df])
        return name2seq
    
    @staticmethod
    def get_spBarcodeDict_keys(genomes_dir):
        genomes = os.listdir(genomes_dir)
        return [genome.split('.')[0] for genome in genomes]
    
    @staticmethod
    def Name2seq(Name2seqDF, name):
        return Name2seqDF[Name2seqDF['read_name'] == name]['barcode'].values[0]
    
    @staticmethod    
    def chrName2spName(chrName):
        return "_".join(chrName.split('_')[0:-1])
    
    def make_SpBarcodeDict(self):
        SpBarcodeDict = {}
        keys = self.get_spBarcodeDict_keys(self.genomes_dir)
        for key in keys:
            SpBarcodeDict[key] = []
        bam = pysam.AlignmentFile(self.bam, 'rb')
        Name2seq = self.make_Name2seqDF()
        for align in tqdm(bam.fetch()):
            sp_name = self.chrName2spName(align.reference_id)
            barcode = self.Name2seq(Name2seq, align.query_name)
            SpBarcodeDict[sp_name] = SpBarcodeDict[sp_name].append(barcode)
        return SpBarcodeDict
            
    def calculate_distibution(self):
        SpBarcodeDict = self.make_SpBarcodeDict()
        barcode_list = np.unique([barcode for species_barcodes in SpBarcodeDict.values() for barcode in species_barcodes])
        M_Statistics_list = []
        for barcode in barcode_list:
            count_list = []
            for key in SpBarcodeDict.keys():
                count = SpBarcodeDict[key].count(barcode)
                count_list.append(count)
            M_Statistics_list.append(max(count_list)/sum(count_list))
        return M_Statistics_list
    
    def Visualize(self, out_path):
        M_Statistics_list = self.calculate_distibution()
        sns.displot(M_Statistics_list, kde=True, x = 'M_statistics')
        plt.savefig(os.path.join(out_path, "M_Statistics_distribution.png"), dpi=600, bbox_inches='tight', format='png')
        
if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(prog='MixedGenomeAlignment')
    parser.add_argument('-b', '--bam')
    parser.add_argument('-f', '--barcode_fq_list')
    parser.add_argument('-g', '--genomes_dir')
    parser.add_argument('-o', '--out_path')
    
    args = parser.parse_args()
    
    MixedGenomeAlignment(args.bam, args.barcode_fq_list, args.genomes_dir).Visualize(args.out_path)