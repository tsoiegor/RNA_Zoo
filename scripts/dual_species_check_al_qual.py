import pysam
from tqdm import tqdm
import matplotlib.pyplot as plt
import argparse

parser = argparse.ArgumentParser(
                    prog='Dual species alignment quality check',
                    desc= 'please provide one positional argument: path to bam file')
parser.add_argument('bam')
args = parser.parse_args()

bam = pysam.AlignmentFile(args.bam)

def read_source_readName(name):
    if 'NS500615' in name:
        return 'Human'
    elif 'HISEQ' in name:
        return 'Mouse'
    else:
        raise NameError
def read_source_chromName(name):
    if 'Human_' in name:
        return 'Human'
    elif 'Mouse_' in name:
        return 'Mouse'
    else:
        raise NameError
def plot_results(match, chrom_mouse_read_human, chrom_human_read_mouse):
    data = [match, chrom_mouse_read_human, chrom_human_read_mouse]
    labels = ['match', 'chrom_mouse_read_human', 'chrom_human_read_mouse']
    colors = sns.color_palette('pastel')

    fig, ax = plt.subplots()

    wedges, texts, autotexts = ax.pie(
        data,
        labels=labels,
        colors=colors,
        startangle=90,
        autopct='%.0f%%',
        wedgeprops=dict(width=0.5),
        labeldistance=2
    )
    ax.legend(wedges, labels, title="Classes", loc="center left", bbox_to_anchor=(1, 0, 0.5, 1))
    plt.show()

match = 0
chrom_mouse_read_human = 0
chrom_human_read_mouse = 0
for alignment in tqdm(bam.fetch()):
    species_chromName = read_source_chromName(alignment.reference_name)
    species_readName = read_source_readName(alignment.query_name)
    if species_chromName == species_readName:
        match += 1
    elif species_chromName == 'Mouse' and species_readName == "Human":
        chrom_mouse_read_human +=1
    elif species_chromName == 'Human' and species_readName == "Mouse":
        chrom_human_read_mouse +=1
    else:
        raise NameError
    
plot_results(match, chrom_mouse_read_human, chrom_human_read_mouse)
