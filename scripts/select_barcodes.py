import pandas as pd
import numpy as np
import os
import argparse
import matplotlib.pyplot as plt

<<<<<<< HEAD


=======
>>>>>>> refs/remotes/origin/main
parser = argparse.ArgumentParser(prog='sample_from_barcodes')

parser.add_argument('-b', '--barcode_file', type=str)
parser.add_argument('-m', '--m_stat', type=str)
parser.add_argument('-t', '--count_thr', type=int)
parser.add_argument('-o', '--out_path', type=str)
parser.add_argument('-M', '--m_thr', type=float)

args = parser.parse_args()
barcodes_total = pd.read_csv(args.barcode_file, header=0, sep=',', )

<<<<<<< HEAD
print('Making barcode2count dict...', flush=True)
barcode2count = dict(zip(barcodes_total['barcodes'], barcodes_total['count']))
print('Done')

if args.m_stat == None:
    count_thr = args.count_thr
    barcodes_total.query('count > @count_thr')['barcodes'].to_csv(os.path.join(args.out_path, 'selected_barcodes.csv'), index=None)
else:
    barcode2M_stat = pd.read_csv(args.m_stat)[['barcode', 'M']]
    count_threshold = args.count_thr
    M_threshold = args.m_thr
    barcode2M_stat['count'] = barcode2M_stat['barcode'].apply(lambda barcode: barcode2count[barcode])

    print('\n*********** STATS: cell loss ***********')

    print(f'''Barcodes before filtration by count: {len(barcodes_total)}, \n
    Barcodes after filtration by count: {len(barcode2M_stat.query('count > @count_threshold'))}, \n
    Impact of filtration by count: {'{}%'.format(np.round(len(barcodes_total.query('count <= @count_threshold'))/len(barcodes_total)*100, 2))}, \n
    Barcodes after filtreation by count and M: {len(barcode2M_stat.query('count > @count_threshold').query('M > @M_threshold'))}, \n
    Impact of filtration by M: {'{}%'.format(np.round(len(barcode2M_stat.query('count > @count_threshold').query('M <= @M_threshold'))/len(barcode2M_stat.query('count > @count_threshold'))*100, 2))}''')

    print('\n*********** STATS: read loss ***********')
    print(f'''
    Reads before filtering: {barcodes_total['count'].to_numpy().sum()}, \n
    Reads after filtering by count: {barcodes_total.query('count > 10000')['count'].to_numpy().sum()}, \n
    Reads lost after filtering by count: {barcodes_total.query('count <= 10000')['count'].to_numpy().sum()}, \n
    Reads lost after filtering by M: {barcode2M_stat.query('count>10000').query('M<=0.75')['count'].to_numpy().sum()}, \n
    Reads left after filtering by count and M: {'{}%'.format(np.round(barcode2M_stat.query('count>10000').query('M>0.75')['count'].to_numpy().sum()/barcode2M_stat.query('count>10000')['count'].to_numpy().sum()*100, 2))}
    ''')
    print('*********** END ***********\n')

    print('\n*********** MAKING PLOT: cell loss ***********')

    MaxM = barcode2M_stat['M'].max()
    MaxCount = barcode2M_stat['count'].max()
    plt.scatter(x=barcode2M_stat.M,  y=barcode2M_stat['count'], s=1000/len(barcode2M_stat))
    plt.hlines(y=count_threshold, xmin=M_threshold, xmax=MaxM, colors='r')
    plt.vlines(x=M_threshold, ymin=count_threshold, ymax=MaxCount, colors='r')
    plt.yscale('log')
    plt.ylabel('count')
    plt.xlabel('M')
    plt.savefig(os.path.join(args.out_path, "scatter_M_Count.png"), dpi=600, bbox_inches='tight', format='png')

    print('*********** DONE ***********\n')

    print('Saving a list of valid barcodes...')


    barcode2M_stat.query('count > @count_threshold').query('M > @M_threshold')['barcode'].to_csv(os.path.join(args.out_path, 'selected_barcodes.csv'), index=None)
=======

barcode2count = barcodes_total[barcodes_total['count']>=args.count_thr]
assert barcode2count['count'].min() >= args.count_thr

barcode2M_stat = pd.read_csv(args.m_stat)[['barcode', 'M']]
count_threshold = args.count_thr
M_threshold = args.m_thr
barcode2count = barcode2count.query('count > @count_threshold')
barcode2M_stat['count'] = barcode2M_stat['barcode'].apply(lambda barcode: barcode2count[barcode2count['barcodes'] == barcode]['count'].item())

print('\n*********** STATS: cell loss ***********')

print(f'''Barcodes before filtration by count: {len(barcodes_total)}, \n
Barcodes after filtration by count: {len(barcode2M_stat.query('count > @count_threshold'))}, \n
Impact of filtration by count: {'{}%'.format(np.round(len(barcodes_total.query('count <= @count_threshold'))/len(barcodes_total)*100, 2))}, \n
Barcodes after filtreation by count and M: {len(barcode2M_stat.query('count > @count_threshold').query('M > @M_threshold'))}, \n
Impact of filtration by M: {'{}%'.format(np.round(len(barcode2M_stat.query('count > @count_threshold').query('M <= @M_threshold'))/len(barcode2M_stat.query('count > @count_threshold'))*100, 2))}''')

print('\n*********** STATS: read loss ***********')
print(f'''
Reads before filtering: {barcodes_total['count'].to_numpy().sum()}, \n
Reads after filtering by count: {barcodes_total.query('count > 10000')['count'].to_numpy().sum()}, \n
Reads lost after filtering by count: {barcodes_total.query('count <= 10000')['count'].to_numpy().sum()}, \n
Reads lost after filtering by M: {barcode2M_stat.query('count>10000').query('M<=0.75')['count'].to_numpy().sum()}, \n
Reads left after filtering by count and M: {'{}%'.format(np.round(barcode2M_stat.query('count>10000').query('M>0.75')['count'].to_numpy().sum()/barcode2M_stat.query('count>10000')['count'].to_numpy().sum()*100, 2))}
''')
print('*********** END ***********\n')

print('\n*********** MAKING PLOT: cell loss ***********')

MaxM = barcode2M_stat['M'].max()
MaxCount = barcode2M_stat['count'].max()
plt.scatter(x=barcode2M_stat.M,  y=barcode2M_stat['count'], s=1000/len(barcode2M_stat))
plt.hlines(y=count_threshold, xmin=M_threshold, xmax=MaxM, colors='r')
plt.vlines(x=M_threshold, ymin=count_threshold, ymax=MaxCount, colors='r')
plt.yscale('log')
plt.ylabel('count')
plt.xlabel('M')
plt.savefig(os.path.join(args.out_path, "scatter_M_Count.png"), dpi=600, bbox_inches='tight', format='png')

print('*********** DONE ***********\n')

print('Saving a list of valid barcodes...')


barcode2M_stat.query('count > @count_threshold').query('M > @M_threshold')['barcode'].to_csv(os.path.join(args.out_path, 'selected_barcodes.csv'), index=None)
>>>>>>> refs/remotes/origin/main
