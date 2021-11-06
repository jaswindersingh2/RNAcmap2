#!/usr/bin/env python

import pandas as pd
import numpy as np
import pickle as pkl
import os
from pathlib import Path
from utils import *
from termcolor import colored
import argparse

base_path = os.path.dirname(os.path.realpath(__file__))
base_path = str(Path(base_path).parent.absolute())

parser = argparse.ArgumentParser()
parser.add_argument('--neff', default='all', type=str, help='Specify Neff values to compare performance.\n Options: no_hit, low, medium, high, no_hit_low_medium, all; default = ''all''\n', metavar='')
parser.add_argument('--dca_method', default='gremlin', type=str, help='Specify dca predictor.\n Options: gremlin, plmc, mfdca, plmdca; default = ''gremlin''\n', metavar='')
parser.add_argument('--RFAM', default=-1, type=int, help='set this to 1 for performance on RFAM mapped RNAs only and to 0 non-RFAM mapped RNAs; default = ''-1''\n', metavar='')
parser.add_argument('--exclude_tRNA', default=0, type=int, help='set this to 1 for performance on RFAM mapped RNAs by excluding tRNAs; default = ''0''\n', metavar='')
args = parser.parse_args()

colors = ['green', 'red', 'yellow', 'blue', 'cyan', 'white']*20

msa_list = ['RNAcmap', 'RNAcmap2_meta'] 

with open(base_path + '/benchmarking/true_labels_pdb', 'rb') as f:
    true_labels_pdb = pkl.load(f)

with open(base_path + '/benchmarking/add_true_labels_pdb2', 'rb') as f:
    add_true_labels_pdb = pkl.load(f)

with open(base_path + '/benchmarking/save_missing_nts_all', 'rb') as f:
    save_missing_nts_all = pkl.load(f)

path_seq = base_path + '/datasets/sequences/'

if args.dca_method == 'gremlin':
	dca_pred = GREMLIN
elif args.dca_method == 'plmc':
	dca_pred = plmc
elif args.dca_method == 'mfdca':
	dca_pred = mfdca_top_L_by_n
elif args.dca_method == 'plmdca':
	dca_pred = plmdca

with open(base_path + '/datasets/high_neff_ids') as f:
	high_neff_ids = f.read().splitlines()
with open(base_path + '/datasets/medium_neff_ids') as f:
	medium_neff_ids = f.read().splitlines()
with open(base_path + '/datasets/low_neff_ids') as f:
	low_neff_ids = f.read().splitlines()
with open(base_path + '/datasets/no_hit_ids') as f:
    no_hit_ids = f.read().splitlines()

with open(base_path + '/datasets/rnacmap2_nr_spotrna_ids') as f:
    rnacmap2_nr_spotrna_ids = f.read().splitlines()

no_hit_rnacmap_meta_ids = ['1ykq_B','3iab_R','3nmu_D','3npn_A','6d12_C','6e8u_B','6wpi_B','5m0h_A']

if args.neff == 'high':
	ids = high_neff_ids
elif args.neff == 'medium':
	ids = medium_neff_ids
elif args.neff == 'low':
	ids = low_neff_ids
elif args.neff == 'no_hit':
	ids = no_hit_ids
elif args.neff == 'no_hit_low_medium':
        ids = no_hit_ids + low_neff_ids + medium_neff_ids
elif args.neff == 'all':
	ids = no_hit_ids + low_neff_ids + medium_neff_ids + high_neff_ids

with open(base_path + '/datasets/rfam_mapped_ids') as f:
    rfam_mapped_ids = f.read().splitlines()
rfam_families = [i.split(' ')[1] for i in rfam_mapped_ids]


if args.exclude_tRNA == 0:
	rfam_mapped_ids = [i.split(' ')[0] for i in rfam_mapped_ids]
elif args.exclude_tRNA == 1:
	rfam_mapped_ids = [i.split(' ')[0] for i in rfam_mapped_ids if i.split(' ')[1]!='RF00005']


if args.RFAM == 1:
	ids = [i for i in ids if i in rfam_mapped_ids]
	msa_list = ['RNAcmap', 'RFAM', 'RNAcmap2_meta'] 
elif args.RFAM == 0:
	ids = [i for i in ids if i not in rfam_mapped_ids]



#print(len(ids))


color_count = 0
save_metrics_all = {}
watson_pairs_dic = {}
wobble_pairs_dic = {}
other_pairs_dic = {}

for n in [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]:

	print('\n\n \t \t \t \t \t Top L/{} contacts \n'.format(n))

	for msa_no, msa_type in enumerate(msa_list[0:]):

		count = 0; save_all_bps = []; save_neff = []
		save_all_metrics = []
		for K,k in enumerate(ids[0:]):
			#print(k)
			with open(path_seq + str(k)) as f:
				temp_1 = pd.read_csv(f, comment='#', delim_whitespace=True, header=None, skiprows=[0]).values
			seq = [j.upper() for j in temp_1[0, 0]]

			try: true_pairs = true_labels_pdb[k]  
			except: true_pairs = add_true_labels_pdb[k]

			watson_pairs_dic[k], wobble_pairs_dic[k], other_pairs_dic[k] = type_pairs(true_pairs, seq)
			true_pairs = watson_pairs_dic[k] + wobble_pairs_dic[k] + other_pairs_dic[k]
			true_pairs = [i for i in true_pairs if abs(i[0]-i[1]) > 3]  # non-local base-pairs


			if msa_type == 'RNAcmap':
				if k not in no_hit_ids:
					pred_pair, dca = dca_pred(k, seq, n, path=base_path + '/predictions/' + args.dca_method + '/rnacmap/')
					with open(base_path + '/predictions/gremlin/' + msa_type.lower() + '/' + str(k) + '.neff') as f:
						temp_2 = pd.read_csv(f, delim_whitespace=True, header=None).values  
					neff = temp_2[0][2]  
				else: pred_pair = []; neff = 0


			elif msa_type == 'RNAcmap2_meta':
				if k in no_hit_rnacmap_meta_ids: # no hit in rnacmap_meta ids
					pred_pair = [];	neff = 0
				elif k in ['2du4_C','2pxb_B','3adb_C','2zni_C','2qus_A','1wz2_C'] + high_neff_ids and k not in ['4qjd_B']: # high neff ids in rnacmap_meta pipeline
						pred_pair, dca = dca_pred(k, seq, n, path=base_path + '/predictions/' + args.dca_method + '/rnacmap_meta/')
						with open(base_path + '/predictions/gremlin/rnacmap_meta' + '/' + str(k) + '.neff') as f:
							temp_2 = pd.read_csv(f, delim_whitespace=True, header=None).values  
						neff = temp_2[0][2]  
				else:
					pred_pair, dca = dca_pred(k, seq, n, path=base_path + '/predictions/' + args.dca_method + '/rnacmap2_meta/')
					with open(base_path + '/predictions/gremlin/' + msa_type.lower() + '/' + str(k) + '.neff') as f:
						temp_2 = pd.read_csv(f, delim_whitespace=True, header=None).values  
					neff = temp_2[0][2]  
					#print(k)

			missing_nts =  save_missing_nts_all[k]

			pred_pair = [i for i in pred_pair if i[0] not in missing_nts and i[1] not in missing_nts]

			pred_correctly = [i for i in pred_pair if i in true_pairs]; #print(pred_correctly)
			pred_wrongly = [i for i in pred_pair if i not in true_pairs]; #print(pred_wrongly)

			tp = 0;fp = 0;fn = 0
			correct_pairs = []
			for i,I in enumerate(pred_pair):
				if I in true_pairs:
				    tp +=1
				    correct_pairs.append(I)
				elif I not in true_pairs:
				    fp += 1
				    # print(tp)
			for i,I in enumerate(true_pairs):
				if I not in pred_pair:
				    fn += 1
			tn = (len(seq))*((len(seq)-1)/2) - tp - fp - fn

			try:
				pre = tp / (tp + fp)
				sen = tp / (tp + fn)
				f1 = 2*((pre*sen)/(pre + sen))
				#with np.errstate(invalid='ignore'):
				mcc = ((tp * tn) - (fp * fn)) / np.sqrt(np.float64((tp + fp) * (tp + fn) * (tn + fn) * (tn + fp)))
			except:
				pre = 0
				sen = 0
				f1 = 0
				mcc = 0; #print(k)
			save_all_bps.append([f1, pre, sen])

			save_all_metrics.append([k, f1, pre, sen, neff, neff/len(seq), len(seq)])
			#print(k, f1)

			save_neff.append(neff)

			count += 1

		all_metrics = np.mean(save_all_bps, axis=0)        

		if msa_no==0: print('\t\t       F1\tPrecision\tSensitivity\tNo. of RNAs\tMedian Neff')		
		print(msa_type + ' '*(len('direct_infernal')-len(msa_type)), end='    ')
		print(colored(" {:.3f} \t {:.3f} \t\t {:.3f}\t\t    {}\t\t  {:.1f}", colors[color_count]).format(all_metrics[0], all_metrics[1], all_metrics[2], count, np.median(save_neff)), end='    ')

		color_count = color_count + 1
		print()

		header_list = ['RNA_id', 'F1_all_bps', 'PR_all_bps', 'SN_all_bps', 'Neff', 'Neff_by_L', 'Seq_Length']

		df = pd.DataFrame(np.array(save_all_metrics), columns=header_list)
		df.to_csv('./figure_3_data/' + msa_type + '_' + args.dca_method + '_top_L_by_' + str(n) + '.csv', index=False)

