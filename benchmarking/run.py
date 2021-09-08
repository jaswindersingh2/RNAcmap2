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
parser.add_argument('--neff', default='low', type=str, help='name of RNA sequence file; default = ''6p2h_A''\n', metavar='')
parser.add_argument('--dca_method', default='gremlin', type=str, help='file consists of list name of RNA sequence files for batch prediction; default = ''datasets/TS1_ids''\n', metavar='')
args = parser.parse_args()

colors = ['green', 'red', 'yellow', 'blue', 'cyan', 'white']*20

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
	dca_pred = mfdca
elif args.dca_method == 'plmdca':
	dca_pred = plmdca

with open(base_path + '/datasets/high_neff_ids') as f:
	high_neff_ids = f.read().splitlines()
with open(base_path + '/datasets/median_neff_ids') as f:
	median_neff_ids = f.read().splitlines()
with open(base_path + '/datasets/low_neff_ids') as f:
	low_neff_ids = f.read().splitlines()


if args.neff == 'high':
	ids = high_neff_ids
elif args.neff == 'median':
	ids = median_neff_ids
elif args.neff == 'low':
	ids = low_neff_ids
elif args.neff == 'all':
	ids = low_neff_ids + median_neff_ids + high_neff_ids

neff_1_ids = ['1ddy_A', '1dk1_B', '1s03_A', '1ykq_B', '2gje_R', '2oiu_P', '2qwy_A', '2zy6_A', '3iab_R', '3nmu_D', '3npn_A', '3x1l_I', '4frn_A', '4k4w_B', '4pjo_1', '4q9q_R', '4r4v_A', '4rzd_A', '4ts0_X', '5bjp_E', '5de5_A', '5gip_G', '5ng6_B', '5ob3_A', '5vof_A', '5wlh_B', '5wti_B', '5xuz_B', '5y85_B', '6d12_C', '6d3p_A', '6e8u_B', '6g7z_B', '6ifn_N', '6iv8_D', '6vff_C', '6wpi_B']

#ids = low_neff_ids + median_neff_ids #+ high_neff_ids

msa_list = ['blastn', 'direct_infernal', 'RNAcmap', 'RNAcmap2', 'RNAcmap_meta', 'RNAcmap2_meta'] 

#print(len(ids))


color_count = 0
save_metrics_all = {}
watson_pairs_dic = {}
wobble_pairs_dic = {}
other_pairs_dic = {}


for msa_no, msa_type in enumerate(msa_list[0:]):

	count = 0; save_all_bps = []; save_neff = []

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

		if msa_type == 'blastn':
			try: 
				pred_pair, dca = dca_pred(k, seq, path=base_path + '/predictions/' + args.dca_method + '/blastn/')   

				with open(base_path + '/predictions/gremlin/' + msa_type.lower() + '/' + str(k) + '.neff') as f:
					temp_2 = pd.read_csv(f, delim_whitespace=True, header=None).values  
				neff = temp_2[0][2]   
			except: pred_pair = []; neff = 0
		elif msa_type == 'direct_infernal':
			try: 
				pred_pair, dca = dca_pred(k, seq, path=base_path + '/predictions/' + args.dca_method + '/direct_infernal/')   

				with open(base_path + '/predictions/gremlin/' + msa_type.lower() + '/' + str(k) + '.neff') as f:
					temp_2 = pd.read_csv(f, delim_whitespace=True, header=None).values  
				neff = temp_2[0][2]   
			except: pred_pair = []; neff = 0
		elif msa_type == 'RFAM':
			pred_pair, dca = dca_pred(k, seq, path=base_path + '/predictions/' + args.dca_method + '/RFAM/')   
			with open(base_path + '/predictions/gremlin/' + msa_type + '/' + str(k) + '.neff') as f:
				temp_2 = pd.read_csv(f, delim_whitespace=True, header=None).values  
			neff = temp_2[0][2]  

		elif msa_type == 'RNAcmap':
			pred_pair, dca = dca_pred(k, seq, path=base_path + '/predictions/' + args.dca_method + '/rnacmap/')
			with open(base_path + '/predictions/gremlin/' + msa_type.lower() + '/' + str(k) + '.neff') as f:
				temp_2 = pd.read_csv(f, delim_whitespace=True, header=None).values  
			neff = temp_2[0][2]  
		elif msa_type == 'RNAcmap2':
			if k in neff_1_ids + low_neff_ids + median_neff_ids:
				pred_pair, dca = dca_pred(k, seq, path=base_path + '/predictions/' + args.dca_method + '/rnacmap2/')

				with open(base_path + '/predictions/gremlin/' + msa_type.lower() + '/' + str(k) + '.neff') as f:
					temp_2 = pd.read_csv(f, delim_whitespace=True, header=None).values  
				neff = temp_2[0][2]  
			elif k in high_neff_ids:
				pred_pair, dca = dca_pred(k, seq, path=base_path + '/predictions/' + args.dca_method + '/rnacmap/')
				with open(base_path + '/predictions/gremlin/rnacmap/' + str(k) + '.neff') as f:
					temp_2 = pd.read_csv(f, delim_whitespace=True, header=None).values  
				neff = temp_2[0][2]  

		elif msa_type == 'RNAcmap_meta':
			if k not in ['5m0h_A']: # 5m0h no hit for rnacmap_meta but msa exits for rnacmap
				pred_pair, dca = dca_pred(k, seq, path=base_path + '/predictions/' + args.dca_method + '/rnacmap_meta/')

				with open(base_path + '/predictions/gremlin/' + msa_type.lower() + '/' + str(k) + '.neff') as f:
					temp_2 = pd.read_csv(f, delim_whitespace=True, header=None).values  
				neff = temp_2[0][2]  
			else: pred_pair = []; neff = 0

		elif msa_type == 'RNAcmap2_meta':
			if k not in ['5m0h_A']: # 5m0h no hit for rnacmap_meta but msa exits for rnacmap
				if k in neff_1_ids + low_neff_ids + median_neff_ids:
					try:
						pred_pair, dca = dca_pred(k, seq, path=base_path + '/predictions/' + args.dca_method + '/rnacmap2_meta/')
						with open(base_path + '/predictions/gremlin/' + msa_type.lower() + '/' + str(k) + '.neff') as f:
							temp_2 = pd.read_csv(f, delim_whitespace=True, header=None).values  
						neff = temp_2[0][2]  
					except:
						pred_pair, dca = dca_pred(k, seq, path=base_path + '/predictions/' + args.dca_method + '/rnacmap_meta/')
						with open(base_path + '/predictions/gremlin/rnacmap_meta/' + str(k) + '.neff') as f:
							temp_2 = pd.read_csv(f, delim_whitespace=True, header=None).values  
						neff = temp_2[0][2]  

				elif k in high_neff_ids:
				    pred_pair, dca = dca_pred(k, seq, path=base_path + '/predictions/' + args.dca_method + '/rnacmap_meta/')

				    with open(base_path + '/predictions/gremlin/rnacmap_meta/' + str(k) + '.neff') as f:
					    temp_2 = pd.read_csv(f, delim_whitespace=True, header=None).values  
				    neff = temp_2[0][2]  
			else: pred_pair = []

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

		save_neff.append(neff)

		count += 1

	all_metrics = np.mean(save_all_bps, axis=0)        

	if msa_no==0: print('\t\t       F1\tPrecision\tSensitivity\tNo. of RNAs\tMedian Neff')		
	print(msa_type + ' '*(len('direct_infernal')-len(msa_type)), end='    ')
	print(colored(" {:.3f} \t {:.3f} \t\t {:.3f}\t\t    {}\t\t  {:.1f}", colors[color_count]).format(all_metrics[0], all_metrics[1], all_metrics[2], count, np.median(save_neff)), end='    ')

	color_count = color_count + 1
	print()

