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
parser.add_argument('--neff', default='all', type=str, help='Specify Neff values to compare performance.\n Options: no_hit, low, medium, high, no_hit_low_medium, all, all_combined; default = ''all''\n', metavar='')
parser.add_argument('--dca_method', default='gremlin', type=str, help='Specify dca predictor.\n Options: gremlin, plmc, mfdca, plmdca; default = ''gremlin''\n', metavar='')
parser.add_argument('--RFAM', default=-1, type=int, help='set this to 1 for performance on RFAM mapped RNAs only and to 0 non-RFAM mapped RNAs; default = ''-1''\n', metavar='')
parser.add_argument('--exclude_tRNA', default=0, type=int, help='set this to 1 for performance on RFAM mapped RNAs by excluding tRNAs; default = ''0''\n', metavar='')
parser.add_argument('--figure', default=0, type=int, help='Specify figure number.\n Options: 2, 4; default = ''0''\n', metavar='')
args = parser.parse_args()

colors = ['green', 'red', 'yellow', 'blue', 'cyan', 'white']*20

msa_list = ['blastn', 'direct_infernal', 'RNAcmap', 'RNAcmap_meta', 'RNAcmap2_meta'] 

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
with open(base_path + '/datasets/medium_neff_ids') as f:
	medium_neff_ids = f.read().splitlines()
with open(base_path + '/datasets/low_neff_ids') as f:
	low_neff_ids = f.read().splitlines()
with open(base_path + '/datasets/no_hit_ids') as f:
    no_hit_ids = f.read().splitlines()

with open(base_path + '/datasets/rnacmap2_nr_spotrna_ids') as f:
    rnacmap2_nr_spotrna_ids = f.read().splitlines()

no_hit_rnacmap_meta_ids = ['1ykq_B','3iab_R','3nmu_D','3npn_A','6d12_C','6e8u_B','6wpi_B','5m0h_A']

high_neff_spotrna = ['3zd4_A','5t83_A','5x2h_B','4xwf_A','4qyz_L','3a3a_A','3g8t_P','4meh_B','4tna_A','4v8n_AV','4ybb_CB','4yvi_C','5ibb_1K','6pmo_A','2czj_B','5ztm_C','4qln_A', '1xbp_9','6ck5_A','5ml7_A','4wzd_2K','5o6u_A','3bo2_B','1u9s_A','3ndb_M','4kqy_A','6dlq_A','6dmd_A','6ufm_B','4yaz_A','4v9k_AW','4wj4_B','4v9i_AY','1kxk_A','3eph_E','6p2h_A', '3diq_A','4lck_C','6n2v_A','6vmy_A','6jxm_B','4v88_A4']  # last id 4v88_A4 is low neff but more than 50k sequences in 1st round of search.

if args.neff == 'high':
	all_ids = [high_neff_ids]
elif args.neff == 'medium':
	all_ids = [medium_neff_ids]
elif args.neff == 'low':
	all_ids = [low_neff_ids]
elif args.neff == 'no_hit':
	all_ids = [no_hit_ids]
elif args.neff == 'no_hit_low_medium':
        all_ids = [no_hit_ids + low_neff_ids + medium_neff_ids]
elif args.neff == 'all':
	all_ids = [no_hit_ids, low_neff_ids, medium_neff_ids, high_neff_ids]
elif args.neff == 'all_combined':
	all_ids = [no_hit_ids + low_neff_ids + medium_neff_ids + high_neff_ids]

with open(base_path + '/datasets/rfam_mapped_ids') as f:
    rfam_mapped_ids = f.read().splitlines()
rfam_families = [i.split(' ')[1] for i in rfam_mapped_ids]


#### unique rfam families ######
#keep_rfam = []
#for i in rfam_families:
#	if i not in keep_rfam:
#		keep_rfam.append(i)
#print(len(rfam_families), len(keep_rfam))


if args.exclude_tRNA == 0:
	rfam_mapped_ids = [i.split(' ')[0] for i in rfam_mapped_ids]
elif args.exclude_tRNA == 1:
	rfam_mapped_ids = [i.split(' ')[0] for i in rfam_mapped_ids if i.split(' ')[1]!='RF00005']


if args.RFAM == 1:
	all_ids = [[i for i in ids if i in rfam_mapped_ids] for ids in all_ids]
	msa_list = ['RNAcmap', 'RFAM', 'RNAcmap2_meta'] 
elif args.RFAM == 0:
	all_ids = [[i for i in ids if i not in rfam_mapped_ids]]

#print(len(all_ids))

color_count = 0
save_metrics_all = {}
watson_pairs_dic = {}
wobble_pairs_dic = {}
other_pairs_dic = {}

for ids in all_ids:

	if len(ids)==21: neff_type = 'no_hit'; print('\n\n \t \t \t \t \t No-hit RNAs\n')
	elif len(ids)==83: neff_type = 'low'; print('\n\n \t \t \t \t \t Low Neff RNAs\n')
	elif len(ids)==31: neff_type = 'medium'; print('\n\n \t \t \t \t \t Medium Neff RNAs\n')
	elif len(ids)==110: neff_type = 'high'; print('\n\n \t \t \t \t \t High Neff RNAs\n')
	elif len(ids)==245: neff_type = 'all'; print('\n\n \t \t \t \t \t All Neff RNAs\n')

	for msa_no, msa_type in enumerate(msa_list[0:]):

		count = 0; save_all_bps = []; save_neff = []; save_all_metrics = []

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
				if k not in no_hit_ids:
					pred_pair, dca = dca_pred(k, seq, path=base_path + '/predictions/' + args.dca_method + '/rnacmap/')
					with open(base_path + '/predictions/gremlin/' + msa_type.lower() + '/' + str(k) + '.neff') as f:
						temp_2 = pd.read_csv(f, delim_whitespace=True, header=None).values  
					neff = temp_2[0][2]  
				else: pred_pair = []; neff = 0


			elif msa_type == 'RNAcmap_meta':
				if k in no_hit_rnacmap_meta_ids: # no hit in rnacmap_meta ids
					pred_pair = [];	neff = 0
				else: 
					pred_pair, dca = dca_pred(k, seq, path=base_path + '/predictions/' + args.dca_method + '/rnacmap_meta/')
					with open(base_path + '/predictions/gremlin/' + msa_type.lower() + '/' + str(k) + '.neff') as f:
						temp_2 = pd.read_csv(f, delim_whitespace=True, header=None).values  
					neff = temp_2[0][2]  

			elif msa_type == 'RNAcmap_meta_spotrna':
				if k in ['1ykq_B', '3nmu_D', '6wpi_B']: # no hit in rnacmap_meta_spotrna ids
					pred_pair = [];	neff = 0
				else: 
					pred_pair, dca = dca_pred(k, seq, path=base_path + '/predictions/' + args.dca_method + '/rnacmap_meta_spotrna/')
					with open(base_path + '/predictions/gremlin/' + msa_type.lower() + '/' + str(k) + '.neff') as f:
						temp_2 = pd.read_csv(f, delim_whitespace=True, header=None).values  
					neff = temp_2[0][2]  
	#			if neff > 50:
	#				print(k)

			elif msa_type == 'RNAcmap2_meta_spotrna':
				if k in ['1ykq_B', '3nmu_D', '6wpi_B']: # no hit in rnacmap_meta_spotrna ids
					pred_pair = [];	neff = 0
				elif k in high_neff_spotrna: 
					pred_pair, dca = dca_pred(k, seq, path=base_path + '/predictions/' + args.dca_method + '/rnacmap_meta_spotrna/')
					with open(base_path + '/predictions/gremlin/rnacmap_meta_spotrna/' + str(k) + '.neff') as f:
						temp_2 = pd.read_csv(f, delim_whitespace=True, header=None).values  
					neff = temp_2[0][2]  
				else: 
					pred_pair, dca = dca_pred(k, seq, path=base_path + '/predictions/' + args.dca_method + '/rnacmap2_meta_spotrna/')
					with open(base_path + '/predictions/gremlin/' + msa_type.lower() + '/' + str(k) + '.neff') as f:
						temp_2 = pd.read_csv(f, delim_whitespace=True, header=None).values  
					neff = temp_2[0][2]  


			elif msa_type == 'RNAcmap2_meta':
				if k in no_hit_rnacmap_meta_ids: # no hit in rnacmap_meta ids
					pred_pair = [];	neff = 0
				elif k in ['2du4_C','2pxb_B','3adb_C','2zni_C','2qus_A','1wz2_C'] + high_neff_ids and k not in ['4qjd_B']: # high neff ids in rnacmap_meta pipeline
						pred_pair, dca = dca_pred(k, seq, path=base_path + '/predictions/' + args.dca_method + '/rnacmap_meta/')
						with open(base_path + '/predictions/gremlin/rnacmap_meta' + '/' + str(k) + '.neff') as f:
							temp_2 = pd.read_csv(f, delim_whitespace=True, header=None).values  
						neff = temp_2[0][2]  
				else:
					pred_pair, dca = dca_pred(k, seq, path=base_path + '/predictions/' + args.dca_method + '/rnacmap2_meta/')
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

			save_neff.append(neff)

			count += 1

		all_metrics = np.mean(save_all_bps, axis=0)        

		if msa_no==0: print('\t\t       F1\tPrecision\tSensitivity\tNo. of RNAs\tMedian Neff')		
		print(msa_type + ' '*(len('direct_infernal')-len(msa_type)), end='    ')
		print(colored(" {:.3f} \t {:.3f} \t\t {:.3f}\t\t    {}\t\t  {:.1f}", colors[color_count]).format(all_metrics[0], all_metrics[1], all_metrics[2], count, np.median(save_neff)), end='    ')

		color_count = color_count + 1

		header_list = ['RNA_id', 'F1_all_bps', 'PR_all_bps', 'SN_all_bps', 'Neff', 'Neff_by_L', 'Seq_Length']
		df = pd.DataFrame(np.array(save_all_metrics), columns=header_list)

		if args.figure == 2: df.to_csv('./figure_2_data/' + msa_type + '_' + neff_type + '.csv', index=False)
		elif args.figure == 4: df.to_csv('./figure_4_data/' + msa_type + '_' + neff_type + '.csv', index=False)

		print()

