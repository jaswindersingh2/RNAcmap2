import pandas as pd
import numpy as np
import pickle as pkl
import os
#from scipy.io import savemat
#import scipy.stats
#from result_analysis import count_pairs, type_pairs, lone_pair


def type_pairs(pairs, sequence):
    sequence = [i.upper() for i in sequence]
    # seq_pairs = [[sequence[i[0]],sequence[i[1]]] for i in pairs]


    AU_pair = []
    GC_pair = []
    GU_pair = []
    other_pairs = []
    for i in pairs:
        if [sequence[i[0]],sequence[i[1]]] in [["A","U"], ["U","A"]]:
            AU_pair.append(i)
        elif [sequence[i[0]],sequence[i[1]]] in [["G","C"], ["C","G"]]:
            GC_pair.append(i)
        elif [sequence[i[0]],sequence[i[1]]] in [["G","U"], ["U","G"]]:
            GU_pair.append(i)
        else:
            other_pairs.append(i)
    watson_pairs_t = AU_pair + GC_pair
    wobble_pairs_t = GU_pair
    other_pairs_t = other_pairs
        # print(watson_pairs_t, wobble_pairs_t, other_pairs_t)
    return watson_pairs_t, wobble_pairs_t, other_pairs_t


def count_pairs(pairs, sequence):
    sequence = [i.upper() for i in sequence]
    seq_pairs = [[sequence[i[0]],sequence[i[1]]] for i in pairs]

    AU_pair = 0
    GC_pair = 0
    GU_pair = 0
    other_pairs = 0
    for i in seq_pairs:
        if i in [["A","U"], ["U","A"]]:
            AU_pair += 1
        elif i in [["G","C"], ["C","G"]]:
            GC_pair += 1
        elif i in [["G","U"], ["U","G"]]:
            GU_pair += 1
        else:
            # print(i)
            other_pairs += 1
    return AU_pair, GC_pair, GU_pair, other_pairs

def replacer(s, newstring, index, nofail=False):
    # raise an error if index is outside of the string
    if not nofail and index not in range(len(s)):
        raise ValueError("index outside given string")

    # if not erroring, but the index is still not in the correct range..
    if index < 0:  # add it to the beginning
        return newstring + s
    if index > len(s):  # add it to the end
        return s + newstring

    # insert the new string between "slices" of the original
    return s[:index] + newstring + s[index + 1:]

def get_pairs(labels):
    pairs = []
    stack = []
    stack_2 = []
    stack_3 = []
    stack_4 = []
    stack_5 = []
    stack_6 = []
    stack_7 = []
    for i, I in enumerate(labels):
        if I == '(':
            stack.append(i)
            #print(I, i)
        elif I == ')':
            pairs.append(sorted([stack[-1], i]))
            #print(I, i)
            del stack[-1]
        elif I == '<':
            stack_2.append(i)
        elif I == '>':
            pairs.append(sorted([stack_2[-1], i]))
            del stack_2[-1]
        elif I == '{':
            stack_3.append(i)
        elif I == '}':
            pairs.append(sorted([stack_3[-1], i]))
            del stack_3[-1]
        elif I == '[':
            stack_4.append(i)
        elif I == ']':
            pairs.append(sorted([stack_4[-1], i]))
            del stack_4[-1]
        elif I == 'A':
            stack_5.append(i)
        elif I == 'a':
            pairs.append(sorted([stack_5[-1], i]))
            del stack_5[-1]
        elif I == 'B':
            stack_6.append(i)
        elif I == 'b':
            pairs.append(sorted([stack_6[-1], i]))
            del stack_6[-1]
        elif I == 'C':
            stack_7.append(i)
        elif I == 'c':
            pairs.append(sorted([stack_7[-1], i]))
            del stack_7[-1]
        elif I == '.' or I == '-':
            continue
        else:
            print(I)
    return pairs

def get_pairs_pseudo_data(labels):
    pairs = []
    stack = []
    stack_2 = []
    stack_3 = []
    stack_4 = []
    stack_5 = []
    stack_6 = []
    stack_7 = []
    for i, I in enumerate(labels):
        if I == '(':
            stack.append(i)
            #print(I, i)
        elif I == ')':
            pairs.append(sorted([stack[-1], i]))
            #print(I, i)
            del stack[-1]
        elif I == '<':
            stack_2.append(i)
        elif I == '>':
            pairs.append(sorted([stack_2[-1], i]))
            del stack_2[-1]
        elif I == '{':
            stack_3.append(i)
        elif I == '}':
            pairs.append(sorted([stack_3[-1], i]))
            del stack_3[-1]
        elif I == '[':
            stack_4.append(i)
        elif I == ']':
            pairs.append(sorted([stack_4[-1], i]))
            del stack_4[-1]
        elif I == 'A':
            stack_5.append(i)
        elif I == 'a':
            pairs.append(sorted([stack_5[-1], i]))
            del stack_5[-1]
        elif I == 'B':
            stack_6.append(i)
        elif I == 'b':
            pairs.append(sorted([stack_6[-1], i]))
            del stack_6[-1]
        elif I == 'C':
            stack_7.append(i)
        elif I == 'c':
            pairs.append(sorted([stack_7[-1], i]))
            del stack_7[-1]
        elif I == ':' or I == '-':
            continue
        else:
            print(I)
    return pairs

def ct_file_output(pairs, seq, id, save_result_path):

    col1 = np.arange(1, len(seq) + 1, 1)
    col2 = np.array([i for i in seq])
    col3 = np.arange(0, len(seq), 1)
    col4 = np.append(np.delete(col1, 0), [0])
    col5 = np.zeros(len(seq), dtype=int)

    for i, I in enumerate(pairs):
        col5[I[0]] = int(I[1]) + 1
        col5[I[1]] = int(I[0]) + 1
    col6 = np.arange(1, len(seq) + 1, 1)
    temp = np.vstack((np.char.mod('%d', col1), col2, np.char.mod('%d', col3), np.char.mod('%d', col4),
                      np.char.mod('%d', col5), np.char.mod('%d', col6))).T
    #os.chdir(save_result_path)
    #print(os.path.join(save_result_path, str(id[0:-1]))+'.spotrna')
    np.savetxt(os.path.join(save_result_path, str(id))+'.ct', (temp), delimiter='\t\t', fmt="%s", header=str(len(seq)) + '\t\t' + str(id) + '\t\t' + 'SPOT-RNA output\n' , comments='')

    return

# ----------------------- find multiplets pairs--------------------------------#
def flatten(x):
    result = []
    for el in x:
        if hasattr(el, "__iter__") and not isinstance(el, str):
            result.extend(flatten(el))
        else:
            result.append(el)
    return result

def multiplets_pairs(pred_pairs):

    pred_pair = [i[:2] for i in pred_pairs]
    temp_list = flatten(pred_pair)
    temp_list.sort()
    new_list = sorted(set(temp_list))
    dup_list = []
    for i in range(len(new_list)):
        if (temp_list.count(new_list[i]) > 1):
            dup_list.append(new_list[i])

    dub_pairs = []
    for e in pred_pair:
        if e[0] in dup_list:
            dub_pairs.append(e)
        elif e[1] in dup_list:
            dub_pairs.append(e)

    temp3 = []
    for i in dup_list:
        temp4 = []
        for k in dub_pairs:
            if i in k:
                temp4.append(k)
        temp3.append(temp4)
        
    return temp3

def multiplets_free_bp(pred_pairs, dca, seq_len):
    L = len(pred_pairs)

    y_pred = np.zeros((seq_len, seq_len))
    for i in dca:
        y_pred[int(i[0]-1), int(i[1]-1)] = i[2]
#        print(i)

    multiplets_bp = multiplets_pairs(pred_pairs)
    save_multiplets = []
    while len(multiplets_bp) > 0:
        remove_pairs = []
        for i in multiplets_bp:
            save_prob = []
            for j in i:
                save_prob.append(y_pred[j[0], j[1]])
            remove_pairs.append(i[save_prob.index(min(save_prob))])
            save_multiplets.append(i[save_prob.index(min(save_prob))])
        pred_pairs = [k for k in pred_pairs if k not in remove_pairs]
        multiplets_bp = multiplets_pairs(pred_pairs)
    save_multiplets = [list(x) for x in set(tuple(x) for x in save_multiplets)]
    assert L == len(pred_pairs)+len(save_multiplets)
    #print(L, len(pred_pairs), save_multiplets)

    return pred_pairs, save_multiplets

########--------------------- parse aliFreeFold output ---------------------------- #########################
def parse_aliFreefold(k, seq):

    with open('/home/jaswinder/mnt/spl/project5/predictions/alifreefold/output_200/' + str(k) + '/rep.db') as f:
        temp = pd.read_csv(f, comment='#', delim_whitespace=True, header=None, skiprows=[0,1]).values    
    labels = [i for i in temp[0][0]]
    assert len(labels) == len(seq)
    pred_pair = get_pairs(labels)

    return pred_pair

########--------------------- parse TurboFold-II output ---------------------------- #########################
def parse_Turbofold(k, seq):

    with open('/home/jaswinder/mnt/spl/project5/predictions/TurboFold2/output_200_MEA/' + str(k) + '/' + str(k) + '.ct') as f:
        temp = pd.read_csv(f, comment='#', header=None, skiprows=[0], nrows=len(seq))
    temp = temp[0].str.split(expand=True,).values
    seq_pred = [i for i in temp[:, 1]]
    label_pred = [sorted(list(i)) for i in zip(temp[:, 0], temp[:, 4])]
    label_pred = [sorted([int(i[0])-1, int(i[1])-1]) for i in label_pred if '0' not in i]
    pred_pair = []
    for i in label_pred:
        if i not in pred_pair:
            pred_pair.append(i)

    with open('/home/jaswinder/mnt/spl/project5/predicted_motifs/Turbofold/' + str(k) + '.dbn') as f:
        temp = pd.read_csv(f, comment='#', delim_whitespace=True, header=None, usecols=[0], skiprows=[0]).values
    seq_pred = [i for i in temp[0,0]]
    labels = [I for i,I in enumerate(temp[1, 0])]
    assert len(seq) == len(seq_pred) == len(labels)
    #pred_pair = get_pairs(labels)

    return pred_pair


########--------------------- parse TurboFold-II-ProbKnot output ---------------------------- #########################
def parse_Turbofold_ProbKnot(k, seq):

    with open('/home/jaswinder/mnt/spl/project5/predictions/TurboFold2/output_100_ProbKnot/' + str(k) + '/' + str(k) + '.ct') as f:
        temp = pd.read_csv(f, comment='#', header=None, skiprows=[0], nrows=len(seq))
    temp = temp[0].str.split(expand=True,).values
    seq_pred = [i for i in temp[:, 1]]
    label_pred = [sorted(list(i)) for i in zip(temp[:, 0], temp[:, 4])]
    label_pred = [sorted([int(i[0])-1, int(i[1])-1]) for i in label_pred if '0' not in i]
    pred_pair = []
    for i in label_pred:
        if i not in pred_pair:
            pred_pair.append(i)

    return pred_pair

########--------------------- MXSCARNA ---------------------------- #########################
def parse_MXSCARNA(k, seq):

    with open('/home/jaswinder/mnt/spl/project5/predictions/MXSCARNA/output_500/' + str(k) + '.fasta') as f:
        temp = pd.read_csv(f, comment='#', delim_whitespace=True, header=None).values

    ref_seq = [i for i in temp[1, 0] if i!='-']
    remove_index = [i for i,I in enumerate(temp[1, 0]) if I=='-']
    dbn_ss_dash = [I for i,I in enumerate(temp[-1,0])]
    pred_pair_dash = get_pairs(dbn_ss_dash)
    for i in pred_pair_dash:
        if i[0] in remove_index or i[1] in remove_index:
            dbn_ss_dash[i[0]] = '.'
            dbn_ss_dash[i[1]] = '.'
    dbn_ss = [I for i,I in enumerate(dbn_ss_dash) if i not in remove_index]
    assert len(dbn_ss) == len(seq)        
    pred_pair = get_pairs(dbn_ss)

    return pred_pair

########--------------------- parse SPARSE output ---------------------------- #########################
def parse_SPARSE(k, seq):

    with open('/home/jaswinder/mnt/spl/project5/predictions/SPARSE/output_200/' + str(k) + '/' + str(k) + '.dbn') as f:
        temp = pd.read_csv(f, comment='#', delim_whitespace=True, header=None, skiprows=[0]).values
    
    seq_pred = [i for i in temp[0,0] if i!='-']
    remove_index = [i for i,I in enumerate(temp[0, 0]) if I=='-']
    dbn_ss_dash = [I for i,I in enumerate(temp[1,0])]
    pred_pair_dash = get_pairs(dbn_ss_dash)
    for i in pred_pair_dash:
        if i[0] in remove_index or i[1] in remove_index:
            dbn_ss_dash[i[0]] = '.'
            dbn_ss_dash[i[1]] = '.'
    dbn_ss = [I for i,I in enumerate(dbn_ss_dash) if i not in remove_index]
    assert len(seq) == len(seq_pred) == len(dbn_ss)
    pred_pair = get_pairs(dbn_ss)

    return pred_pair

######## --------------------- parse PETfold output ---------------------------- #########################
def parse_PETfold(k, seq):
    with open('/home/jaswinder/mnt/spl/project5/predictions/PETfold/output_1k/' + str(k) + '.dbn') as f:
        temp = pd.read_csv(f, comment='#', delim_whitespace=True, header=None, skiprows=[0]).values

    seq_pred = [i for i in temp[0,0] if i!='-']
    remove_index = [i for i,I in enumerate(temp[0, 0]) if I=='-']
    labels = [I for i,I in enumerate(temp[1, 0]) if i not in remove_index]
    assert len(seq) == len(seq_pred) == len(labels)
    pred_pair = get_pairs(labels)

    return pred_pair

######## --------------------- parse RNAalifold output ---------------------------- #########################
def parse_RNAalifold(k, seq):
    with open('/home/jaswinder/mnt/spl/project5/predictions/RNAalifold/output_100_MEA/' + str(k) + '.dbn') as f:
        temp = pd.read_csv(f, comment='#', delim_whitespace=True, header=None).values

    seq_pred = [i for i in temp[0,0] if i!='-']
    remove_index = [i for i,I in enumerate(temp[0, 0]) if I=='-']
    labels = [I for i,I in enumerate(temp[1, 0]) if i not in remove_index]
    assert len(seq) == len(seq_pred) == len(labels)
    pred_pair = get_pairs(labels)

    return pred_pair

######## --------------------- parse Centroid_alifold-MEA output ---------------------------- #########################
def parse_Centroid_alifold_MEA(k, seq):
    with open('/home/jaswinder/mnt/spl/project5/predictions/Centroid_alifold/output_500_MEA/' + str(k) + '.dbn') as f:
        temp = pd.read_csv(f, comment='#', delim_whitespace=True, header=None, usecols=[0], skiprows=[0]).values
    seq_pred = [i for i in temp[0,0]]
    labels = [I for i,I in enumerate(temp[1, 0])]
    assert len(seq) == len(seq_pred) == len(labels)
    pred_pair = get_pairs(labels)

    return pred_pair

######## --------------------- parse Centroid_alifold-NC output ---------------------------- #########################
def parse_Centroid_alifold_NC(k, seq):
    with open('/home/jaswinder/mnt/spl/project5/predictions/Centroid_alifold/output_1k_NC/' + str(k) + '.dbn') as f:
        temp = pd.read_csv(f, comment='#', delim_whitespace=True, header=None, usecols=[0], skiprows=[0]).values
    seq_pred = [i for i in temp[0,0]]
    labels = [I for i,I in enumerate(temp[1, 0])]
    assert len(seq) == len(seq_pred) == len(labels)
    pred_pair = get_pairs(labels)

    return pred_pair

######## --------------------- parse Centroid_alifold output ---------------------------- #########################
def CentroidaliFold_prob(k, seq):
    with open('/home/jaswinder/mnt/spl/project5/predictions/Centroid_alifold/output_100_prob/' + str(k) + '.prob') as f:
        temp = pd.read_csv(f, comment='#', header=None)

    y_pred = np.zeros((len(seq), len(seq)))
    for i in temp[:][0]:
        a = i.split(' ')
        for j in a[2:-1]:
            k = j.split(':')
            y_pred[int(a[0]) - 1, int(k[0]) - 1] = float(k[1])
    #print(y_pred)
    tri_inds = np.triu_indices(y_pred.shape[0], k=1)

    out_pred = y_pred[tri_inds]
    outputs = out_pred[:, None]
    seq_pairs = [[tri_inds[0][j], tri_inds[1][j], ''.join([seq[tri_inds[0][j]], seq[tri_inds[1][j]]])] for j in
                 range(tri_inds[0].shape[0])]

    outputs_T = np.greater_equal(outputs, 0.186)
    pred_pairs = [i for I, i in enumerate(seq_pairs) if outputs_T[I]]
    pred_pairs = [i[:2] for i in pred_pairs]

    return pred_pairs

######## --------------------- parse RNAalifold probability output ---------------------------- #########################
def RNAalifold_prob(id, seq):
    with open('/home/jaswinder/mnt/spl/project5/predictions/RNAalifold/output_100_prob/' + str(id) + '.prob') as f:
        temp = pd.read_csv(f, comment='#', header=None).values
    y_pred = np.zeros((len(seq), len(seq)))
    for i in temp[:,0]:
        a = i.split(' ')
        y_pred[int(a[0]) - 1, int(a[1]) - 1] = float(a[2])

    tri_inds = np.triu_indices(y_pred.shape[0], k=1)

    out_pred = y_pred[tri_inds]
    outputs = out_pred[:, None]
    seq_pairs = [[tri_inds[0][j], tri_inds[1][j], ''.join([seq[tri_inds[0][j]], seq[tri_inds[1][j]]])] for j in
                 range(tri_inds[0].shape[0])]

    outputs_T = np.greater_equal(outputs, 0.516)
    pred_pairs = [i for I, i in enumerate(seq_pairs) if outputs_T[I]]
    pred_pairs = [i[:2] for i in pred_pairs]

    return pred_pairs

######## --------------------- parse CentroidFold (MEA) output ---------------------------- #########################
def CentroidFold_mea(k, seq):
    #print(k)
    with open('/home/jaswinder/mnt/spl/project5/predictions_single/CentroidFold/' + str(k) + '.dbn_mea') as f:
        temp = pd.read_csv(f, comment='#', delim_whitespace=True, header=None, usecols=[0], skiprows=[0]).values
    seq_pred = [i for i in temp[0,0]]
    labels = [I for i,I in enumerate(temp[1, 0])]
    assert len(seq) == len(seq_pred) == len(labels)
    pred_pairs = get_pairs(labels)

    return pred_pairs

######## --------------------- parse CentroidFold (MFE) output ---------------------------- #########################
def CentroidFold_mfe(k, seq):
    with open('/home/jaswinder/mnt/spl/project5/predictions_single/CentroidFold/' + str(k) + '.dbn_mfe') as f:
        temp = pd.read_csv(f, comment='#', delim_whitespace=True, header=None, usecols=[0], skiprows=[0]).values
    seq_pred = [i for i in temp[0,0]]
    labels = [I for i,I in enumerate(temp[1, 0])]
    assert len(seq) == len(seq_pred) == len(labels)
    pred_pairs = get_pairs(labels)

    return pred_pairs

########--------------------- parse GREMLIN output ---------------------------- #########################
def GREMLIN(k, seq, path):

	###### read GREMLIN output #######
	with open(path + str(k) + '.dca') as f:
		temp = pd.read_csv(f, comment='#', delim_whitespace=True, header=None, skiprows=[0], usecols=[0,1,2]).values
		
	###### read down-weight dca value for |i-j| < 4 by 0.01 #######
	dca = []
	for i in temp:
		if abs(i[0]-i[1]) < 4.0:
			dca.append([i[0],i[1], 0.01*i[2]])	
		else: 
			dca.append([i[0],i[1],i[2]])	
	dca = np.array(dca)
	dca = np.flipud(dca[dca[:,2].argsort()])

	#### consider top L/4 dca values ######
	dca = dca[0:int(len(seq)/4)]

	#### get base-pair index values starting from 0 ########
	pred_pair = [[int(i[0]), int(i[1])] for i in dca]

	return pred_pair, dca

########--------------------- parse plmc output ---------------------------- #########################
def plmc(k, seq, path):
		
############ read gremlin output ##############################
	with open(path + str(k) + '.dca') as f:
		temp = pd.read_csv(f, comment='#', delim_whitespace=True, header=None, usecols=[0,2,5]).values

	###### read down-weight dca value for |i-j| < 4 by 0.01 #######
	dca = []
	for i in temp:
		if abs(i[0]-i[1]) < 4.0:
			dca.append([i[0],i[1], 0.01*i[2]])	
		else: 
			dca.append([i[0],i[1],i[2]])	

	dca = np.array(dca)
	dca = np.flipud(dca[dca[:,2].argsort()])

	#### consider top L/4 dca values ######
	dca = dca[0:int(len(seq)/4)]

	#### get base-pair index values starting from 0 ########
	pred_pair = [[int(i[0]-1), int(i[1]-1)] for i in dca]

	return pred_pair, dca

########--------------------- parse mfdca output output ---------------------------- #########################
def mfdca(k, seq, path):

	with open(path + '/MFDCA_output_' + k + '/MFDCA_apc_fn_scores_' + k + '.txt') as f:
		temp = pd.read_csv(f, comment='#', delim_whitespace=True, header=None, skiprows=[0,1,2,3,4,5,6,7,8,9,10], usecols=[0,1,2]).values

	###### read down-weight dca value for |i-j| < 4 by 0.01 #######
	dca = []
	for i in temp:
		if abs(i[0]-i[1]) < 4.0:
			dca.append([i[0],i[1], 0.01*i[2]])	
		else: 
			dca.append([i[0],i[1],i[2]])	

	dca = np.array(dca)
	dca = np.flipud(dca[dca[:,2].argsort()])

	#### consider top L/4 dca values ######
	dca = dca[0:int(len(seq)/4)]

	#### get base-pair index values starting from 0 ########
	pred_pair = [[int(i[0]-1), int(i[1]-1)] for i in dca]

	return pred_pair, dca

########--------------------- parse plmdca output output ---------------------------- #########################
def plmdca(k, seq, path):
		
############ read gremlin output ##############################

	with open(path + '/PLMDCA_output_' + k + '/PLMDCA_apc_fn_scores_' + k + '.txt') as f:
		temp = pd.read_csv(f, comment='#', delim_whitespace=True, header=None, skiprows=[0,1,2,3,4,5,6,7,8,9,10,11], usecols=[0,1,2]).values

	###### read down-weight dca value for |i-j| < 4 by 0.01 #######
	dca = []
	for i in temp:
		if abs(i[0]-i[1]) < 4.0:
			dca.append([i[0],i[1], 0.01*i[2]])	
		else: 
			dca.append([i[0],i[1],i[2]])	

	dca = np.array(dca)
	dca = np.flipud(dca[dca[:,2].argsort()])

	#### consider top L/4 dca values ######
	dca = dca[0:int(len(seq)/4)]

	#### get base-pair index values starting from 0 ########
	pred_pair = [[int(i[0]-1), int(i[1]-1)] for i in dca]

	return pred_pair, dca

######## --------------------- parse RNAalifold output ---------------------------- #########################
def rnaalifold(k, seq, path):

    with open(path + str(k) + '.dbn') as f:
        temp = pd.read_csv(f, comment='#', delim_whitespace=True, header=None).values

    #print(temp)

    seq_pred = [i for i in temp[0,0] if i!='-']
    remove_index = [i for i,I in enumerate(temp[0, 0]) if I=='-']
    labels = [I for i,I in enumerate(temp[1, 0]) if i not in remove_index]

    if k!='3ivn_A':
        assert len(seq) == len(seq_pred) == len(labels)

    pred_pair = get_pairs(labels)

    return pred_pair

######## --------------------- parse Centroid_alifold-MEA output ---------------------------- #########################
def centroidalifold(k, seq, path):

    with open(path + str(k) + '.dbn') as f:
        temp = pd.read_csv(f, comment='#', delim_whitespace=True, header=None, usecols=[0], skiprows=[0]).values

    seq_pred = [i for i in temp[0,0]]
    labels = [I for i,I in enumerate(temp[1, 0])]
    if k!='3ivn_A':
        assert len(seq) == len(seq_pred) == len(labels)
    pred_pair = get_pairs(labels)

    return pred_pair

######## --------------------- parse R-scape output ---------------------------- #########################
def rscape(k, seq, path):
	
	with open(path + str(k) + '_1.sorted.cov') as f:
		temp = pd.read_csv(f, comment='#', header=None, delim_whitespace=True, skiprows=[0,1,2,3,4,5]).values

	pred_contacts = []
	for i in temp:
		pred_contacts.append([int(i[1]-1), int(i[2]-1)])

	pred_contacts = [i for i in pred_contacts[0:int(len(seq)/4)]]

	return pred_contacts

######## --------------------- parse CaCoFold output ---------------------------- #########################
def cacofold(k, seq, path):
	
	with open(path + str(k) + '.bps') as f:
		temp = pd.read_csv(f, comment='#', header=None, delim_whitespace=True).values

	pred_pair = [[pair[0]-1,pair[1]-1] if pair[0]<pair[1] else [pair[1]-1,pair[0]-1] for pair in temp]

	return pred_pair


######## --------------------- parse RNAfold output ---------------------------- #########################
def RNAfold(k, seq, path):
    with open(path + str(k) + '.dbn') as f:
        temp = pd.read_csv(f, comment='#', delim_whitespace=True, header=None, usecols=[0], skiprows=[0,3,4,5]).values
    seq_pred = [i for i in temp[0,0]]
    labels = [I for i,I in enumerate(temp[1, 0])]
#    print(''.join(seq)); print(''.join(seq_pred)); print(''.join(labels))
    if k!='3ivn_A':
        assert len(seq) == len(seq_pred) == len(labels)
    pred_pairs = get_pairs(labels)

    return pred_pairs


########--------------------- parse GREMLIN output ---------------------------- #########################
def GREMLIN_ct(k, seq):

	with open('/home/jaswinder/mnt/spl/project5/predictions/GREMLIN/' + str(k) + '.ct') as f:
		temp = pd.read_csv(f, comment='#', header=None, skiprows=[0], nrows=len(seq))
	temp = temp[0].str.split(expand=True,).values
	seq_pred = [i for i in temp[:, 1]]
	assert len(seq_pred) == len(seq)
	label_pred = [sorted(list(i)) for i in zip(temp[:, 0], temp[:, 4])]
	label_pred = [sorted([int(i[0])-1, int(i[1])-1]) for i in label_pred if '0' not in i]
	pred_pair = []
	for i in label_pred:
		if i not in pred_pair:
		    pred_pair.append(i)

	return pred_pair

########--------------------- parse GREMLIN output ---------------------------- #########################
def template_prediction(k, seq):

	with open('/home/jaswinder/mnt/spl/project5/predictions/template_prediction/' + str(k) + '.ct') as f:
		temp = pd.read_csv(f, comment='#', header=None, skiprows=[0], nrows=len(seq))
	temp = temp[0].str.split(expand=True,).values
	seq_pred = [i for i in temp[:, 1]]
	assert len(seq_pred) == len(seq)
	label_pred = [sorted(list(i)) for i in zip(temp[:, 0], temp[:, 4])]
	label_pred = [sorted([int(i[0])-1, int(i[1])-1]) for i in label_pred if '0' not in i]
	pred_pair = []
	for i in label_pred:
		if i not in pred_pair:
		    pred_pair.append(i)

	return pred_pair

######## --------------------- parse mxfold output ---------------------------- #########################
def mxfold(k, seq):
    with open('/home/jaswinder/mnt/spl/project5/predictions_single/mxfold/' + str(k) + '.dbn') as f:
        temp = pd.read_csv(f, comment='#', delim_whitespace=True, header=None, usecols=[0], skiprows=[0,2]).values
    seq_pred = [i for i in temp[0,0]]
    labels = [I for i,I in enumerate(temp[1, 0])]
    assert len(seq) == len(seq_pred) == len(labels)
    pred_pairs = get_pairs(labels)

    return pred_pairs

######## --------------------- parse mxfold2 output ---------------------------- #########################
def mxfold2(k, seq):
    with open('/home/jaswinder/mnt/spl/project5/predictions_single/mxfold2/' + str(k) + '.dbn') as f:
        temp = pd.read_csv(f, comment='#', delim_whitespace=True, header=None, usecols=[0], skiprows=[0]).values
    seq_pred = [i for i in temp[0,0]]
    labels = [I for i,I in enumerate(temp[1, 0])]
    assert len(seq) == len(seq_pred) == len(labels)
    pred_pairs = get_pairs(labels)

    return pred_pairs

######## --------------------- parse CONTRAFold output ---------------------------- #########################
def CONTRAFold(k, seq):
    with open('/home/jaswinder/mnt/spl/project5/predictions_single/CONTRAFold/' + str(k) + '.dbn.nc') as f:
        temp = pd.read_csv(f, comment='#', delim_whitespace=True, header=None, usecols=[0], skiprows=[0,2]).values
    seq_pred = [i for i in temp[0,0]]
    labels = [I for i,I in enumerate(temp[1, 0])]
    assert len(seq) == len(seq_pred) == len(labels)
    pred_pairs = get_pairs(labels)

    return pred_pairs

######## --------------------- parse ContextFold output ---------------------------- #########################
def ContextFold(k, seq):
    #print(k)
    with open('/home/jaswinder/mnt/spl/project5/predictions_single/ContextFold/' + str(k) + '.pred') as f:
        temp = pd.read_csv(f, comment='#', delim_whitespace=True, header=None, skiprows=[0]).values
    seq_pred = [i for i in temp[0,0]]
    labels = [I for i,I in enumerate(temp[1, 0])]
    assert len(seq) == len(seq_pred) == len(labels)
    pred_pairs = get_pairs(labels)

    return pred_pairs

######## --------------------- parse Knotty output ---------------------------- #########################
def Knotty(k, seq):
    #print(k)
    with open('/home/jaswinder/mnt/spl/project5/predictions_single/Knotty/' + str(k) + '.dbn') as f:
        temp = pd.read_csv(f, comment='#', delim_whitespace=True, header=None).values
    seq_pred = [i for i in temp[0,0]]
    labels = [I for i,I in enumerate(temp[1, 0])]
    assert len(seq) == len(seq_pred) == len(labels)
    pred_pairs = get_pairs(labels)

    return pred_pairs

######## --------------------- parse IPknot output ---------------------------- #########################
def IPknot(k, seq):
    with open('/home/jaswinder/mnt/spl/project5/predictions_single/IPknot/' + str(k) + '.dbn') as f:
        temp = pd.read_csv(f, comment='#', delim_whitespace=True, header=None, usecols=[0], skiprows=[0]).values
    seq_pred = [i for i in temp[0,0]]
    labels = [I for i,I in enumerate(temp[1, 0])]
    assert len(seq) == len(seq_pred) == len(labels)
    pred_pairs = get_pairs(labels)

    return pred_pairs


######## --------------------- parse ProbKnot output ---------------------------- #########################
def ProbKnot(k, seq):
    with open('/home/jaswinder/mnt/spl/project5/predictions_single/ProbKnot/' + str(k) + '.dbn') as f:
        temp = pd.read_csv(f, comment='#', delim_whitespace=True, header=None, skiprows=[0]).values
    seq_pred = [i for i in temp[0,0]]
    labels = [I for i,I in enumerate(temp[1, 0])]
    assert len(seq) == len(seq_pred) == len(labels)
    pred_pairs = get_pairs(labels)

    return pred_pairs

######## --------------------- parse RNAstructure output ---------------------------- #########################
def RNAstructure(k, seq):
    with open('/home/jaswinder/mnt/spl/project5/predictions_single/RNAstructure/' + str(k) + '.dbn') as f:
        temp = pd.read_csv(f, comment='#', delim_whitespace=True, header=None, skiprows=[0]).values
    seq_pred = [i for i in temp[0,0]]
    labels = [I for i,I in enumerate(temp[1, 0])]
    assert len(seq) == len(seq_pred) == len(labels)
    pred_pairs = get_pairs(labels)

    return pred_pairs

######## --------------------- parse RNAshapes MEA output ---------------------------- #########################
def RNAshapes(k, seq):
    with open('/home/jaswinder/mnt/spl/project5/predictions_single/RNAshapes_MEA/' + str(k)) as f:
        temp = pd.read_csv(f, comment='#', delim_whitespace=True, header=None, skiprows=[0]).values
    seq_pred = [i for i in temp[0,0]]
    labels = [I for i,I in enumerate(temp[1, 0])]
    assert len(seq) == len(seq_pred) == len(labels)
    pred_pairs = get_pairs(labels)

    return pred_pairs

######## --------------------- parse pkiss output ---------------------------- #########################
def pkiss(k, seq):
    with open('/home/jaswinder/mnt/spl/project5/predictions_single/pkiss/' + str(k)) as f:
        temp = pd.read_csv(f, comment='#', delim_whitespace=True, header=None, skiprows=[0]).values
    seq_pred = [i for i in temp[0,0]]
    labels = [I for i,I in enumerate(temp[1, 0])]
    assert len(seq) == len(seq_pred) == len(labels)
    pred_pairs = get_pairs(labels)

    return pred_pairs

######## --------------------- parse E2Efold output ---------------------------- #########################
def E2Efold(k, seq):
    with open('/home/jaswinder/mnt/spl/project5/predictions_single/E2Efold/' + str(k) + '.dbn') as f:
        temp = pd.read_csv(f, comment='#', delim_whitespace=True, header=None, skiprows=[0]).values
    seq_pred = [i for i in temp[0,0]]
    labels = [I for i,I in enumerate(temp[1, 0])]
    assert len(seq) == len(seq_pred) == len(labels)
    pred_pairs = get_pairs(labels)

    return pred_pairs

######## --------------------- parse 2dRNA output ---------------------------- #########################
def two_dRNA(k, seq):
    with open('/home/jaswinder/mnt/spl/project5/predictions_single/2dRNA/' + str(k)) as f:
        temp = pd.read_csv(f, comment='#', delim_whitespace=True, header=None, skiprows=[0]).values
    seq_pred = [i for i in temp[0,0]]
    labels = [I for i,I in enumerate(temp[1, 0])]
    assert len(seq) == len(seq_pred) == len(labels)
    pred_pairs = get_pairs(labels)

    return pred_pairs

######## --------------------- parse CycleFold output ---------------------------- #########################
def CycleFold(k, seq):
    with open('/home/jaswinder/mnt/spl/project5/predictions_single/CycleFold/' + str(k) + '.dbn') as f:
        temp = pd.read_csv(f, comment='#', delim_whitespace=True, header=None, skiprows=[0]).values
    seq_pred = [i for i in temp[0,0]]
    labels = [I for i,I in enumerate(temp[1, 0])]
    assert len(seq) == len(seq_pred) == len(labels)
    pred_pairs = get_pairs(labels)

    return pred_pairs

######## --------------------- parse LinearFold output ---------------------------- #########################
def LinearFold(k, seq):
    with open('/home/jaswinder/mnt/spl/project5/predictions_single/LinearFold/' + str(k) + '.dbn') as f:
        temp = pd.read_csv(f, comment='#', delim_whitespace=True, header=None).values
    seq_pred = [i for i in temp[0,0]]
    labels = [I for i,I in enumerate(temp[1, 0])]
    assert len(seq) == len(seq_pred) == len(labels)
    pred_pairs = get_pairs(labels)

    return pred_pairs

######## --------------------- parse LinearPartition output ---------------------------- #########################
def LinearPartition(I, seq):
	with open('/home/jaswinder/mnt/spl/project5/predictions_single/LinearPartition/' + I + '.prob', 'r') as f:
		prob = pd.read_csv(f, delimiter=None, delim_whitespace=True, header=None).values
	y_pred =  np.zeros((len(seq), len(seq)))
	for i in prob:
		y_pred[int(i[0])-1, int(i[1])-1] = i[2]

	tri_inds = np.triu_indices(y_pred.shape[0], k=1)

	out_pred = y_pred[tri_inds]
	outputs = out_pred[:, None]
	seq_pairs = [[tri_inds[0][j], tri_inds[1][j], ''.join([seq[tri_inds[0][j]], seq[tri_inds[1][j]]])] for j in
		         range(tri_inds[0].shape[0])]

	outputs_T = np.greater_equal(outputs, 0.198)
	pred_pairs = [i for I, i in enumerate(seq_pairs) if outputs_T[I]]
	pred_pairs = [i[:2] for i in pred_pairs]

	return pred_pairs, y_pred

######## --------------------- parse LinearPartition-V output ---------------------------- #########################
def LinearPartition_V(I, seq):
	with open('/home/jaswinder/mnt/spl/project5/predictions_single/LinearPartition_V/' + I + '.prob', 'r') as f:
		prob = pd.read_csv(f, delimiter=None, delim_whitespace=True, header=None).values
	y_pred =  np.zeros((len(seq), len(seq)))
	for i in prob:
		y_pred[int(i[0])-1, int(i[1])-1] = i[2]

	tri_inds = np.triu_indices(y_pred.shape[0], k=1)

	out_pred = y_pred[tri_inds]
	outputs = out_pred[:, None]
	seq_pairs = [[tri_inds[0][j], tri_inds[1][j], ''.join([seq[tri_inds[0][j]], seq[tri_inds[1][j]]])] for j in
		         range(tri_inds[0].shape[0])]

	outputs_T = np.greater_equal(outputs, 0.198)
	pred_pairs = [i for I, i in enumerate(seq_pairs) if outputs_T[I]]
	pred_pairs = [i[:2] for i in pred_pairs]

	return pred_pairs, y_pred

############ load base-pair prob form SPOT-RNA ##############################
def parse_spot_rna(k, seq):
    y_pred = np.loadtxt('/home/jaswinder/mnt/spl/project5/predictions_single/SPOT-RNA/' + k + '.prob', delimiter='\t')
    tri_inds = np.triu_indices(y_pred.shape[0], k=1)

    out_pred = y_pred[tri_inds]
    outputs = out_pred[:, None]
    seq_pairs = [[tri_inds[0][j], tri_inds[1][j], ''.join([seq[tri_inds[0][j]], seq[tri_inds[1][j]]])] for j in
                 range(tri_inds[0].shape[0])]

    outputs_T = np.greater_equal(outputs, 0.335)
    pred_pairs = [i for I, i in enumerate(seq_pairs) if outputs_T[I]]
    pred_pairs = [i[:2] for i in pred_pairs]

    return pred_pairs, y_pred

#########------------- ensemble of spot-rna-profile output ---------------------#########
def spotrna2(k, seq, path):

    threshold=0.795

    y_pred = np.loadtxt(path + k + '_outputs/' + k + '.prob')

    tri_inds = np.triu_indices(y_pred.shape[0], k=0)

    out_pred = y_pred[tri_inds]
    seq_pairs = [[tri_inds[0][j], tri_inds[1][j], ''.join([seq[tri_inds[0][j]], seq[tri_inds[1][j]]])] for j in range(tri_inds[0].shape[0])]

    outputs = out_pred[:, None]

    outputs_T = np.greater_equal(outputs, threshold)
    pred_pair = [i[0:2] for I, i in enumerate(seq_pairs) if outputs_T[I]]

    return pred_pair, y_pred

#########------------- ensemble of spot-rna-profile-dt output ---------------------#########
def parse_spotrna_profile_dt(k, seq, threshold=0.675):

    y_pred = np.loadtxt('/mnt/ssd/Documents/rna_contact_profile/cmp_others/SPOT-RNA/output_dt/model_all/' + k + '.prob')
    tri_inds = np.triu_indices(y_pred.shape[0], k=0)

    out_pred = y_pred[tri_inds]
 #   out_true = output[tri_inds]
    seq_pairs = [[tri_inds[0][j], tri_inds[1][j], ''.join([seq[tri_inds[0][j]], seq[tri_inds[1][j]]])] for j in range(tri_inds[0].shape[0])]

#    labels = out_true[:, None]
    outputs = out_pred[:, None]

    outputs_T = np.greater_equal(outputs, threshold)
    pred_pair = [i[0:2] for I, i in enumerate(seq_pairs) if outputs_T[I]]

    return pred_pair

#########------------- ensemble of spot-rna3 output ---------------------#########
def parse_spotrna3(k, seq, threshold=0.8600):

    #y_pred = final_test_output_1[k][0]
    y_pred = final_test2_output_1[k][0]
#    print(y_pred[1])
    tri_inds = np.triu_indices(y_pred.shape[0], k=0)

    out_pred = y_pred[tri_inds]
 #   out_true = output[tri_inds]
    seq_pairs = [[tri_inds[0][j], tri_inds[1][j], ''.join([seq[tri_inds[0][j]], seq[tri_inds[1][j]]])] for j in range(tri_inds[0].shape[0])]

#    labels = out_true[:, None]
    outputs = out_pred[:, None]

    outputs_T = np.greater_equal(outputs, threshold)
    pred_pair = [i[0:2] for I, i in enumerate(seq_pairs) if outputs_T[I]]

    return pred_pair

######## --------------------- parse base-pair probability CentroidaliFold output ---------------------------- #########################
def CentroidaliFold_bp_prob(id, seq):
    with open('/home/jaswinder/mnt/spl/project5/predictions/Centroid_alifold/output_100_prob/' + str(id) + '.prob') as f:
        temp = pd.read_csv(f, comment='#', header=None)

    output_pred = np.zeros((len(seq), len(seq)))
    #print(output_pred.shape)
    for i in temp[:][0]:
        a = i.split(' ')
        for j in a[2:-1]:
            k = j.split(':')
            output_pred[int(a[0]) - 1, int(k[0]) - 1] = float(k[1])
    return output_pred

######## --------------------- parse base-pair probability RNAalifold output ---------------------------- #########################
def RNAalifold_bp_prob(id, seq):
    with open('/home/jaswinder/mnt/spl/project5/predictions/RNAalifold/output_100_prob/' + str(id) + '.prob') as f:
        temp = pd.read_csv(f, comment='#', header=None).values
    #print(temp.shape)
    output_pred = np.zeros((len(seq), len(seq)))
    for i in temp[:,0]:
        a = i.split(' ')
        #print(a)
        output_pred[int(a[0]) - 1, int(a[1]) - 1] = float(a[2])**2
    #print(output_pred)
    return output_pred

######## --------------------- parse base-pair probability LinearPartition output ---------------------------- #########################
def LinearPartition_bp_prob(I, seq):
	with open('/home/jaswinder/mnt/spl/project5/predictions_single/LinearPartition/' + I + '.prob', 'r') as f:
		prob = pd.read_csv(f, delimiter=None, delim_whitespace=True, header=None).values
	bp_prob_lp =  np.zeros((len(seq), len(seq)))
	for i in prob:
		bp_prob_lp[int(i[0])-1, int(i[1])-1] = i[2]
#	bp_prob_lp = bp_prob_lp + np.transpose(bp_prob_lp)

	return bp_prob_lp

######## --------------------- parse base-pair probability LinearPartition-V output ---------------------------- #########################
def LinearPartition_V_bp_prob(I, seq):
	with open('/home/jaswinder/mnt/spl/project5/predictions_single/LinearPartition_V/' + I + '.prob', 'r') as f:
		prob = pd.read_csv(f, delimiter=None, delim_whitespace=True, header=None).values
	bp_prob_lp =  np.zeros((len(seq), len(seq)))
	for i in prob:
		bp_prob_lp[int(i[0])-1, int(i[1])-1] = i[2]
#	bp_prob_lp = bp_prob_lp + np.transpose(bp_prob_lp)

	return bp_prob_lp

######## --------------------- parse base-pair probability CONTRAfold output ---------------------------- #########################
def CONTRAfold_bp_prob(id,seq):
    with open('/home/jaswinder/mnt/spl/project5/predictions_single/CONTRAFold/' + str(id) + '.prob') as f:
        temp = pd.read_csv(f, comment='#', header=None)

    output_pred = np.zeros((len(seq), len(seq)))
    for i in temp[:][0]:
        a = i.split(' ')
        for j in a[2:]:
            k = j.split(':')
            output_pred[int(a[0]) - 1, int(k[0]) - 1] = float(k[1])
        # print(j)
    return output_pred

######## --------------------- parse base-pair probability CONTRAfold (MFE) output ---------------------------- #########################
def CentroidFold_mfe_bp_prob(id,seq):
    with open('/home/jaswinder/mnt/spl/project5/predictions_single/CentroidFold/' + str(id) + '.prob_mfe') as f:
        temp = pd.read_csv(f, comment='#', header=None)
    #print(len(seq))
    output_pred = np.zeros((len(seq), len(seq)))
    for i in temp[:][0]:
        a = i.split(' ')
        for j in a[2:-1]:
            k = j.split(':')
            output_pred[int(a[0]) - 1, int(k[0]) - 1] = float(k[1])
    return output_pred

######## --------------------- parse base-pair probability CONTRAfold (MEA) output ---------------------------- #########################
def CentroidFold_mea_bp_prob(id,seq):
    with open('/home/jaswinder/mnt/spl/project5/predictions_single/CentroidFold/' + str(id) + '.prob_mea') as f:
        temp = pd.read_csv(f, comment='#', header=None)

    output_pred = np.zeros((len(seq), len(seq)))
    for i in temp[:][0]:
        a = i.split(' ')
        for j in a[2:-1]:
            k = j.split(':')
            output_pred[int(a[0]) - 1, int(k[0]) - 1] = float(k[1])
    return output_pred

######## --------------------- parse base-pair probability RNAfold output ---------------------------- #########################
def RNAfold_bp_prob(id, seq):
    with open('/home/jaswinder/mnt/spl/project5/predictions_single/RNAfold/' + str(id) + '.prob') as f:
        temp = pd.read_csv(f, comment='#', header=None).values
    #print(temp.shape)
    output_pred = np.zeros((len(seq), len(seq)))
    for i in temp[:,0]:
        a = i.split(' ')
        #print(a)
        output_pred[int(a[0]) - 1, int(a[1]) - 1] = float(a[2])**2
    #print(output_pred)
    return output_pred

############ load base-pair prob form SPOT-RNA ##############################
def spot_rna_bp_prob(k, seq):
    y_pred = np.loadtxt('/mnt/ssd/Documents/project5/programs/SPOT-RNA/outputs/' + k + '.prob', delimiter='\t')

    return y_pred

#########------------- ensemble of spot-rna-profile output ---------------------#########
def spotrna_profile_bp_prob_dt(k, seq):

    y_pred = np.loadtxt('/mnt/ssd/Documents/rna_contact_profile/cmp_others/SPOT-RNA/output_dt/model_all/' + k + '.prob')

    return y_pred

#########------------- ensemble of spot-rna-profile output ---------------------#########
def spotrna_profile_bp_prob(k, seq):

    y_pred = np.loadtxt('/mnt/ssd/Documents/rna_contact_profile/cmp_others/SPOT-RNA/output_tl/model_all/' + k + '.prob')

    return y_pred
