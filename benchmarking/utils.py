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
    np.savetxt(os.path.join(save_result_path, str(id))+'.ct', (temp), delimiter='\t\t', fmt="%s", header=str(len(seq)) + '\t\t' + str(id) + '\n', comments='')

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

	#### consider top L/3 dca values ######
	dca = dca[0:int(len(seq)/3)]

	#### get base-pair index values starting from 0 ########
	pred_pair = [[int(i[0]), int(i[1])] for i in dca]

	return pred_pair, dca

########--------------------- parse plmc output ---------------------------- #########################
def plmc(k, seq, path):
		
    ######### read PLMC output #############
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

	#### consider top L/3 dca values ######
	dca = dca[0:int(len(seq)/3)]

	#### get base-pair index values starting from 0 ########
	pred_pair = [[int(i[0]-1), int(i[1]-1)] for i in dca]

	return pred_pair, dca

########--------------------- parse mfdca output output ---------------------------- #########################
def mfdca(k, seq, path):

    ######### read mfDCA output #############
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

	#### consider top L/3 dca values ######
	dca = dca[0:int(len(seq)/3)]

	#### get base-pair index values starting from 0 ########
	pred_pair = [[int(i[0]-1), int(i[1]-1)] for i in dca]

	return pred_pair, dca

########--------------------- parse mfdca output output ---------------------------- #########################
def mfdca_top_L_by_n(k, seq, n, path):

    ######### read mfDCA output #############
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

	#### consider top L/n dca values ######
	dca = dca[0:int(len(seq)/n)]

	#### get base-pair index values starting from 0 ########
	pred_pair = [[int(i[0]-1), int(i[1]-1)] for i in dca]

	return pred_pair, dca


########--------------------- parse plmdca output output ---------------------------- #########################
def plmdca(k, seq, path):
		
    ######### read plmDCA output #############
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

	#### consider top L/3 dca values ######
	dca = dca[0:int(len(seq)/3)]

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

    return pred_pair, None


######## --------------------- parse CaCoFold output ---------------------------- #########################
def cacofold(k, seq, path):
	
	with open(path + str(k) + '.bps') as f:
		temp = pd.read_csv(f, comment='#', header=None, delim_whitespace=True).values

	pred_pair = [[pair[0]-1,pair[1]-1] if pair[0]<pair[1] else [pair[1]-1,pair[0]-1] for pair in temp]

	return pred_pair, None


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

######## --------------------- parse SPOT-RNA output ---------------------------- #########################
def spotrna(k, seq, path):
    with open(path + str(k) + '.dbn') as f:
        temp = pd.read_csv(f, comment='#', delim_whitespace=True, header=None, skiprows=[0]).values
    seq_pred = [i for i in temp[0,0]]
    labels = [I for i,I in enumerate(temp[1, 0])]
    assert len(seq) == len(seq_pred) == len(labels)
    pred_pairs = get_pairs(labels)

    return pred_pairs


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

#########------------- get base-pairs from SPOT-RNA2 output ---------------------#########
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



