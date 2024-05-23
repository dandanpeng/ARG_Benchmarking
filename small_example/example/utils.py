import tsinfer
import tsdate
import tskit
import numpy as np

import pandas as pd
import os

import time
import datetime
import arg_needle_lib as anl
from asmc.preparedecoding import *

os.environ['OPENBLAS_NUM_THREADS'] = '1'

def tsfunc(data, sel_pos, N, u):
  data = data.astype(int)

  #data[:, 0] = data[:, 0] - 1
  with open(r.new_temp + '/data.npy', 'wb') as f:
      np.save(f, data)

  with tsinfer.SampleData(sequence_length = r.len_hap + 10) as sample_data:
    for i in range(0,len(data[:,0])):
        sample_data.add_site(data[i,0], data[i, 1:], ["A", "T"], ancestral_allele = 0)
    sample_data.finalise()
  
  start_time = time.time()
  inferred_ts = tsinfer.infer(sample_data)
  inferred_ts = inferred_ts.simplify(keep_unary = False)
  
  dated_ts = tsdate.date(inferred_ts, Ne=N, mutation_rate = u)
  print("-----%s seconds -------" % (time.time() - start_time))

  #dated_ts.dump(r.new_temp + "/dated_tree")
  
  #dated_tree = tskit.load(r.new_temp + "/dated_tree")
  middle_tree = dated_ts.at(sel_pos)
  middle_tree = middle_tree.split_polytomies() 
  
  return(middle_tree.as_newick())
  
  # if (n+1) % 2 == 0: # if even, return trees on either side
  #   for i in range((n+1)/2,(n+1)/2+1):
  #     tr = inferred_ts.at(i)
  #     middle_tree.append(tr.newick())
  # else: # return tree at 
  #   tr = inferred_ts.at((n+2)/2)
  #   middle_tree.append(tr.newick())
  # 
  # return middle_tree
  
  #trees = []
  #for i in inferred_ts.trees():
  #  trees.append(i.newick())
  #return trees

def aw_input(file, pos, len_hap, temp_dir):
    
    #os.chdir(temp_dir)
    df = pd.read_csv(file, skiprows = 1, header = None)
    
    df_expand = df.iloc[:,0].str.split('', expand = True)
    df_expand = df_expand.drop(columns = [0, df_expand.shape[1] - 1])
    
    df_t = df_expand.T
    
    df_t = df_t.replace('0', 'A')
    df_t = df_t.replace('1', 'T')
    
    hap = df_t[df_t.columns].apply(lambda x: ''.join(x), axis = 1)
    
    #pos = pd.read_csv(file, nrows = 1, header = None, sep = ' ')
    #pos = pos.dropna(axis = 'columns')
    
    #pos = list((pos.iloc[0, :] * len_hap).astype(int))

    #print([pos.index(item) for item, count in collections.Counter(pos).items() if count > 1])
    '''
    for item, count in collections.Counter(pos).items():
        if count > 1:
            idx = pos.index(item) 
            for j in range(1, count):           
                pos[idx + j] = item + j
    '''            
    names = []
    for i in range(df_t.shape[1]):
        names.append('n' + str(i))
    os.chdir(temp_dir)
    with open("aw.sites", "w") as f:
        f.write("NAMES" + "\t" + "\t".join(names) + "\n")
        f.write("\t".join(["REGION", "chr", "1", str(len_hap)]) + "\n")
        for i in range(len(pos)):
            f.write(str(pos[i]) + "\t" + hap.iloc[i] + "\n")
    os.chdir("../")
    

def an_newick(argn_file, newick_file):
    arg = anl.deserialize_arg(argn_file)
    arg.populate_children_and_roots()
    an_newick = anl.arg_to_newick(arg, verbose = True)
    sep_newick = an_newick.split("\n")
    pos_array = np.zeros((len(sep_newick),2))
    for i in range(len(sep_newick)):
        if(sep_newick[i] != ''):
            pos_array[i, 0] = float(sep_newick[i].split("):")[0].split('[')[1].split(',')[0])
            pos_array[i, 1] = float(sep_newick[i].split("):")[0].split('[')[1].split(',')[1])
    tree_index =  np.where( (pos_array[:,0] <= 10e4) & (pos_array[:,1] >= 10e4))[0][0]
    tree_newick = sep_newick[tree_index].split("):")[1]

    newick_file = open(newick_file, "w")
    newick_file.write(tree_newick+';')
    newick_file.close()
    
def create_decoding_quantities(demo_file, disc_file, freq_file):
    dq = prepare_decoding(
            demography = demo_file,
            discretization = disc_file,
            frequencies = freq_file,
            samples = 300,
            mutation_rate = 2e-8,
            #recombination_rate = 2.5e-8,
    )
    dq.save_decoding_quantities(r.new_temp + "/custom")

def get_singer_tree(singer_tree_seq, sel_pos):
    tree_seq = tskit.load(singer_tree_seq)
    sel_tree = tree_seq.at(sel_pos)
    sel_tree = sel_tree.split_polytomies()
    return(sel_tree.as_newick())

