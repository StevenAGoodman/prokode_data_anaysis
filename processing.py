import sys
# inputs
feature_string = sys.argv[1] # "@@@@@@@@100001000000" #
data_organism = sys.argv[2] # "Escherichia coli" # 
genome_length = int(sys.argv[3]) #   # 
prokode_dir = "/workspaces/prokode" # sys.argv[4] # 

import pandas as pd
import numpy as np
import json
import statistics
import cvxpy
import re
import os
sys.path.append(prokode_dir)
from scipy.optimize import curve_fit
from src.maths import *

# be sure to run run.py before running this file
network_loc = prokode_dir + "/src/network.json"
data_dir = prokode_dir + "/beta_training/GEO_expression_data/" + data_organism

cell_volume = 1e-15 # L

# global jazz
num_data = 10
sensor_normal_dist = 10
basal_rate = 3
decay_rate = np.log(2)/300
Np = 6000
Kd_p = 0.1
Nns = 4600000

def protein_amnts_from_mRNA_amnts(mRNA_amnts):
     protein_amnts = mRNA_amnts
     return protein_amnts

def function_for_timepoint(data_t0, data_t1, protein_data_t0, N_rnap, N_ribo, dt, gene_key, tf_key):
     # MOST CRITICAL PART TO DETERMINE!!!
     protein_data_t1 = [None] * len(gene_key)

     for i in range(len(gene_key)):
          gene = gene_key[i]
          gene_mRNA_t0 = data_t0[i]

          gene_info_dict = search_network_json(network_loc, network_key, gene)

          if gene_info_dict == None:
               print(f"\t\t\t\tFailed: {gene}\tgene from data not found in genome")
               continue
          
          regulators_dict = gene_info_dict["regulators"]
          Kd_rnap_gene = regulators_dict["polymerase"]

          overall_mRNA_change_rate = (data_t1[i] - data_t0[i]) / dt
          coefficient_arr, beta_all, N_rnap, N_ribo, R_max_trans, Kd_ribo, protein_decay_rate_farr, growth_rate_fid, beta_function_fid, translation_rate_fid = beta_from_overall_mRNA(gene, gene_mRNA_t0, overall_mRNA_change_rate, protein_data_t0, gene_info_dict, regulators_dict, gene_key, tf_key, genome_length, Kd_rnap_gene, N_rnap, N_ribo, feature_string)

          # update protein amnts
          overall_protein_change_rate = translation_rate(gene, gene_mRNA_t0, protein_data_t0, N_ribo, gene_info_dict, R_max_trans, Kd_ribo, translation_rate_fid) - protein_decay_rate(gene, protein_data_t0, gene_info_dict, gene_key, growth_rate_fid, protein_decay_rate_farr) * protein_data_t0[i]
          protein_data_t1[i] = protein_data_t0[i] + overall_protein_change_rate * dt

     return coefficient_arr, beta_all, protein_data_t1, N_rnap, N_ribo, beta_function_fid

def fit_function(coefficient_matrix, beta_all_arr, feature_id):
     # clean matrices
     empty_columns = []
     for i in range(len(list(zip(*coefficient_matrix)))):
          column = list(zip(*coefficient_matrix))[i]
          if column == [0] * len(coefficient_matrix):
               empty_columns.append(i)
          else:
               continue

     # curve fitting parameters
     if feature_id == "@":
          def beta_function(P_arr, *beta_arr):
               assert len(P_arr) == len(beta_arr)
               return np.array(P_arr) @ np.array(beta_arr).T
          func_to_fit = beta_function
     elif feature_id == "A":
          def beta_function(P_arr, *beta_arr):
               assert len(P_arr) == len(beta_arr)
               return np.log10(np.array(P_arr)) @ np.log10(np.array(beta_arr)).T
          func_to_fit = beta_function
     elif feature_id == "B":
          def beta_function(P_arr, *beta_arr):
               assert len(P_arr) == len(beta_arr)
               return np.array(P_arr) @ np.log10(np.array(beta_arr)).T
          func_to_fit = beta_function
     elif feature_id == "C":
          def beta_function(P_arr, *beta_arr):
               assert len(P_arr) == len(beta_arr)
               return np.array(P_arr) @ np.log2(np.array(beta_arr)).T
          func_to_fit = beta_function

     n_features = len(coefficient_matrix[0])
     p = [1] * n_features
     beta_arr, covar_uncertainty = curve_fit(func_to_fit, coefficient_matrix, beta_all_arr, p0 = p)

     for index in empty_columns:
          beta_arr[index] = None

     return beta_arr, covar_uncertainty

network_key = create_network_key(network_loc)
data_file_list = os.listdir(data_dir)

for data_file in data_file_list:
     # get data
     data_df = pd.read_csv(f"{data_dir}/{data_file}", index_col=0)
     data_df = data_df.multiply(1 / 6.02e23).multiply(1 / cell_volume)
     gene_key = list(data_df.index) # may have to subtracted header
     tf_key = [ tf for val in json.load(open(network_loc, 'r')).values() for tf in val["regulators"].keys() ]
     tf_key = list(set(tf_key)) # filter repeats

     # prep stuffs
     # beta_all (1, n_samples) coefficient_matrix (n_samples, n_coefficients)
     beta_all_arr, coefficient_matrix = ([], [])
     results_matrix = []

     # lil preprocessing
     groups = [ [] for i in range(100) ]
     for i in range(len(data_df.axes[1])):
          colname = data_df.columns[i]
          group_n = re.search("\\| (\\d+)", colname)
          if group_n == None:
               group_n = "0"
          else:
               data_df = data_df.rename(columns = {colname: colname[:colname.find("|")]})
               group_n = group_n.group(1)
          groups[int(group_n)].append(i)
     groups = [x for x in groups if x != []]

     # do the stuff
     for i in range(len(groups)):
          print(f"\tgroup {i}")
          group_data_indecies = sorted(groups[i])

          protein_data_matrix = [list(data_df.iloc[:,0])]
          N_rnap = 2200
          N_ribo = 3000

          col_names = []

          for l in range(len(group_data_indecies[:-1])):
               print(f"\t\tdata pair {l}")
               n = group_data_indecies[l]
               m = group_data_indecies[l+1]
               print(n)
               print(m)
               # input: pair of time points
               data_t0 = list(data_df.iloc[:,n])
               data_t1 =  list(data_df.iloc[:,m])
               print(data_df.columns[m])
               print(data_df.columns[n])
               dt = float(data_df.columns[m]) - float(data_df.columns[n])

               protein_data_t0 = protein_data_matrix[group_data_indecies.index(n)]

               coefficient_arr, beta_all, protein_data_t1, N_rnap, N_ribo, beta_function_fid = function_for_timepoint(data_t0, data_t1, protein_data_t0, N_rnap, N_ribo, dt, gene_key, tf_key)
               coefficient_matrix.append(coefficient_arr)
               beta_all_arr.append(beta_all)

               col_names.append(f"{data_df.columns[n]} to {data_df.columns[m]}")
               protein_data_matrix.append(protein_data_t1)
               print(f"\t\t\tbeta all: {beta_all}")

# overall function fitting for all data points in an organisms
beta_arr, covar = fit_function(coefficient_matrix, beta_all_arr, beta_function_fid)

# export to csv
results_matrix = np.array(results_matrix)
beta_names = tf_key
results_df = pd.DataFrame(results_matrix.T, index = beta_names, columns = col_names)
results_df.to_csv(f"{prokode_dir}/beta_training/results/{data_file[data_file.find('beta_training')+13:data_file.find('.py')]}_group{i}.csv")
