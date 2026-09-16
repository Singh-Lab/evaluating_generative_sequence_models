import glob 
import os
import sys 
import numpy as np 
import pandas as pd 
import subprocess
from config import DATA_PATH, HOME_PATH 


if __name__ == "__main__":
	input_name = sys.argv[1]

	with open(f"./fastas/{input_name}.fasta", "r") as f1:
		lines = f1.readlines()
	f1.close()


	lines = list(map(lambda x : x.replace("\n", "").replace(">", ""), lines))

	headers = lines[::2]
	seqs = lines[1::2]

	header_to_seq = dict(zip(headers, seqs))

	protein_mpnn_res_path = DATA_PATH + f"/protein_MPNN/model_outputs/{input_name}/"


	results = []

	alphabet = "ACDEFGHIKLMNPQRSTVWYX"

	alphabet_index = {aa : i for (i, aa) in enumerate([*alphabet])}


	for header, seq in zip(headers, seqs):
		identifier = header.lower()

		npz_path = protein_mpnn_res_path + identifier + ".npz"
		if not os.path.exists(npz_path):
			print("missing", identifier)
			continue 

		data = np.load(npz_path)

		data_len = np.shape(data)[0]

		if len(seq) != np.shape(data)[0]:

			d1 = {
				"identifier" : identifier,
				"header" : "output_data_unequal_length",
				"seq" : seq,
				"log_likelihood": "output_data_unequal_length",
				"average_log_likelihood" : "output_data_unequal_length",
				"average_aa_probability" : "output_data_unequal_length",
				"seq_len" : seq_len,
				"data_len" : data_len
			}
			results.append(d1)
			continue
		
		seq_probs = []
		for pos, wt_aa in enumerate([*seq]):
			aa_index = alphabet_index[wt_aa]
			aa_prob = data[pos][aa_index]
			seq_probs.append(aa_prob)

		logP = np.sum(np.log(seq_probs))
		mean_prob = np.mean(seq_probs)

		seq_len = len(seq) 
		d1 = {
			"identifier" : identifier,
			"header" : header,
			"seq" : seq,
			"log_likelihood": logP,
			"average_log_likelihood" : logP / seq_len,
			"average_aa_probability" : mean_prob,
			"seq_len" : seq_len,
			"data_len" : data_len
		}

		results.append(d1)

	dfo = pd.DataFrame(results)
	dfo.to_csv(f"./proteinmpnn_likelihoods/{input_name}_proteinMPNN.csv")






		# pass

