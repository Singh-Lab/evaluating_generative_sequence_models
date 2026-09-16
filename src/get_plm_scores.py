from transformers import (
	EsmForMaskedLM,
	BertTokenizer,
	EsmTokenizer,
	pipeline,
	BertForMaskedLM
)
import numpy as np
import pandas as pd 
import copy
import torch
import sys 

from config import HUGGINGFACE_CACHE_PATH

def softmax_1d(x):
	# Subtract max for numerical stability
	e_x = np.exp(x - np.max(x))
	return e_x / e_x.sum()

if __name__ == "__main__":

	input_name = sys.argv[1]

	with open(f"./fastas/{input_name}.fasta", "r") as f1:
		lines = f1.readlines()
	f1.close()
	lines = list(map(lambda x : x.replace("\n", "").replace(">", ""), lines))
	headers = lines[::2]
	seqs = lines[1::2]
	header_to_seq = dict(zip(headers, seqs))


	########################################################

	key_to_model_string = {
		"esm_1b" : "facebook/esm1b_t33_650M_UR50S",
		"esm_1v" : "facebook/esm1v_t33_650M_UR90S_1",
		"esm_2_15B" : "facebook/esm2_t48_15B_UR50D",
		"esm_2_3B" :  "facebook/esm2_t36_3B_UR50D",
		"esm_2_650M" : "facebook/esm2_t33_650M_UR50D",
		"esm_2_150M" : "facebook/esm2_t30_150M_UR50D",
		"esm_2_35M" :  "facebook/esm2_t12_35M_UR50D",
		"esm_2_8M" :   "facebook/esm2_t6_8M_UR50D",
	}

	model_key = sys.argv[2]
	
	device = "cuda" if torch.cuda.is_available() else "cpu"
	model_name = key_to_model_string[model_key]
	local_files_only = True

	tokenizer = EsmTokenizer.from_pretrained(
		model_name,
		do_lower_case=False,
		cache_dir=HUGGINGFACE_CACHE_PATH,
		local_files_only=local_files_only
	)

	vocab = tokenizer.get_vocab()

	model = EsmForMaskedLM.from_pretrained(
		model_name,
		cache_dir=HUGGINGFACE_CACHE_PATH,
		local_files_only=local_files_only
	)

	model.eval()
	model.to(device)

	results = []


	for counter, (header, sequence) in enumerate(header_to_seq.items()):

		if counter % 100 == 0 : 
			print("counter ", counter)

		seq_len = len(sequence)

		wt_seq = copy.deepcopy([*sequence])
		
		tokenized_input = tokenizer(
			wt_seq,
			return_tensors="pt",
			is_split_into_words=True
		).to(device)

		# will be used to map back to logits
		


		output = model(**tokenized_input, output_attentions=False, output_hidden_states = True)

		# this is the embedding of the seqeunce : Size  = (seq_len + 2, plm_dim )
		final_hidden = output.hidden_states[-1].cpu().detach().numpy()[0]
		# CLS-pooled embedding is the special token at the beginning of seq
		cls_pooled_embedding = final_hidden[0, :] # Size = (plm_dim)
		# EOS-pooled embedding is the special token at the end of seq
		eos_pooled_embedding = final_hidden[-1, :] # Size = (plm_dim)

		mean_pooled_embedding = np.mean(final_hidden[1:-1, ], axis = 0) # Size = (plm_dim) - removing the EOS and CLS tokens in the mean pooling

		logits = output.logits.cpu().detach().numpy() # will be shape (1, seq_len + 2, vocab_size)

		wt_probs = []


		for idx, aa in enumerate(wt_seq):
			# using idx+1 because of the CLS token is at position 0.
			# the first M token will be at position 1
			# only computing the scores over the AAs, not special tokens
			aa_distribution = logits[0][idx + 1] # will be size (vocab_size)
			probs = softmax_1d(aa_distribution)
			aa_token_idx = vocab[aa]
			aa_prob = probs[aa_token_idx]
			wt_probs.append(aa_prob) # probability of the amoino acid at the position


		wt_probs = np.asarray(wt_probs) # size (seq_len)
		average_prob = np.sum(wt_probs) / seq_len
		likelihood = np.prod(wt_probs)
		log_likelihood = np.sum(np.log(wt_probs))


		r1 = {
			"name" : header,
			"likelihood" : likelihood,
			"log_likelihood" : log_likelihood,
			"avg_log_likelihood" : log_likelihood / seq_len,
			"average_prob" : average_prob,
			"length" : seq_len
		}

		results.append(r1)
		
		# print(counter, header, average_prob, log_likelihood, log_likelihood / seq_len)

	out_fn = f"./plm_outputs/{input_name}_{model_key}.csv"
	dfo = pd.DataFrame(results)
	dfo.to_csv(out_fn)
	print("done")



