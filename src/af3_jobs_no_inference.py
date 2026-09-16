
import sys 
import json
import os
import glob 
from config import DATA_PATH, AF3_SHARED_PATH
from a_get_esmfold_structures import read_fasta

def chunk_list(lst, chunk_size=20):
	return [
		lst[i:i + chunk_size]
		for i in range(0, len(lst), chunk_size)
	]





if __name__ == "__main__":
	input_name = sys.argv[1]

	if not os.path.exists(DATA_PATH + f"/af3/input_json/{input_name}"):
		os.mkdir(DATA_PATH + f"/af3/input_json/{input_name}")

	if not os.path.exists(DATA_PATH + f"/af3/aln_data_outputs/{input_name}"):
		os.mkdir(DATA_PATH + f"/af3/aln_data_outputs/{input_name}")



	with open(f"./fastas/{input_name}.fasta", "r") as f1:
		lines = f1.readlines()
	f1.close()


	lines = list(map(lambda x : x.replace("\n", "").replace(">", ""), lines))

	headers = lines[::2]
	seqs = lines[1::2]


	for seq_name, seq in zip(headers, seqs):
		af3_json = {
			"name": seq_name,
			"modelSeeds": [0],
			"sequences": [
				{
					"protein": {
						"id": "A",
						"sequence" : seq
					}
				}
			],
			"dialect": "alphafold3",
			"version": 1
		}
		
		data_fn = DATA_PATH + f"/af3/input_json/{input_name}/{seq_name}.json"

		with open(data_fn, "w") as f1:
			json.dump(af3_json, f1)
		f1.close()


	header_to_seq = dict(zip(headers, seqs))

	num_seqs_per_job = 3

	
	cwd = os.getcwd()


	needed_idents = set()
	seen_idents = set()

	completed_files = glob.glob(DATA_PATH + f"/af3/aln_data_outputs/{input_name}/*")

	for fn in completed_files:
		seen_idents.add(fn.split("/")[-1])

	for h in headers:
		if h.lower() not in seen_idents:
			needed_idents.add(h)


	# for s in seen_idents:
	# 	if "91535b_gan_selected_mdh_round3" in s.lower():
	# 		print(s)
	# 		print("yoooooo")

		

	# print(len(seen_idents), len(needed_idents), len(headers))

	# print(needed_idents)

	

	constant_headers = [
		"#SBATCH --cpus-per-task=8",
		"#SBATCH --mem=256G",
		"#SBATCH --nodes=1" ,
		"#SBATCH --ntasks=1" ,
		"#SBATCH --time=04:30:00",
	]

	print(cwd)
	
	MAX_JOBS = 50

	list_chunks = chunk_list(list(needed_idents), chunk_size = num_seqs_per_job)
	print(len(list_chunks), len(needed_idents)) 
	quit()
	
	for i in range(len(list_chunks)):

		job_fn = cwd + f"/jobs/{input_name}_job_{i}.sh"
		
		if i > MAX_JOBS:
			quit()

		with open(job_fn, "w") as f1:
			f1.write("#!/bin/sh\n")
			f1.write(f"#SBATCH --output={cwd}/outputs/slurm-%j.out\n")

			for header in constant_headers:
				f1.write(header + "\n")

			f1.write("\nmodule load alphafold/3.0.1\n\n")

			seqs_needed = list_chunks[i]

			for seq_name in seqs_needed:
				data_fn = DATA_PATH + f"/af3/input_json/{input_name}/{seq_name}.json"

				af_args = [
					"run_alphafold.py", 
					"--db_dir", 
					AF3_SHARED_PATH,
					"--model_dir",
					AF3_SHARED_PATH,
					"--json_path",
					data_fn,
					"--output_dir",
					f"{DATA_PATH}/af3/aln_data_outputs/{input_name}",
					"--norun_inference"
				]

				af_args_str = " ".join(af_args)

				f1.write(af_args_str + "\n\n\n")

		f1.close()

		os.system("sbatch " + job_fn)






