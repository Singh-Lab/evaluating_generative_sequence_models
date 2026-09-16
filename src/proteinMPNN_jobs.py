import glob 
import os
import sys 
import subprocess
from config import DATA_PATH, HOME_PATH 

if __name__ == "__main__":
	input_name = sys.argv[1]

	if not os.path.exists(DATA_PATH + f"/protein_MPNN/json_inputs/{input_name}"):
		os.mkdir(DATA_PATH + f"/protein_MPNN/json_inputs/{input_name}")

	if not os.path.exists(DATA_PATH + f"/protein_MPNN/model_outputs/{input_name}"):
		os.mkdir(DATA_PATH + f"/protein_MPNN/model_outputs/{input_name}")


	with open(f"./fastas/{input_name}.fasta", "r") as f1:
		lines = f1.readlines()
	f1.close()


	lines = list(map(lambda x : x.replace("\n", "").replace(">", ""), lines))

	headers = lines[::2]
	seqs = lines[1::2]

	header_to_seq = dict(zip(headers, seqs))

	print(len(headers), len(seqs))

	pdb_path = DATA_PATH + f"/af3/consolidated_structures/{input_name}/"
	
	all_pdb_files = glob.glob(pdb_path + "*.pdb")
	all_pdb_files = set(map(lambda x : x.split("/")[-1].replace(".pdb", ""), all_pdb_files))

	print(len(all_pdb_files))

	header_to_pdb_fn = {}

	for header, seq in zip(headers, seqs):
		ident = header.lower()
		pdb_fn = pdb_path + header.lower() + ".pdb"

		if ident in all_pdb_files:
			all_pdb_files.remove(ident)

		if not os.path.exists(pdb_fn):
			print("missing", header, pdb_fn)

			continue 

		header_to_pdb_fn[header] = pdb_fn

	if len(header_to_pdb_fn) == len(headers):
		print("found all pdb files, continuing")

	else:
		print("missing certain pdb files, fix before continuing ")
		quit()
	
	for header, pdb_fn in header_to_pdb_fn.items():
		sequence = header_to_seq[header]
		seq_len = len(sequence)
		designable_res = "A1-A" + str(seq_len) 

		json_flags = [
			"--pdb_dir",
			"./tmp_pdb_dir/",
			"--designable_res",
			designable_res
		]

		os.system("cp " + pdb_fn + " ./tmp_pdb_dir/")

		# creates proteinmpnn_res_specs.json
		command_1 = ["python", HOME_PATH + "/proteinMPNN/proteinmpnn/run/generate_json.py"] + json_flags

		subprocess.run(
			command_1,
			check=True
		)

		print("FINISHED COMMAND 1 ", header)

		out_folder = DATA_PATH + f"/protein_MPNN/model_outputs/{input_name}"
		proteinmpnn_flags = [
			"--dump_probs",
			"--pdb_dir", "./tmp_pdb_dir/", 
			"--backbone_noise", "0.0",
			"--num_seq_per_target","5",
			"--batch_size", "1", # likely doesn't need changing
			"--sampling_temp", "0.1", # specify what sampling temperatures to use
			"--out_folder", out_folder,
			"--design_specs_json", "proteinmpnn_res_specs.json",
			"--model_name", "v_48_020" # specify which model to use

		]

		command_2 = ["python", HOME_PATH + "/proteinMPNN/proteinmpnn/run/run_protein_mpnn.py"] + proteinmpnn_flags
	
		subprocess.run(
			command_2,
			check=True
		)

		os.system("rm ./tmp_pdb_dir/*.pdb")

		print("FINISHED COMMAND 2 ", header)

	# module load anaconda3/2025.6
	# conda activate /scratch/gpfs/MONA/jf9645/conda_envs/mpnn_cu12.4
	# python /scratch/gpfs/MONA/jf9645/proteinMPNN/proteinmpnn/run/generate_json.py @json.flags
	# python /scratch/gpfs/MONA/jf9645/proteinMPNN/proteinmpnn/run/run_protein_mpnn.py @proteinmpnn.flags
	# python /scratch/gpfs/MONA/jf9645/proteinMPNN/proteinmpnn/run/helper_scripts/other_tools/view_probs.py -input mpnn_out/p53_human.npz


