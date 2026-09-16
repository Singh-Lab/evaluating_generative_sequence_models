import glob 
import os
import sys 
import subprocess
import json 
from config import DATA_PATH, HOME_PATH 
import pandas as pd 

if __name__ == "__main__":
	input_name = sys.argv[1]
	af3_dir = DATA_PATH + f"/af3/structure_outputs/{input_name}/"
	print(af3_dir)
	with open(f"./fastas/{input_name}.fasta", "r") as f1:
		lines = f1.readlines()
	f1.close()


	lines = list(map(lambda x : x.replace("\n", "").replace(">", ""), lines))

	headers = lines[::2]
	seqs = lines[1::2]

	counter = 0

	results = []
	for header in headers:
		header = header.lower()

		counter += 1

		seq_struct_path = DATA_PATH + f"/af3/structure_outputs/{input_name}/{header}/"
		seq_struct_results = glob.glob(seq_struct_path + "*")



		if not os.path.exists(seq_struct_path + f"{header}_summary_confidences.json"):
			print("header doesn't exist", header)
			esults.append([header, "error"])
			continue 

		try:

			with open(seq_struct_path + f"{header}_summary_confidences.json", "r") as f1:
				data = json.load(f1)
			f1.close()

			results.append([header, data["ptm"]])

		except Exception as e:
			print("error with josn reading", header)
			print(e)

	dfo = pd.DataFrame(results, columns = ["header", "pTM"])
	dfo.to_csv(f"./ptm_results/{input_name}.csv")


