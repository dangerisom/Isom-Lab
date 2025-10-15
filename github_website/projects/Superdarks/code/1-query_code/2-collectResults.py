# Copyright (c) 2025 Isom Lab
# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.
# If a copy of the MPL was not distributed with this file, You can obtain one at
# https://mozilla.org/MPL/2.0/.

result_path = "/projectnb/isomlab/dan/in-and-out/output/alphaFold_uniprot_v4/"
result_directory = "2025.09.27.7tmps-foldseek_e0p01/"
result_analysis_tag = "-analysis_foldseek_e0p01"
save_directory = result_path + result_directory
min_tm_score = 0.5 # 0.45
min_residues = 200 # 200

# result_path = "/projectnb/isomlab/dan/in-and-out/output/alphaFold_uniprot_v4/2023.11.12.7tmps/"
# result_directory = "2025.08.20.7tmps-rhod_glp1/"
# save_directory = result_path + result_directory
# min_tm_score = 0.5 # 0.45
# min_residues = 200 # 200

# Dependencies...
#################################################################################################
from os import listdir, mkdir, sep
from os.path import exists
from decimal import Decimal, ROUND_HALF_UP

# Classes...
#################################################################################################
class Result:

	"""
	Class:
	------
	Result

	Description:
	------------
	This class **parses and stores structural alignment results**  
	from a space-separated result file line.

	================================================================================
	📌 Attributes
	================================================================================

	result_file_line : str
		The **original file line** from the results file.

	data : list
		A list of **parsed values** extracted from `result_file_line`.

	query_pdb_file_path : str
		Path to the **query PDB file** used in structural comparison.

	target_pdb_file_path : str
		Path to the **target (mobile) PDB file** that was aligned to the query.

	tm_score : float
		The **Template Modeling (TM) score**, **normalized** by the target structure.  
			📌 **Higher TM-score (~1.0) → Better structural match.**  
			📌 **Lower TM-score (~0.0) → Poor structural match.**  

	rmsd : float
		The **Root Mean Square Deviation (RMSD)** between aligned C-alpha atoms.  
			📌 **Lower RMSD → More similar structures.**  
			📌 **Higher RMSD → More deviation in atomic positions.**  

	query_length : int
		The **number of residues** in the **query structure**.

	n_align : int
		The **number of aligned residues** between the **query and target structures**.

	================================================================================
	📌 Methods
	================================================================================

	__repr__():
		Returns a **CSV-formatted string representation** of the alignment result.  
		The string includes:
			✅ **TM-score**  
			✅ **RMSD**  
			✅ **Query length**  
			✅ **Aligned residues count**  
			✅ **Query PDB file path**  
			✅ **Target PDB file path**  

	"""

	def __init__(self, result_file_line):

		self.result_file_line = result_file_line
		self.data = result_file_line.strip().split()
		self.query_pdb_file_path = self.data[0]
		self.target_pdb_file_path = self.data[1]
		self.tm_score_target = Decimal(str(float(self.data[2]))).quantize(Decimal('0.01'), rounding=ROUND_HALF_UP) 	# TM score normalized by the target structure (i.e. TM1)
		self.tm_score_query = Decimal(str(float(self.data[3]))).quantize(Decimal('0.01'), rounding=ROUND_HALF_UP) 	# TM score normalized by the query structure (i.e., TM2)
		self.tm_score = 0.0 																						# highest of the two scores TM1 and TM2 
		self.rmsd = float(self.data[4])
		self.query_length = int(self.data[8])
		self.n_align = int(self.data[10])

	def __repr__(self):

		return_string = "%f,%f,%i,%i,%s,%s\n" % (
													self.tm_score,
													self.rmsd,
													self.query_length,
													self.n_align,
													self.query_pdb_file_path,
													self.target_pdb_file_path

												)
		return return_string

#################################################################################################
#
# Main code...
#
#################################################################################################

"""
Description:
------------
This script **processes structural alignment results**, filters hits based on  
**TM-score and residue count**, and **saves the filtered results** for further analysis.

================================================================================
📌 Variables & Parameters
================================================================================

result_path : str
	Base directory where **alignment results** are stored.

result_directory : str
	Subdirectory containing **specific results** for the analysis.

save_directory : str
	Output directory where **filtered results, statistics, and histograms**  
	are saved.

min_tm_score : float (default=0.45)
	Minimum **TM-score** threshold for filtering alignment hits.

min_residues : int (default=150)
	Minimum **residue count** for filtering alignment hits.

================================================================================
📌 Class: Result
================================================================================
This class **parses and stores structural alignment results**  
from a **space-separated results file line**.

### Attributes:
	✅ `result_file_line` → Original result file line.  
	✅ `query_pdb_file_path` → Path to **query PDB file**.  
	✅ `target_pdb_file_path` → Path to **target (mobile) PDB file**.  
	✅ `tm_score` → **Template Modeling Score** (higher = better match).  
	✅ `rmsd` → **Root Mean Square Deviation** (lower = better match).  
	✅ `query_length` → Length of **query structure**.  
	✅ `n_align` → Number of **aligned residues**.  

### Methods:
	✅ `__repr__()` → Returns **CSV-formatted string representation** of the result.

================================================================================
📌 Process Overview
================================================================================

1️⃣ **Create an Output Directory**
	- Ensures that `save_directory` exists before writing results.

2️⃣ **Collect & Filter Hits**
	- Iterates over **alignment result files** (`tmalign-results.txt`).
	- Applies **two filters**:
		✅ **Filter 1** → TM-score ≥ `min_tm_score`.  
		✅ **Filter 2** → Query residue count ≥ `min_residues`.  
	- Stores **valid hits** in a dictionary.

3️⃣ **Save Filtered Results**
	- Writes **filtered alignment results** to a **CSV file**.

4️⃣ **Generate Summary Statistics**
	- Writes the total **number of screened models** and **filtered hits**.

5️⃣ **Generate TM-Score Histogram**
	- Categorizes TM-scores into **bins** (e.g., 0.45–0.5, 0.5–0.55, etc.).
	- Saves histogram **as a CSV file**.

================================================================================
🔹 Summary
================================================================================
This script **processes, filters, and organizes** structural alignment results.  
It applies **TM-score and residue count filters**, saves **filtered results**,  
computes **summary statistics**, and generates a **TM-score histogram**.
"""

# Make a directory for the analysis results...
#################################################################################################
result_directories = listdir(result_path + result_directory)
save_directory = save_directory + result_directory[:-1] + result_analysis_tag + "/"
if not exists(save_directory):
	mkdir(save_directory)

# Collect hit results...
# Filter 1: min_tm_score
# Filter 2: min_residues
#################################################################################################
i, j, hits = 0, 0, {}
for d1 in result_directories:
	for f1 in listdir(result_path + result_directory + d1):
		if "tmalign-results.txt" in f1:
			print(result_path + result_directory + d1 + sep + f1)
			print("Pair-wise models screened:", i)
			f = open(result_path + result_directory + d1 + sep + f1, 'r')
			k, lines = 0, f.readlines()
			f.close()
			while k < len(lines) - 1:
				if "#PDBchain1" in lines[k]:
					r = Result(lines[k+1])
					# Filter 1: min_tm_score
					# Now adjusted to match the highest score between TM1 or TM2
					if r.tm_score_target >= r.tm_score_query:
						r.tm_score = r.tm_score_target
					else:
						r.tm_score = r.tm_score_query
					if r.tm_score >= min_tm_score:
						# Filter 2: min_residues
						if r.query_length >= min_residues:
							j += 1
							# print(j, r.query_pdb_file_path, r.tm_score, r.query_length, r.n_align)
							hits[(r.tm_score, r.query_pdb_file_path)] = r
						# else:
						# 	print(j, r.query_pdb_file_path, r.tm_score, r.query_length, r.n_align)
					i += 1
				k += 1

print("Total pair-wise models screened:", i)

# Write csv results to file...
#################################################################################################
f = open(save_directory + result_directory[:-1] + "_results.csv", "w")
f.write("TM score,RMSD,Query length,N align,Query PDB path,Target PDB path\n")
keys = sorted(hits)
for key in keys:
	r = hits[key]
	f.write(str(r))
f.close()

# Write result stats to file...
#################################################################################################
f = open(save_directory + result_directory[:-1] + "_results_stats.txt", "w")
f.write("Total pair-wise models screened: " + str(i) + "\n")
f.write("Number of initial hits: " + str(len(hits)) + "\n")
f.close()

# Write tm_score histogram to file...
#################################################################################################
bins = {
	
		(0.45, 0.5): 0,
		(0.5, 0.55): 0,
		(0.55, 0.60): 0,
		(0.60, 0.65): 0,
		(0.65, 0.70): 0,
		(0.70, 0.75): 0,
		(0.75, 0.80): 0,
		(0.80, 0.85): 0,
		(0.85, 0.90): 0,
		(0.90, 0.95): 0,
		(0.95, 1.0): 0,

		}
keys, bin_keys = sorted(hits), sorted(bins)
for key in keys:
	r = hits[key]
	if r.tm_score >= 0.95:
		bins[(0.95, 1.0)] += 1
	for bin_key in bin_keys[:-1]:
		if r.tm_score >= bin_key[0] and r.tm_score < bin_key[1]:
			bins[bin_key] += 1
f = open(save_directory + result_directory[:-1] + "_results_histogram.csv", "w")
f.write("Bin,Occurrences\n")
for bin_key in bin_keys:
	bin_string = str(bin_key[0]) + "-" + str(bin_key[1])
	f.write(bin_string + "," + str(bins[bin_key]) + "\n")
f.close()

















