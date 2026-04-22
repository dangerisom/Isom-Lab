# Copyright (c) 2025 Isom Lab
# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.
# If a copy of the MPL was not distributed with this file, You can obtain one at
# https://mozilla.org/MPL/2.0/.

# Filtered result collection...
################################################################################################
hit_file_path = "/projectnb/isomlab/dan/in-and-out/output/alphaFold_uniprot_v4/2025.09.27.7tmps-foldseek_e0p01/2025.09.27.7tmps-foldseek_e0p01-analysis_foldseek_e0p01/"
hit_file_name = "2025.09.27.7tmps-foldseek_e0p01_results.csv"
aligned_hits_file_path = "/projectnb/isomlab/dan/in-and-out/output/alphaFold_uniprot_v4/2025.09.27.7tmps-foldseek_e0p01/2025.09.27.7tmps-foldseek_e0p01-analysis_foldseek_e0p01/2025.09.27.7tmps-foldseek_e0p01_results-unfiltered_hits/"

# hit_file_path = "/projectnb/isomlab/dan/in-and-out/output/alphaFold_uniprot_v4/2025.09.07.7tmps-frizzled_vs_human_proteome/2025.09.07.7tmps-frizzled_vs_human_proteome-analysis_final/"
# hit_file_name = "2025.09.07.7tmps-frizzled_vs_human_proteome_results.csv"
# aligned_hits_file_path = "/projectnb/isomlab/dan/in-and-out/output/alphaFold_uniprot_v4/2025.09.07.7tmps-frizzled_vs_human_proteome/2025.09.07.7tmps-frizzled_vs_human_proteome-analysis_final/2025.09.07.7tmps-frizzled_vs_human_proteome_results-unfiltered_hits/"

# hit_file_path = "/projectnb/isomlab/dan/in-and-out/output/alphaFold_uniprot_v4/2025.08.20.7tmp_sensitivity_test/frizzled/2025.08.20.7tmps-a2a_frizzled-analysis/"
# hit_file_name = "2025.08.20.7tmps-a2a_frizzled_results.csv"
# aligned_hits_file_path = "/projectnb/isomlab/dan/in-and-out/output/alphaFold_uniprot_v4/2025.08.20.7tmp_sensitivity_test/frizzled/2025.08.20.7tmps-a2a_frizzled-analysis/2025.08.20.7tmps-a2a_frizzled_results-unfiltered_hits/"



# Calculation details...
#################################################################################################
where_query_pdb_files_are_located = "/projectnb/isomlab/dan/in-and-out/input/alphaFold_v4/alphaFoldQueryStructures-Triton/rhod/"
superdark_script = "/home/dgi5/scripts/proteomes-214M-v4/1-job-code/1-query/4-parseHits.py"
number_of_calcs_per_job = 50 # Takes about 10 minutes per 500 jobs.
swell_hull = 0.0
by_surface = 1
job_resubmissions = []

# Choose UM cluster...
# triton or pegasus
cluster = "triton"

# UM cluster information...
# Pegasus queues: debug or general
# Triton queues: short or normal
queue = "normal" 

"""
=========================================================================================================
#  Script: **Filtered Structural Hit Collection**
#  Description: 
#  This script processes **aligned structural hits** by filtering and organizing them for downstream analysis.  
#  It **splits large datasets into manageable chunks**, submits computational jobs to a **high-performance 
#  cluster (Triton)**, and automates the job monitoring and cleanup processes.
=========================================================================================================

📌 **CONFIGURATION PARAMETERS**
---------------------------------------------------------------------------------------------------------
🔹 **File Paths**
    - `hit_file_path` → Path to the **alignment results directory**.
    - `hit_file_name` → Name of the **CSV file** containing UniProt IDs.
    - `aligned_hits_file_path` → Path to **unfiltered aligned hits**.

🔹 **Cluster & Job Management**
    - `cluster` → Specifies the cluster (**Triton** or **Pegasus**).
    - `queue` → Job submission queue (**short** or **normal**).
    - `number_of_calcs_per_job` → Controls the **batch size** of jobs.
    - `superdark_script` → The Python script (`4-parseHits.py`) that processes filtered results.
    - `swell_hull` → Expansion factor for hull-based filtering.
    - `by_surface` → Boolean flag to determine **surface-based** filtering.

🔹 **Job Resubmission**
    - `job_resubmissions` → Allows **resubmission** of specific job numbers (empty if none).

=========================================================================================================
📌 **DEPENDENCIES**
---------------------------------------------------------------------------------------------------------
- `os`, `shutil`, `subprocess`, `time` → Used for **file handling**, **system commands**,  
  **job submission**, and **job monitoring**.

=========================================================================================================
📌 **EXECUTION FLOW**
---------------------------------------------------------------------------------------------------------
🔹 **Step 1: Retrieve and Organize Aligned Hits**
    - Reads **PDB structure files** from the aligned hits directory.
    - Splits them into **subgroups** based on `number_of_calcs_per_job`.

🔹 **Step 2: Create Necessary Directories**
    - **Creates or refreshes** three directories:
        1️⃣ **`input_files_directory/`** → Stores job input files.
        2️⃣ **`save_directory/`** → Stores **filtered structural hits**.
        3️⃣ **`log_directory/`** → Stores **job logs**.

🔹 **Step 3: Generate Job Submission Files**
    - Converts **each subgroup** into a **comma-separated list** of PDB files.
    - Writes **input files** for job submission.

🔹 **Step 4: Submit Jobs to the Cluster**
    - **Uses `bsub`** to submit jobs to the **batch queue**.
    - Assigns:
        ✅ **Job name** → `"j3-<job_number>"`.  
        ✅ **Standard output/error logs**.  

🔹 **Step 5: Monitor Jobs & Handle Resubmissions**
    - If `job_resubmissions` is **not empty**, only selected jobs are resubmitted.

=========================================================================================================
📌 **SUMMARY**
---------------------------------------------------------------------------------------------------------
✅ Automates **filtered collection** of **aligned protein structures**.  
✅ Runs **batch jobs** on **Triton or Pegasus cluster**.  
✅ Implements **job monitoring & error handling**.  
✅ Organizes **filtered hits** into structured directories.  
✅ Cleans up **temporary calculation files** after completion.  

=========================================================================================================
"""


# Dependencies...
from os import mkdir, listdir, remove
from os.path import exists
from shutil import rmtree
from subprocess import run
from time import sleep

##################################################################################################
# Evenly divide the list of aligned hit PDB files to be processed...
##################################################################################################

# Get the pdb file paths...
file_names = listdir(aligned_hits_file_path)
file_paths = []
for file_name in file_names:
	if ".pdb" in file_name:
		file_paths.append(aligned_hits_file_path + file_name)

# Evenly the split the list of aligned hit PDB file paths...
myList = file_paths 
N = int(len(myList)/number_of_calcs_per_job)
if not N:
	N = 1
pdb_groups = [myList[(i*len(myList))//N:((i+1)*len(myList))//N] for i in range(N)]

if not job_resubmissions:

	# Make a directory for the job input files..
	input_files_directory = hit_file_path + hit_file_name.split(".csv")[0] + "_input_files_aligned_hits/"
	if input_files_directory != "/":
		if not exists(input_files_directory):
			mkdir(input_files_directory)
		else:
			rmtree(input_files_directory)
			mkdir(input_files_directory)

	# Make a save directory...
	save_directory = hit_file_path + hit_file_name.split(".csv")[0] + "-filtered_hits/"
	if save_directory != "/":
		if not exists(save_directory):
			mkdir(save_directory)
		else:
			rmtree(save_directory)
			mkdir(save_directory)

	# Make a log directory...
	log_directory = hit_file_path + hit_file_name.split(".csv")[0] + "-filtered_hits_logs/"
	if log_directory != "/":
		if not exists(log_directory):
			mkdir(log_directory)
		else:
			rmtree(log_directory)
			mkdir(log_directory)

# Convert each sublist of pdb_groups into a string of comma-delimited pdb_groups paths...
job_number, pdb_input_file_paths = 1, []
for pdb_group in pdb_groups:

	if not pdb_group:
		continue

	pdb_group_string = ""
	for file_path in pdb_group:
		pdb_group_string += file_path + "\n"
	f = open(input_files_directory + str(job_number) + "_input_file_aligned_hits.txt", "w")
	f.write(pdb_group_string)
	f.close()
	pdb_input_file_paths.append(input_files_directory + str(job_number) + "_input_file_aligned_hits.txt")
	job_number += 1

# Farm the job submissions...
################################################################################################################
job_number = 1
for pdb_input_file_path in pdb_input_file_paths:

	if job_resubmissions:
		if job_number not in job_resubmissions:
			job_number += 1
			continue

	i_log_directory = log_directory + "job-"+str(job_number)+".log"
	i_error_directory = log_directory + "job-"+str(job_number)+"-error.log"

	run(["bsub", "-q", queue, "-J", "j3-" + str(job_number), "-o", i_log_directory, "-e", i_error_directory, "python", superdark_script, "-a", where_query_pdb_files_are_located, "-b", pdb_input_file_path, "-c", save_directory, "-d", str(swell_hull), "-e", str(by_surface)])	

	job_number += 1
