# Copyright (c) 2025 Isom Lab
# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.
# If a copy of the MPL was not distributed with this file, You can obtain one at
# https://mozilla.org/MPL/2.0/.

"""
=========================================================================================================
#  Script: **Pairwise Structural Alignment with TM-align**
#  Description: 
#  This script **performs large-scale pairwise structure alignment** of a **query structure**  
#  against **214,528,851 AlphaFold Database predictions**. It distributes computations  
#  across **1000 compute nodes** on the **Triton HPC cluster**.

=========================================================================================================

📌 **CONFIGURATION PARAMETERS**
---------------------------------------------------------------------------------------------------------
🔹 **Cluster Information**
    - `cluster` → Specifies the cluster (**Triton** or **Pegasus**).
    - `queue` → Job submission queue (**short** or **normal**).

🔹 **Save Directories**
    - `save_directory` → Root directory for storing results.
    - `calculation_name` → Subdirectory for this specific calculation.
    - `calculation_save_path` → Full path where results are stored.

🔹 **File Locations**
    - `where_query_pdb_files_are_located` → Directory containing **query PDB structures**.
    - `where_target_cif_set_files_are_located` → Location of **AlphaFold structure submission sets**.
    - `superdark_script` → Path to the **alignment script** (`1-superDark.py`).

🔹 **Job Resubmission**
    - `job_resubmission_file` → If resubmitting failed jobs, this file contains **job numbers**.
    - `job_resubmission_numbers` → List of **job IDs** to be resubmitted.

=========================================================================================================
📌 **DEPENDENCIES**
---------------------------------------------------------------------------------------------------------
- `os`, `shutil`, `subprocess` → Used for **file handling**, **job submission**,  
  **directory creation**, and **monitoring job progress**.

=========================================================================================================
📌 **EXECUTION FLOW**
---------------------------------------------------------------------------------------------------------
🔹 **Step 1: Create Necessary Directories**
    - **Creates or refreshes** directories for:
        1️⃣ **Calculation results**  
        2️⃣ **Log files**  
        3️⃣ **Job submission inputs**  

🔹 **Step 2: Collect Job Submission Sets**
    - If **not resubmitting**, retrieves **all** 1000 job set files.
    - If **resubmitting**, collects **only the failed job set files**.

🔹 **Step 3: Prepare Job Submission Files**
    - Generates **query structure paths** from the **query directory**.
    - Reads **target structure paths** from the **submission sets**.
    - Creates **combined query-target file lists** for TM-align processing.

🔹 **Step 4: Submit Jobs to the Cluster**
    - **Creates job submission files** (each containing query-target pairs).
    - Uses `bsub` to submit jobs to the **batch queue**.
    - Assigns:
        ✅ **Job name** → `"j1-<global_job_number>-<job_number>"`  
        ✅ **Standard output/error logs**  

🔹 **Step 5: Monitor Job Progress & Handle Resubmissions**
    - If `job_resubmission_file` exists, only selected jobs are resubmitted.

=========================================================================================================
📌 **SUMMARY**
---------------------------------------------------------------------------------------------------------
✅ Automates **large-scale TM-align** pairwise structure comparisons.  
✅ Runs **distributed computing** over **1000 nodes**.  
✅ Implements **job monitoring & error handling**.  
✅ Organizes **alignment results** into structured directories.  
✅ Cleans up **temporary calculation files** after completion.  

=========================================================================================================
"""



# Cluster information...
##################################################################################
cluster = "triton"
queue = "normal" # normal or short

# Dependencies...
##################################################################################
from os import mkdir, listdir
from os.path import exists
from shutil import rmtree
import subprocess

# Save locations...
##################################################################################
save_directory = "/projectnb/isomlab/dan/in-and-out/output/alphaFold_uniprot_v4/"
# calculation_name = "2023.11.12.7tmps/"
# calculation_name = "2025.08.20.7tmps-rhod_mglu/"
calculation_name = "2025.09.27.7tmps-foldseek_e0p01/"

calculation_save_path = save_directory + calculation_name
if calculation_name and calculation_name != "/":
	# Create (or delete and create anew) the directory for program output...
	if not exists(calculation_save_path):
		mkdir(calculation_save_path)
	else:
		rmtree(calculation_save_path)
		mkdir(calculation_save_path)

# File locations...
##################################################################################
where_query_pdb_files_are_located = "/projectnb/isomlab/dan/in-and-out/input/alphaFold_v4/alphaFoldQueryStructures-Triton/rhod/"
# where_target_cif_set_files_are_located = "/scratch/alphafold4/submission_sets/"
where_target_cif_set_files_are_located = "/projectnb/isomlab/dan/in-and-out/output/alphaFold_uniprot_v4/2025.09.27.7tmps-foldseek_e0p01_models/submission_sets/"
superdark_script = "/home/dgi5/scripts/proteomes-214M-v4/1-job-code/1-query/1-superDark.py"

# This file is necessary for partial job resubmissions...
# If not a job resubmission, set to None.
##################################################################################
job_resubmission_file = "" # "/home/dgi5/scripts/proteomes-214M-v4/0-job-management/6_job_resubmissions.txt"
job_resubmission_numbers = []
if job_resubmission_file:
	f = open(job_resubmission_file, "r")
	for line in f:
		job_number = line.strip()
		job_resubmission_numbers.append(job_number)
	f.close()

# Collect the job submission sets files...
##################################################################################
set_file_paths = []
if not job_resubmission_file:

	# Collect the 1,000 job submission sets files...
	##################################################################################
	job_submission_files = sorted(listdir(where_target_cif_set_files_are_located))
	for job_submission_file in job_submission_files:
		set_file_paths.append(where_target_cif_set_files_are_located + job_submission_file)

else:

	# Collect the specific job resubmission sets files...
	##################################################################################
	job_submission_files = sorted(listdir(where_target_cif_set_files_are_located))
	for job_submission_file in job_submission_files:
		set_file_number = job_submission_file.split("-")[1].split(".txt")[0]
		if set_file_number in job_resubmission_numbers:
			set_file_paths.append(where_target_cif_set_files_are_located + job_submission_file)

# Submit the jobs...
######################################################################################################

# Comment out if doing a full calculation
# Condition responding to reviewers
# Limits analysis analysis to ~5% of 2.14M, so ~10M
job_numbers = {n: n for n in range(1, 31)}

for set_file_path in set_file_paths:

	# Create save locations...
	##################################################################################################
	global_job_number = set_file_path.split("set-")[1].split(".")[0]
	# Comment out if doing a full calculation
	if int(global_job_number) not in job_numbers:
		continue
	query_label = "job-" + global_job_number
	logResultsDirectory = calculation_save_path + query_label + "/"

	# Create (or delete and create anew) the directory for program output...
	##################################################################################################
	if not exists(logResultsDirectory):
		mkdir(logResultsDirectory)
	else:
		rmtree(logResultsDirectory)
		mkdir(logResultsDirectory)

	# Set the target paths...
	##################################################################################################
	target_paths = []
	f = open(set_file_path, 'r')
	for line in f:
		target_paths.append(line.strip())
	f.close()

	# Set the query paths...
	##################################################################################################
	query_paths = []
	file_names = sorted(listdir(where_query_pdb_files_are_located))
	for file_name in file_names:
		# if "-ssco.pdb" in file_name:
		# 	source_pdb_file_name = file_name
		# 	query_paths.append(where_query_pdb_files_are_located + source_pdb_file_name)
		if ".pdb" in file_name:
			if file_name.endswith("-ssco.pdb"):
				continue
			if file_name.endswith("-s.pdb"):
				continue
			if file_name.endswith("-h.pdb"):
				continue
			source_pdb_file_name = file_name
			query_paths.append(where_query_pdb_files_are_located + source_pdb_file_name)

	# Set the combined paths...
	##################################################################################################
	combined_paths = []
	for query_path in query_paths:
		for target_path in target_paths:
			combined_paths.append(target_path + "*" + query_path)

	# Write the job submission file. Each line contains two models to align...
	##################################################################################################
	job_number = 1
	f_job_path = logResultsDirectory + "job-"+str(job_number)+"-tmalign-submission.txt"
	f_job = open(f_job_path, "w")
	for combined_path in combined_paths:
		f_job.write(combined_path + "\n")
	f_job.close()

	# Create log and error file for the job submission file...
	##################################################################################################
	log_directory = logResultsDirectory + "job-"+str(job_number)+"-tmalign.log"
	error_directory = logResultsDirectory + "job-"+str(job_number)+"-tmalign-error.log"

	# Submit the job...
	##################################################################################################
	print("submission-" + str(global_job_number) + "-" + str(job_number))
	subprocess.run(["bsub", "-q", queue, "-J", "j1-" + str(global_job_number) + "-" + str(job_number), "-o", log_directory, "-e", error_directory, "python", superdark_script, "-a", f_job_path, "-b", cluster, "-c", str(job_number), "-d", logResultsDirectory])

