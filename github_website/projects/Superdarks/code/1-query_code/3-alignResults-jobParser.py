# Copyright (c) 2025 Isom Lab
# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.
# If a copy of the MPL was not distributed with this file, You can obtain one at
# https://mozilla.org/MPL/2.0/.


# Specific UniProt collection...
#################################################################################################
# hit_file_path = "/projectnb/isomlab/dan/in-and-out/output/alphaFold_uniprot_v4/2023.11.12.7tmps/2023.11.12.7tmps-analysis_final/"
# hit_file_name = "2023.11.12.7tmps_results.csv"

hit_file_path = "/projectnb/isomlab/dan/in-and-out/output/alphaFold_uniprot_v4/2025.09.27.7tmps-foldseek_e0p01/2025.09.27.7tmps-foldseek_e0p01-analysis_foldseek_e0p01/"
hit_file_name = "2025.09.27.7tmps-foldseek_e0p01_results.csv"

# hit_file_path = "/projectnb/isomlab/dan/in-and-out/output/alphaFold_uniprot_v4/2025.08.20.7tmp_sensitivity_test/mglu/2025.08.20.7tmps-rhod_mglu-analysis_test/"
# hit_file_name = "2025.08.20.7tmps-rhod_mglu_results.csv"

# Calculation details...
#################################################################################################
cluster = "triton"
number_of_calcs_per_job = 100 # Takes about 4 minutes per 500 jobs..
superdark_script = "/home/dgi5/scripts/proteomes-214M-v4/1-job-code/1-query/3-alignResults.py"

# UM cluster information...
# Pegasus queues: debug or general
# Triton queues: short or normal
queue = "normal" 

# Time to sleep in seconds between each check for running job(s)...
sleep_wait = 30 

"""
=========================================================================================================
#  Script: **Batch TM-align Structural Alignment**
#  Description: 
#  This script automates the **pairwise structural alignment** of query structures using **TM-align**
#  against a set of **UniProt-derived protein structures**. It distributes the calculations across
#  a **high-performance computing cluster (Triton)**, monitors job execution, logs results, and
#  performs cleanup operations.
=========================================================================================================

📌 **CONFIGURATION PARAMETERS**
---------------------------------------------------------------------------------------------------------
🔹 **File Paths**
    - `hit_file_path` → Directory containing the list of target structures for alignment.
    - `hit_file_name` → Name of the CSV file with UniProt IDs of structures to be aligned.

🔹 **Cluster Configuration**
    - `cluster` → Specifies the cluster (**Triton**).
    - `queue` → Specifies the job queue (**normal**).
    - `number_of_calcs_per_job` → Controls how many calculations are submitted per job.
    - `superdark_script` → The Python script (`3-alignResults.py`) that processes the alignments.

🔹 **Job Monitoring**
    - `sleep_wait` → Defines the waiting time (seconds) between job status checks.

=========================================================================================================
📌 **DEPENDENCIES**
---------------------------------------------------------------------------------------------------------
- `os`, `shutil`, `subprocess`, `time` → Used for file handling, system commands, and cluster job management.

=========================================================================================================
📌 **CLASSES**
---------------------------------------------------------------------------------------------------------
🔹 **`Result`**
    - Represents a **pairwise alignment result**.
    - **Attributes:**
        - `tm_score`: **Template modeling score** (higher = better match).
        - `rmsd`: **Root Mean Square Deviation** (lower = better match).
        - `query_length`: Length of the **query** structure.
        - `n_align`: Number of **aligned residues**.
        - `query_pdb_file_path`: File path for the **query** structure.
        - `target_pdb_file_path`: File path for the **target** structure.
    - **Methods:**
        - `__repr__()`: Returns a **CSV-formatted** string of the alignment result.

=========================================================================================================
📌 **EXECUTION FLOW**
---------------------------------------------------------------------------------------------------------
🔹 **Step 1: Load and Parse Structural Matches**
    - Reads **target structures** from the CSV file (`hit_file_name`).
    - Creates a list of **`Result`** objects.

🔹 **Step 2: Split Targets into Job Groups**
    - Divides the list into **evenly-sized** groups.
    - Ensures each **job processes `number_of_calcs_per_job` structures**.

🔹 **Step 3: Create Job Submission Directories**
    - Creates **three directories**:
        1️⃣ **`input_files_directory/`** → Stores job input files.
        2️⃣ **`save_directory/`** → Stores alignment results.
        3️⃣ **`log_directory/`** → Stores job logs.

🔹 **Step 4: Generate Job Submission Files**
    - Writes the **pairwise structure pairs** into job-specific input files.
    - Prepares the **input file paths** for batch submission.

🔹 **Step 5: Submit Jobs to the Cluster**
    - **Uses `bsub`** to submit each job to the Triton **batch queue**.
    - Assigns:
        ✅ Job name → `"j3-<job_number>"`.
        ✅ Standard output/error logs.

🔹 **Step 6: Monitor Job Completion**
    - Uses **`bjobs`** to check for running jobs.
    - **Waits (`sleep_wait` seconds) before rechecking.**

🔹 **Step 7: Cleanup**
    - **Deletes intermediate rotation matrix files** (`.txt`) generated by TM-align.

=========================================================================================================
📌 **SUMMARY**
---------------------------------------------------------------------------------------------------------
✅ Automates large-scale **protein structure alignment**.  
✅ Runs **distributed jobs** on the **Triton cluster**.  
✅ Implements **job monitoring & error handling**.  
✅ Cleans up **temporary calculation files** after completion.  

=========================================================================================================
"""



# Dependencies...
#################################################################################################
from os import mkdir, listdir, remove
from os.path import exists
from shutil import rmtree
from subprocess import run
from time import sleep

# Classes...
#################################################################################################
class Result:

	def __init__(self, result_file_line):

		self.result_file_line = result_file_line
		self.data = result_file_line.strip().split(",")
		self.tm_score = float(self.data[0]) # TM score normalized by the target structure
		self.rmsd = float(self.data[1])
		self.query_length = int(self.data[2])
		self.n_align = int(self.data[3])
		self.query_pdb_file_path = self.data[4]
		self.target_pdb_file_path = self.data[5]

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

# Get each initial structure similarity match...
##################################################################################################
result_objects = []
f = open(hit_file_path + hit_file_name, "r")
for line in f:
	try:
		result_object = Result(line)
		result_objects.append(result_object)
	except:
		pass
f.close()

# Evenly the split matches...
##################################################################################################
myList = result_objects
N = int(len(myList)/number_of_calcs_per_job)
if not N:
	N = 1
result_groups = [myList[(i*len(myList))//N:((i+1)*len(myList))//N] for i in range(N)]

# Make a directory for the job input files...
##################################################################################################
input_files_directory = hit_file_path + hit_file_name.split(".csv")[0] + "_input_files/"
if input_files_directory != "/":
	if not exists(input_files_directory):
		mkdir(input_files_directory)
	else:
		rmtree(input_files_directory)
		mkdir(input_files_directory)

# Make a save directory...
##################################################################################################
save_directory = hit_file_path + hit_file_name.split(".csv")[0] + "-unfiltered_hits/"
if save_directory != "/":
	if not exists(save_directory):
		mkdir(save_directory)
	else:
		rmtree(save_directory)
		mkdir(save_directory)

# Make a log directory...
##################################################################################################
log_directory = hit_file_path + hit_file_name.split(".csv")[0] + "-unfiltered_hits_logs/"
if log_directory != "/":
	if not exists(log_directory):
		mkdir(log_directory)
	else:
		rmtree(log_directory)
		mkdir(log_directory)

# Write initial match groups to job submission input files...
##################################################################################################
job_number, result_object_file_paths = 1, []
for result_group in result_groups:

	if not result_group:
		continue

	result_object_group_string = ""
	for result_object in result_group:
		result_object_group_string += str(result_object)
	f = open(input_files_directory + str(job_number) + "_input_file.txt", "w")
	f.write(result_object_group_string)
	f.close()

	# Save the input file path for submission...
	result_object_file_paths.append(input_files_directory + str(job_number) + "_input_file.txt")
	job_number += 1		

# Farm the job submissions...
################################################################################################################
job_number = 1
for result_object_file_path in result_object_file_paths:

	i_log_directory = log_directory + "job-"+str(job_number)+".log"
	i_error_directory = log_directory + "job-"+str(job_number)+"-error.log"

	run(["bsub", "-q", queue, "-J", "j3-" + str(job_number), "-o", i_log_directory, "-e", i_error_directory, "python", superdark_script, "-f", result_object_file_path, "-c", cluster, "-s", save_directory])	

	job_number += 1

# Monitor the jobs and delete intermediate, TMAlign-generated .txt rotation matrix files...
################################################################################################################
jobs_running = 1
while jobs_running:
	t = run(["bjobs"], capture_output=True)
	stdout = t.stdout.decode()
	if "j3-" in stdout:
		print("still running...3-alignResults-jobParser.py")
		sleep(sleep_wait) 
		continue
	else:
		jobs_running = 0

		# Remove temporary mobile files...
		print("Removing temporary calculation files....")
		filenames = listdir(save_directory)
		for filename in filenames:
			for phrase in [".txt"]:
				if phrase in filename:
					print(filename)
					remove(save_directory + filename)
					break