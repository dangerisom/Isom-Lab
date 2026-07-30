# Copyright (c) 2025 Isom Lab
# Released under the MIT License; see LICENSE at the repository root.

f_job_path = ""
cluster = ""
job_number = ""
logResultsDirectory = ""

import sys, getopt
# Get command line input, output, and fail file arguments...
argv = sys.argv[1:]
# Try to collect the program options and arguments...
try:
	opts, args = getopt.getopt(argv,"a:b:c:d:",["f_job_path=", "cluster=", "job_number=", "logResultsDirectory="])
except getopt.GetoptError:
	print("Argument error....exiting submission.")
	sys.exit(2)

# Update parameters from the submitted options...
for opt, arg in opts:

	if opt in ("-a", "--f_job_path"):
		f_job_path = arg
	if opt in ("-b", "--cluster"):
		cluster = arg
	if opt in ("-c", "--job_number"):
		job_number = arg
	if opt in ("-d", "--logResultsDirectory"):
		logResultsDirectory = arg

# Set the location of the necessary binary locations on the cluster
if cluster == "pegasus":
	tmAlignBin = "/nethome/dgi5/compiledSoftware/tmalign/TMalign"
else:
	tmAlignBin = "/home/dgi5/compiledSoftware/tmalign/TMalign"

# Python dependencies...
import subprocess

##################################################################################################
# Structural alignment functions...
##################################################################################################

def runTMalign(mobile_path, query_path, f_output):

	"""
	Function:
	---------
	runTMalign()

	Description:
	------------
	This function **executes TM-align**, aligning a **mobile (target) structure**  
	to a **query structure**, and saves the **alignment output** to a file.

	================================================================================
	📌 Arguments
	================================================================================

	mobile_path : str
		Path to the **mobile (target) structure file** (PDB or CIF format)  
		which will be aligned to the query structure.

	query_path : str
		Path to the **query structure file** (PDB or CIF format)  
		that serves as the **reference** for alignment.

	f_output : file object
		A **file object** where the TM-align output will be **saved**.  
			📌 The file must be **opened in write mode** before calling this function.  

	================================================================================
	📌 Process Overview
	================================================================================

	1️⃣ **Execute TM-align**
		- Runs `tmAlignBin` to **align `mobile_path` to `query_path`**  
		  using the **default TM-align scoring function**.

	2️⃣ **Save Output**
		- Redirects TM-align's **alignment output** to `f_output`.

	================================================================================
	📌 TM-align Output Format (`-outfmt 2`)
	================================================================================
	This function **uses TM-align's output format 2**, which provides:
		✅ **TM-score** (template modeling score).  
		✅ **RMSD** (root mean square deviation).  
		✅ **Sequence alignment information**.  
		✅ **Rotation matrix** for structural superposition.  

	================================================================================
	🔹 Summary
	================================================================================
	This function **runs TM-align**, aligns a **mobile structure** to a **query structure**,  
	and writes the alignment results to a **specified output file**.
	"""

	subprocess.run([tmAlignBin, mobile_path, query_path, "-outfmt", "2"], stdout=f_output)

#################################################################################################
#
# Main code...
#
#################################################################################################

"""
Description:
------------
This script **submits structural alignment jobs** using **TM-align**,  
aligning **mobile (target) structures** to **query structures**  
on a **cluster computing environment**.

================================================================================
📌 Command-Line Arguments
================================================================================

f_job_path : str
	Path to the **job file**, which contains **pairs of structures**  
	(mobile and query) separated by `*`.

cluster : str
	Specifies the **cluster name** where the job is being executed.  
		✅ `"pegasus"` → Uses TM-align located at `/nethome/dgi5/compiledSoftware/tmalign/TMalign`  
		✅ **Other clusters** → Uses TM-align at `/home/dgi5/compiledSoftware/tmalign/TMalign`  

job_number : str
	Unique **job identifier**, used to label the **alignment output file**.

logResultsDirectory : str
	Directory where the **TM-align results file** will be saved.

================================================================================
📌 Process Overview
================================================================================

1️⃣ **Parse Command-Line Arguments**
	- Uses `getopt` to retrieve:
		✅ `f_job_path` (input file path)  
		✅ `cluster` (computing cluster name)  
		✅ `job_number` (job identifier)  
		✅ `logResultsDirectory` (where results will be saved)  

2️⃣ **Set TM-align Binary Path**
	- Selects the **TM-align binary location** based on the cluster.

3️⃣ **Read the Job File**
	- Opens `f_job_path` and reads **pairs of structure files**  
	  (mobile and query paths separated by `*`).

4️⃣ **Run TM-align for Each Structure Pair**
	- Calls `runTMalign()` to:
		✅ Align the **mobile structure** to the **query structure**.  
		✅ Save the alignment results in **logResultsDirectory**.  

5️⃣ **Save the Alignment Results**
	- Writes the **alignment output** to `job-<job_number>-tmalign-results.txt`.

================================================================================
📌 Function: runTMalign()
================================================================================
This function **executes TM-align**, aligning a **mobile structure**  
to a **query structure**, and saves the **alignment output** to a file.

### Arguments:
	✅ `mobile_path` → Path to **mobile (target) structure file**.  
	✅ `query_path` → Path to **query structure file**.  
	✅ `f_output` → File object to save **TM-align output**.  

### Process:
	1️⃣ **Runs TM-align** with `-outfmt 2`.  
	2️⃣ **Redirects TM-align output** to `f_output`.  

================================================================================
🔹 Summary
================================================================================
This script **executes TM-align jobs** on a **cluster**,  
aligning **multiple protein structures** and saving the results  
in a **log directory** for further analysis.
"""


##################################################################################################
# Align the structures 
##################################################################################################
	
f = open(f_job_path, "r")
lines = f.readlines()
f.close()
f_output_path = logResultsDirectory + "job-"+str(job_number)+"-tmalign-results.txt"
f_output = open(f_output_path, "w")
for line in lines:
	line = line.strip()
	mobile_path = line.split("*")[0]
	query_path = line.split("*")[1]
	runTMalign(mobile_path, query_path, f_output)
f_output.close()


