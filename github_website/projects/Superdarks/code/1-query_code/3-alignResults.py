# Copyright (c) 2025 Isom Lab
# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.
# If a copy of the MPL was not distributed with this file, You can obtain one at
# https://mozilla.org/MPL/2.0/.

result_object_file_path = ""
cluster = ""
save_directory = ""

import sys, getopt
# Get command line input, output, and fail file arguments...
argv = sys.argv[1:]
# Try to collect the program options and arguments...
try:
	opts, args = getopt.getopt(argv,"f:c:s:",["result_object_file_path=", "cluster=", "save_directory="])
except getopt.GetoptError:
	print("Argument error....exiting submission.")
	sys.exit(2)

# Update parameters from the submitted options...
for opt, arg in opts:

	if opt in ("-f", "--result_object_file_path"):
		result_object_file_path = arg
	if opt in ("-c", "--cluster"):
		cluster = arg
	if opt in ("-s", "--save_directory"):
		save_directory = arg

# Set the location of the necessary binary locations on the cluster
if cluster == "pegasus":
	tmAlignBin = "/nethome/dgi5/compiledSoftware/tmalign/TMalign"
else:
	tmAlignBin = "/home/dgi5/compiledSoftware/tmalign/TMalign"

# Dependencies...
#################################################################################################
from os import sep
from pdbFile import PDBfile as PDBfile
from pdbFile_cif import PDBfile as PDBfile_cif
from subprocess import run

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
	from a **CSV-formatted result file line**.

	================================================================================
	📌 Attributes
	================================================================================

	result_file_line : str
		The **original file line** from the results CSV file.

	data : list
		A list of **parsed values** extracted from `result_file_line`.

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

	target_pdb_file_path : str
		Path to the **target PDB file** used in structural comparison.

	query_pdb_file_path : str
		Path to the **query PDB file** used in structural comparison.

	================================================================================
	📌 Class Overview
	================================================================================
	This class **reads a CSV result file line**, extracts alignment data,  
	and stores **structural comparison metrics** (TM-score, RMSD, alignment length)."""

	def __init__(self, result_file_line):

		# File line in .csv format...
		self.result_file_line = result_file_line
		self.data = result_file_line.strip().split(",")
		self.tm_score = float(self.data[0]) # TM score normalized by the target structure
		self.rmsd = float(self.data[1])
		self.query_length = int(self.data[2]) 
		self.n_align = int(self.data[3]) 
		self.target_pdb_file_path = self.data[4]
		self.query_pdb_file_path = self.data[5]

# Functions...
#################################################################################################
def useMatrixToAlign(score_string, mobile_path, query_path, save_directory):

	"""
	Arguments:
	----------
	score_string : str
		A string representing the **alignment score** (e.g., TM-score or RMSD),  
		used to label the output **aligned structure file**.

	mobile_path : str
		Path to the **mobile structure file** (PDB or CIF format)  
		which will be **aligned** to the query structure.

	query_path : str
		Path to the **query structure file** (PDB or CIF format)  
		that serves as the **reference** for alignment.

	save_directory : str
		Directory where the **aligned structure file** and rotation  
		matrix file will be saved.

	Returns:
	--------
	None  
		This function **saves the transformed mobile structure**  
		as an **aligned PDB file** in `save_directory`.

	================================================================================
	📌 Function Overview
	================================================================================
	This function **aligns a mobile structure to a query structure**  
	using **TM-align** and applies the **rotation matrix** to  
	transform the coordinates of the **mobile structure**.

	================================================================================
	🔹 Main Steps
	================================================================================

	1️⃣ **Extract File Names**
		- Extracts **base names** from `mobile_path` and `query_path`,  
		  removing `.pdb`, `.cif`, or `-ssco.pdb` extensions.

	2️⃣ **Run TM-align to Compute the Rotation Matrix**
		- Calls `tmAlignBin` to **align the mobile structure**  
		  to the query and **save the transformation matrix**.

	3️⃣ **Extract the Rotation Matrix from TM-align Output**
		- Reads the **rotation matrix file** generated by TM-align.  
		- Parses **translation (t1, t2, t3)** and **rotation matrix (u11 - u33)** values.

	4️⃣ **Load the Mobile Structure**
		- Opens the **untransformed mobile structure** (PDB or CIF).  
		- Extracts **atomic coordinates** for transformation.

	5️⃣ **Apply the Rotation Matrix**
		- Uses the **computed transformation matrix** to update the **(x, y, z) coordinates**  
		  of the **mobile structure atoms**.

	6️⃣ **Save the Transformed Structure**
		- Writes the **transformed mobile structure** to a new **PDB file**,  
		  labeled with the **alignment score**.

	================================================================================
	🔹 Summary
	================================================================================
	This function **aligns a mobile protein structure to a query structure**  
	using **TM-align** and applies a **rotation matrix** to transform its coordinates.  
	The aligned structure is **saved as a PDB file** in the specified `save_directory`.
	"""

	if ".cif" in mobile_path:
		mobile_name = mobile_path.split("/")[-1].split(".cif")[0]
	else:
		mobile_name = mobile_path.split("/")[-1].split(".pdb")[0]

	if "-ssco.pdb" in query_path:
		query_name = query_path.split("/")[-1].split("-ssco.pdb")[0]
	elif ".cif" in query_path:
		query_name = query_path.split("/")[-1].split(".cif")[0]
	else:
		query_name = query_path.split("/")[-1].split(".pdb")[0]

	matrix_file_name = mobile_name + "--" + query_name + ".txt"
	run([tmAlignBin, mobile_path, query_path, "-m", save_directory + matrix_file_name])

	# Open the rotation matrix file...
	f = open(save_directory + matrix_file_name, "r")
	l = f.readlines()
	r1s, r2s, r3s = l[2], l[3], l[4]
	r1, r2, r3 = r1s.split()[1:], r2s.split()[1:], r3s.split()[1:]
	f.close()

	#  -------- Rotation matrix to rotate Chain_1 to Chain_2 ------
	#  m          t(m)         u(m,1)         u(m,2)         u(m,3)
	#  1     59.6536152118  -0.0731286234   0.6894102883   0.7206702844
	#  2     10.3013940839  -0.6162468656   0.5369168133  -0.5761598182
	#  3     -9.3468354017  -0.7841504990  -0.4862445783   0.3855829419
	t1, t2, t3 = float(r1[0]), float(r2[0]), float(r3[0])
	u11,u21, u31 = float(r1[1]), float(r2[1]), float(r3[1])
	u12, u22, u32 = float(r1[2]), float(r2[2]), float(r3[2])
	u13, u23, u33 = float(r1[3]), float(r2[3]), float(r3[3])

	# Using the rotation matrix from the -ssco.pdb query and mobile alignment...
	# ... transform the coordinates of the mobile structure into alignment...
	###########################################################################

	# Open the untransformed mobile structure....
	pdb_file_path_mobile = sep.join(mobile_path.split(sep)[:-1]) + sep
	pdb_file_name_mobile = mobile_path.split(sep)[-1]
	if ".cif" in pdb_file_name_mobile:
		mobile_protein = PDBfile_cif(pdb_file_path_mobile, pdb_file_name_mobile)
	else:
		mobile_protein = PDBfile(pdb_file_path_mobile, pdb_file_name_mobile)
	mobile_atom_dict = mobile_protein.atoms
	mobile_atom_keys = sorted(mobile_atom_dict)
	mobileAtoms = []
	for mobile_atom_key in mobile_atom_keys:
		mobileAtoms.append(mobile_atom_dict[mobile_atom_key])

	# Transform the coordinates of the mobile structure to the query structure...
	# Label the transformed mobile structure as chain A...
	for atom in mobileAtoms:
		x, y, z = atom.x, atom.y, atom.z
		rx = t1 + u11*x + u12*y + u13*z
		ry = t2 + u21*x + u22*y + u23*z
		rz = t3 + u31*x + u32*y + u33*z
		atom.x, atom.y, atom.z = rx, ry, rz
		atom.chain_identifier = "A"

	# Write the transformed mobile structure to file...
	pdb_file_name = mobile_name + "--" + query_name + score_string + ".pdb"
	f = open(save_directory + pdb_file_name, "w")
	for atom in mobileAtoms:
		if ".cif" in pdb_file_name_mobile:
			f.write(atom.get_pdb_format())
		else:
			f.write(str(atom))
	f.close()

#################################################################################################
#
# Main code...
#
#################################################################################################

"""
Description:
------------
This script **aligns multiple mobile (target) structures** to their respective **query structures**  
using **rotation matrices** computed from TM-align.

================================================================================
📌 Process Overview
================================================================================

1️⃣ **Open the Result File**
	- Reads `result_object_file_path`, which contains **alignment results**  
	  in CSV format.

2️⃣ **Parse Each Alignment Result**
	- Iterates through each **line in the results file**.
	- Creates a **`Result` object** to extract:
		✅ **TM-score**  
		✅ **Target (mobile) structure file path**  
		✅ **Query structure file path**  

3️⃣ **Generate an Alignment Score Label**
	- Constructs `score_string` from the **TM-score** to label the  
	  aligned structure files.

4️⃣ **Align Each Mobile Structure**
	- Calls `useMatrixToAlign()` to:
		✅ Apply the **rotation matrix** from TM-align.  
		✅ Transform the **mobile structure coordinates**.  
		✅ Save the **aligned structure** in `save_directory`.

5️⃣ **Close the Result File**
	- Ensures the results file is properly closed after processing.

================================================================================
🔹 Summary
================================================================================
This script **iterates over a set of alignment results**,  
aligns each **mobile structure to its query**, and saves  
the **transformed structures** in `save_directory`.
"""

# The terminology mobile and target structure are synonymous...
f = open(result_object_file_path, "r")
for line in f:

	r = Result(line)
	score_string = "--" + str(r.tm_score)
	useMatrixToAlign(score_string, r.target_pdb_file_path, r.query_pdb_file_path, save_directory)
	
f.close()
