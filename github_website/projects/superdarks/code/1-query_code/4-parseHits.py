# Copyright (c) 2025 Isom Lab
# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.
# If a copy of the MPL was not distributed with this file, You can obtain one at
# https://mozilla.org/MPL/2.0/.

where_query_pdb_files_are_located = ""
pdb_input_file_path = ""
save_directory = ""
swell_hull = 0.0
by_surface = 0

##################################################################################################
# If python is run from command line with pdbDirectoryPath option...
# This approach is used when submitting bsub python jobs on the clusters...
##################################################################################################

# Enables pdbDirectoryPath to be assigned by command land, if desired... 
import sys, getopt
jobNumber=""
# Get command line input, output, and fail file arguments...
argv = sys.argv[1:]
# Try to collect the program options and arguments...
try:
	opts, args = getopt.getopt(argv,"a:b:c:d:e:", ["where_query_pdb_files_are_located=", "pdb_input_file_path=", "save_directory=", "swell_hull=", "by_surface="])
except getopt.GetoptError:
	print("Argument error....exiting submission.")
	sys.exit(2)

# Update parameters from the submitted options...
for opt, arg in opts:

	if opt in ("-a", "--where_query_pdb_files_are_located"):
		where_query_pdb_files_are_located = arg
	if opt in ("-b", "--pdb_input_file_path"):
		pdb_input_file_path = arg
	if opt in ("-c", "--save_directory"):
		save_directory = arg
	if opt in ("-d", "--swell_hull"):
		swell_hull = float(arg)
	if opt in ("-e", "--by_surface"):
		by_surface = float(arg)

# Dependencies...
##################################################################################################
from os import sep
from subprocess import run
from pdbFile import PDBfile, PseudoAtom
from convexHull3D_2_0 import convexHull3D
from compGeometry import circumSphere, distance, Vertex4D

# PDB functions...
##################################################################################################

def open_pdb(pdb_file_path, pdb_file_name, pdb_format="pdb", zip_status=0):

	"""Arguments:
	----------
	pdb_file_path : str
		Path to the **PDB file**, specifying the directory where the file is located.

	pdb_file_name : str
		Name of the **PDB file**, including the filename without the path.

	pdb_format : str, optional (default="pdb")
		Specifies the **format of the input file**. Supported formats:
			✅ `"pdb"` → **Protein Data Bank (PDB) format**  
			✅ `"mmCIF"` → **Macromolecular Crystallographic Information File (mmCIF) format**  

	zip_status : int, optional (default=0)
		Indicates whether the PDB file is **compressed**:
			✅ `0` → File is **not zipped**  
			✅ `1` → File is **compressed**  

	Returns:
	--------
	object : 
		A **PDBfile object**, which represents the parsed **protein structure**  
		and contains **residue and atomic information**.

	================================================================================
	📌 Function Overview
	================================================================================
	This function **loads and parses a PDB or mmCIF file**, creating a **PDBfile object**  
	that can be used for **structural analysis and geometric computations**.

	================================================================================
	🔹 Main Steps
	================================================================================

	1️⃣ **Determine File Format**
		- If `pdb_format == "mmCIF"`, imports and initializes **`pdbFile_cif.PDBfile`**.  
		- Otherwise, imports and initializes **`pdbFile.PDBfile`** for standard PDB files.

	2️⃣ **Create a PDBfile Object**
		- Initializes a **PDBfile instance**, which:
			✅ Reads the **atomic structure** from the file.  
			✅ Parses **residue and coordinate data**.  
			✅ Stores **topology information** for downstream analysis.  

	3️⃣ **Return the PDBfile Object**
		- The function returns a **PDBfile object**, which contains all parsed  
		  **protein structure data**, allowing for **further analysis**.

	================================================================================
	🔹 Summary
	================================================================================
	This function **loads a protein structure file (PDB or mmCIF)** and returns  
	a **PDBfile object**, enabling further **topological and structural analysis**.
	"""

	# Create an instance of a PDBfile (really protein) object...
	##############################################################################################
	if pdb_format == "mmCIF":
		import pdbFile_cif
		pdb_file_object = pdbFile_cif.PDBfile(pdb_file_path, pdb_file_name, zip_status=zip_status)
	else:
		import pdbFile
		pdb_file_object = pdbFile.PDBfile(pdb_file_path, pdb_file_name, zip_status=zip_status)	
	return pdb_file_object

def get_ca_atoms(residue_dict, origin=""):

	"""Arguments:
	----------
	residue_dict : dict
		A dictionary of **residues**, where:
			🔹 **Keys** → Residue identifiers (e.g., `(residue_number, chain_id, insertion_code)`).  
			🔹 **Values** → Residue objects containing **atomic data**.

	origin : str, optional (default="")
		A string label assigned to the **C-alpha (`CA`) atoms**, used for **tracking  
		their source** in downstream analysis.  
			✅ If provided, assigns `origin` to `caAtom.symbol`.  
			✅ If empty (`""`), no label is added.  

	Returns:
	--------
	dict : 
		A dictionary of **C-alpha (`CA`) atoms** from the input residues, where:
			✅ **Keys** → Residue identifiers from `residue_dict`.  
			✅ **Values** → Corresponding **C-alpha atoms (`CA`)**.  

	================================================================================
	📌 Function Overview
	================================================================================
	This function **extracts the C-alpha (`CA`) atoms** from a dictionary of residues.  
	It ensures that only residues **containing a C-alpha atom** are included,  
	allowing for **structural and geometric calculations**.

	================================================================================
	🔹 Main Steps
	================================================================================

	1️⃣ **Sort Residue Keys**
		- Ensures that residues are processed **in sequential order**.

	2️⃣ **Extract C-alpha Atoms**
		- Iterates through the **sorted residue dictionary**:
			✅ **Attempts to retrieve the `CA` atom** for each residue.  
			✅ If found, stores it in `caAtoms`.  
			✅ If an `origin` string is provided, assigns it to `caAtom.symbol`.  
			✅ If no `CA` atom is found, the residue is **skipped**.  

	3️⃣ **Return Dictionary of C-alpha Atoms**
		- The function returns a **filtered dictionary** containing **only C-alpha atoms**,  
		  which are essential for **structural alignment and analysis**.

	================================================================================
	🔹 Summary
	================================================================================
	This function **isolates the C-alpha (`CA`) atoms** from a set of residues,  
	ensuring that only **backbone-representative atoms** are retained for  
	further **geometric processing and topology calculations**.
	"""

	residue_keys = sorted(residue_dict)
	caAtoms = {}
	i = 0
	while i < len(residue_keys):
		key = residue_keys[i]
		try:
			# Collect the C alpha atom of the side chain...
			caAtom = residue_dict[key].atoms["CA"]
			if caAtom:
				if origin:
					caAtom.symbol = origin
				caAtoms[key] = caAtom 
		except:
			# No C alpha atom for side chain...
			pass
		i += 1
	return caAtoms

##################################################################################################
# Computational geometry functions...
##################################################################################################

def trim(convex_hull, pdb_file_path, pdb_file_name, pdb_format="pdb", zip_status=0):

	"""Arguments:
	----------
	convex_hull : object
		A **ConvexHull3D** object representing the **convex hull of the query structure**.  
		This is used to determine which residues from the **hit structure** are **inside the hull**.

	pdb_file_path : str
		Path to the **hit PDB file**, specifying the directory where the file is located.

	pdb_file_name : str
		Name of the **hit PDB file**, including the filename without the path.

	pdb_format : str, optional (default="pdb")
		Specifies the **format of the input file**. Supported formats:
			✅ "pdb" → **Protein Data Bank (PDB) format**
			✅ "cif" → **Crystallographic Information File (CIF) format**

	zip_status : int, optional (default=0)
		Indicates whether the PDB file is **zipped**:
			✅ `0` → File is **not zipped**  
			✅ `1` → File is **compressed**  

	Returns:
	--------
	tuple : 
		A tuple containing:
			✅ `residuesInHull` → A dictionary of **hit residues inside the convex hull**  
				(matching the query structure's topology).  
			✅ `residuesAll` → A dictionary of **all hit structure residues**,  
				including those outside the convex hull.  

	================================================================================
	📌 Function Overview
	================================================================================
	This function **trims a hit protein structure** by filtering **only those residues**  
	that **fall inside the convex hull** of the query structure.  
	It ensures that only **topologically relevant residues** are retained for analysis.

	================================================================================
	🔹 Main Steps
	================================================================================

	1️⃣ **Load the Hit Structure**
		- Opens the **PDB/CIF file** and extracts **residue and atomic information**.
		- Collects **all residue keys** from the **hit structure**.

	2️⃣ **Filter Residues Using Convex Hull**
		- Iterates through **each residue** in the **hit structure**:
			✅ **Retrieves the C-alpha atom** of the residue.  
			✅ **Tests if the atom is inside the convex hull** using `test_if_this_vertex_is_in_the_hull()`.  
			✅ **If inside**, the residue is stored in `residuesInHull`.  
			✅ **All residues are stored in `residuesAll`** (including those outside the hull).  

	3️⃣ **Return Trimmed and Unfiltered Residues**
		- Returns:
			✅ **Filtered residues (`residuesInHull`)** → Only residues **inside the convex hull**.  
			✅ **Unfiltered residues (`residuesAll`)** → The complete set of hit residues.  

	================================================================================
	🔹 Summary
	================================================================================
	This function **removes extraneous hit structure residues** by applying **convex hull constraints**.  
	It ensures **alignment to the query topology** by keeping **only residues inside the convex hull**.
	"""

	# 1) Open the hit PDB file...
	##############################################################################################
	protein = open_pdb(pdb_file_path, pdb_file_name, pdb_format=pdb_format, zip_status=zip_status)
	residueKeys, mids = sorted(protein.residues), []

	# 2) Identify the subset of hit residues inside the convex hull of the query structure...
	##############################################################################################
	residueKeys, residuesInHull, residuesAll = sorted(protein.residues), {}, {}
	for residueKey in residueKeys:

		# Try to get the c-alpha atom for the hit residue...
		try:
			ca = protein.residues[residueKey].atoms["CA"]
		except:
			continue

		# Test if the c-alpha atom of the hit residue is inside the convex hull...
		result = convex_hull.test_if_this_vertex_is_in_the_hull(ca.v)
		if result:
			residuesInHull[residueKey] = protein.residues[residueKey]
		residuesAll[residueKey] = protein.residues[residueKey]

	return (residuesInHull, residuesAll)

def get_vertex_4d(atom):

	"""Arguments:
	----------
	atom : object
		A **PDBfile.atom** object representing a single atomic coordinate 
		from a protein structure.

	Returns:
	--------
	Vertex4D : object
		A **4D vertex representation** of the given atom, where:
			✅ (x, y, z) → Original **3D coordinates** of the atom.  
			✅ `u = x² + y² + z²` → The computed **fourth coordinate** for geometric processing.  
			✅ `vertex.data` → Stores the **original atomic data**.  
			✅ `vertex.id` → Stores the **atom serial number**.  

	================================================================================
	📌 Function Overview
	================================================================================
	This function **converts a 3D atomic coordinate** into a **4D vertex representation**  
	by adding an extra coordinate, `u = x² + y² + z²`. This transformation is commonly  
	used in **geometric analysis**, including **convex hull calculations**.

	================================================================================
	🔹 Main Steps
	================================================================================

	1️⃣ **Extract Atomic Coordinates**
		- Retrieves **x, y, z** coordinates from the given **PDBfile.atom** object.

	2️⃣ **Compute 4D Coordinate (`u`)**
		- Calculates `u` using the formula:  
		  \t📌 `u = x² + y² + z²`
		- This mapping enables **higher-dimensional geometric operations**.

	3️⃣ **Create a `Vertex4D` Object**
		- Initializes a **`Vertex4D` object** with the calculated (x, y, z, u) coordinates.
		- Stores:
			✅ **Original atom data** (`vertex.data = atom`).  
			✅ **Atom serial number** (`vertex.id = atom.atom_serial`).  

	4️⃣ **Return the `Vertex4D` Object**
		- The function returns a **4D vertex**, useful for **topological and geometric analyses**.

	================================================================================
	🔹 Summary
	================================================================================
	This function **transforms protein atomic coordinates into 4D space** for  
	**structural and geometric computations**, including **surface modeling**  
	and **fold topology analysis**.
	"""

	x, y, z= atom.x, atom.y, atom.z
	u = atom.x**2 + atom.y**2 + atom.z**2
	vertex = Vertex4D((x,y,z,u),setC=1)
	vertex.data = atom
	vertex.id = atom.atom_serial
	return vertex

def calculate_query_surface(pdb_file_path, pdb_file_name, pdb_format="pdb", zip_status=0):

	"""Arguments:
	----------
	pdb_file_path : str
		Path to the **query PDB file**, specifying the directory where the file is located.

	pdb_file_name : str
		Name of the **query PDB file**, including the filename without the path.

	pdb_format : str, optional (default="pdb")
		Specifies the **format of the input file**. Supported formats:
			✅ "pdb" → **Protein Data Bank (PDB) format**
			✅ "cif" → **Crystallographic Information File (CIF) format**

	zip_status : int, optional (default=0)
		Indicates whether the PDB file is **zipped**:
			✅ `0` → File is **not zipped**  
			✅ `1` → File is **compressed**  

	Returns:
	--------
	tuple : 
		A tuple containing:
			✅ `query_surface` → A **pHinder-calculated surface representation**  
				of the query structure, including its **facet geometry**.
			✅ `query_ca_atoms` → A list of **C-alpha atoms** extracted from  
				the query structure.
			✅ `query_residues_all` → A dictionary of **query residues**,  
				mapped by residue keys.

	================================================================================
	📌 Function Overview
	================================================================================
	This function **calculates the 3D molecular surface** of a query protein structure.  
	It uses **convex hull-based geometry** to create a **pHinder-calculated surface representation**,  
	which can be used for **structural comparison and filtering**.

	================================================================================
	🔹 Main Steps
	================================================================================

	1️⃣ **Load the Query Structure**
		- Opens the **PDB/CIF file** and extracts **residue and atomic information**.
		- Collects **C-alpha (CA) atoms** from the structure.

	2️⃣ **Convert Query Residues to 4D Vertices**
		- Maps **C-alpha atoms** into **4D space** for geometric analysis.

	3️⃣ **Calculate the Query Surface**
		- Uses **pHinder's `calculateSurface` function** to compute a **convex hull-based  
		  surface representation** of the query structure.
		- Uses a **circumsphere radius of 8.0 Å** to ensure **higher permissiveness** in surface detection.

	4️⃣ **Write Surface Representation to PyMOL**
		- Saves the **calculated surface** as a **PyMOL visualization script**  
		  for further analysis.

	5️⃣ **Return Surface and Query Residues**
		- Returns the **computed surface representation**, **C-alpha atoms**,  
		  and **query residue dictionary**.

	================================================================================
	🔹 Summary
	================================================================================
	This function **generates a 3D molecular surface** for a query structure, allowing  
	for **topological analysis, structural comparisons, and surface-based filtering**.  
	It is particularly useful for **aligning protein structures** using geometric constraints.
	"""

	# Try surface
	from compGeometry import Vertex4D
	from pHinderSurface import calculateSurface, inLocalSurface, writeSurface_as_pymol_script

	# 1) Calculate the convex hull of the query structure...
	##############################################################################################
	query_protein = open_pdb(pdb_file_path, pdb_file_name, pdb_format=pdb_format, zip_status=zip_status)
	query_residues_all = query_protein.residues
	query_ca_atoms = get_ca_atoms(query_residues_all)

	query_vertices = [] 
	for residue_key in query_protein.residues:
		try:
			# ca_vertex = query_vertices.residues[residue_key].atoms["CA"].v
			ca_vertex_4d = get_vertex_4d(query_protein.residues[residue_key].atoms["CA"])
			query_vertices.append(ca_vertex_4d)
		except:
			pass

	# Normally circumsphere is 6.5 for protein surface in other pHinder calculations...
	# I am using 8.0 here to be more permissive...
	# However, being more permissive results in surfaces with fewer facets...
	# Therefore, the allowSmallSurfaces optional argument must be turned on... 
	query_surface = calculateSurface(query_vertices, 8.0, highResolutionSurface=0, allowSmallSurfaces=1)

	writeSurface_as_pymol_script(query_surface[1], pdb_file_name.split(".pdb")[0], pdb_file_path)

	return(query_surface, query_ca_atoms, query_residues_all)

def trim_using_surface(query_surface, core_cutoff, margin_cutoff, pdb_file_path, pdb_file_name, pdb_format="pdb", zip_status=0):

	"""Arguments:
	----------
	query_surface : object
		A **pHinder-calculated quad surface** representation of the query structure,  
		used for **spatial filtering** of hit residues.

	core_cutoff : float
		The **distance threshold (in Ångströms)** defining the **core region**  
		of the query surface. Residues **inside this threshold** are considered **deeply embedded**.

	margin_cutoff : float
		The **distance threshold (in Ångströms)** defining the **margin region**  
		of the query surface. Residues **outside both core and margin thresholds** are excluded.

	pdb_file_path : str
		Path to the **hit PDB file**, specifying the directory where the file is located.

	pdb_file_name : str
		Name of the **hit PDB file**, including the filename without the path.

	pdb_format : str, optional (default="pdb")
		Specifies the **format of the input file**. Supported formats:
			✅ "pdb" → **Protein Data Bank (PDB) format**
			✅ "cif" → **Crystallographic Information File (CIF) format**

	zip_status : int, optional (default=0)
		Indicates whether the PDB file is **zipped**:
			✅ `0` → File is **not zipped**  
			✅ `1` → File is **compressed**  

	Returns:
	--------
	tuple : 
		A tuple containing:
			✅ `hit_residues_in_hull` → A dictionary of **filtered hit residues** that pass  
				the **query surface-based trimming criteria**.
			✅ `hit_residues_all` → A dictionary containing **all unfiltered residues**  
				from the **hit structure**.

	================================================================================
	📌 Function Overview
	================================================================================
	This function **trims hit structure residues** based on their **proximity to the query surface**.  
	Using **core and margin cutoffs**, it removes **residues that do not align well** with the query topology.

	================================================================================
	🔹 Main Steps
	================================================================================

	1️⃣ **Load the Hit Structure**
		- Opens the **hit PDB/CIF file** and extracts **residue information**.
		- Collects **C-alpha (CA) atoms** for spatial analysis.

	2️⃣ **Convert Hit Residues to 4D Vertices**
		- Maps **C-alpha atoms** into **4D space** for better geometric classification.

	3️⃣ **Classify Residues Using Query Surface**
		- Uses **`inLocalSurface`** from `pHinderSurface` to determine **surface-based classifications**.
			✅ **Classification [1] = 1** → Residue is **outside** the core/margin → REMOVE  
			✅ **Classification [1] ≠ 1** → Residue is **inside** the core/margin → KEEP  

	4️⃣ **Filter Hit Residues**
		- Retains **only residues classified as inside the query surface**.
		- Stores **filtered hit residues** in `hit_residues_in_hull`.

	5️⃣ **Return Trimmed and Unfiltered Residue Dictionaries**
		- Returns **two dictionaries**:
			✅ **Filtered residues (`hit_residues_in_hull`)** → Only residues **inside** the query surface.  
			✅ **Unfiltered residues (`hit_residues_all`)** → The complete set of hit residues.

	================================================================================
	🔹 Summary
	================================================================================
	This function **removes extraneous hit structure residues** by applying **query surface constraints**.  
	It ensures **alignment to the query topology** by keeping **only residues within core/margin thresholds**.
	"""

	# 1) Open the hit PDB file...
	##############################################################################################
	hit_protein = open_pdb(pdb_file_path, pdb_file_name, pdb_format=pdb_format, zip_status=zip_status)
	hit_residues_all = hit_protein.residues
	hit_ca_atoms = get_ca_atoms(hit_residues_all)

	hit_vertices = [] 
	for residue_key in hit_protein.residues:
		try:
			ca_vertex_4d = get_vertex_4d(hit_protein.residues[residue_key].atoms["CA"])
			hit_vertices.append(ca_vertex_4d)
		except:
			pass

	from pHinderSurface import inLocalSurface
	classifications, hit_residues_in_hull = inLocalSurface(hit_vertices, query_surface, core_cutoff, margin_cutoff), {}
	for classification in classifications:
		if classification[1] != 1:
			v_4d =  classification[0]
			hit_residues_in_hull[v_4d.data.residue.key] = v_4d.data.residue

	return (hit_residues_in_hull, hit_residues_all)

def more_trim_using_surface(query_surface, core_cutoff, margin_cutoff, hit_residues_in_surface):

	"""Arguments:
	----------
	query_surface : object
		The **surface representation** of the query structure, used to determine 
		whether residues in the hit structure fall within the defined **core** or **margin** regions.

	core_cutoff : float
		The **distance threshold (in Ångströms)** defining the **core region** 
		of the query surface. Residues within this threshold are considered **deeply embedded**.

	margin_cutoff : float
		The **distance threshold (in Ångströms)** defining the **margin region** 
		of the query surface. Residues falling outside both **core** and **margin** regions 
		are excluded.

	hit_residues_in_surface : dict
		A dictionary of **hit structure residues** mapped to their **C-alpha atomic positions**.  
			🔹 **Keys** → Residue identifiers  
			🔹 **Values** → Residue objects containing **atomic data**  

	Returns:
	--------
	dict : 
		A **filtered dictionary of hit residues** that remain **after surface-based trimming**, where:
			✅ **Residues inside the core/margin** are **retained**.  
			✅ **Residues outside the surface-defined thresholds** are **removed**.  

	================================================================================
	📌 Function Overview
	================================================================================
	This function **trims hit structure residues** based on their **proximity to the query surface**.  
	Using **core and margin cutoffs**, it removes **residues that do not align well** with the query topology.

	================================================================================
	🔹 Main Steps
	================================================================================

	1️⃣ **Extract C-alpha Atom Coordinates**
		- Collects **C-alpha (CA) atom vertices** from `hit_residues_in_surface`.

	2️⃣ **Classify Residues Using Query Surface**
		- Uses the **`inLocalSurface`** function to determine **surface-based classifications**.
			✅ **Classification [1] = 1** → Residue is **outside** the core/margin → REMOVE  
			✅ **Classification [1] ≠ 1** → Residue is **inside** the core/margin → KEEP  

	3️⃣ **Filter Hit Residues**
		- Retains **only residues classified as inside the query surface**.
		- Stores **filtered hit residues** in a new dictionary.

	4️⃣ **Return Trimmed Residue Dictionary**
		- The function returns a **cleaned dictionary of hit residues** that fit within 
		  the **query structure’s surface constraints**.

	================================================================================
	🔹 Summary
	================================================================================
	This function **removes extraneous hit structure residues** by applying **query surface constraints**.  
	It ensures **alignment to the query topology** by keeping **only residues within core/margin thresholds**.
	"""

	hit_vertices = []
	for r_key in hit_residues_in_surface:
		# ca_vertex = hit_residues_in_surface[r_key].atoms['CA'].v
		ca_vertex = hit_residues_in_surface[r_key].v
		hit_vertices.append(ca_vertex)

	from pHinderSurface import inLocalSurface
	classifications, hit_residues_in_surface = inLocalSurface(hit_vertices, query_surface, core_cutoff, margin_cutoff), {}
	for classification in classifications:
		if classification[1] != 1:
			v_4d =  classification[0]
			hit_residues_in_surface[v_4d.data.residue.key] = v_4d.data.residue.atoms['CA']

	return hit_residues_in_surface

def calculate_query_ssco_hull(pdb_file_path, pdb_file_name, pdb_format="pdb", zip_status=0, swell_hull=0.0):


	"""Arguments:
	----------
	pdb_file_path : str
		Path to the **query PDB file**, specifying the directory where the file is located.

	pdb_file_name : str
		Name of the **query PDB file**, including the filename without the path.

	pdb_format : str, optional (default="pdb")
		Specifies the **format of the input file**. Supported formats:
			✅ "pdb" → **Protein Data Bank (PDB) format**
			✅ "cif" → **Crystallographic Information File (CIF) format**

	zip_status : int, optional (default=0)
		Indicates whether the PDB file is **zipped**:
			✅ `0` → File is **not zipped**  
			✅ `1` → File is **compressed**  

	swell_hull : float, optional (default=0.0)
		Defines the **expansion factor (in Ångströms)** for swelling the convex hull.  
			📌 If `0.0`, no expansion is applied.  
			📌 If `> 0.0`, the convex hull is **expanded outward** by the given distance.  

	Returns:
	--------
	tuple : 
		A tuple containing:
			✅ `queryHull` → The **convex hull** of the query structure.  
			✅ `query_ca_atoms` → A list of **C-alpha atoms** from residues inside the convex hull.  
			✅ `query_residues_all` → A dictionary of **residues inside the convex hull**, sorted sequentially.  

	================================================================================
	📌 Function Overview
	================================================================================
	This function **computes the convex hull of a query protein structure** from a PDB file.  
	It identifies **residues inside the convex hull**, collects **C-alpha atoms**, and optionally 
	**expands the hull** for further structural analysis.

	================================================================================
	🔹 Main Steps
	================================================================================

	1️⃣ **Load the Query Structure**
		- Opens the **PDB/CIF file** and extracts **atomic information**.
		- Collects **C-alpha atom coordinates**.

	2️⃣ **Compute the Convex Hull**
		- Constructs a **3D convex hull** around the query structure.
		- Identifies residues **inside** and **on the surface** of the hull.

	3️⃣ **Sort Residues by Sequence Order**
		- Ensures residues inside the hull are **sorted sequentially**.

	4️⃣ **Extract C-alpha Atoms**
		- Collects **C-alpha atoms** from residues **inside the hull**.

	5️⃣ **(Optional) Swell the Convex Hull**
		- If `swell_hull > 0`, expands the hull **outward** by the given distance.
		- Saves the **swollen hull coordinates** to a text file for visualization.

	6️⃣ **Return the Hull and Residues**
		- Returns the **query convex hull**, **C-alpha atoms**, and **sorted residue dictionary**.

	================================================================================
	🔹 Summary
	================================================================================
	This function generates a **convex hull representation** of a protein structure 
	for use in **topological analysis, structural comparisons, and fold detection**.  
	It allows optional **hull expansion** for refined **surface-based** analysis.
	"""

	# 1) Calculate the convex hull of the query structure...
	##############################################################################################
	hullProtein = open_pdb(pdb_file_path, pdb_file_name, pdb_format=pdb_format, zip_status=zip_status)
	hullVertices = []
	for residue_key in hullProtein.residues:
		try:
			ca_vertex = hullProtein.residues[residue_key].atoms["CA"].v
			hullVertices.append(ca_vertex)
		except:
			pass
	queryHull = convexHull3D(hullVertices)
	
	# 2) Identify the subset of residues inside the convex hull...
	##############################################################################################
	v_in = queryHull.get_the_vertices_inside_of_the_hull()
	query_residues_all = {}
	for v in v_in:
		ca_atom = v.data
		residue = ca_atom.residue
		query_residues_all[residue.key] = residue
	v_of = queryHull.get_the_vertices_of_the_hull()
	for v in v_of:
		ca_atom = v.data
		residue = ca_atom.residue
		query_residues_all[residue.key] = residue	

	# 3) Sort the query residues inside the convex hull so they are in sequential order...
	##############################################################################################
	keys, tmp = sorted(query_residues_all), {}
	for key in keys:
		tmp[key] = query_residues_all[key]
	query_residues_all = tmp

	# 4) Collect the c-alpha atoms of the query residues inside the convex hull...
	##############################################################################################
	query_ca_atoms = get_ca_atoms(query_residues_all, origin="")

	# 5) If selected, using the facet normal vectors, swell the hull by the factor swell_hull...
	##############################################################################################
	if swell_hull:
		print("Swelling the convex hull by...", swell_hull, "Angstroms...")
		swollenHull = queryHull.modulateHull(multiplier=swell_hull) 
		# Write the swollen convex hull coordinates for visualization...
		surfaceFacetsString = ""
		for f in swollenHull.getFacets():
			t = f
			format = "%15.6f%15.6f%15.6f "
			surfaceFacetsString += format % (t.v1.x, t.v1.y, t.v1.z) + \
									format % (t.v2.x, t.v2.y, t.v2.z) + \
									format % (t.v3.x, t.v3.y, t.v3.z) + "\n"
		f = open("cHull-swollen-" + str(swell_hull) + ".txt", "w")
		f.write(surfaceFacetsString)
		f.close()
		queryHull = swollenHull

	return (queryHull, query_ca_atoms, query_residues_all)

def renumber(residue_dict):


	"""Arguments:
	----------
	residue_dict : dict
		A dictionary of residues, where:
			🔹 Keys are tuples representing residue identifiers: **(original_residue_number, chain_id, insertion_code)**
			🔹 Values are residue objects containing atomic information.

	Returns:
	--------
	dict : 
		A **renumbered dictionary of residues**, where:
			✅ Residue numbers are **sequentially assigned** starting from 1.
			✅ Residue keys are updated accordingly.
			✅ Atom sequence numbers and charges are updated to reflect the new numbering.

	================================================================================
	📌 Function Overview
	================================================================================
	This function **renumbers protein residues** by sorting them **sequentially** and 
	assigning **new residue numbers starting from 1**. It ensures **consistent numbering** 
	for downstream structural analysis.

	================================================================================
	🔹 Main Steps
	================================================================================

	1️⃣ **Sort Residues by Original Sequence Number**
		- Extracts residue keys and **sorts them** based on the **original residue number**.

	2️⃣ **Renumber Residues Sequentially**
		- Assigns new **residue numbers starting from 1**.
		- Creates a **new dictionary** with updated keys and values.

	3️⃣ **Update Atomic Information**
		- Iterates over each **renumbered residue**:
			✅ Updates **residue sequence numbers** in atomic data.  
			✅ Maintains **original charge information** from the old numbering.  

	4️⃣ **Return Renumbered Residue Dictionary**
		- The function returns a **fully renumbered dictionary**, ensuring that 
		  **residues and their atoms are consistently indexed**.

	================================================================================
	🔹 Summary
	================================================================================
	This function **standardizes protein residue numbering** for analysis by ensuring 
	**sequential ordering** of residues and updating atomic information accordingly.
	"""

	# 1) Sort the residues in residue_dict by residue number... 
	##############################################################################################
	residue_keys = sorted(residue_dict)

	# 2) Renumber the residues in residue_dict starting at 1...
	##############################################################################################
	i, res_num_1, sorted_residue_dict_renumbered = 1, residue_keys[0][0] - 1, {}
	for residue_key in residue_keys:
		new_residue_num = i
		new_residue_key = (i, residue_key[1], residue_key[2])
		sorted_residue_dict_renumbered[new_residue_key] = residue_dict[residue_key]
		for atom in sorted_residue_dict_renumbered[new_residue_key].atoms:
			sorted_residue_dict_renumbered[new_residue_key].atoms[atom].residue_sequence_number = new_residue_num
			sorted_residue_dict_renumbered[new_residue_key].atoms[atom].charge = residue_key[0]
		i += 1
	return sorted_residue_dict_renumbered

def get_4D_vertices(atom_list):

	"""Arguments:
	----------
	atom_list : list
		A list of PDBfile.atom objects, each representing an atomic coordinate 
		from a protein structure.

	Returns:
	--------
	list : 
		A list of 4D vertex objects (`Vertex4D`), where each vertex represents 
		the atom in a transformed **four-dimensional space (x, y, z, u)**.

	================================================================================
	📌 Function Overview
	================================================================================
	This function **converts a list of atomic coordinates** into a **4D vertex representation** 
	by adding an extra coordinate, `u = x² + y² + z²`. The function helps in **geometric 
	and topological analysis** of protein structures.

	================================================================================
	🔹 Main Steps
	================================================================================

	1️⃣ **Extract Atomic Coordinates**
		- Iterates through `atom_list` and extracts **x, y, z** coordinates from each atom.

	2️⃣ **Compute 4D Coordinate (`u`)**
		- Calculates `u` using the formula:  
		  \t📌 `u = x² + y² + z²`
		- This transformation maps **3D points into a 4D space**, aiding in geometric 
		  calculations such as **convex hull construction**.

	3️⃣ **Create `Vertex4D` Objects**
		- Each atom is **wrapped in a `Vertex4D` object**, storing:
			✅ (x, y, z, u) coordinates  
			✅ **Reference to original atomic data**  
			✅ **Atom serial ID**  

	4️⃣ **Return the List of 4D Vertices**
		- The function returns a list of `Vertex4D` objects, which can be used 
		  in **topological analysis**, such as **fold recognition** and **convex hull computations**.

	================================================================================
	🔹 Summary
	================================================================================
	This function **transforms protein atomic coordinates into 4D space** for efficient 
	geometric and structural analysis. The resulting vertices can be used in **fold 
	topology calculations**, **surface mapping**, and **convex hull generation**.
	"""

	# 1) Create a list of 4D vertices from a list of PDBfile.atom objects...
	##############################################################################################
	vertices_4D = []
	for atom in atom_list:
		x, y, z = atom.x, atom.y, atom.z
		u = atom.x**2 + atom.y**2 + atom.z**2
		vertex = Vertex4D((x,y,z,u),setC=1)
		vertex.data = atom
		vertex.id = atom.atom_serial
		vertices_4D.append(vertex)
	return vertices_4D

def calculate_difference_topology(
	query_residues_all, 
	hit_residues_trimmed_subset, 
	core_cutoff=0, 
	margin_cutoff=0, 
	query_surface=None):

	"""Arguments:
	----------
	query_residues_all : dict
		A dictionary of **query structure residues**, where:
			🔹 **Keys** → Residue identifiers.  
			🔹 **Values** → Residue objects containing atomic and sequence information.  

	hit_residues_trimmed_subset : dict
		A dictionary of **hit structure residues that have been trimmed**, where:
			🔹 **Keys** → Residue identifiers.  
			🔹 **Values** → Residue objects containing atomic and sequence information.  

	core_cutoff : int, optional (default=0)
		Defines the **core cutoff distance (in Ångströms)** when refining  
		the hit structure to remove **excess noise** in the topology.

	margin_cutoff : int, optional (default=0)
		Defines the **margin cutoff distance (in Ångströms)** when refining  
		the hit structure to remove **unnecessary residues** from the comparison.

	query_surface : object, optional (default=None)
		A **surface representation of the query structure**, used for  
		further refinement of the hit structure when it is excessively  
		populated with atoms.

	Returns:
	--------
	tuple : 
		A tuple containing:
			✅ `coverage_dict` → A dictionary mapping **query residues to their matched hit residues**,  
				indicating **which residues are covered** in the hit structure.  
			✅ `coverage_gaps` → A dictionary containing **gaps in query structure coverage**,  
				including their **size and position**.  
			✅ `coverage_gaps_cummulative` → A **cumulative measure** of all **coverage gaps**  
				relative to the query structure length.  
			✅ `coverage` → The **fraction of the query structure** covered by the hit structure.  
			✅ `n_terminal_percent_truncation` → The **percentage of truncation**  
				at the **N-terminal** of the query structure.  
			✅ `c_terminal_percent_truncation` → The **percentage of truncation**  
				at the **C-terminal** of the query structure.  
			✅ `hit_coverage_dict` → A dictionary mapping **hit residues that contribute  
				to query coverage**, indicating how well the hit structure aligns  
				with the query topology.  

	================================================================================
	📌 Function Overview
	================================================================================
	This function **compares the topology of a query structure** against a hit structure,  
	determining **structural coverage, truncation, and gaps**. It identifies **missing residues,  
	matches, and deviations** in the hit structure relative to the query.

	================================================================================
	🔹 Main Steps
	================================================================================

	1️⃣ **Initialize Coverage Dictionaries**
		- Creates dictionaries to track:
			✅ **Coverage of the query structure**.  
			✅ **First-come-first-serve (FCFS) topology matching**.  

	2️⃣ **Refine the Hit Structure (If Needed)**
		- If the hit structure is **overly dense**, applies **surface-based trimming**  
		  to remove **excess residues**.

	3️⃣ **Determine Residue Matches**
		- Loops through **query and hit residues**, assigning **best-matching residues**  
		  based on **C-alpha distance filtering**.

	4️⃣ **Filter Noise & Identify Gaps**
		- Removes **singleton (isolated) contacts** that do not form continuous  
		  topology matches.  
		- Identifies **coverage gaps** in the **query structure**.  

	5️⃣ **Classify Coverage Runs**
		- Segments **sequential contiguous topology coverage** into **runs**.  
		- Removes **low-diversity coverage islands**.  
		- Ensures that **runs are contiguous in the hit structure**.  

	6️⃣ **Fill Gaps in the Hit Structure**
		- Iteratively **fills small missing coverage gaps** to improve  
		  topology continuity.

	7️⃣ **Compute Topology Coverage Metrics**
		- Calculates:
			✅ **Overall fold coverage** of the query structure.  
			✅ **N-terminal and C-terminal truncation scores**.  
			✅ **Cumulative coverage gaps**.  

	================================================================================
	🔹 Summary
	================================================================================
	This function **measures the similarity of two protein structures**  
	by evaluating their **topological coverage, truncation, and gaps**.  
	It is designed for **protein fold recognition, topology validation,  
	and structural comparison analysis**.
	"""
	
	# The filter cut_off distance for a contact between a query and hit residue in the fold topology....
	cut_off = 6.0

	# Step 1: Initialize the dictionary for quantifying coverage of the query structure...
	coverage_dict = {}
	for renumbered_query_residue in sorted(query_residues_all):
		coverage_dict[renumbered_query_residue] = False

	# Step 2: Initialize the dictionary for quantifying coverage of the query structure...
	# fcfs = first come first serve
	fcfs_coverage_dict = {}
	for renumbered_query_residue in sorted(query_residues_all):
		fcfs_coverage_dict[renumbered_query_residue] = False

	# Step 3) Quantify coverage of the query structure by the hit structure...
	#########################################################################################################################################
	keys_q, keys_h = sorted(query_residues_all), sorted(hit_residues_trimmed_subset)

	# Step 3.a) If the hit structure is messy and floods the zone with atoms, de-noise by further trimming hit_residues_trimmed_subset...
	if len(keys_h)/len(keys_q) > 1 and query_surface:
		cut_off = 4.0
		hit_residues_in_surface = more_trim_using_surface(query_surface, core_cutoff, margin_cutoff, hit_residues_trimmed_subset)
		hit_residues_trimmed_subset = hit_residues_in_surface
		keys_h = sorted(hit_residues_trimmed_subset)

	# Step 3.b)
	# Loop over the c-alpha atoms of the query residues and apply the first-come-first-serve algorithm to calculate fold topology coverage... 
	# When searching the query topology from first to last residue,
	# against the hit topology from first to last residue, 
	# missing coverage can be attributed to one of these causes:
	#
	# 1 -- a local deviation in hit fold
	# 2 -- a large deviation in hit fold, or internal gap
	# 3 -- an N- or C-terminal truncation
	# 4 -- a large deviation in hit topology (i.e., different helical arrangements through a presumed bilayer)
	#
	# Predictions missing 1-4 are typically rank 1. Otherwise predications associated with 1-4 are typically rank 2 or 3. 	
	i, j, first_come_first_serve = 0, 0, {}
	while i < len(keys_q):

		# Loop over the c-alpha atoms of the hit residues...
		k, ca_q, min_d, min_q, min_h = 0, query_residues_all[keys_q[i]].atoms["CA"], 1000000, None, None
		while k < len(keys_h):
			
			# Get the hit c-alpha atom...
			# ca_h = hit_residues_trimmed_subset[keys_h[k]].atoms["CA"]
			ca_h = hit_residues_trimmed_subset[keys_h[k]]

			# Apply the topology match distance filter...
			# try:
			d = distance(ca_q, ca_h)
			if d < cut_off:
				if keys_q[i] not in first_come_first_serve:
					fcfs_coverage_dict[keys_q[i]] = keys_h[k]
			# except:
			# 	print(type(ca_q), type(ca_h))

			# Increment the inner loop...
			k += 1
		# Increment the outer loop...
		i += 1

	# Step 3.c) Initialize the dictionary for quantifying coverage of the query topology by the hit structure...
	hit_coverage_dict = {}
	for hit_residue in sorted(hit_residues_trimmed_subset):
		hit_coverage_dict[hit_residue] = False

	# Step 3.d) Update the hit coverage dictionary...
	coverage_dict = fcfs_coverage_dict
	keys = sorted(coverage_dict)
	for key in keys:
		query_r, hit_r = key, coverage_dict[key]
		if hit_r:
			hit_coverage_dict[hit_r] = True

	# Step 4) Iteratively remove singleton noise in coverage gaps...
	# Singleton noise refers to spurious contacts between the query and hit structures that do not produce runs of contiguous topology matching...
	coverage_revision = True
	while coverage_revision:
		coverage_revision = False
		i, coverage_dict_keys = 2, list(coverage_dict)
		while i < len(coverage_dict) - 2:

			# Evaluate the query residue +/- two adjacent residue registers...
			r0, r1, r2, r3, r4  = coverage_dict_keys[i-2], coverage_dict_keys[i-1], coverage_dict_keys[i], coverage_dict_keys[i+1], coverage_dict_keys[i+2]
			if coverage_dict[r2]:

				if not coverage_dict[r1] and not coverage_dict[r3]:
					coverage_dict[r2] = False
					coverage_revision = True
				if not coverage_dict[r0] and coverage_dict[r1] and not coverage_dict[r3]:
					coverage_dict[r2] = False
					coverage_dict[r1] = False
					coverage_revision = True
				if not coverage_dict[r1] and coverage_dict[r3] and not coverage_dict[r4]:
					coverage_dict[r2] = False
					coverage_dict[r3] = False
					coverage_revision = True
				if coverage_dict[r1] and coverage_dict[r3] and not coverage_dict[r0] and not coverage_dict[r4]:
					coverage_dict[r1] = False
					coverage_dict[r2] = False
					coverage_dict[r3] = False
					coverage_revision = True

			i += 1

	# Step 5: Calculate the sequential contigous topology coverage runs...
	runs, current_run = [], {}
	for query_r in coverage_dict:
		# Empty coverage site or end of run...
		if not coverage_dict[query_r]:
			if current_run:
				runs.append(current_run)
				current_run = {}
		else:
			current_run[query_r] = coverage_dict[query_r]
	# Get the last run...
	if current_run:
		runs.append(current_run)

	# Step 6: Calculate the residue diversity of each sequential contigous topology coverage run...
	run_diversities = []
	for run in runs:

		hit_diversity = {}
		for query_r in run:
			hit_r = run[query_r]
			hit_diversity[hit_r] = None

		# Remove the low diversity coverage island...
		if len(hit_diversity) <= 5:
			for query_r in run:
				coverage_dict[query_r] = False

		# Runs are calculated using the query topology...
		# Using the query run information, ensure an equivalent contiguous topology coverage run in the hit coverage dictionary...
		else:

			run_fill = {}
			for run_key in run:
				run_fill[run[run_key]] = None
			run_fill_keys = sorted(run_fill)

			hit_r_keys = sorted(hit_residues_trimmed_subset)
			first_run_index = hit_r_keys.index(run_fill_keys[0])
			last_run_index = hit_r_keys.index(run_fill_keys[-1])

			while first_run_index != last_run_index:
				hit_r_key = hit_r_keys[first_run_index]
				hit_coverage_dict[hit_r_key] = True
				first_run_index += 1
			hit_r_key = hit_r_keys[first_run_index]
			hit_coverage_dict[hit_r_key] = True

	# Step 7) The opposite of Step 4: Iteratively fill singleton noise coverage gaps...
	# Singleton noise refers to spurious MISSING contacts between the query and hit structures that break runs of contiguous topology matching...
	coverage_revision = True
	while coverage_revision:
		coverage_revision = False
		i, hit_coverage_dict_keys = 2, list(hit_coverage_dict)
		while i < len(hit_coverage_dict) - 2:

			# Evaluate the query residue +/- two adjacent residue registers...
			r0, r1, r2, r3, r4  = hit_coverage_dict_keys[i-2], hit_coverage_dict_keys[i-1], hit_coverage_dict_keys[i], hit_coverage_dict_keys[i+1], hit_coverage_dict_keys[i+2]
			if not hit_coverage_dict[r2]:

				if hit_coverage_dict[r1] and hit_coverage_dict[r3]:
					hit_coverage_dict[r2] = True
					coverage_revision = True
				if hit_coverage_dict[r0] and not hit_coverage_dict[r1] and hit_coverage_dict[r3]:
					hit_coverage_dict[r2] = True
					hit_coverage_dict[r1] = True
					coverage_revision = True
				if hit_coverage_dict[r1] and not hit_coverage_dict[r3] and hit_coverage_dict[r4]:
					hit_coverage_dict[r2] = True
					hit_coverage_dict[r3] = True
					coverage_revision = True
				if not hit_coverage_dict[r1] and not hit_coverage_dict[r3] and hit_coverage_dict[r0] and hit_coverage_dict[r4]:
					hit_coverage_dict[r1] = True
					hit_coverage_dict[r2] = True
					hit_coverage_dict[r3] = True
					coverage_revision = True

			i += 1

	# Step 8) Calculate topological coverage of the query structure by the hit structure....
	coverage_count, coverage_map_2d_from_n = 0, []
	for renumbered_query_residue in coverage_dict:
		if coverage_dict[renumbered_query_residue]:
			coverage_count += 1
		coverage_map_2d_from_n.append((renumbered_query_residue[0], coverage_count))
	coverage = coverage_count/len(coverage_dict)

	# Step 9) Calculate N-terminal truncation score...
	n_terminal_percent_truncation = 0
	for x in coverage_map_2d_from_n:
		if x[1] == 15: # XXXX 10/2 changed from 10
			n_terminal_percent_truncation = x[0]/len(coverage_dict)
			break

	# Step 9.a) Invert the 2D coverage map....
	coverage_count, coverage_map_2d_from_c = 0, []
	coverage_dict_keys = list(coverage_dict)
	coverage_dict_keys.reverse()
	for coverage_dict_key in coverage_dict_keys:
		if coverage_dict[coverage_dict_key]:
			coverage_count += 1
		coverage_map_2d_from_c.append((coverage_dict_key[0], coverage_count))

	# Step 9.b) Calculate C-terminal truncation score...
	c_terminal_percent_truncation = 0
	for x in coverage_map_2d_from_c:
		if x[1] == 15: # XXXX 10/2 changed from 10
			c_terminal_percent_truncation = (len(coverage_dict) - x[0])/len(coverage_dict)
			break

	# Step 10) Calculate the coverage gaps...
	coverage_gaps, coverage_gaps_cummulative, collection_gap = {}, 0, {}
	for renumbered_query_residue in coverage_dict:
		if not coverage_dict[renumbered_query_residue]:
			collection_gap[renumbered_query_residue] = None
		else:
			if collection_gap:
				collection_gap_keys = list(collection_gap)
				first_gap_residue = collection_gap_keys[0]
				last_gap_residue = collection_gap_keys[-1]
				coverage_gap = last_gap_residue[0] - first_gap_residue[0] + 1
				if coverage_gap > 0.1*len(coverage_dict):
					coverage_gaps_cummulative += coverage_gap
					coverage_gaps[(coverage_gap, first_gap_residue, last_gap_residue)] = collection_gap
				collection_gap = {}

	# Step 11) Get the last coverage gap...
	if collection_gap:
		collection_gap_keys = list(collection_gap)
		first_gap_residue = collection_gap_keys[0]
		last_gap_residue = collection_gap_keys[-1]
		coverage_gap = last_gap_residue[0] - first_gap_residue[0] + 1
		if coverage_gap > 0.1*len(coverage_dict):
			coverage_gaps_cummulative += coverage_gap
			coverage_gaps[(coverage_gap, first_gap_residue, last_gap_residue)] = collection_gap

	# Step 12) Calculate cummulative coverage gaps...
	coverage_gaps_cummulative = coverage_gaps_cummulative/len(coverage_dict)

	return (coverage_dict, coverage_gaps, coverage_gaps_cummulative, coverage, n_terminal_percent_truncation, c_terminal_percent_truncation, hit_coverage_dict)

##################################################################################################
# Score the TM-score hits as full or partial structural matches...
##################################################################################################

f = open(pdb_input_file_path, "r")
pdb_input_file_paths = []
for line in f:
	pdb_input_file_paths.append(line.strip())
f.close()

"""This script scores structural predictions by comparing a query protein structure against predicted hit structures. 
It evaluates structural coverage, topology differences, and classification ranks based on topology alignment and truncation.

================================================================================
📌 Main Steps
================================================================================

1️⃣ **Initialize Data Structures**
	- `query_hulls` and `query_surfaces` store convex hulls and surface representations of query structures.
	- Iterates through multiple **PDB input files**.

2️⃣ **Compute Query Structure Representation**
	- If **convex hull-based scoring** is used:
		- Computes or retrieves the **3D convex hull** of the query protein (`query_hulls`).
		- Renumbers residues for consistency.
	- If **surface-based scoring** is used:
		- Computes or retrieves **query surface representation** (`query_surfaces`).
		- Renumbers residues accordingly.

3️⃣ **Trim Hit Structures**
	- Depending on whether **convex hull** or **surface** filtering is used:
		- Filters hit residues that **fall within the convex hull** of the query.
		- Or trims hit residues using **surface-based filtering**.

4️⃣ **Extract C-alpha Atom (CA) Structural Representation**
	- Extracts **C-alpha atoms** (`CA`) from the trimmed hit structure.
	- Computes **SSCO (secondary structure contacts)** from `CA` atoms.

5️⃣ **Compute Topology Score & Structural Coverage**
	- Calls `calculate_difference_topology()`, which determines:
		- **Residue coverage** of the hit structure on the query.
		- **Coverage gaps**, cumulative coverage, and truncation at **N-terminal** and **C-terminal** regions.

6️⃣ **Parse Structural Metrics**
	- Extracts **TM-score** and **UniProt ID** from the file name.
	- TM-score is a common **structural similarity metric** for protein folds.

7️⃣ **Generate PDB Outputs**
	- **Generates PDB files** for both **trimmed** and **full** versions of the hit structure.
	- The **hit residues aligned to the query** are stored in different files.

8️⃣ **Classify Structural Matches**
	- Assigns **ranks (r1, r2, r3)** based on coverage quality:
		- **r1 (High-quality match)** → Minimal truncation, high coverage, small gaps.
		- **r2 (Intermediate match)** → Moderate truncation or gaps, but still structurally related.
		- **r3 (Low-quality match)** → Large truncations, severe coverage gaps.
	- Classification is based on:
		- **Coverage percentage**
		- **N-/C-terminal truncation**
		- **Presence of coverage gaps**
		- **Cumulative coverage gap threshold (10-15%)**

9️⃣ **Save Ranked PDB Files**
	- Saves **aligned** and **trimmed** PDB structures into different files:
		- **Rank 1 (best matches)**
		- **Rank 2 (moderate matches)**
		- **Rank 3 (weak matches)**
	- Files are named based on **rank, TM-score, coverage, and query/hit names**.

================================================================================
🔹 Summary
================================================================================
This script compares predicted protein structures against a known query structure 
using topological alignment and filtering methods (convex hull/surface). 
It assigns scores and ranks based on similarity, truncation, and structural integrity, 
then saves PDB files of aligned structures.
"""

############################################################################################
# Score each structure prediction...
############################################################################################
from ssco import *
query_hulls, query_surfaces = {}, {}
for pdb_input_file_path in pdb_input_file_paths:

	if not pdb_input_file_path:
		continue

	pdb_input_file_name = pdb_input_file_path.split(sep)[-1]
	pdb_input_file_path = sep.join(pdb_input_file_path.split(sep)[:-1]) + sep

	##########################################################################################
	# If it has not been calculated, calculate the 3D convex hull of the query structure...
	##########################################################################################

	query_ssco_pdb_file_path = where_query_pdb_files_are_located
	query_ssco_pdb_file_name = pdb_input_file_name.split("--")[1] + "-ssco.pdb"
	from os.path import exists
	if not exists(query_ssco_pdb_file_path + query_ssco_pdb_file_name):
		grab_ssco(query_ssco_pdb_file_path, 
				query_ssco_pdb_file_path, 
				pdb_input_file_name.split("--")[1] + ".pdb", 
				write_helix=True, write_sheet=True, write_combined=True,
				return_mode="write")
	# calculate_query_ssco_hull(query_ssco_pdb_file_path, query_ssco_pdb_file_name, swell_hull=swell_hull)
		
	# If hull...
	if not by_surface:
		if query_ssco_pdb_file_path + query_ssco_pdb_file_name not in query_hulls:

			query_hulls[query_ssco_pdb_file_path + query_ssco_pdb_file_name] = calculate_query_ssco_hull(query_ssco_pdb_file_path, query_ssco_pdb_file_name, swell_hull=swell_hull)
			query_residues_all = query_hulls[query_ssco_pdb_file_path + query_ssco_pdb_file_name][2]
			change_to_renumbered = list(query_hulls[query_ssco_pdb_file_path + query_ssco_pdb_file_name])
			change_to_renumbered[2] = renumber(query_residues_all)
			query_hulls[query_ssco_pdb_file_path + query_ssco_pdb_file_name] = tuple(change_to_renumbered)

		queryHull = query_hulls[query_ssco_pdb_file_path + query_ssco_pdb_file_name][0]
		query_residues_all = query_hulls[query_ssco_pdb_file_path + query_ssco_pdb_file_name][2]

	# If surface...
	else:
		if query_ssco_pdb_file_path + query_ssco_pdb_file_name not in query_surfaces:

			query_surfaces[query_ssco_pdb_file_path + query_ssco_pdb_file_name] = calculate_query_surface(query_ssco_pdb_file_path, query_ssco_pdb_file_name, pdb_format="pdb", zip_status=0)
			query_residues_all = query_surfaces[query_ssco_pdb_file_path + query_ssco_pdb_file_name][2]
			change_to_renumbered = list(query_surfaces[query_ssco_pdb_file_path + query_ssco_pdb_file_name])
			change_to_renumbered[2] = renumber(query_residues_all)
			query_surfaces[query_ssco_pdb_file_path + query_ssco_pdb_file_name] = tuple(change_to_renumbered)

		query_surface = query_surfaces[query_ssco_pdb_file_path + query_ssco_pdb_file_name][0][1]
		query_residues_all = query_surfaces[query_ssco_pdb_file_path + query_ssco_pdb_file_name][2]

	##########################################################################################
	# Open and trim the hit structure using the convex hull of the query structure...
	##########################################################################################

	if not by_surface:
		# If hull...	
		hit_result = trim(queryHull, pdb_input_file_path, pdb_input_file_name, pdb_format="pdb", zip_status=0)
		hit_residues_in_hull = hit_result[0]
		hit_residues_all = hit_result[1]

	else:
		# If surface...
		hit_result = trim_using_surface(query_surface, 8, 8, pdb_input_file_path, pdb_input_file_name, pdb_format="pdb", zip_status=0)
		hit_residues_in_surface = hit_result[0]
		hit_residues_all = hit_result[1]

	##########################################################################################
	# Grab the SSCO of the trimmed hit structure...
	##########################################################################################
	
	if not by_surface:
		# If hull...
		hit_cas = []
		for residue_key in hit_residues_in_hull:
			hit_cas.append(hit_residues_in_hull[residue_key].atoms["CA"])
		# ca_ssco = grab_ssco(hit_cas)
		# hit_residues_in_hull = ca_ssco
		ca_ssco = grab_ssco(ca_atoms=hit_cas, return_mode="atoms")
		hit_residues_in_hull = ca_ssco

	else:
		# If surface...
		hit_cas = []
		for residue_key in hit_residues_in_surface:
			hit_cas.append(hit_residues_in_surface[residue_key].atoms["CA"])
		# ca_ssco = grab_ssco(hit_cas)
		# hit_residues_in_surface = ca_ssco
		ca_ssco = grab_ssco(ca_atoms=hit_cas, return_mode="atoms")
		hit_residues_in_surface = ca_ssco

	##########################################################################################
	# Calculate topology score and fractional coverage...
	##########################################################################################
	
	# Scores...
	if not by_surface:
		# If hull...
		topology_score = calculate_difference_topology(query_residues_all, hit_residues_in_hull)
	else:
		# If surface...
		topology_score = calculate_difference_topology(query_residues_all, hit_residues_in_surface, core_cutoff=4, margin_cutoff=4, query_surface=query_surface)

	coverage_dict = topology_score[0]
	coverage_gaps = topology_score[1]
	coverage_gaps_cummulative = topology_score[2]
	coverage = topology_score[3]
	n_terminal_percent_truncation = topology_score[4]
	c_terminal_percent_truncation = topology_score[5]
	hit_coverage_dict = topology_score[6]

	# On cluster...
	# Parse the tm_score, rmsd, and uniprot ID from the hit file name...
	###########################################################
	tm_score = float(pdb_input_file_name.split("--")[-1].split(".pdb")[0])
	hit_file_name = pdb_input_file_name.split("--")[0]
	query_file_name = pdb_input_file_name.split("--")[1]

	##########################################################################################
	# Create PDB strings for the trimmed and full hit structures...
	##########################################################################################

	# Calculated query topology coverage and display using the pHinder convex hull...
	if not by_surface:

		# If hull...
		trimmed_pdb = ""
		for r in hit_residues_in_hull:
			trimmed_pdb += str(hit_residues_in_hull[r])
		full_pdb = ""
		for r in hit_residues_all:
			full_pdb += str(hit_residues_all[r])

	# Calculated query topology coverage and display using the pHinder surface...
	else:

		trimmed_pdb = ""
		for r in hit_coverage_dict:
			if hit_coverage_dict[r]: # hit residue is coincident on query topology...
				trimmed_pdb += str(hit_residues_in_surface[r].residue)
		full_pdb = ""
		for r in hit_residues_all:
			full_pdb += str(hit_residues_all[r])

	##########################################################################################
	# Classify score: the ultimate heuristic...
	# The takehome: tm_score in insufficient for robustly characterizing fold matches...
	##########################################################################################

	print("Calculating score for:", hit_file_name)
	classification, rank_1_match, rank_2_match, rank_3_match, rank_flag = "", {}, {}, {}, "--OK"
	pre_key = [coverage, tm_score, hit_file_name, query_file_name]

	if n_terminal_percent_truncation >= 0.15:

		if coverage >= 0.5:

			classification = "r2"
			key = tuple(pre_key + [classification,])
			rank_2_match[key] = (trimmed_pdb, full_pdb)
			rank_flag = "--NTC"

		else:

			classification = "r3"
			key = tuple(pre_key + [classification,])
			rank_3_match[key] = (trimmed_pdb, full_pdb)
			rank_flag = "--LC"

	elif c_terminal_percent_truncation >= 0.15:

		if coverage >= 0.5:

			classification = "r2"
			key = tuple(pre_key + [classification,])
			rank_2_match[key] = (trimmed_pdb, full_pdb)
			rank_flag = "--CTC"

		else:

			classification = "r3"
			key = tuple(pre_key + [classification,])
			rank_3_match[key] = (trimmed_pdb, full_pdb)
			rank_flag = "--LC"

	elif coverage_gaps:

		# A cummulative coverage gap > 10% is unacceptable for rank 1
		if coverage_gaps_cummulative > 0.15:

			if coverage_gaps_cummulative > 0.5:

				classification = "r3"
				key = tuple(pre_key + [classification,])
				rank_3_match[key] = (trimmed_pdb, full_pdb)
				rank_flag = "--LC"

			else:

				classification = "r2"
				key = tuple(pre_key + [classification,])
				rank_2_match[key] = (trimmed_pdb, full_pdb)
				rank_flag = "--ICG"

		# New...was just an else for r1.
		# elif tm_score < 0.6 and coverage >= 0.8: # Least stringent path to r1
		elif coverage > 0.7:

			classification = "r1"
			key = tuple(pre_key + [classification,])
			rank_1_match[key] = (trimmed_pdb, full_pdb)

		else:

			classification = "r2"
			key = tuple(pre_key + [classification,])
			rank_2_match[key] = (trimmed_pdb, full_pdb)
			rank_flag = "--ICG"	

	elif coverage < 0.7:

		if coverage >= 0.5:

			classification = "r2"
			key = tuple(pre_key + [classification,])
			rank_2_match[key] = (trimmed_pdb, full_pdb)
			rank_flag = "--LC"

		else:

			classification = "r3"
			key = tuple(pre_key + [classification,])
			rank_3_match[key] = (trimmed_pdb, full_pdb)
			rank_flag = "--LC"

	else: # Most stringent path to r1

		classification = "r1"
		key = tuple(pre_key + [classification,])
		rank_1_match[key] = (trimmed_pdb, full_pdb)		

	##################################################################################################
	# Copy aligned PDB files...
	##################################################################################################

	if rank_1_match:

		scores = sorted(rank_1_match)
		scores.reverse()
		for score in scores:

			# Full...
			file_name = "1--" + str(round(score[0],2)) + "-" + str(round(score[1],2)) + "--" + score[2] + "--" + score[3] + rank_flag + "-full.pdb" 
			f = open(save_directory + file_name, "w")
			f.write(rank_1_match[score][1])
			f.close()

			# Trimmed...
			file_name = "1--" + str(round(score[0],2)) + "-" + str(round(score[1],2)) + "--" + score[2] + "--" + score[3] + rank_flag + "-trim.pdb" 
			f = open(save_directory + file_name, "w")
			f.write(rank_1_match[score][0])
			f.close()

	if rank_2_match:

		scores = sorted(rank_2_match)
		scores.reverse()
		for score in scores:

			# Full...
			file_name = "2--" + str(round(score[0],2)) + "-" + str(round(score[1],2)) + "--" + score[2] + "--" + score[3] + rank_flag + "-full.pdb" 
			f = open(save_directory + file_name, "w")
			f.write(rank_2_match[score][1])
			f.close()

			# Trimmed...
			file_name = "2--" + str(round(score[0],2)) + "-" + str(round(score[1],2)) + "--" + score[2] + "--" + score[3] + rank_flag + "-trim.pdb" 
			f = open(save_directory + file_name, "w")
			f.write(rank_2_match[score][0])
			f.close()

	if rank_3_match:

		scores = sorted(rank_3_match)
		scores.reverse()
		for score in scores:

			# Full...
			file_name = "3--" + str(round(score[0],2)) + "-" + str(round(score[1],2)) + "--" + score[2] + "--" + score[3] + rank_flag + "-full.pdb" 
			f = open(save_directory + file_name, "w")
			f.write(rank_3_match[score][1])
			f.close()

			# Trimmed...
			file_name = "3--" + str(round(score[0],2)) + "-" + str(round(score[1],2)) + "--" + score[2] + "--" + score[3] + rank_flag + "-trim.pdb" 
			f = open(save_directory + file_name, "w")
			f.write(rank_3_match[score][0])
			f.close()