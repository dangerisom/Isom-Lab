# Import dependencies here for enhanced performance.
####################################################
from os import sep
import gzip
from compGeometry import Vertex
from decimal import *

singleLetter = {"ALA":"A", "ASP":"D", "ASN":"N", "ARG":"R", "CYS":"C", "GLY":"G", "GLU":"E", "GLN":"Q", "HIS":"H", "ILE":"I", "LEU":"L", "LYS":"K", "MET":"M", "PRO":"P", "PHE":"F", "SER":"S", "THR":"T", "TYR":"Y", "TRP":"W", "VAL":"V"}

NONPOLAR  = ["ALA","ILE","LEU","MET","PHE","PRO","VAL","GLY"]
POLAR     = ["ASN","GLN","SER","THR","TYR","TRP"]
IONIZABLE = ["ASP", "GLU", "HIS", "CYS", "LYS", "ARG", "MLY", "MLZ", "M3L"]
ALL_SIDECHAINS = NONPOLAR + POLAR + IONIZABLE
ACTIVE    = ["GLY","ARG","ASP","GLU","HIS","LYS","CYS","SER","TYR","THR","PHE","TRP"]

def assess_cif_fields(cif_options):

	cif_options_all = {
						"_atom_site.group_PDB": False,
						"_atom_site.id":False,
						"_atom_site.type_symbol":False,
						"_atom_site.label_atom_id":False, 
						"_atom_site.label_alt_id":False, 
						"_atom_site.label_comp_id":False, 
						"_atom_site.label_asym_id":False, 
						"_atom_site.label_entity_id":False,
						"_atom_site.label_seq_id":False, 
						"_atom_site.pdbx_PDB_ins_code":False, 
						"_atom_site.Cartn_x":False, 
						"_atom_site.Cartn_y":False, 
						"_atom_site.Cartn_z":False, 
						"_atom_site.occupancy":False, 
						"_atom_site.B_iso_or_equiv":False, 
						"_atom_site.Cartn_x_esd":False, 
						"_atom_site.Cartn_y_esd":False, 
						"_atom_site.Cartn_z_esd":False, 
						"_atom_site.occupancy_esd":False, 
						"_atom_site.B_iso_or_equiv_esd":False, 
						"_atom_site.pdbx_formal_charge":False, 
						"_atom_site.auth_seq_id":False, 
						"_atom_site.auth_comp_id":False, 
						"_atom_site.auth_asym_id":False, 
						"_atom_site.auth_atom_id":False, 
						"_atom_site.pdbx_PDB_model_num":False}

	for cif_option in cif_options.split("\n"):
		if cif_option in cif_options_all:
			cif_options_all[cif_option] = True

	return cif_options_all

def parse_cif_format(pdb_line):
	pass

class PDBfile:

	def __init__(self, pdbFilePath="", pdbFileName="", pdbFileAsString="", zip_status=0):

		self.zip = zip_status
		self.pdbFilePath = pdbFilePath
		self.pdbFileName = pdbFileName
		self.pdbCode = "NONE"
		if self.pdbFileName:
			self.pdbCode = pdbFileName.split(".cif")[0]
		self.chain_selected_in_file_name, self.chains = 0, []
		if "-chains." in self.pdbCode:
			self.chain_selected_in_file_name = 1
			chain_string = self.pdbCode.split("-chains.")[1]
			for chain in chain_string.split("."):
				if chain:
					self.chains.append(chain)

		self.pdbFileLines = []
		self.atoms = {}
		self.hetatoms = {}
		self.residues, self.res_bychain = {}, {}
		self.het_residues, self.het_bychain = {}, {}

		self.cif_options_all = {
								"_atom_site.group_PDB": False,
								"_atom_site.id":False,
								"_atom_site.type_symbol":False,
								"_atom_site.label_atom_id":False, 
								"_atom_site.label_alt_id":False, 
								"_atom_site.label_comp_id":False, 
								"_atom_site.label_asym_id":False, 
								"_atom_site.label_entity_id":False,
								"_atom_site.label_seq_id":False, 
								"_atom_site.pdbx_PDB_ins_code":False, 
								"_atom_site.Cartn_x":False, 
								"_atom_site.Cartn_y":False, 
								"_atom_site.Cartn_z":False, 
								"_atom_site.occupancy":False, 
								"_atom_site.B_iso_or_equiv":False, 
								"_atom_site.Cartn_x_esd":False, 
								"_atom_site.Cartn_y_esd":False, 
								"_atom_site.Cartn_z_esd":False, 
								"_atom_site.occupancy_esd":False, 
								"_atom_site.B_iso_or_equiv_esd":False, 
								"_atom_site.pdbx_formal_charge":False, 
								"_atom_site.auth_seq_id":False, 
								"_atom_site.auth_comp_id":False, 
								"_atom_site.auth_asym_id":False, 
								"_atom_site.auth_atom_id":False, 
								"_atom_site.pdbx_PDB_model_num":False}

		# Collect the line in the PDB file
		if self.zip and self.pdbFilePath and self.pdbFileName:
			f = gzip.GzipFile(self.pdbFilePath + self.pdbFileName,"rb")
			for line in f:
				line = line.decode()
				self.pdbFileLines.append(line)
				possible_cif_option = line.strip()
				if possible_cif_option in self.cif_options_all:
					self.cif_options_all[possible_cif_option] = True
		elif pdbFileAsString:
			for line in pdbFileAsString.split("\n"):
				possible_cif_option = line
				if possible_cif_option in self.cif_options_all:
					self.cif_options_all[possible_cif_option] = True
				self.pdbFileLines.append(line + "\n")
		elif self.pdbFilePath and self.pdbFileName:
			f = open(self.pdbFilePath + self.pdbFileName, "r")
			for line in f:
				self.pdbFileLines.append(line)
				possible_cif_option = line.strip()
				if possible_cif_option in self.cif_options_all:
					self.cif_options_all[possible_cif_option] = True

		# Loop over the PDB file lines to parse ATOM and HETATM records
		for line in self.pdbFileLines:

			# ATOM record
			if line[0:4] == "ATOM":

				# Get and save the atom information
				atom = Atom(line, self.cif_options_all)
				self.atoms[atom.atom_serial] = atom

				# Update residues with atom
				if atom.residue_key in self.residues:
					self.residues[atom.residue_key].addatom(atom)
					#self.res_bychain[atom.chain_identifier][atom.residue_key] = self.residues[atom.residue_key]
				else:
					residue = Residue(atom)
					self.residues[residue.key] = residue
					if atom.chain_identifier not in self.res_bychain:
						self.res_bychain[atom.chain_identifier] = {}
						self.res_bychain[atom.chain_identifier][atom.residue_key] = residue
					else:
						self.res_bychain[atom.chain_identifier][atom.residue_key] = residue

				# Update chains identifiers
				if not self.chain_selected_in_file_name and atom.chain_identifier not in self.chains:
					self.chains.append(atom.chain_identifier)

			# HETATM record
			if line[0:6] == "HETATM":

				# Get and save the atom information
				atom = Atom(line, self.cif_options_all)
				self.atoms[atom.atom_serial] = atom

				# Update residues with atom
				if atom.residue_key in self.residues:
					self.het_residues[atom.residue_key].addatom(atom)
					#self.het_bychain[atom.chain_identifier][atom.residue_key] = self.het_residues[atom.residue_key]
				else:
					residue = Residue(atom)
					self.het_residues[residue.key] = residue
					if atom.chain_identifier not in self.het_bychain:
						self.het_bychain[atom.chain_identifier] = {}
						self.het_bychain[atom.chain_identifier][atom.residue_key] = residue
					else:
						self.het_bychain[atom.chain_identifier][atom.residue_key] = residue

				# Update chains identifiers
				if not self.chain_selected_in_file_name and atom.chain_identifier not in self.chains:
					self.chains.append(atom.chain_identifier)

		# FOR EACH SIDECHAIN, LINK SUB (ATOMS) TO SUPER (RESIDUE).
		# SET THE SHORT-HAND REPRESENTATIONS OF EACH SIDECHAIN.
		##########################################################
		if self.residues:
			for residue in self.residues.values():
				residue.setAltConfBackboneAtoms()
				residue.setsub()
				residue.set_sidechain_representations()

		# CREATE CHAIN INSTANCES FOR THE RESIDUES IN THE PDB FILE.
		#########################################################
		self.res_chains = {}
		for chain in self.res_bychain.keys():
			self.res_chains.update({chain:Chain(self.pdbFileName,
												self.pdbCode,
												chain,
												self.res_bychain[chain],
												peptide=1,
												het_residues=\
												self.het_residues)})

		# CREATE CHAIN INSTANCES FOR THE HETERO ATOMS IN THE PDB FILE.
		##############################################################
		self.het_chains = {}
		for chain in self.het_bychain:
			self.het_chains.update({chain:Chain(self.pdbFileName,
												self.pdbCode,
												chain,
												self.het_bychain[chain],
												peptide=0)})

		# # SET THE NUMBER OF PROTEIN CHAINS IN THE PDB FILE.
		# ###################################################
		# self.numchains = len(self.res_bychain)

		# #SET THE NUMBER OF LIGANDS IN THE PDB FILE.
		# ###########################################
		# for chain in self.het_bychain:
		# 	self.numligands += len(self.het_bychain[chain])

	def res_atom_hash(self):
		string = ""
		keys = self.atoms.keys()
		keys.sort()
		for key in keys:
			atom = self.atoms[key]
			string += "%-6s|%s" % (str(atom.atom_serial),str(atom))
		return string
			
	def het_atom_hash(self):
		string = ""
		keys = self.hetatoms.keys()
		keys.sort()
		for key in keys:
			atom = self.hetatoms[key]
			string += "%-6s|%s" % (str(atom.atom_serial),str(atom))
		return string

	def convert_to_pdb_format(self):

		pdb_formatted_string = ""
		for atom in self.atoms:
			pdb_formatted_string += self.atoms[atom].get_pdb_format()
		for atom in self.hetatoms:
			pdb_formatted_string += self.atoms[atom].get_pdb_format()
		return pdb_formatted_string

	# def general_information(self):
	# 	string = ""
	# 	string += "%-6s%s\n" % ("NCHN",str(self.numchains))
	# 	string += "%-6s%s\n" % ("NLIG",str(self.numligands))
	# 	string += "%-6s%s\n" % ("LLIG",str(self.num_atoms_of_largest_ligand()))
	# 	string += self.compound
	# 	return string

	# def num_atoms_of_largest_ligand(self):
	# 	l = 0
	# 	for chain in self.het_bychain.values():
	# 		for res in chain.values():
	# 			if len(res.atoms) > l:
	# 				l = len(res.atoms)
	# 	return l
			
class Chain:

	def __init__(self, pdbfilename, pdbcode, chain_identifier, resdict,peptide=1, het_residues=None):

		# GENERAL INFORMATION FOR THE CHAIN.
		####################################
		self.pdbFileName      = pdbfilename
		self.pdbCode          = pdbcode
		self.chain_identifier = chain_identifier
		self.residues         = resdict
		self.length           = len(self.residues)
		self.gap              = None
		
		# INDICATES WHETHER THE CHAIN IS A PEPTIDE CHAIN OR HETERO CHAIN.
		#################################################################
		self.peptide = peptide
		
		# ALL OF THE HETERO RESIDUES, I.E., LIGANDS, THAT ARE  ASSOCIATED
		# WITH THE PROTEIN, NOT JUST THIS CHAIN.
		##################################################################
		self.het_residues = het_residues

		# CHECK IF THERE ARE GAPS IN THE PROTEIN SEQUENCE.
		##################################################
		if self.peptide:
			self.gap = self.maxgap()

		# SIDECHAINS GROUPED BY TYPE AND COMBINED TYPES.
		################################################
		self.nonpolar_sc           = {}
		self.nonpolar_polar_sc     = {}
		self.nonpolar_ionizable_sc = {}
		self.polar_sc              = {}
		self.polar_ionizable_sc    = {}
		self.ionizable_sc          = {}
		self.active_sc             = {}

		# SORT SIDECHAINS BY TYPE.
		##########################
		for r in self.residues:
			if self.residues[r].name in NONPOLAR:
				self.nonpolar_sc.update({r:self.residues[r]})
			if self.residues[r].name in NONPOLAR+POLAR:
				self.nonpolar_polar_sc.update({r:self.residues[r]})
			if self.residues[r].name in NONPOLAR+IONIZABLE:
				self.nonpolar_ionizable_sc.update({r:self.residues[r]})
			if self.residues[r].name in POLAR:
				self.polar_sc.update({r:self.residues[r]})
			if self.residues[r].name in POLAR+IONIZABLE:
				self.polar_ionizable_sc.update({r:self.residues[r]})
			if self.residues[r].name in IONIZABLE:
				self.ionizable_sc.update({r:self.residues[r]})
			if self.residues[r].name in ACTIVE:
				self.active_sc.update({r:self.residues[r]})
		
	def maxgap(self):
		keys = sorted(self.residues)
		gap = 1
		resnum = self.residues[keys[0]].num
		for key in keys[1:]:
			current_gap = key[0]-resnum
			if current_gap > gap:
				gap = current_gap
			resnum = key[0]
		return gap

	def get_c_alpha_atoms(self,type="all",return_Vertex=0,return_Atom=1):
		atoms = []
		if type == "all":
			for res in self.residues.values():
				try: atoms.append(res.atoms['CA'])
				except: continue
		elif type == "nonpolar":
			for res in self.nonpolar_sc.values():
				try: atoms.append(res.atoms['CA'])
				except: continue
		elif type == "nonpolar-polar":
			for res in self.nonpolar_polar_sc.values():
				try: atoms.append(res.atoms['CA'])
				except: continue
		elif type == "nonpolar-ionizable":
			for res in self.nonpolar_ionizable_sc.values():
				try: atoms.append(res.atoms['CA'])
				except: continue
		elif type == "polar":
			for res in self.polar_sc.values():
				try: atoms.append(res.atoms['CA'])
				except: continue
		elif type == "polar-ionizable":
			for res in self.polar_ionizable_sc.values():
				try: atoms.append(res.atoms['CA'])
				except: continue
		elif type == "ionizable":
			for res in self.ionizable_sc.values():
				try: atoms.append(res.atoms['CA'])
				except: continue
		elif type == "active":
			for res in self.active_sc.values():
				try: atoms.append(res.atoms['CA'])
				except: continue

		# Return a list of Atom or Vertex objects.
				##########################################
		if return_Vertex:
			vertices = []
			for atom in atoms:
				vertex = atom.v
				vertex.data = atom.Atom_LO
				vertices.append(vertex)
			return vertices
		else:
			return atoms

	def get_reduced_sidechain_atoms(self,type="all"):
		atoms = []
		if type == "all":
			for res in self.residues.values():
				if res.scom: atoms.append(res.scom)
		elif type == "nonpolar":
			for res in self.nonpolar_sc.values():
				if res.scom: atoms.append(res.scom)
		elif type == "nonpolar-polar":
			for res in self.nonpolar_polar_sc.values():
				if res.scom: atoms.append(res.scom)
		elif type == "nonpolar-ionizable":
			for res in self.nonpolar_ionizable_sc.values():
				if res.scom: atoms.append(res.scom)
		elif type == "polar":
			for res in self.polar_sc.values():
				if res.scom: atoms.append(res.scom)
		elif type == "polar-ionizable":
			for res in self.polar_ionizable_sc.values():
				if res.scom: atoms.append(res.scom)
		elif type == "ionizable":
			for res in self.ionizable_sc.values():
				if res.scom: atoms.append(res.scom)
		elif type == "active":
			for res in self.active_sc.values():
				if res.scom: atoms.append(res.scom)
		return atoms

	def get_terminal_sidechain_atoms(self, type="all", return_Vertex=0, return_Atom=1):
		
		atoms = []
		if type == "all":
			for res in self.residues.values():
				if res.ter: atoms.append(res.ter)
		elif type == "nonpolar":
			for res in self.nonpolar_sc.values():
				if res.ter: atoms.append(res.ter)
		elif type == "nonpolar-polar":
			for res in self.nonpolar_polar_sc.values():
				if res.ter: atoms.append(res.ter)
		elif type == "nonpolar-ionizable":
			for res in self.nonpolar_ionizable_sc.values():
				if res.ter: atoms.append(res.ter)
		elif type == "polar":
			for res in self.polar_sc.values():
				if res.ter: atoms.append(res.ter)
		elif type == "polar-ionizable":
			for res in self.polar_ionizable_sc.values():
				if res.ter: atoms.append(res.ter)
		elif type == "ionizable":
			for res in self.ionizable_sc.values():
				if res.ter: atoms.append(res.ter)
		elif type == "active":
			for res in self.active_sc.values():
				if res.ter: atoms.append(res.ter)

		# Return a list of Atom or Vertex objects.
		##########################################
		if return_Vertex:
			vertices = []
			for atom in atoms:
				vertex = atom.v
				vertex.data = atom.Atom_LO
				vertices.append(vertex)
			return vertices
		else:
			return atoms

	def get_backbone_atoms(self,type="all"):
		atoms = []
		if type == "all":
			for res in self.residues.values():
				atoms += res.get_backbone_atoms()
		elif type == "nonpolar":
			for res in self.nonpolar_sc.values():
				atoms.append(res.get_backbone_atoms())
		elif type == "nonpolar-polar":
			for res in self.nonpolar_polar_sc.values():
				atoms.append(res.get_backbone_atoms())
		elif type == "nonpolar-ionizable":
			for res in self.nonpolar_ionizable_sc.values():
				atoms.append(res.get_backbone_atoms())
		elif type == "polar":
			for res in self.polar_sc.values():
				atoms.append(res.get_backbone_atoms())
		elif type == "polar-ionizable":
			for res in self.polar_ionizable_sc.values():
				atoms.append(res.get_backbone_atoms())
		elif type == "ionizable":
			for res in self.ionizable_sc.values():
				atoms.append(res.get_backbone_atoms())
		elif type == "active":
			for res in self.active_sc.values():
				atoms.append(res.get_backbone_atoms())
		return atoms

	def get_sidechain_atoms(self,type="all"):
		atoms = []
		if type == "all":
			for res in self.residues.values():
				atoms += res.get_sidechain_atoms()
		elif type == "nonpolar":
			for res in self.nonpolar_sc.values():
				atoms.append(res.get_sidechain_atoms())
		elif type == "nonpolar-polar":
			for res in self.nonpolar_polar_sc.values():
				atoms.append(res.get_sidechain_atoms())
		elif type == "nonpolar-ionizable":
			for res in self.nonpolar_ionizable_sc.values():
				atoms.append(res.get_sidechain_atoms())
		elif type == "polar":
			for res in self.polar_sc.values():
				atoms.append(res.get_sidechain_atoms())
		elif type == "polar-ionizable":
			for res in self.polar_ionizable_sc.values():
				atoms.append(res.get_sidechain_atoms())
		elif type == "ionizable":
			for res in self.ionizable_sc.values():
				atoms.append(res.get_sidechain_atoms())
		elif type == "active":
			for res in self.active_sc.values():
				atoms.append(res.get_sidechain_atoms())
		return atoms

	def general_information(self):
		string = ""
		string += "%-6s%s\n" % ("LCHN",str(self.length))
		string += "%-6s%s\n" % ("GAPS",str(self.gap))
		return string

class Residue:

	def __init__(self,Atom_instance):

			# 2013.05.23
			# This likely not a comprehensive list of backbone atom names. 
			# For example, today I had to add the atom names D and DA that
			# originated in PDB files containing neutron diffraction data.
			# The unwanted behavior is that unrecognized backbone atoms 
			# will be classified as sidechain atoms. This, in turn, has
			# affects the calculation of minimum sidechain depth because
			# the sidechain, which very well may be buried in the protein
			# core, may also have backbone atoms that are not. This leads,
			# in relativey rare cases, to mis-classification.
			###############################################################
			self.nonSidechainAtomNames = ['N','CA','C','O','D','H','DA']

			self.num  = Atom_instance.residue_sequence_number
			self.name = Atom_instance.residue_name
			self.chn  = Atom_instance.chain_identifier
			self.key  = (self.num,self.name,self.chn)
			self.branched = 0
			self.missingSidechainAtoms = 0

			# IF THE RESIDUES HAS MORE THAN ONE CONFORMATION,
			# ALL OF THE ATOMS OF ALL OF THE CONFORMATIONS ARE
			# INCLUDED IN THE self.atoms ATTRIBUTE.
			################################################## 
			self.hasAlternateConformation = 0
			self.alt = {}
			self.atoms = {}
			self.atomList = [Atom_instance,] # In order encountered in PDB file
			atom_name = Atom_instance.atom_name.split()[0]
			if Atom_instance.atom_alternate_location:
				self.alt.update({Atom_instance.atom_alternate_location:
								{atom_name:Atom_instance}})
			self.atoms.update({atom_name:Atom_instance})

			# ALTERNATIVE ATOM REPRESENTATIONS.
			###################################
			self.ter  = None
			self.bcom = None
			self.scom = None

			# AN ATTRIBUTE FOR ASSOCIATING ADDITIONAL DATA WITH THE CLASS.
			##############################################################
			self.xtra_info = None

	def addatom(self,Atom_instance):

		# Uses that last alternate conformation of each residue atom encountered in the PDB file.
		# It would be more robust to use all conformations, but I have yet to find a need to justify the trouble.
		#########################################################################################################
		atom_name = Atom_instance.atom_name.split()[0]		
		self.atoms.update({atom_name:Atom_instance})
		self.atomList.append(Atom_instance)
		if Atom_instance.atom_alternate_location:
			if Atom_instance.atom_alternate_location in self.alt:
				self.alt[Atom_instance.atom_alternate_location].update({atom_name:Atom_instance})
			else:
				self.alt.update({Atom_instance.atom_alternate_location:{atom_name:Atom_instance}})

	def change_chain_identifier(self, chain):

		for atom in self.atoms:
			atom = self.atoms[atom]
			atom.chain_identifier = chain

		for atom in self.alt:
			atom = self.alt[atom]
			atom.chain_identifier = chain

		for atom in self.atomList:
			atom.chain_identifier = chain

		self.chn  = chain
		self.key  = (self.num,self.name,self.chn)

	def setAltConfBackboneAtoms(self):
		
		if self.hasAlternateConformation:
			# Update each alternate conformation with its backbone atoms.
						#############################################################
			for confKey in self.alt.keys():
				self.alt[confKey].update(self.atoms)
		
	def setsub(self):

		# Updated 2013.07.23.
		# Link each atom to is residue object.
		######################################
		for atom in self.atoms.values():
			atom.residue = self

		# Also link atoms in alternate conformations.
		#############################################
		for conf in self.alt.values():
				for atom in conf.values():
						atom.residue = self

	def set_sidechain_representations(self):
		self.ter = self.get_terminal_sidechain_atom()
		if not self.ter:
			self.ter = self.get_reduced_sidechain_atom()
			if not self.ter:
				try: self.ter = self.atoms["CA"]
				except: self.ter = self.get_reduced_backbone_atom()
		self.bcom = self.get_reduced_backbone_atom()
		self.scom = self.get_reduced_sidechain_atom()

	def get_atoms(self):
		atoms = []
		if self.alt:
			for conf in self.alt.values():
				for atom in conf.values():
					atoms.append(atom.v)
		else:
			for atom in self.atoms.keys():
				atoms.append(self.atoms[atom].v)
		return atoms        

	def get_backbone_atoms(self):
		atoms = []
		if self.alt:
			for conf in self.alt.values():
				for atom in conf.values():
					if atom.atom_name in self.nonSidechainAtomNames:
						atoms.append(atom)
		else:
			for atom in self.atoms.keys():
				if atom in self.nonSidechainAtomNames:
					atoms.append(self.atoms[atom])
		return atoms

	def get_polar_atoms(self):
		atoms = []
		if self.alt:
			for conf in self.alt.values():
				for atom in conf.values():
					if atom[0] in ["O","N","S"]:
						atoms.append(atom.v)
		else:
			for atom in self.atoms.keys():
				if atom[0] in ["O","N","S"]:
					atoms.append(self.atoms[atom].v)
		return atoms

	def get_sidechain_atoms(self):
		atoms = []
		# v_107 !!!! systemic major issues with self.alt that need to be fixed.
		if 0:#self.alt:
			
			for conf in self.alt.values():
				for atom in conf.values():
					if atom.atom_name not in self.nonSidechainAtomNames:
											if atom not in atoms:
												atoms.append(atom)
					elif atom.residue_name == "GLY":
						try: atoms.append(conf['CA'])
						except: continue
						if not atoms:
							try: atoms.append(conf['CA'])
							except: pass
		else:
			
			for atom in self.atoms.keys():
				if atom not in self.nonSidechainAtomNames:
					atoms.append(self.atoms[atom])
				elif self.atoms[atom].residue_name == "GLY":
					try: atoms.append(self.atoms['CA'])
					except: continue
					if not atoms:
						try: atoms.append(self.atoms['CA'])
						except: pass                  

			#print(atoms)
			return atoms

	def get_terminal_sidechain_atom(self):


		# This function is only relevant for amino acid sidechains.
		###########################################################
		if self.name not in NONPOLAR+POLAR+IONIZABLE:
			return None

		terminal_atoms = {"GLY":"CA", "ALA":"CB", "ILE":"CD1", "MET":"CE", "CYS":"SG",
							"SER":"OG", "LYS":"NZ", "MLY":"NZ", "MLZ":"NZ", "M3L":"NZ"}
		
#        multiple_terminal_atoms = {"THR":("OG1","CG2"),
#                                    "LEU":("CD1","CD2"),
#                                    "PHE":("CG","CD1","CD2","CE1","CE2","CZ"),
#                                    "TYR":("CG","CD1","CD2","CE1","CE2","CZ","OH"),
#                                    "PRO":("CB","CG","CD"),
#                                    "VAL":("CG1","CG2"),
#                                    "ASN":("OD1","ND2"),
#                                    "GLN":("OE1","NE2"),
#                                    "TRP":("CG","CD1","CD2","NE1","CE2","CE3","CZ2","CZ3","CH2"),
#                                    "ARG":("NH1","NH2"),
#                                    "ASP":("OD1","OD2"),
#                                    "GLU":("OE1","OE2"),
#                                    "HIS":("ND1","NE2")}

		# The last atom in the keyed value tuple gets set as the PSA name...
		multiple_terminal_atoms = {"THR":("CG2", "OG1"),
									"LEU":("CD1", "CD2"),
									"PHE":("CZ",),
									"TYR":("OH",),
									"PRO":("CB", "CG", "CD"),
									"VAL":("CG1", "CG2"),
									"ASN":("OD1", "ND2"),
									"GLN":("OE1", "NE2"),
									"TRP":("CZ3", "CH2"),
									"ARG":("NH1", "NH2"),
									"ASP":("OD1", "OD2"),
									"GLU":("OE1", "OE2"),
									"HIS":("ND1", "NE2")}

		# Return the terminal side chain atom.
		######################################
		if self.name in terminal_atoms:
			
			try:
				
				# Create a PSA instance of the terminal side chain atom.
				########################################################
				terAtom = self.atoms[terminal_atoms[self.name]]
				ter = PseudoAtom()
				ter.x = terAtom.x 
				ter.y = terAtom.y
				ter.z = terAtom.z
				ter.atom_serial = terAtom.atom_serial + 1
				ter.residue_sequence_number = terAtom.residue_sequence_number
				ter.atom_name = terAtom.atom_name
				ter.residue = self
				ter.residue_name = ter.residue.name
				ter.chain_identifier = terAtom.chain_identifier
				ter.reinitialize()
				self.ter = ter
				return ter
			
			except:
				
				#print("Warning: Cannot calculate terminal side chain atom because there are no side chain atoms.", self.name, self.num, self.chn)
				self.missingSidechainAtoms = 1
				return None

		# Return the centroid representation of the terminal side chain atom.
		#####################################################################
		else:
			
			self.branched = 1
			n, x, y, z = 0., 0., 0., 0.
			atomSerial = 0
			atomName = "TER"
			residueSequenceNumber = 0
			chain_identifier = ""
			for key in multiple_terminal_atoms[self.name]:
				
				try:

					# Represent the terminal atom of the side chain as the centroid.
					################################################################
					atom = self.atoms[key]
					n += 1 
					x += atom.x
					y += atom.y
					z += atom.z
					if atom.atom_serial > atomSerial:
						atomSerial = atom.atom_serial
					atomName = atom.atom_name
					residueSequenceNumber = atom.residue_sequence_number
					chain_identifier = atom.chain_identifier
			
				except:
					
					#print("Warning: Cannot calculate terminal side chain atom because there are no side chain atoms.", self.name, self.num, self.chn)
					self.missingSidechainAtoms = 1
					return None

			# Using the centroid coordinates create a PSA instance representing the terminal side chain atom.
			#################################################################################################
			ter = PseudoAtom()
			ter.atom_name = atomName 
			ter.residue = self
			ter.atom_serial = atomSerial + 1
			ter.residue_sequence_number = residueSequenceNumber
			ter.residue_name = ter.residue.name
			ter.chain_identifier = chain_identifier
			ter.x = x/n; ter.y = y/n; ter.z = z/n
			ter.reinitialize()
			self.ter = ter
			return ter
			
	def get_preterminal_sidechain_atom(self):

		# IF THE RESIDUE HAS AN ALTERNATE CONFORMATION, BY DEFAULT
		# REPRESENT IT WITH ITS C-ALPHA ATOM.
				###########################################################
		if self.alt:
			try: return self.atoms["CA"]
			except: return None

		# THIS FUNCTION IS ONLY RELEVANT FOR AMINO ACID SIDECHAINS.
				###########################################################
		if self.name not in NONPOLAR+POLAR+IONIZABLE:
			return None

		preterminal_atoms = {"ALA":"CA", "ARG":"CZ", "ASN":"CG", "ASP":"CG", "CYS":"CB",
								"GLN":"CD", "GLU":"CD", "HIS":"CG", "ILE":"CG1", "LEU":"CG",
								"LYS":"CE", "MET":"SD", "PHE":"CB", "PRO":"CA", "SER":"CB",
								"THR":"CB", "TRP":"CB", "TYR":"CZ", "VAL":"CB"}

		if self.name in preterminal_atoms:
			try:
				ter = self.atoms[preterminal_atoms[self.name]]
				ter.residue = self
				ter.reinitialize()
				return ter
			except:
				return None

	def get_reduced_backbone_atom(self):

		# IF THE RESIDUE HAS AN ALTERNATE CONFORMATION, BY DEFAULT
		# REPRESENT IT WITH ITS C-ALPHA ATOM.
		###########################################################
		if self.alt:
			try: return self.atoms["CA"].v
			except: return None

		# THIS FUNCTION IS ONLY RELEVANT FOR AMINO ACID SIDECHAINS.
		###########################################################
		if self.name not in NONPOLAR+POLAR+IONIZABLE:
			return None

		"""Represent the backbone atoms of the residue by their center
			of mass.
		"""
		atoms = self.get_backbone_atoms()
		n, x, y, z = 0., 0., 0., 0.
		for atom in atoms:
			n += 1; x += atom.x; y += atom.y; z += atom.z
		com = Atom()
		com.record_name = "ATOM  "
		keys = sorted(self.atoms)
		com.atom_serial = self.atoms[keys[0]].atom_serial
		com.atom_name = "COM"
		com.residue_name = self.name
		com.chain_identifier = self.chn
		com.residue_sequence_number = self.num
		com.residue_key = self.key
		# There are backbone atoms for the sidechain.
		#############################################
		if n:
			com.x = x/n; com.y = y/n; com.z = z/n
		# There are no backbone atoms for the sidechain.
		################################################
		else:
			com.x = 0.; com.y = 0.; com.z = 0.
		com.residue = Residue(com)
		com.residue.name = self.name
		com.residue.num = self.num
		com.residue.chn = self.chn
		com.reinitialize()
		com.residue = self
		return com.v

	def get_reduced_sidechain_atom(self):

		# IF THE RESIDUE HAS AN ALTERNATE CONFORMATION, BY DEFAULT
		# REPRESENT IT WITH ITS C-ALPHA ATOM.
		###########################################################
		if self.alt:
			try: return self.atoms["CA"].v
			except: return None

		# THIS FUNCTION IS ONLY RELEVANT FOR AMINO ACID SIDECHAINS.
		###########################################################
		if self.name not in NONPOLAR+POLAR+IONIZABLE:
			return None

		"""Represent the sidechain atoms of the residue by their center
			of mass.
		"""
		atoms = self.get_sidechain_atoms()
		if not atoms:
			return None
		n, x, y, z = 0., 0., 0., 0.
		for atom in atoms:
			n += 1; x += atom.x; y += atom.y; z += atom.z
		com = Atom()
		com.record_name = "ATOM  "
		keys = sorted(self.atoms)
		com.atom_serial = self.atoms[keys[0]].atom_serial
		com.atom_name = "COM"
		com.residue_name = self.name
		com.chain_identifier = self.chn
		com.residue_sequence_number = self.num
		com.residue_key = self.key
		# There are sidechain atoms.
		#############################################
		if n:
			com.x = x/n; com.y = y/n; com.z = z/n
		# There are no sidechain atoms.
		################################################
		else:
			com.x = 0.; com.y = 0.; com.z = 0.
		com.residue = Residue(com)
		com.residue.name = self.name
		com.residue.num = self.num
		com.residue.chn = self.chn
		com.reinitialize()
		com.residue = self
		return com

	def __repr__(self):
		string = ""
		atoms = {}
		for atom in self.atoms.values():
			atoms[atom.atom_serial] = atom
		sortedKeys = sorted(atoms.keys())
		for sortedKey in sortedKeys:
			string += str(atoms[sortedKey])
		return string

class Atom:

	def __init__(self, pdb_file_line="", cif_options={}):

		# From PDB
		self.format = "cif"
		self.pdb_file_line = pdb_file_line
		self.cif_options = cif_options
		self.pdb_code = " "
		self.record_name = " "
		self.atom_serial = 0
		self.atom_name = " "
		self.atom_alternate_location = " "
		self.residue_name = ""
		self.chain_identifier = " "
		self.entity_id = " "
		self.residue_sequence_number = 0
		self.residue_insertion_code = " "
		self.x = 0.
		self.y = 0.
		self.z = 0.
		self.occupancy = 0.
		self.temperature_factor = 0.
		self.sigma_x = 0.
		self.sigma_y = 0.
		self.sigma_z = 0.
		self.sigma_occupancy = 0.
		self.sigma_temperature_factor = 0.
		self.segment_identifier = " "
		self.symbol = " "
		self.charge = 0

		# Isom lab classifiers
		self.atom_key = " "
		self.residue_key = " "
		self.residue = None
		self.exposed, self.extended = 0, 0
		self.margin = 0
		self.core = 0

		# CIF information
		self.cif_dictionary = {}

		# CIF output format field widths
		self.cif_field_widths = {}

		# Parse the ATOM information
		if self.pdb_file_line and cif_options:

			atom_details, cif_option_keys = self.pdb_file_line.strip().split(), list(cif_options)
			i = 0
			while i < len(atom_details):

				atom_detail = atom_details[i]
				cif_option = cif_option_keys[i]

				if cif_option == "_atom_site.group_PDB":
					try:
						self.record_name = "%-6s" % atom_detail
						self.cif_dictionary["_atom_site.group_PDB"] = self.record_name
					except:
						pass
				if cif_option == "_atom_site.id":
					try:
						self.atom_serial = int(atom_detail)
						self.cif_dictionary["_atom_site.id"] = self.atom_serial
					except:
						pass
				if cif_option == "_atom_site.type_symbol":
					try:
						self.symbol = atom_detail
						self.cif_dictionary["_atom_site.type_symbol"] = self.symbol
					except:
						pass
				if cif_option == "_atom_site.label_atom_id":
					try:
						self.atom_name = atom_detail
						self.cif_dictionary["_atom_site.label_atom_id"] = self.atom_name
					except:
						pass
				if cif_option == "_atom_site.label_alt_id":
					try:
						self.atom_alternate_location = atom_detail
						self.cif_dictionary["_atom_site.label_alt_id"] = self.atom_alternate_location
					except:
						pass
				if cif_option == "_atom_site.label_comp_id":
					try:
						self.residue_name = atom_detail
						self.cif_dictionary["_atom_site.label_comp_id"] = self.residue_name
					except:
						pass
				if cif_option == "_atom_site.label_asym_id":
					try:
						self.chain_identifier = atom_detail
						self.cif_dictionary["_atom_site.label_asym_id"] = self.chain_identifier
					except:
						pass
				if cif_option == "_atom_site.label_entity_id":
					try:
						self.entity_id = atom_detail
						self.cif_dictionary["_atom_site.label_entity_id"] = self.entity_id
					except:
						pass
				if cif_option == "_atom_site.label_seq_id":
					try:
						self.residue_sequence_number = int(atom_detail)
						self.cif_dictionary["_atom_site.label_seq_id"] = self.residue_sequence_number
					except:
						pass
				if cif_option == "_atom_site.pdbx_PDB_ins_code":
					try:
						self.residue_insertion_code = atom_detail
						self.cif_dictionary["_atom_site.pdbx_PDB_ins_code"] = self.residue_insertion_code
					except:
						pass
				if cif_option == "_atom_site.Cartn_x":
					try:
						self.x = float(atom_detail)
						self.cif_dictionary["_atom_site.Cartn_x"] = self.x
					except:
						pass
				if cif_option == "_atom_site.Cartn_y":
					try:
						self.y = float(atom_detail)
						self.cif_dictionary["_atom_site.Cartn_y"] = self.y
					except:
						pass
				if cif_option == "_atom_site.Cartn_z":
					try:
						self.z = float(atom_detail)
						self.cif_dictionary["_atom_site.Cartn_z"] = self.z
					except:
						pass
				if cif_option == "_atom_site.occupancy":
					try:
						self.occupancy = float(atom_detail)
						self.cif_dictionary["_atom_site.occupancy"] = self.occupancy
					except:
						pass
				if cif_option == "_atom_site.B_iso_or_equiv":
					try:
						self.temperature_factor = float(atom_detail)
						self.cif_dictionary["_atom_site.B_iso_or_equiv"] = self.temperature_factor
					except:
						pass
				if cif_option == "_atom_site.Cartn_x_esd":
					try:
						self.sigma_x = float(atom_detail)
						self.cif_dictionary["_atom_site.Cartn_x_esd"] = self.sigma_x
					except:
						pass
				if cif_option == "_atom_site.Cartn_y_esd":
					try:
						self.sigma_y = float(atom_detail)
						self.cif_dictionary["_atom_site.Cartn_y_esd"] = self.sigma_y
					except:
						pass
				if cif_option == "_atom_site.Cartn_z_esd":
					try:
						self.sigma_z = float(atom_detail)
						self.cif_dictionary["_atom_site.Cartn_z_esd"] = self.sigma_z
					except:
						pass
				if cif_option == "_atom_site.occupancy_esd":
					try:
						self.sigma_occupancy = float(atom_detail)
						self.cif_dictionary["_atom_site.occupancy_esd"] = self.sigma_occupancy
					except:
						pass
				if cif_option == "_atom_site.B_iso_or_equiv_esd":
					try:
						self.sigma_temperature_factor = float(atom_detail)
						self.cif_dictionary["_atom_site.B_iso_or_equiv_esd"] = self.sigma_temperature_factor
					except:
						pass
				if cif_option == "_atom_site.pdbx_formal_charge":
					try:
						self.charge = float(atom_detail)
						self.cif_dictionary["_atom_site.pdbx_formal_charge"] = self.charge
					except:
						pass
				# if cif_option == "_atom_site.auth_seq_id":
				# 	self.record_name = atom_detail
				# if cif_option == "_atom_site.auth_comp_id":
				# 	self.record_name = atom_detail
				# if cif_option == "_atom_site.auth_asym_id":
				# 	self.record_name = atom_detail
				# if cif_option == "_atom_site.auth_atom_id":
				# 	self.record_name = atom_detail
				# if cif_option == "_atom_site.pdbx_PDB_model_num":
				# 	self.record_name = int(atom_detail)

				i += 1

			# Set key values
			self.residue_key = (self.residue_sequence_number, self.residue_name, self.chain_identifier)
			self.atom_key = (self.residue_sequence_number, self.residue_name, self.chain_identifier, self.atom_name)

			# Create vertex object for computational geometry applications
			self.v = Vertex((self.x,self.y,self.z), data=self, unique_id=self.atom_serial)

	def reinitialize(self):

		if self.residue:
			self.residue.num  = self.residue_sequence_number
			self.residue.name = self.residue_name
			self.residue.chn  = self.chain_identifier
		self.residue_key = (self.residue_sequence_number, self.residue_name, self.chain_identifier)
		self.atom_key = (self.residue_sequence_number, self.residue_name, self.chain_identifier, self.atom_name)
		self.v = Vertex((self.x,self.y,self.z), data=self, unique_id=self.atom_serial)

	# def copy(self):

	# 	return Atom(str(self))

	# def get_cif_options_string(self):

	# 	cif_format = "#\n"
	# 	cif_format += "loop_\n"
	# 	for cif_option in self.cif_options:
	# 		if cif_option in self.cif_dictionary:
	# 			cif_format += cif_option + "\n"
	# 	return cif_format

	def get_pdb_format(self):

		import pdbFile
		atom = pdbFile.Atom()
		atom.record_name = self.record_name
		atom.atom_serial = self.atom_serial
		atom.atom_name = self.atom_name
		if self.atom_alternate_location == ".":
			self.atom_alternate_location = ""
		atom.alternate_location = self.atom_alternate_location
		atom.residue_name = self.residue_name
		atom.chain_identifier = self.chain_identifier
		atom.residue_sequence_number = self.residue_sequence_number
		if self.residue_insertion_code == "?":
			self.residue_insertion_code = ""
		atom.residue_insertion_code = self.residue_insertion_code 
		atom.x = self.x
		atom.y = self.y
		atom.z = self.z
		atom.occupancy = self.occupancy
		atom.temperature_factor = self.temperature_factor
		atom.segment_identifier = self.segment_identifier
		atom.symbol = self.symbol 
		atom.charge = self.charge
		return str(atom)

	def __repr__(self):

		atom_format = ""
		for cif_option in self.cif_options:
			if cif_option in self.cif_dictionary:
				atom_format += str(self.cif_dictionary[cif_option]) + " "
		return atom_format + "\n"	

class Atom_LO:

	def __init__(self):
		self.atom_serial = None
		self.residueKey = None
		self.repr = ""

	def __repr__(self):
		return self.repr

class PseudoAtom:

	def __init__(self):

		self.format = "pdb"
		self.pdbfileline = ""
		self.record_name = "ATOM  "
		self.atom_serial = 1
		self.atom_name = " O  " #one space to the left, two to the right
		self.atom_alternate_location = ""
		self.residue_name = "PSA"
		self.chain_identifier = "A"
		self.residue_sequence_number= 1
		self.residue_insertion_code = ""
		self.x = 0.
		self.y = 0.
		self.z = 0.
		self.occupancy = 0.
		self.temperature_factor = 0.
		self.segment_identifier = ""
		self.symbol = ""
		self.charge = 0
		self.atom_key = (self.residue_sequence_number, self.residue_name, self.chain_identifier, self.atom_name)
		self.residue_key = (self.residue_sequence_number, self.residue_name, self.chain_identifier)
		self.residue = Residue(self)
		self.exposed, self.extended = 0, 0
		self.margin = 0
		self.core = 0
		self.pdbCode = ""

		# INSTANTIATE A LOW-OVERHEAD ATOM INSTANCE.
		###########################################
		self.Atom_LO = Atom_LO()
		self.Atom_LO.atom_serial = self.atom_serial
		self.Atom_LO.residueKey = self.residue_key
		self.Atom_LO.repr = self.__repr__()

		self.v = Vertex((self.x,self.y,self.z), data=self, unique_id=self.atom_serial)

	def reinitialize(self):
		
		# Reset the residue key information.
		####################################
		if self.residue:
			self.residue.num  = self.residue_sequence_number
			self.residue.name = self.residue_name
			self.residue.chn  = self.chain_identifier
		self.residue_key = (self.residue_sequence_number, self.residue_name, self.chain_identifier)
		self.atom_key = (self.residue_sequence_number, self.residue_name, self.chain_identifier, self.atom_name)

		# REINITIALIZE ATOM_LO.
		#######################
		self.Atom_LO = Atom_LO()
		self.Atom_LO.atom_serial = self.atom_serial
		self.Atom_LO.residueKey = self.residue_key
		self.Atom_LO.repr = self.__repr__()

		# REINITIALIZE V.
		#################
		self.v = Vertex((self.x,self.y,self.z), data=self, unique_id=self.atom_serial)
							
	def copy(self):
		return Atom(str(self))

	def __repr__(self, forceHex=0):
		
		if self.chain_identifier == "NULL":
			chain_identifier = " "
		else:
			chain_identifier = self.chain_identifier

		# # Robust representation of atom_serial.
		# #######################################
		# stringAtomSerial = ""
		# if self.hex_atom_serial and forceHex:
		# 	stringAtomSerial = "%5s" % hex(self.atom_serial)[2:]
		# else:
		# 	stringAtomSerial = "%5i" % self.atom_serial
		# #if self.atom_serial <= 99999:
		# #    stringAtomSerial = "%5i" % self.atom_serial
		# #else:
		# #    stringAtomSerial = "%5s" % hex(self.atom_serial)[2:]

		# # Robust representation of residue_sequence_number.
		# ####################################################   
		# stringResidueNumber = ""
		# if self.hex_residue_sequence_number and forceHex:
		# 	stringResidueNumber = "%3s" % hex(self.residue_sequence_number)[2:]
		# else:
		# 	stringResidueNumber = "%3i" % self.residue_sequence_number
		# #if self.residue_sequence_number <= 999:
		# #    stringResidueNumber = "%03i" % self.residue_sequence_number
		# #else:
		# #    stringResidueNumber = "%3s" % hex(self.residue_sequence_number)[2:]

		# fixed
		# Format the atom data in standard PDB file format.
		###################################################
		pdbformat  = "%6s%5s %-4s%1s%3s %-1s%4s%1s"
		# Check for 2-character chain identifier.
		# If so, adjust chain and residue number fields accordingly.
		############################################################
		if len(chain_identifier) > 1:
			pdbformat  = "%6s%5s %-4s%1s%3s %-2s%3s%1s"
		pdbformat += "   %8.3f%8.3f%8.3f%6.2f%6.2f      %-4s%2s%2s\n"

		return pdbformat % (self.record_name,
							self.atom_serial,
							self.atom_name,
							self.atom_alternate_location,
							self.residue_name,
							self.chain_identifier,
							self.residue_sequence_number,
							self.residue_insertion_code,
							self.x,
							self.y,
							self.z,
							self.occupancy,
							self.temperature_factor,
							self.segment_identifier,
							self.symbol,
							self.charge)
