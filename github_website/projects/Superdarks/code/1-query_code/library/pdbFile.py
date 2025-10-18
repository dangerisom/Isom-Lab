# Copyright (c) 2025 Isom Lab
# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.
# If a copy of the MPL was not distributed with this file, You can obtain one at
# https://mozilla.org/MPL/2.0/.

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

class PDBfile:

	def __init__(self, pdbFilePath="", pdbFileName="", pdbFileAsString="", twoCharacterChain=0, zip_status=0):

		# Open the PDB file, or set of PDB files.
		#########################################
		self.zip = zip_status
		self.pdbFilePath = pdbFilePath
		self.pdbFileName = pdbFileName
		self.pdbCode = pdbFileName.split(".pdb")[0]
		self.chain_selected_in_file_name, self.chains = 0, []
		if "-chains." in self.pdbCode:
			self.chain_selected_in_file_name = 1
			chain_string = self.pdbCode.split("-chains.")[1]
			for chain in chain_string.split("."):
				if chain:
					self.chains.append(chain)

		self.atoms = {}
		self.hetatoms = {}
		self.residues, self.res_bychain = {}, {}
		self.het_residues, self.het_bychain = {}, {}

		self.isNmrStructure = 0

		# PDB files that need to be re-written.
		#######################################
		self.rewritePdbFiles = []

		# Check PDB files for NULL chains. Rename chain A.
		##################################################
		makeChainA = 0
		if self.zip:
			pdbfile = gzip.GzipFile(self.pdbFilePath + self.pdbFileName,"rb").readlines()
			for line in pdbfile:
				line = line.decode() # Have to do because of Gzipfile byte issue in Python3
				if line[0:4] == "ATOM" and line[21:22] == " ":
					makeChainA = 1
					self.rewritePdbFiles.append(self.pdbFilePath + self.pdbFileName)
					print("PDB file has NULL chain.  Rename as chain A....")
					break
		else:
			# pdbfile = open(self.pdbFilePath + self.pdbFileName, 'U').readlines()
			pdbfile = open(self.pdbFilePath + self.pdbFileName, 'r').readlines() # Change for Python 3.11
			for line in pdbfile:
				if line[0:4] == "ATOM" and line[21:22] == " ":
					makeChainA = 1
					self.rewritePdbFiles.append(self.pdbFilePath + self.pdbFileName)
					print("PDB file has NULL chain.  Rename as chain A....")
					break

		# Check PDB files for insertion codes.
		######################################
		if self.zip:
			pdbfile = gzip.GzipFile(self.pdbFilePath + self.pdbFileName,"rb").readlines()
			for line in pdbfile:
				line = line.decode() # Have to do because of Gzipfile byte issue in Python3
				if line[0:4] == "ATOM" and line[26:27] != " ":
					self.rewritePdbFiles.append(self.pdbFilePath + self.pdbFileName)
					print("PDB file has residue insertion codes.  Must renumber....")
					break
		else:
			# pdbfile = open(self.pdbFilePath + self.pdbFileName, 'U').readlines()
			pdbfile = open(self.pdbFilePath + self.pdbFileName, 'r').readlines() # Change for Python 3.11
			for line in pdbfile:
				if line[0:4] == "ATOM" and line[26:27] != " ":
					self.rewritePdbFiles.append(self.pdbFilePath + self.pdbFileName)
					print("PDB file has residue insertion codes.  Must renumber....")
					break
					
		# If necessary, overwrite PDB files to remove insertion codes.
		##############################################################
		if self.rewritePdbFiles:
			i = 0
			fixedPdbFileLines = ""
			if self.zip:
				pdbfile = gzip.GzipFile(self.pdbFilePath + self.pdbFileName,"rb").readlines()
				for line in pdbfile:
					line = line.decode() # Have to do because of Gzipfile byte issue in Python3
					if line[0:4] == "ATOM":
						atom = Atom(line)
						if makeChainA:
							atom.chain_identifier = "A"
						if atom.atom_name == " N":
							i += 1
						atom.residue_sequence_number = i
						atom.residue_insertion_code = ""
						atom.reinitialize()
						fixedPdbFileLines += str(atom)
					else:
						fixedPdbFileLines += line  
				# Rewrite the PDB file.
				#######################
				pdbFileRewrite = gzip.GzipFile(self.pdbFilePath + self.pdbFileName,"wb")
				pdbFileRewrite.write(fixedPdbFileLines)
				pdbFileRewrite.close()
			else:
				pdbfile = open(self.pdbFilePath + self.pdbFileName, 'r').readlines()
				for line in pdbfile:
					if line[0:4] == "ATOM":
						atom = Atom(line)
						#print("1", atom)
						if makeChainA:
							atom.chain_identifier = "A"
						if atom.atom_name == " N":
							i += 1
						atom.residue_sequence_number = i
						atom.residue_insertion_code = ""
						atom.reinitialize()
						#print("2", atom)
						fixedPdbFileLines += str(atom)
					else:
						fixedPdbFileLines += line
				# Rewrite the PDB file.
				#######################
				pdbFileRewrite = open(self.pdbFilePath + self.pdbFileName,"w")
				pdbFileRewrite.write(fixedPdbFileLines)
				pdbFileRewrite.close()

		# Check PDB files for hexidecimal atom serial numbers or residue sequence numbers.
		##################################################################################
		self.rewritePdbFiles = []
		if self.zip:
			pdbfile = gzip.GzipFile(self.pdbFilePath + self.pdbFileName,"rb").readlines()
			for line in pdbfile:
				line = line.decode() # Have to do because of Gzipfile byte issue in Python3
				if line[0:4] == "ATOM" or line[0:6] == "HETATM":
					atomInstance = Atom(line)
					#print(atomInstance.hex_atom_serial, atomInstance.hex_residue_sequence_number)
					if atomInstance.hex_atom_serial or atomInstance.hex_residue_sequence_number:
						self.rewritePdbFiles.append(self.pdbFilePath + self.pdbFileName)
						print("PDB file has hexidecimal atom serial or residue numbers.  Must renumber....")
						print("By default, a single-character chain identifier, 5-digit atom serial, and 4-digit residue number is assumed.")
						break
		else:
			# pdbfile = open(self.pdbFilePath + self.pdbFileName, 'U').readlines()
			pdbfile = open(self.pdbFilePath + self.pdbFileName, 'r').readlines() # Change for Python 3.11
			for line in pdbfile:
				if line[0:4] == "ATOM" or line[0:6] == "HETATM":
					atomInstance = Atom(line)
					if atomInstance.hex_atom_serial or atomInstance.hex_residue_sequence_number:
						self.rewritePdbFiles.append(self.pdbFilePath + self.pdbFileName)
						print("PDB file has hexidecimal atom serial or residue numbers.  Must renumber....")
						print("By default, a single-character chain identifier, 5-digit atom serial, and 4-digit residue number is assumed.")
						break

		# If necessary, overwrite PDB files to remove insertion codes.
		##############################################################
		if self.rewritePdbFiles:
			i = 0
			fixedPdbFileLines = ""
			if self.zip:
				pdbfile = gzip.GzipFile(self.pdbFilePath + self.pdbFileName,"rb").readlines()
				for line in pdbfile:
					line = line.decode() # Have to do because of Gzipfile byte issue in Python3
					if line[0:4] == "ATOM" or line[0:6] == "HETATM":
						atom = Atom(line)
						# The default behavior of invoking an Atom instance is to remove the hexidecimal value.
						fixedPdbFileLines += str(atom)
					else:
						fixedPdbFileLines += line  
				# Rewrite the PDB file.
				#######################
				pdbFileRewrite = gzip.GzipFile(self.pdbFilePath + self.pdbFileName,"wb")
				pdbFileRewrite.write(fixedPdbFileLines)
				pdbFileRewrite.close()
			else:
				pdbfile = open(self.pdbFilePath + self.pdbFileName, 'U').readlines()
				for line in pdbfile:
					if line[0:4] == "ATOM" or line[0:6] == "HETATM":
						atom = Atom(line)
						# The default behavior of invoking an Atom instance is to remove the hexidecimal value.
						fixedPdbFileLines += str(atom)
					else:
						fixedPdbFileLines += line
				# Rewrite the PDB file.
				#######################
				pdbFileRewrite = open(self.pdbFilePath + self.pdbFileName,"w")
				pdbFileRewrite.write(fixedPdbFileLines)
				pdbFileRewrite.close()

		# Check PDB files for incorrectly formatted heteroatom record names.
		####################################################################
		self.rewritePdbFilesHeteroAtoms = []
		if self.zip:
			pdbfile = gzip.GzipFile(self.pdbFilePath + self.pdbFileName,"rb").readlines()
			for line in pdbfile:
				line = line.decode() # Have to do because of Gzipfile byte issue in Python3
				if line[0:4] == "ATOM" and line[17:20] not in ALL_SIDECHAINS:
					self.rewritePdbFilesHeteroAtoms.append(self.pdbFilePath + self.pdbFileName)
					print("PDB file has incorrectly formatted heteroatoms.  Must reformat....")
					break
		else:
			# pdbfile = open(self.pdbFilePath + self.pdbFileName, 'U').readlines()
			pdbfile = open(self.pdbFilePath + self.pdbFileName, 'r').readlines() # Change for Python 3.11
			for line in pdbfile:
				if line[0:4] == "ATOM" and line[17:20] not in ALL_SIDECHAINS:
					self.rewritePdbFilesHeteroAtoms.append(self.pdbFilePath + self.pdbFileName)
					print("PDB file has incorrectly formatted heteroatoms.  Must reformat....")
					break

		# If necessary, overwrite PDB files to properly format heteroatoms record names.
		###############################################################################
		if self.rewritePdbFilesHeteroAtoms:
			fixedPdbFileLines = ""
			if self.zip:
				pdbfile = gzip.GzipFile(self.pdbFilePath + self.pdbFileName,"rb").readlines()
				for line in pdbfile:
					line = line.decode() # Have to do because of Gzipfile byte issue in Python3
					if line[0:4] == "ATOM" and line[17:20] not in ALL_SIDECHAINS:
						atom = Atom(line)
						atom.record_name = "HETATM"
						atom.reinitialize()
						fixedPdbFileLines += str(atom)
					else:
						fixedPdbFileLines += line  
				# Rewrite the PDB file.
				#######################
				pdbFileRewrite = gzip.GzipFile(self.pdbFilePath + self.pdbFileName,"wb")
				pdbFileRewrite.write(fixedPdbFileLines)
				pdbFileRewrite.close()
			else:
				# pdbfile = open(self.pdbFilePath + self.pdbFileName, 'U').readlines()
				pdbfile = open(self.pdbFilePath + self.pdbFileName, 'r').readlines() # Change for Python 3.11
				for line in pdbfile:
					if line[0:4] == "ATOM" and line[17:20] not in ALL_SIDECHAINS:
						atom = Atom(line)
						atom.record_name = "HETATM"
						atom.reinitialize()
						fixedPdbFileLines += str(atom)
					else:
						fixedPdbFileLines += line
				# Rewrite the PDB file.
				#######################
				pdbFileRewrite = open(self.pdbFilePath + self.pdbFileName,"w")
				pdbFileRewrite.write(fixedPdbFileLines)
				pdbFileRewrite.close()

		# Open the set of one or more PDB files using the absolute file paths listed in pdbFilePathList.
		################################################################################################
		self.pdbfile = []
		if self.zip:
			self.pdbfile += gzip.GzipFile(self.pdbFilePath + self.pdbFileName,"rb").readlines()
		else:
			# self.pdbfile += open(self.pdbFilePath + self.pdbFileName, 'U').readlines()
			self.pdbfile += open(self.pdbFilePath + self.pdbFileName, 'r').readlines() # Change for Python 3.11

		# self.numchains  = 0
		# self.numligands = 0

		# # THE COMPOUND INFORMATION FOR THE FILE.
		# ########################################
		# self.compound = ""

		# # ALL OF THE ATOMS AND HETERO ATOMS IN THE PDB FILE.
		# ####################################################
		# self.atoms    = {}
		# self.hetatoms = {}

		# # MAX AND MIN COORDINATES.
		# ##########################
		# self.maxx = -1000000.
		# self.maxy = -1000000.
		# self.maxz = -1000000.
		# self.minx = 1000000.
		# self.miny = 1000000.
		# self.minz = 1000000.

		# # # # COLLECT CHAIN IDENTIFIERS.
		# # # ############################
		# # # self.chains = []

		# # # ORGANIZED RESIDUES.
		# # #####################
		# # self.residues, self.res_bychain = {}, {}
		# # self.het_residues, self.het_bychain = {}, {}

		# # # WATER MOLECULES.
		# # ##################
		# # self.water, self.wat_bychain = {}, {}

		# # # Identify hexidecimal integers in atom serial and residue numbers. 
		# # # Such occurrences indicate a truly legal PDB file.  All other files, 
		# # # specifically those with more than 999 residues must be legalized
		# # # to accommodate 2-character chain identifiers.
		# # ####################################################################
		# # #legalize = 1
		# # #for line in self.pdbfile:
		# # #
		# # #    # "3e8" is the hexidecimal value for 1000.
		# # #    ##########################################
		# # #    if twoCharacterChain:
		# # #        legalize = 0
		# # #    elif "3e8" in line:
		# # #        print("PDB file has > 99,999 atom serial numbers or > 999 residues sequence numbers and is legal.")
		# # #        legalize = 0
		# # #        break
		
		i, resnum, previousAtom = 0, 0, Atom()
		for line in self.pdbfile:

			if self.zip:
				line = line.decode() # Have to do because of Gzipfile byte issue in Python3
			
			# Identify NMR structures
			#########################
			if "EXPDTA" in line:
				self.isNmrStructure = 1
			
			# Only use first model in NMR structure files.
			##############################################
			if "ENDMDL" in line and self.isNmrStructure:
				break

			# Mark the end of one chain and the start of the next.
			# Also, reset residue counter.
			######################################################
			if line[0:3] == "TER":
				i = 0

			# if line[0:6] == "COMPND":
			# 	self.compound += line
			
			# Smart atom parse that adaptively deals with 1 and 2-character chain names, 
			# atom serial numbers that exceed 99,999 entries, and residue numbers that exceed 999.
			#######################################################################################
			if line[0:4] == "ATOM":

				# There may be more than 999 residues, but the residue number in not in hexidecimal format.
				# The slice of the chain and residue_sequence_number indices have been changed to accommodate
				# this discrepancy and "legalize" the entry.  Further output the PDB file line will
				# be in the legal formation (i.e. hexidecimal format).
				##############################################################################################
				#if legalize:
				#    atom = Atom(line, legalize=1)
				#else:
				#    atom = Atom(line)
				atom = Atom(line, twoCharacterChain=twoCharacterChain)
				if atom.chain_identifier == "NULL":
					atom.chain_identifier = "A"
					atom.reinitialize()

				# Catch the transition to > 999 residues. A rare occurrence that needs to be legalized.
				#######################################################################################
				#if not legalize and atom.residue_sequence_number == 0 and previousAtom.residue_sequence_number == 999:
				#    legalize = 1

				# Save previous atom information in the case it is needed for legalization status.
				##################################################################################
				previousAtom = atom

				# Update chains identifiers.
				############################
				if not self.chain_selected_in_file_name and atom.chain_identifier not in self.chains:
					self.chains.append(atom.chain_identifier)

				# # Survey for max and min of each dimension.
				# # This is for all the atoms in the file, not individual protein chains.
				# #######################################################################
				# if atom.x > self.maxx:
				# 	self.maxx = atom.x
				# if atom.y > self.maxy:
				# 	self.maxy = atom.y
				# if atom.z > self.maxz:
				# 	self.maxz = atom.z
				# if atom.x < self.minx:
				# 	self.minx = atom.x
				# if atom.y < self.miny:
				# 	self.miny = atom.y
				# if atom.z < self.minz:
				# 	self.minz = atom.z
				
				self.atoms.update({atom.atom_serial:atom})

				# Residues.
				######################################################################################
				if atom.residue_key in self.residues:
					chain   = atom.chain_identifier
					residue = self.residues[atom.residue_key]
					self.residues[atom.residue_key].addatom(atom)
					residue.addatom(atom)
					self.res_bychain[chain].update({atom.residue_key:residue})
				else:
					chain   = atom.chain_identifier
					residue = Residue(atom)
					self.residues.update({residue.key:residue})
					if atom.chain_identifier not in self.res_bychain:
						self.res_bychain.update({chain:
												{atom.residue_key:residue}})
					else:
						self.res_bychain[atom.chain_identifier].update({atom.residue_key:residue})

				# Update chains identifiers
				if not self.chain_selected_in_file_name and atom.chain_identifier not in self.chains:
					self.chains.append(atom.chain_identifier)

			if line[0:6] == "HETATM":

				atom = Atom(line, twoCharacterChain=twoCharacterChain)
				self.hetatoms.update({atom.atom_serial:atom})

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

				# # Water.
				# #######################################################################################
				# if atom.residue_name == "HOH":
				# 	self.hetatoms.update({atom.atom_serial:atom})
				# 	if atom.residue_key in self.water:
				# 		chain   = atom.chain_identifier
				# 		residue = self.water[atom.residue_key]
				# 		self.water[atom.residue_key].addatom(atom)
				# 		self.wat_bychain[chain].update({atom.residue_key:residue})
				# 	else:
				# 		chain   = atom.chain_identifier
				# 		residue = Residue(atom)
				# 		self.water.update({residue.key:residue})
				# 		if atom.chain_identifier not in self.wat_bychain:
				# 			self.wat_bychain.update({chain:
				# 										{atom.residue_key:residue}})
				# 		else:
				# 			self.wat_bychain[atom.chain_identifier].update({atom.residue_key:residue})
				# # Ligands.
				# #######################################################################################
				# else:
				# 	self.hetatoms.update({atom.atom_serial:atom})
				# 	if atom.residue_key in self.het_residues:
				# 		chain   = atom.chain_identifier
				# 		residue = self.het_residues[atom.residue_key]
				# 		self.het_residues[atom.residue_key].addatom(atom)
				# 		self.het_bychain[chain].update({atom.residue_key:residue})
				# 	else:
				# 		chain   = atom.chain_identifier
				# 		residue = Residue(atom)
				# 		self.het_residues.update({residue.key:residue})
				# 		if atom.chain_identifier not in self.het_bychain:
				# 			self.het_bychain.update({atom.chain_identifier:
				# 										{atom.residue_key: residue}})
				# 		else:
				# 			self.het_bychain[atom.chain_identifier].update({atom.residue_key:residue})


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

	def __init__(self,pdbfilename,pdbcode,chain_identifier,
					resdict,peptide=1,het_residues=None):

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
			if Atom_instance.alternate_location:
				self.alt.update({Atom_instance.alternate_location:
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
		if Atom_instance.alternate_location:
			if Atom_instance.alternate_location in self.alt:
				self.alt[Atom_instance.alternate_location].update({atom_name:Atom_instance})
			else:
				self.alt.update({Atom_instance.alternate_location:{atom_name:Atom_instance}})

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

	def __init__(self,pdbfileline='', twoCharacterChain=0, legalize=0):

		self.format = "pdb"
		self.pdbfileline = pdbfileline
		self.record_name = ""
		self.atom_serial, self.hex_atom_serial = 0, 0
		self.atom_name = ""
		self.alternate_location = ""
		self.residue_name = ""
		self.chain_identifier = ""
		self.residue_sequence_number, self.hex_residue_sequence_number = 0, 0
		self.residue_insertion_code = ""
		self.x = 0.
		self.y = 0.
		self.z = 0.
		self.occupancy = 0.
		self.temperature_factor = 0.
		self.segment_identifier = ""
		self.symbol = ""
		self.charge = 0
		self.atom_key = ""
		self.residue_key = ""
		self.residue = None
		self.exposed, self.extended = 0, 0
		self.margin = 0
		self.core = 0
		self.pdbCode = ""

		if self.pdbfileline:
			
			line = self.pdbfileline.split('\n')[0]
			self.record_name = line[0:6]  
			try:
				# < 99,999 atoms in the PDB file.
				#################################
				self.atom_serial = int(line[6:11]) 
			except:
				# > 99,999 atoms in the PDB file.
				# After 99,999 atoms numbering is hexidecimal.
				##############################################
				self.atom_serial = int(line[6:11], 16)
				self.hex_atom_serial = 1

			self.atom_name = line[12:16].split()[0]
			
			# Hydrogens need four character spaces.
			# All other atoms only need at most three character spaces.
			# In the case of non-hydrogen atoms, add a space preceding the atom name.
			# This preserves the proper PDB format when and if the Atom() instance is written to file.
			##########################################################################################
			if len(self.atom_name) != 4:
				self.atom_name = " " + self.atom_name
					
			self.alternate_location = line[16:17]
			if self.alternate_location == " ":
				self.alternate_location = ""
			self.residue_name = line[17:20] 
			#Remove any white-space in the name of the residue.
			temp_residue_name = ""
			for char in self.residue_name:
				if char != " ":
					temp_residue_name += char
			self.residue_name = temp_residue_name 

			# Accept 1 and 2-character chains by argument.
			##############################################
			if twoCharacterChain:
				# 2-character chain.
				#############################
				self.chain_identifier = line[21:23]
			else:
				# 1-character chain.
				##########################
				self.chain_identifier = line[21:22]

			# Accept 1 and 2-character chains adaptively thru legalize.
			###########################################################
			#if legalize:
			#    # Illegal 1-character chain.
			#    #############################
			#    self.chain_identifier = line[21:22]
			#else:
			#    # Legal 2-character chain.
			#    ##########################
			#    self.chain_identifier = line[21:23]

			if self.chain_identifier == " ":
				self.chain_identifier = "NULL"
			else:
				# Remove whitespace if necessary.
				#################################
				self.chain_identifier = self.chain_identifier.split()[0]

			# Accepts 3 and 4-character residue numbers by argument.
			########################################################
			if twoCharacterChain:
				try:
					# int() with base 10.
					#####################
					self.residue_sequence_number = int(line[23:26])
				except:
					# residue is numbered using a hexidecimal system.
					#################################################
					self.residue_sequence_number = int(line[23:26], 16)
					self.hex_residue_sequence_number = 1
			else:
				try:
					self.residue_sequence_number = int(line[22:26])
				except:
					# > 999 residues in the PDB file.
					# After 999 residues numbering is hexidecimal.
					##############################################
					self.residue_sequence_number = int(line[22:26], 16)
					self.hex_residue_sequence_number = 1
					#print("Here.....Hex", line[22:26], self.residue_sequence_number, self.hex_residue_sequence_number)

			# Accepts 3 and 4-character residue numbers adaptively thru legalize.
			#####################################################################
			#if legalize:
			#    try:
			#        # int() with base 10.
			#        #####################
			#        self.residue_sequence_number = int(line[22:26])
			#    except:
			#        # residue is numbered using a hexidecimal system.
			#        #################################################
			#        self.residue_sequence_number = int(line[22:26], 16)
			#else:
			#    try:
			#        self.residue_sequence_number = int(line[23:26])
			#    except:
			#        # > 999 residues in the PDB file.
			#        # After 999 residues numbering is hexidecimal.
			#        ##############################################
			#        self.residue_sequence_number = int(line[23:26], 16)

			self.residue_insertion_code = line[26:27]
			self.x = float(line[30:38]) 
			self.y = float(line[38:46])
			self.z = float(line[46:54])
			if line[54:60]:
				self.occupancy = float(line[54:60])
			if line[60:66]:
				self.temperature_factor = float(line[60:66])
			# XXXX 2019.04.10 shifted the next three lines left...test
			self.segment_identifier = line[72:76] 
			self.symbol = line[76:78]
			try:
				self.charge = int(line[78:79])
			except:
				self.charge = 0

			self.residue_key = (self.residue_sequence_number, self.residue_name, self.chain_identifier)
			self.atom_key = (self.residue_sequence_number, self.residue_name, self.chain_identifier, self.atom_name)

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

		# Copy the current state of the instance.
		#########################################
		return Atom(str(self))

	def __repr__(self, forceHex=0):
				
		###############################################
		# If the conditional argument forceHex is true,
		# the string respresentation of the atom serial 
		# number and residue number will be written in 
		# hexidecimal format.
		###############################################

		if self.chain_identifier == "NULL":
			chain_identifier = " "
		else:
			chain_identifier = self.chain_identifier
				
		# Robust representation of atom_serial.
		#######################################
		stringAtomSerial = ""
		if self.hex_atom_serial and forceHex:
			stringAtomSerial = "%5s" % hex(self.atom_serial)[2:]
		else:
			stringAtomSerial = "%5i" % self.atom_serial
		#if self.atom_serial <= 99999:
		#    stringAtomSerial = "%5i" % self.atom_serial
		#else:
		#    stringAtomSerial = "%5s" % hex(self.atom_serial)[2:]

		# Robust representation of residue_sequence_number.
		####################################################   
		stringResidueNumber = ""
		if self.hex_residue_sequence_number and forceHex:
			stringResidueNumber = "%3s" % hex(self.residue_sequence_number)[2:]
		else:
			stringResidueNumber = "%3i" % self.residue_sequence_number
		#if self.residue_sequence_number <= 999:
		#    stringResidueNumber = "%03i" % self.residue_sequence_number
		#else:
		#    stringResidueNumber = "%3s" % hex(self.residue_sequence_number)[2:]

		# fixed
		# Format the atom data in standard PDB file format.
		###################################################
		pdbformat  = "%6s%5s %-4s%1s%3s %-1s%4s%1s"
		# Check for 2-character chain identifier.
		# If so, adjust chain and residue number fields accordingly.
		############################################################
		if len(chain_identifier) > 1:
			pdbformat  = "%6s%5s %-4s%1s%3s %-2s%3s%1s"
		pdbformat += "   %8.3f%8.3f%8.3f%6.2f%6.2f      %-4s%2s%2i\n"

		return pdbformat % (self.record_name,
							stringAtomSerial,
							self.atom_name,
							self.alternate_location,
							self.residue_name,
							chain_identifier,
							stringResidueNumber,
							self.residue_insertion_code,
							self.x,
							self.y,
							self.z,
							self.occupancy,
							self.temperature_factor,
							self.segment_identifier,
							self.symbol,
							self.charge)

class Atom_LO:

	def __init__(self):
		self.atom_serial = None
		# Key is number, name, chain.
		#############################
		self.residueKey = None
		self.repr = ""

	def __repr__(self):
		return self.repr

class PseudoAtom:

	def __init__(self):

		self.format = "pdb"
		self.pdbfileline = ""
		self.record_name = "ATOM  "
		self.atom_serial, self.hex_atom_serial = 1, 0
		self.atom_name = " O  " #one space to the left, two to the right
		self.alternate_location = ""
		self.residue_name = "PSA"
		self.chain_identifier = "A"
		self.residue_sequence_number, self.hex_residue_sequence_number = 1, 0
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

		# Robust representation of atom_serial.
		#######################################
		stringAtomSerial = ""
		if self.hex_atom_serial and forceHex:
			stringAtomSerial = "%5s" % hex(self.atom_serial)[2:]
		else:
			stringAtomSerial = "%5i" % self.atom_serial
		#if self.atom_serial <= 99999:
		#    stringAtomSerial = "%5i" % self.atom_serial
		#else:
		#    stringAtomSerial = "%5s" % hex(self.atom_serial)[2:]

		# Robust representation of residue_sequence_number.
		####################################################   
		stringResidueNumber = ""
		if self.hex_residue_sequence_number and forceHex:
			stringResidueNumber = "%3s" % hex(self.residue_sequence_number)[2:]
		else:
			stringResidueNumber = "%3i" % self.residue_sequence_number
		#if self.residue_sequence_number <= 999:
		#    stringResidueNumber = "%03i" % self.residue_sequence_number
		#else:
		#    stringResidueNumber = "%3s" % hex(self.residue_sequence_number)[2:]

		# fixed
		# Format the atom data in standard PDB file format.
		###################################################
		pdbformat  = "%6s%5s %-4s%1s%3s %-1s%4s%1s"
		# Check for 2-character chain identifier.
		# If so, adjust chain and residue number fields accordingly.
		############################################################
		if len(chain_identifier) > 1:
			pdbformat  = "%6s%5s %-4s%1s%3s %-2s%3s%1s"
		pdbformat += "   %8.3f%8.3f%8.3f%6.2f%6.2f      %-4s%2s%2i\n"

		return pdbformat % (self.record_name,
							stringAtomSerial,
							self.atom_name,
							self.alternate_location,
							self.residue_name,
							chain_identifier,
							stringResidueNumber,
							self.residue_insertion_code,
							self.x,
							self.y,
							self.z,
							self.occupancy,
							self.temperature_factor,
							self.segment_identifier,
							self.symbol,
							self.charge)
