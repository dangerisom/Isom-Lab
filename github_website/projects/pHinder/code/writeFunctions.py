# Copyright (c) 2025 Isom Lab
# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.
# If a copy of the MPL was not distributed with this file, You can obtain one at
# https://mozilla.org/MPL/2.0/.

import gzip
import pdbFile, pdbFile_cif

# Have to write PDB formatted files instead of mmCIF formatted files,
# because mmCIF doesn't allow me to output CONECT data for the network edges.

def writeNetworkToFile(networks, triangulation, networkName, outputDirectory, pdbCode, chain=None, include_residues=1, zipit=0):

    # Rank order networks by size...
    ###################################################################################
    network_keys = list(networks)
    ranked_network_keys = []
    for network_key in network_keys:
        ranked_network_keys.append((len(network_key), network_key))
    ranked_network_keys = sorted(ranked_network_keys)
    ranked_network_keys.reverse() 

    # Write PDB and CONECT information order by decreasing network sizes...
    ###################################################################################
    residue_networks, network_connect = {}, {}
    for ranked_network_key in ranked_network_keys:
        network = networks[ranked_network_key[1]]
        residue_network_nodes, connect = {}, {}
        for node in sorted(network):
            # Need to include atom_serial in case residue key is the same from different PDB contributions...
            node_residue_key = tuple([triangulation[node].s1.data.atom_serial] + list(triangulation[node].s1.data.residue.key))
            node_residue = triangulation[node].s1.data.residue
            residue_network_nodes[node_residue_key] = node_residue
            connect[node_residue_key] = {node_residue:None}
            for s2 in triangulation[node].s2s:
                for neighbor_node in s2:
                    # Need to include atom_serial in case residue key is the same from different PDB contributions...
                    neighbor_residue_key = tuple([triangulation[neighbor_node].s1.data.atom_serial] + list(triangulation[neighbor_node].s1.data.residue.key))
                    neighbor_residue = triangulation[neighbor_node].s1.data.residue
                    residue_network_nodes[neighbor_residue_key] = neighbor_residue
                    connect[node_residue_key][neighbor_residue] = None
        network_connect[ranked_network_key[1]] = connect
        residue_networks[ranked_network_key[1]] = residue_network_nodes

    # Write full amino acid side chain PDB version with network and CONECT information...
    ####################################################################################
    #chain_string_i, chain_string = 0, 4*"ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789"
    
    # Write ATOM information
    if include_residues:
        if zipit:
            if chain:
                outFile = outputDirectory + pdbCode + "-" + networkName + "ResidueNetwork-" + chain + ".pdb.gz"
                outputTopoPDB = gzip.open(outFile,"wt")
            else:
                outFile = outputDirectory + pdbCode + "-" + networkName + "ResidueNetwork.pdb.gz"
                outputTopoPDB = gzip.open(outFile,"wt")
        else:
            if chain:
                outFile = outputDirectory + pdbCode + "-" + networkName + "ResidueNetwork-" + chain + ".pdb"
                outputTopoPDB = open(outFile,"w")
            else:
                outFile = outputDirectory + pdbCode + "-" + networkName + "ResidueNetwork.pdb"
                outputTopoPDB = open(outFile,"w")

        network_number, atoms_string = 1, ""
        for residue_network in residue_networks:
            keys = sorted(residue_networks[residue_network])
            for key in keys:
                residue = residue_networks[residue_network][key]
                for atom_key in residue.atoms:
                    atom = residue.atoms[atom_key]
                    if atom.format == "pdb":
                        atom = pdbFile.Atom(str(atom))
                        atom.segment_identifier = network_number
                        atoms_string += str(atom)
                    else:
                        atom = pdbFile_cif.Atom(str(atom), atom.cif_options)
                        atom.segment_identifier = network_number
                        atoms_string += atom.get_pdb_format()
            network_number += 1
        outputTopoPDB.write(atoms_string)

    # Write CONECT information
    conectString, network_number = "", 1
    for network in network_connect:
        connect = network_connect[network]
        conectString += "\nREMARK BEGIN NETWORK " + str(network_number) + " : " + str(len(connect)) + " nodes\n"
        for network_node in sorted(connect):
            conectString += "%-6s" % "CONECT"
            network_nodes = connect[network_node]
            for node in network_nodes:
                tsc_atom = node.get_terminal_sidechain_atom()
                if tsc_atom:
                    atom_serial = tsc_atom.atom_serial - 1 # The -1 corrects for how the tscAtom is instantiated... 
                if not tsc_atom:
                    atom_key = list(node.atoms)[0]
                    tsc_atom = node.atoms[atom_key]
                    atom_serial = tsc_atom.atom_serial # tscAtom is PSA...
                conectString += "%5i" % atom_serial
            conectString += "\n" 
        conectString += "REMARK END NETWORK " + str(network_number) + " : " + str(len(connect)) + " nodes\n"
        # chain_string_i += 1
        network_number += 1
    if include_residues:
        outputTopoPDB.write(conectString)
        outputTopoPDB.close()

    # Write terminal side chain atom PDB version with network and CONECT information...
    ###################################################################################
    if zipit:
        if chain:
            outFile = outputDirectory + pdbCode + "-" + networkName + "AtomNetwork-" + chain + ".pdb.gz"
            outputTopoPDB = gzip.open(outFile,"wt")
        else:
            outFile = outputDirectory + pdbCode + "-" + networkName + "AtomNetwork.pdb.gz"
            outputTopoPDB = gzip.open(outFile,"wt")
    else:
        if chain:
            outFile = outputDirectory + pdbCode + "-" + networkName + "AtomNetwork-" + chain + ".pdb"
            outputTopoPDB = open(outFile,"w")
        else:
            outFile = outputDirectory + pdbCode + "-" + networkName + "AtomNetwork.pdb"
            outputTopoPDB = open(outFile,"w")

    # chain_string_i, chain_string = 0, 4*"ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789"
    network_number, atoms_string = 1, ""
    for residue_network in residue_networks:
        keys = sorted(residue_networks[residue_network])
        for key in keys:
            residue = residue_networks[residue_network][key]
            tsc_atom = residue.get_terminal_sidechain_atom()
            if tsc_atom:
                tsc_atom.atom_serial = tsc_atom.atom_serial - 1 # The -1 corrects for how the tscAtom is instantiated... 
            if not tsc_atom:
                atom_key = list(residue.atoms)[0]
                tsc_atom = residue.atoms[atom_key]

            if tsc_atom.format == "pdb":
                tsc_atom = pdbFile.Atom(str(tsc_atom))
                tsc_atom.segment_identifier = network_number
                atoms_string += str(tsc_atom)
            else:
                tsc_atom = pdbFile_cif.Atom(str(tsc_atom), tsc_atom.cif_options)
                tsc_atom.segment_identifier = network_number
                atoms_string += tsc_atom.get_pdb_format()
        network_number += 1

    outputTopoPDB.write(atoms_string)
    outputTopoPDB.write(conectString)
    outputTopoPDB.close()

