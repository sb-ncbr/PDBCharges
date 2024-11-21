from biotite.structure.io import pdbx as biotite_mmCIF


def check_num_of_chains(mmCIF_file: str) -> bool:
    """
    The maximum number of chains that moleculekit can work with is 52.
    The function returns True if the number of chains is less than 52.
    Otherwise, it returns False.

    :param mmCIF_file: mmCIF file with the structure that should be checked for the number of chains
    :return: True or False
    """
    biotite_mmCIF_file = biotite_mmCIF.CIFFile.read(mmCIF_file)
    num_of_unique_label_asym_ids = len(set(biotite_mmCIF_file.block["atom_site"]["label_asym_id"].as_array()))
    num_of_unique_auth_asym_ids = len(set(biotite_mmCIF_file.block["atom_site"]["auth_asym_id"].as_array()))
    if num_of_unique_label_asym_ids > 52 or num_of_unique_auth_asym_ids > 52:
        return False
    return True

def change_chain_names(mmCIF_file: str,
                       output_mmCIF_file: str,
                       label_asym_ids_map: dict = None,
                       auth_asym_ids_map: dict = None) -> (dict, dict):
    """
    This function is used to change the names of chains in an mmCIF file to single letter names.
    The function also returns two dictionaries that can be used as arguments
    label_asym_ids_map and auth_sym_ids_map in a second call of this function
    to rename the chains back to their original names.

    :param mmCIF_file: mmCIF file with the structure in which the chain names should be changed
    :param output_mmCIF_file: mmCIF file to store the structure which changed chain names
    :param label_asym_ids_map: A dictionary defining how label_asym_id values should be renamed.
    :param auth_asym_ids_map: A dictionary defining how auth_asym_id values should be renamed.
    :return: two dictionaries that can be used to return chain names to their original state
    """

    biotite_mmCIF_file = biotite_mmCIF.CIFFile.read(mmCIF_file)
    label_asym_ids = biotite_mmCIF_file.block["atom_site"]["label_asym_id"].as_array()
    auth_asym_ids = biotite_mmCIF_file.block["atom_site"]["auth_asym_id"].as_array()

    if label_asym_ids_map is None and auth_asym_ids_map is None:
        allowed_chain_names = "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ"
        label_asym_ids_map = {key: allowed_chain_names[index] for index, key in enumerate(sorted(set(label_asym_ids)))}
        auth_asym_ids_map = {key: allowed_chain_names[index] for index, key in enumerate(sorted(set(auth_asym_ids)))}

    biotite_mmCIF_file.block["atom_site"]["label_asym_id"] = [label_asym_ids_map[key] for key in label_asym_ids]
    try:
        biotite_mmCIF_file.block["struct_conn"]["ptnr1_label_asym_id"] = [label_asym_ids_map[key] for key in biotite_mmCIF_file.block["struct_conn"]["ptnr1_label_asym_id"].as_array()]
        biotite_mmCIF_file.block["struct_conn"]["ptnr2_label_asym_id"] = [label_asym_ids_map[key] for key in biotite_mmCIF_file.block["struct_conn"]["ptnr2_label_asym_id"].as_array()]
    except KeyError:
        pass


    biotite_mmCIF_file.block["atom_site"]["auth_asym_id"] = [auth_asym_ids_map[key] for key in auth_asym_ids]
    try:

        biotite_mmCIF_file.block["struct_conn"]["ptnr1_auth_asym_id"] = [auth_asym_ids_map[key] for key in biotite_mmCIF_file.block["struct_conn"]["ptnr1_auth_asym_id"].as_array()]
        biotite_mmCIF_file.block["struct_conn"]["ptnr2_auth_asym_id"] = [auth_asym_ids_map[key] for key in biotite_mmCIF_file.block["struct_conn"]["ptnr2_auth_asym_id"].as_array()]
    except KeyError:
        pass



    biotite_mmCIF_file.write(output_mmCIF_file)

    # these three lines can be probably removed, after biotite 1.0.2 will be released
    mmcif_string = open(output_mmCIF_file).read()
    repaired_mmcif_string = mmcif_string.replace("\n# ", "\n# \n")
    open(output_mmCIF_file, "w").write(repaired_mmcif_string)

    return {value:key for key, value in label_asym_ids_map.items()}, {value:key for key, value in auth_asym_ids_map.items()}
    


