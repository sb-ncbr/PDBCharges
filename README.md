# PDBCharges
## How to run PDBCharges in Docker container
    # build container
    docker build -t local/pdbcharges .

    # download the Chemical Component Dictionary
    curl -s https://files.wwpdb.org/pub/pdb/data/monomers/components-pub.sdf.gz | gunzip -c > components-pub.sdf

    # create folder to store results
    mkdir results
    
    # run calculation on the 1alf.pdb file from the examples folder
    docker run --rm --name PDBcharges \
        -v components-pub.sdf:/opt/components-pub.sdf \
        -v examples:/opt/PDBCharges/examples \
        -v results:/opt/PDBCharges/results \
        local/pdbcharges \
        calculate_charges_workflow.py --CCD_file /opt/components-pub.sdf --PDB_file examples/1alf.pdb --data_dir results

