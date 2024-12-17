# PDBCharges
## How to run PDBCharges in Docker container
    # build container
    docker build -t local/pdbcharges .

    # download the Chemical Component Dictionary
    curl -s https://files.wwpdb.org/pub/pdb/data/monomers/components-pub.sdf.gz | gunzip -c > components-pub.sdf

    # create folder to store results (the folder must be empty)
    mkdir results
    
    # run calculation on the 1alf.pdb file from the examples folder
    docker run --rm --name PDBcharges \
        -v ./components-pub.sdf:/opt/components-pub.sdf \
        -v ./examples:/opt/PDBCharges/examples \
        -v ./results:/opt/PDBCharges/results \
        local/pdbcharges \
        calculate_charges_workflow.py --CCD_file /opt/components-pub.sdf --PDB_file examples/1alf.pdb --data_dir results

    # or you can calculate the charges on another structure from PDB database, e.g.:
    curl -O https://files.rcsb.org/download/6wlv.pdb

    docker run --rm --name PDBcharges \
        -v ./components-pub.sdf:/opt/components-pub.sdf \
        -v ./6wlv.pdb:/opt/PDBCharges/6wlv.pdb \
        -v ./results:/opt/PDBCharges/results \
        local/pdbcharges \
        calculate_charges_workflow.py --CCD_file /opt/components-pub.sdf --PDB_file 6wlv.pdb --data_dir results
