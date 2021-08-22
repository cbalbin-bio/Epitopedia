# Epitopedia

## Getting started

The quickest way to start using Epitopedia is by downloading the docker container which contains all the dependencies preinstalled

```shell
docker pull epitopedia
```

Epitopedia requires the [PDB in mmCIF](https://www.wwpdb.org/ftp/pdb-ftp-sites) format, Epitopedia DB and EPI-SEQ DB. Epitopedia DB and EPI-SEQ DB can be downloaded [here](https://fiudit-my.sharepoint.com/:u:/g/personal/cbalbin_fiu_edu/EWW8XKxSx09CvWC2mhzp8_sBntrnXX9mju4SbA0_ygUFMA?e=zKcZvU)


To download the PDB DB
```shell
rsync -rlpt -v -z --delete --port=33444 \
rsync.rcsb.org::ftp_data/structures/divided/mmCIF/ ./mmCIF
```

OR
<br>

To download the only the PDB files present in Epitopedia DB(EPI-PDB) you can supply the pdb_id_list.txt to resync
```shell
rsync -rlpt -v -z --delete --port=33444 --include-from=/path/to/pdb_id_list.txt rsync.rcsb.org::ftp_data/structures/divided/mmCIF/ ./mmCIF
```


To run Epitopedia provide the paths to to the various directories as shown below.

The data directory should contain Epitopedia DB (epitopedia.sqlite3) and EPI-SEQ (EPI-SEQ.fasta*) which can be downloaded [here](https://fiudit-my.sharepoint.com/:u:/g/personal/cbalbin_fiu_edu/EWW8XKxSx09CvWC2mhzp8_sBntrnXX9mju4SbA0_ygUFMA?e=zKcZvU).

The mmcif directory should point to the sharded PDB directory in mmCIF format as downloaded above.

The output directory is where the output files will be written

```shell
docker run --rm -it -p 5000:5000 \
-v /Path/to/Data/Dir/:/app/data \
-v /Path/to/PDBDB/Dir/:/app/mmcif \
-v /Path/to/Output/Dir/:/app/output \
epitopedia run_epitopedia.py 6VXX_A
```

Epitopedia will output the following files

File Name | Description
------------ | -------------
Content from cell 1 | Content from cell 2
Content in the first column | Content in the second column




Epitopedia can run on multiple input structures to represent a conformational ensemble. To do so, simply provide a list of structures in the format PDBID_CHAINID as shown below.
```shell
run_epitopedia.py 6VXX_A 6XR8_A
```

Epitopedia defaults to a span length of 5, surface accesbility cutoff of 20% surface accesbility span legnth of 3, and no taxa filter, but these parameters can be set using the follow flags.

Flag | Description
------------ | -------------
--span | Minimum span length for a hit to progress
--rasa | Cutoff for relative accessible surface area
--rasa_span | Minimum consecutive accesssible residues to consider a hit a SeqBMM
--taxid_filter | Taxid filter; example to filter out all Coronaviridae --taxid_filter 11118


## Epitopedia database generation

Epitopedia uses IEDB and PDB to generate Epitopedia DB, which is used in the molecular mimicry search.

Generation of the database takes some time (~12 hours). Thus a pregenerated database is provided above.

To create the database, download [IEDB](https://www.iedb.org/downloader.php?file_name=doc/iedb_public.sql.gz) and a mmCIF version of PDB

Point the container to the approriate path's for the IEDB, PDB DB (mmCIF format) and a data directory where the databases will be written to


```shell
docker run --rm -it \
-v /Path/To/iedb_public.sql:/app/iedb \
-v /Path/to/mmCIF/Dir/:/app/mmcif \
-v /Path/to/Data/Dir/:/app/data \
mimicrypipeline generate_database.py
```

