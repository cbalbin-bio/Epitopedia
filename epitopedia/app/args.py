import argparse

parser = argparse.ArgumentParser(description="Epitopedia")
parser.add_argument(
    "--PDB-IDS",
    type=str,
    nargs="+",
    help="List of PDB_IDS formatted as PDBID_CHAIN. e.g. 6xr8_A",
)
parser.add_argument(
    "--span", type=int, default=5, help="Minimum span length to consider a mimic"
)
parser.add_argument(
    "--rasa", type=float, default=0.20, help="Relative accessible surface area cutoff"
)
parser.add_argument(
    "--rasa-span",
    type=int,
    default=3,
    help="Minimum span length for surface accessibility filter",
)
parser.add_argument(
    "--taxid-filter",
    type=str,
    nargs="+",
    help="Filter all taxaony at or below the level described by this taxid",
)
parser.add_argument(
    "--rmsd", type=float, help="Max RMSD to consider match a structural mimic"
)
parser.add_argument("--view", type=str, help="View results from a previous run")
parser.add_argument("--port", type=int, default=5000, help="Port used by webserver")
parser.add_argument(
    "--use-afdb", action="store_true", help="Include AFDB in database generation"
)
parser.add_argument(
    "--gplddt",
    type=float,
    default=0.70,
    help="Minimum global pLDDT score a structure predicted by alphafold must have to be considered",
)
parser.add_argument(
    "--lplddt",
    type=float,
    default=0.90,
    help="Minimum average local pLDDT score a region predicted by alphafold must have to be considered",
)
parser.add_argument("--headless", action="store_true", help="Do not start up webserver")


args = parser.parse_args()
