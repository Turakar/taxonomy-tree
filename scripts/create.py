import time

from taxonomy_tree.build import make_taxonomy_tree

print("Creating taxonomy database for NCBI...")
prior = time.time()
make_taxonomy_tree("taxonomy_ncbi.db")
print(f"NCBI: {time.time() - prior} s")

print("Creating taxonomy database for NCBI+GTDB...")
prior = time.time()
make_taxonomy_tree("taxonomy_ncbi_and_gtdb.db", include_gtdb=True)
print(f"NCBI+GTDB: {time.time() - prior} s")
