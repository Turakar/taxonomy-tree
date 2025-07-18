import time

from taxonomy_tree import Taxonomy

print("Creating taxonomy database for NCBI...")
taxonomy = Taxonomy("taxonomy.db")
prior = time.time()
taxonomy.create_db()
print(f"NCBI: {time.time() - prior} s")

del taxonomy

print("Creating taxonomy database for NCBI+GTDB...")
taxonomy = Taxonomy("taxonomy_ncbi_and_gtdb.db")
prior = time.time()
taxonomy.create_db(add_gtdb_taxonomy=True)
print(f"NCBI+GTDB: {time.time() - prior} s")
