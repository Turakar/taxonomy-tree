import itertools
import re
import subprocess
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterator

import duckdb
import polars as pl


@dataclass(frozen=True)
class TaxonomyEntry:
    identifier: str
    rank: str
    name: str | None = None

    def __str__(self) -> str:
        if self.name is None:
            return f"{self.rank}: {self.identifier}"
        else:
            return f"{self.rank}: {self.identifier} ({self.name})"

    def __repr__(self) -> str:
        return f'TaxonomyEntry("{str(self)}")'


class Taxonomy:
    def __init__(self, db_path: str | Path):
        self.db_path = Path(db_path)
        self._db_: duckdb.DuckDBPyConnection | None = None

    def create_db(self, add_gtdb_taxonomy: bool = False) -> None:
        print("Creating taxonomy database")
        # Create a temporary directory to store the downloaded files
        with (
            tempfile.TemporaryDirectory(prefix="taxonomy-") as tmpdirname,
            duckdb.connect(str(self.db_path)) as db,
        ):
            # Taxonomy files are available under (including a readme file):
            # https://ftp.ncbi.nih.gov/pub/taxonomy/taxdmp.zip: nodes.dmp and names.dmp
            # Assembly summary files are available under:
            # Four master files reporting data for either GenBank or RefSeq genome assemblies
            # are available under https://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/
            # assembly_summary_genbank.txt            - current GenBank genome assemblies
            # assembly_summary_genbank_historical.txt - replaced and suppressed GenBank genome
            #                                           assemblies
            # assembly_summary_refseq.txt             - current RefSeq genome assemblies
            # assembly_summary_refseq_historical.txt  - replaced and suppressed RefSeq genome
            #                                           assemblies

            # Download the files
            taxdmp_zip = str(Path(tmpdirname) / "taxdmp.zip")
            subprocess.run(
                ["curl", "-o", taxdmp_zip, "https://ftp.ncbi.nih.gov/pub/taxonomy/taxdmp.zip"],
                check=True,
            )
            subprocess.run(["unzip", "-o", taxdmp_zip], cwd=tmpdirname, check=True)

            # Load the names
            print("Loading names")
            names_df = (
                pl.scan_csv(
                    f"{tmpdirname}/names.dmp",
                    separator="|",
                    has_header=False,
                    new_columns=["tax_id", "name_txt", "unique_name", "name_class"],
                    quote_char=None,
                )
                .select(
                    pl.col("tax_id").cast(pl.Utf8).str.strip_chars("\t"),
                    pl.col("name_txt").cast(pl.Utf8).str.strip_chars("\t"),
                    pl.col("name_class").cast(pl.Utf8).str.strip_chars("\t"),
                )
                .filter(pl.col("name_class") == "scientific name")
                .select(
                    pl.col("tax_id").alias("taxid"), pl.col("name_txt").alias("scientific_name")
                )
                .collect()
            )
            db.execute("CREATE TABLE names AS SELECT * FROM names_df")
            db.execute("CREATE INDEX idx_names_taxid ON names (taxid)")
            print(f"Loaded {len(names_df)} names")

            # Load the nodes
            print("Loading nodes")
            nodes_df = (
                pl.scan_csv(  # noqa: F841
                    f"{tmpdirname}/nodes.dmp",
                    separator="|",
                    has_header=False,
                    new_columns=[
                        "tax_id",
                        "parent_tax_id",
                        "rank",
                        "embl_code",
                        "division_id",
                        "inherited_div_flag",
                        "genetic_code_id",
                        "inherited_GC_flag",
                        "mitochondrial_genetic_code_id",
                        "inherited_MGC_flag",
                        "GenBank_hidden_flag",
                        "hidden_subtree_root_flag",
                        "comments",
                    ],
                    quote_char=None,
                )
                .select(
                    pl.col("tax_id").cast(pl.Utf8).str.strip_chars("\t").alias("taxid"),
                    pl.col("parent_tax_id")
                    .cast(pl.Utf8)
                    .str.strip_chars("\t")
                    .alias("parent_taxid"),
                    pl.col("rank").cast(pl.Utf8).str.strip_chars("\t"),
                )
                .collect()
            )
            db.execute("CREATE TABLE nodes AS SELECT * FROM nodes_df")
            db.execute("CREATE INDEX idx_nodes_taxid ON nodes (taxid)")
            db.execute("CREATE INDEX idx_nodes_parent_taxid ON nodes (parent_taxid)")
            print(f"Loaded {len(nodes_df)} nodes")

            # Load the assembly summary
            print("Loading assembly summary")

            def scan_assembly_summary(url: str) -> pl.LazyFrame:
                filename = url.split("/")[-1]
                local_path = Path(tmpdirname) / filename
                subprocess.run(["curl", "-o", local_path, url], check=True)
                return pl.scan_csv(
                    local_path,
                    new_columns=[
                        "assembly_accession",
                        "bioproject",
                        "biosample",
                        "wgs_master",
                        "refseq_category",
                        "taxid",
                    ],
                    separator="\t",
                    has_header=False,
                    skip_lines=2,
                    quote_char=None,
                ).select(
                    pl.col("assembly_accession"),
                    pl.col("taxid").cast(pl.Utf8),
                )

            assemblies_df = pl.concat(  # noqa: F841
                (
                    scan_assembly_summary(
                        "https://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_genbank.txt"
                    ),
                    scan_assembly_summary(
                        "https://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_genbank_historical.txt"
                    ),
                    scan_assembly_summary(
                        "https://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_refseq.txt"
                    ),
                    scan_assembly_summary(
                        "https://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_refseq_historical.txt"
                    ),
                ),
                how="vertical",
            ).collect()
            db.execute(
                "CREATE TABLE assemblies AS SELECT * FROM assemblies_df",
            )
            db.execute("CREATE INDEX idx_assemblies_taxid ON assemblies (taxid)")
            db.execute("CREATE INDEX idx_assemblies_accession ON assemblies (assembly_accession)")
            print(f"Loaded {len(assemblies_df)} assemblies")

            # Optionally load prokaryotic taxonomy from the Genome Taxonomy Database (GTDB)
            if add_gtdb_taxonomy:
                self._add_gtdb_taxonomy(db, tmpdirname)

            print("Done")

    def _add_gtdb_taxonomy(
        self,
        db: duckdb.DuckDBPyConnection,
        tmpdirname: str,
        gtdb_release: str = "R10-RS226",
        db_prefix: str = "gtdb",
    ) -> None:
        """Add taxonomy for prokaryotes from the Genome Taxonomy Database (GTDB).
        GTDB taxonomy is based on genome trees inferred using FastTree from an aligned concatenated set of 120 single
        copy marker proteins for Bacteria, and with IQ-TREE from a concatenated set of 53 marker proteins for Archaea.

        Assembly accessions in GTDB are prefixed with "GB_GCA_" for GenBank and "RS_GCF_" for RefSeq. We keep this
        naming scheme to make them distinct from the NCBI assembly accessions in the assemblies table.

        The GTDB taxonomy is added to the nodes and names tables, using a prefix of "gtdb:" for the taxids. Taxids are
        taxonomic names in GTDB prefixed with a marker for the rank, so they are not numeric like NCBI taxids, e.g.
        "gtdb:d__Bacteria", "gtdb:c__Aenigmatarchaeia".

        :param db: The DuckDB connection to the database.
        :param tmpdirname: The temporary directory where the GTDB taxonomy files will be downloaded.
        :param gtdb_release: The GTDB release to use, e.g. "R10-RS226". Given versions are expected to follow the
        versioning scheme defined at https://gtdb.ecogenomic.org/faq#what-is-the-gtdb-versioning-scheme. The version
        can be found at the RELEASE_NOTES.txt of the GTDB release, e.g. at
        https://data.ace.uq.edu.au/public/gtdb/data/releases/release226/226.0/RELEASE_NOTES.txt.
        :param db_prefix: The prefix to use for the GTDB taxids in the database, by default "gtdb".
        :return:
        """
        # check version
        gtdb_version = re.match(r"^R(\d+)-RS(?P<refseq_release>\d+)$", gtdb_release)
        if gtdb_version is None:
            raise ValueError(
                f"GTDB release {gtdb_release} does not match expected versioning scheme defined at "
                f"https://gtdb.ecogenomic.org/faq#what-is-the-gtdb-versioning-scheme."
            )
        rs = gtdb_version.group(
            "refseq_release"
        )  # the RefSeq release is the relevant version number for the URL

        # create a root node for the GTDB taxonomy and anchor it to the root of the NCBI taxonomy
        gtdb_root = f"{db_prefix}:root"
        ncbi_tax_root_node = "1"  # NCBI taxonomy defines the root node with taxid 1
        db.execute(
            "INSERT INTO nodes (taxid, parent_taxid, rank) VALUES (?, ?, ?)",
            (gtdb_root, ncbi_tax_root_node, "no rank"),
        )
        db.execute(
            "INSERT INTO names (taxid, scientific_name) VALUES (?, ?)",
            (gtdb_root, gtdb_root),
        )

        # set name of GTDB root node
        db.execute(
            f"INSERT INTO names (taxid, scientific_name) VALUES ('{db_prefix}:root', '{db_prefix}:root')"
        )

        prefixes = {
            "Bacteria": "bac120",  # named after the 120 marker proteins for Bacteria
            "Archaea": "ar53",  # named after the 53 marker proteins for Archaea
        }
        for domain in [
            "Bacteria",
            "Archaea",
        ]:  # GTDB provides separate taxonomy files for Bacteria and Archaea
            # Download the GTDB taxonomy files
            prefix = prefixes[domain]
            url = f"https://data.ace.uq.edu.au/public/gtdb/data/releases/release{rs}/{rs}.0/{prefix}_taxonomy_r{rs}.tsv"
            gtdb_taxonomy_file = Path(tmpdirname) / f"gtdb_{domain}_taxonomy.tsv"
            subprocess.run(
                ["curl", "-o", gtdb_taxonomy_file, url],
                check=True,
            )

            # Load the GTDB taxonomy
            print(f"Loading GTDB taxonomy for {domain}")
            gtdb_df = (
                pl.scan_csv(
                    gtdb_taxonomy_file,
                    separator="\t",
                    has_header=False,
                    new_columns=["assembly_id", "lineage"],
                )
                .select(
                    pl.col("assembly_id").cast(pl.Utf8),
                    pl.col("lineage").cast(pl.Utf8),
                )
                .collect()
            )

            self._gtdb_add_nodes_and_names(
                domain=domain,
                gtdb_df=gtdb_df,
                db_prefix=db_prefix,
                db=db,
            )

            self._gtdb_add_assemblies(
                domain=domain,
                gtdb_df=gtdb_df,
                db_prefix=db_prefix,
                db=db,
            )

    @staticmethod
    def _gtdb_add_nodes_and_names(
        domain: str,
        gtdb_df: pl.DataFrame,
        db_prefix: str = "gtdb",
        db: duckdb.DuckDBPyConnection | None = None,
    ):
        # parse lineage to nodes
        ranks = ["domain", "phylum", "class", "order", "family", "group", "species"]
        tax_table = (
            gtdb_df.unique("lineage")
            .select(
                pl.col("lineage")
                .str.split(";")
                .list.to_struct(n_field_strategy="max_width", fields=ranks)
                .alias("taxonomy")
            )
            .unnest("taxonomy")
        )

        # update nodes table
        for parent_rank, rank in itertools.pairwise(["root", *ranks]):
            nodes_table_update = (  # noqa: F841
                tax_table.select(  # insert root column and map everything to root
                    pl.lit("root").alias("root"), pl.col("*")
                )
                .select(
                    pl.concat_str(
                        [
                            pl.lit(f"{db_prefix}:"),
                            pl.col(rank),
                        ]
                    ).alias("taxid"),
                    pl.concat_str(
                        [
                            pl.lit(f"{db_prefix}:"),
                            pl.col(parent_rank),
                        ]
                    ).alias("parent_taxid"),
                    pl.lit(rank).alias("rank"),
                )
                .unique()
            )

            db.execute(
                "INSERT INTO nodes (taxid, parent_taxid, rank) SELECT * FROM nodes_table_update"
            )

        # update names table
        names_table_update = (  # noqa: F841
            tax_table.unpivot()
            .select(
                pl.concat_str(
                    [
                        pl.lit(f"{db_prefix}:"),
                        pl.col("value"),
                    ]
                ).alias("taxid"),
                pl.col("value").str.slice(3).alias("scientific_name"),
            )
            .unique()
        )
        db.execute("INSERT INTO names (taxid, scientific_name) SELECT * FROM names_table_update")

        print(f"Loaded taxonomic tree nodes from {len(gtdb_df)} {domain} entries")

    @staticmethod
    def _gtdb_add_assemblies(
        domain: str,
        gtdb_df: pl.DataFrame,
        db_prefix: str = "gtdb",
        db: duckdb.DuckDBPyConnection | None = None,
    ):
        print(f"Add assemblies from GTDB for domain {domain}")

        assemblies_with_taxa = gtdb_df.select(  # noqa: F841
            # GTDB assembly accessions are prefixed (with RS_ or GB_ depending on source (RefSeq/GenBank)) NCBI assembly
            # accessions. By slicing the first three characters, we get the NCBI assembly accession.
            pl.col("assembly_id").str.slice(3).alias("assembly_id"),
            pl.concat_str(
                [
                    pl.lit(f"{db_prefix}:"),
                    pl.col("lineage").str.split(";").list.last(),
                ]
            ).alias("species"),
        )

        # upsert assemblies in GTDB to the assemblies table
        # overwrite existing entries with the same assembly accession
        db.execute("""
            BEGIN TRANSACTION;
            DELETE FROM assemblies WHERE assembly_accession IN (SELECT assembly_id FROM assemblies_with_taxa);
            INSERT INTO assemblies (assembly_accession, taxid) SELECT * FROM assemblies_with_taxa;
            COMMIT;
        """)
        print(f"Loaded assemblies for {domain} with {len(gtdb_df)} entries")

    @property
    def _db(self) -> duckdb.DuckDBPyConnection:
        if self._db_ is None:
            self._db_ = duckdb.connect(str(self.db_path), read_only=True)
        return self._db_

    def contains(self, identifier: str) -> bool:
        if self.get_identifier_type(identifier) == "assembly_accession":
            assembly_fetchone = self._db.execute(
                "SELECT taxid FROM assemblies WHERE assembly_accession = ?",
                (identifier,),
            ).fetchone()
            return assembly_fetchone is not None
        else:
            taxid_fetchone = self._db.execute(
                "SELECT taxid FROM nodes WHERE taxid = ?",
                (identifier,),
            ).fetchone()
            return taxid_fetchone is not None

    def get_identifier_type(self, identifier: str) -> str:
        if identifier.isdigit():
            return "taxid"
        if identifier.startswith("gtdb:"):  # GTDB identifiers are prefixed with "gtdb:"
            return "taxid"
        if identifier.startswith("GCA_") or identifier.startswith("GCF_"):
            return "assembly_accession"
        raise ValueError(f"Identifier {identifier} is not a valid taxid or assembly accession.")

    def get_assembly_taxid(self, identifier: str) -> str:
        if self.get_identifier_type(identifier) != "assembly_accession":
            raise ValueError(f"Identifier {identifier} is not a valid assembly accession.")
        taxid_fetchone = self._db.execute(
            "SELECT taxid FROM assemblies WHERE assembly_accession = ?",
            (identifier,),
        ).fetchone()
        if taxid_fetchone is None:
            raise ValueError(f"Assembly accession {identifier} not found in the database.")
        return taxid_fetchone[0]

    def find_scientific_name(self, identifier: str) -> str:
        if self.get_identifier_type(identifier) == "assembly_accession":
            identifier = self.get_assembly_taxid(identifier)

        # Now we can get the scientific name
        scientific_name_fetchone = self._db.execute(
            "SELECT scientific_name FROM names WHERE taxid = ?",
            (identifier,),
        ).fetchone()
        if scientific_name_fetchone is None:
            raise ValueError(f"Taxid {identifier} not found in the database.")
        return scientific_name_fetchone[0]

    def find_lineage(
        self,
        identifier: str,
        stop_rank: str | None = None,
        with_names: bool = False,
        skip: set[str] | None = None,
    ) -> Iterator[TaxonomyEntry]:
        if skip is None:
            skip = set()
        if identifier in skip:
            return

        if self.get_identifier_type(identifier) == "assembly_accession":
            # If the identifier is an assembly accession, we need to get the taxid first
            taxid = self.get_assembly_taxid(identifier)
            yield TaxonomyEntry(
                identifier=identifier,
                rank="assembly",
            )
            if taxid in skip:
                return
        else:
            taxid = identifier

        # Fetch the node
        node_fetchone = self._db.execute(
            "SELECT parent_taxid, rank FROM nodes WHERE taxid = ?",
            (taxid,),
        ).fetchone()
        if node_fetchone is None:
            raise ValueError(f"Taxid {taxid} not found in the database.")
        parent_taxid, rank = node_fetchone
        if with_names:
            name = self.find_scientific_name(taxid)
        else:
            name = None
        yield TaxonomyEntry(
            identifier=taxid,
            rank=rank,
            name=name,
        )

        # Stop criteria
        if (stop_rank is not None and rank == stop_rank) or parent_taxid == taxid:
            return

        # Recursively fetch the parent node
        for parent_entry in self.find_lineage(
            parent_taxid,
            stop_rank=stop_rank,
            with_names=with_names,
            skip=skip | {identifier},
        ):
            yield parent_entry

    def find_children(
        self,
        identifier: str,
        stop_rank: str | None = None,
        with_names: bool = False,
        skip: set[str] | None = None,
    ) -> Iterator[TaxonomyEntry]:
        if skip is None:
            skip = set()

        if self.get_identifier_type(identifier) == "assembly_accession":
            # Assemblies do not have children
            return

        # Get children
        children_fetchall = self._db.execute(
            "SELECT taxid, rank FROM nodes WHERE parent_taxid = ?",
            (identifier,),
        ).fetchall()
        for child_taxid, child_rank in children_fetchall:
            if child_taxid in skip:
                continue
            if with_names:
                child_name = self.find_scientific_name(child_taxid)
            else:
                child_name = None
            yield TaxonomyEntry(
                identifier=child_taxid,
                rank=child_rank,
                name=child_name,
            )
            # Stop criteria
            if stop_rank is not None and child_rank == stop_rank:
                continue
            # Recursively fetch the children
            child_children = self.find_children(
                child_taxid,
                stop_rank=stop_rank,
                with_names=with_names,
                skip=skip,
            )
            for child_child in child_children:
                yield child_child
        assemblies_fetchall = self._db.execute(
            "SELECT assembly_accession FROM assemblies WHERE taxid = ?",
            (identifier,),
        ).fetchall()
        for (assembly_accession,) in assemblies_fetchall:
            if assembly_accession in skip:
                continue
            yield TaxonomyEntry(
                identifier=assembly_accession,
                rank="assembly",
            )

    def depth_first_search(
        self, identifier: str, with_names: bool = False
    ) -> Iterator[TaxonomyEntry]:
        lineage = self.find_lineage(identifier, with_names=with_names)
        visited = set()
        for entry in lineage:
            if entry.identifier in visited:
                continue
            yield entry
            visited.add(entry.identifier)
            children = self.find_children(
                entry.identifier, with_names=with_names, skip=visited.copy()
            )
            for child in children:
                if child.identifier in visited:
                    continue
                visited.add(child.identifier)
                yield child

    def find_nearest_neighbor(self, identifier: str, targets: set[str]) -> str | None:
        for entry in self.depth_first_search(identifier):
            if entry.identifier in targets:
                return entry.identifier
        return None

    def __getstate__(self) -> Any:
        return {"db_path": str(self.db_path)}

    def __setstate__(self, state: Any) -> None:
        self.db_path = Path(state["db_path"])
        self._db_ = None


if __name__ == "__main__":
    import pickle

    taxonomy = Taxonomy("taxonomy.db")
    if not taxonomy.db_path.exists():
        taxonomy.create_db()

    print(f"9606: {taxonomy.find_scientific_name('9606')}")

    taxonomy = pickle.loads(pickle.dumps(taxonomy))

    print(f"GCA_000001405.29: {taxonomy.find_scientific_name('GCA_000001405.29')}")
    print(
        f"and its lineage: {list(taxonomy.find_lineage('GCA_000001405.29', with_names=True, stop_rank='class'))}"
    )
    print(f"Children of 4897: {list(taxonomy.find_children('4897', with_names=True))}")
    children_314146 = list(taxonomy.find_children("314146", stop_rank="species"))
    print("10090 child of 314146:", "10090" in {c.identifier for c in children_314146})
    print(
        f"Nearest neighbor of 9606 in (4897, 10090): {taxonomy.find_nearest_neighbor('9606', {'4897', '10090'})}"
    )
