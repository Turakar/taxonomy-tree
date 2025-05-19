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

    def create_db(self) -> None:
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

            print("Done")

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
