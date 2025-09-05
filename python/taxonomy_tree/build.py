import itertools
import re
import subprocess
import tempfile
from pathlib import Path

import polars as pl

from . import taxonomy_tree as _rust


def make_taxonomy_tree(output_path: str | Path, include_gtdb: bool = False) -> None:
    with tempfile.TemporaryDirectory(prefix="taxonomy-") as tmpdir:
        builder = _rust.TaxonomyTreeBuilder()
        _load_ncbi(tmpdir, builder)
        if include_gtdb:
            _load_gtdb(tmpdir, builder)
        print("Finalizing and writing taxonomy tree")
        builder.write(str(output_path))
        print(f"Wrote taxonomy tree to {output_path}")


def _load_ncbi(tmpdir: str, builder: _rust.TaxonomyTreeBuilder) -> None:
    print("Creating taxonomy database")
    # Create a temporary directory to store the downloaded files
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
    taxdmp_zip = str(Path(tmpdir) / "taxdmp.zip")
    subprocess.run(
        ["curl", "-o", taxdmp_zip, "https://ftp.ncbi.nih.gov/pub/taxonomy/taxdmp.zip"],
        check=True,
    )
    subprocess.run(["unzip", "-o", taxdmp_zip], cwd=tmpdir, check=True)

    # Load taxonomy nodes (excluding assemblies) with names
    print("Loading taxonomy nodes")
    names_df = (
        pl.scan_csv(
            f"{tmpdir}/names.dmp",
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
        .select(pl.col("tax_id").alias("taxid"), pl.col("name_txt").alias("scientific_name"))
    )
    nodes_df = (
        pl.scan_csv(  # noqa: F841
            f"{tmpdir}/nodes.dmp",
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
            pl.col("parent_tax_id").cast(pl.Utf8).str.strip_chars("\t").alias("parent_taxid"),
            pl.col("rank").cast(pl.Utf8).str.strip_chars("\t"),
        )
        .join(
            names_df,
            left_on="taxid",
            right_on="taxid",
            how="left",
        )
        .collect()
    )
    for row in nodes_df.iter_rows(named=True):
        builder.insert(
            identifier=row["taxid"],
            name=row["scientific_name"],
            rank=row["rank"],
            # root node points to itself as parent
            parent=row["parent_taxid"] if row["parent_taxid"] != row["taxid"] else None,
            overwrite=False,
        )
    print(f"Loaded {len(nodes_df)} internal nodes")

    # Load assemblies and link to taxids
    print("Loading assembly summary")

    def scan_assembly_summary(url: str) -> pl.LazyFrame:
        filename = url.split("/")[-1]
        local_path = Path(tmpdir) / filename
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

    assemblies_df = (
        pl.concat(  # noqa: F841
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
        )
        .unique("assembly_accession")
        .filter(pl.col("taxid").is_in(nodes_df["taxid"]))
        .collect()
    )
    for row in assemblies_df.iter_rows(named=True):
        builder.insert(
            identifier=row["assembly_accession"],
            name=row["assembly_accession"],
            rank="assembly",
            parent=row["taxid"],
            overwrite=False,
        )
    print(f"Loaded {len(assemblies_df)} assemblies")


def _load_gtdb(
    tmpdir: str,
    builder: _rust.TaxonomyTreeBuilder,
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
    builder.insert(gtdb_root, gtdb_root, "cellular root", ncbi_tax_root_node, overwrite=False)

    prefixes = {
        "Bacteria": "bac120",  # named after the 120 marker proteins for Bacteria
        "Archaea": "ar53",  # named after the 53 marker proteins for Archaea
    }
    for domain in [
        "Bacteria",
        "Archaea",
    ]:  # GTDB provides separate taxonomy files for Bacteria and Archaea
        print(f"Loading GTDB taxonomy for {domain}")
        # Download the GTDB taxonomy files
        prefix = prefixes[domain]
        url = f"https://data.ace.uq.edu.au/public/gtdb/data/releases/release{rs}/{rs}.0/{prefix}_taxonomy_r{rs}.tsv"
        gtdb_taxonomy_file = Path(tmpdir) / f"gtdb_{domain}_taxonomy.tsv"
        subprocess.run(
            ["curl", "-o", gtdb_taxonomy_file, url],
            check=True,
        )

        # Load the GTDB taxonomy
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
        _gtdb_add_internal_nodes(
            builder=builder,
            domain=domain,
            gtdb_df=gtdb_df,
            db_prefix=db_prefix,
        )
        _gtdb_add_assemblies(
            builder=builder,
            domain=domain,
            gtdb_df=gtdb_df,
            db_prefix=db_prefix,
        )


def _gtdb_add_internal_nodes(
    builder: _rust.TaxonomyTreeBuilder,
    domain: str,
    gtdb_df: pl.DataFrame,
    db_prefix: str,
):
    # Parse lineage to nodes
    ranks = ["domain", "phylum", "class", "order", "family", "group", "species"]
    tax_table = (
        gtdb_df.lazy()
        .unique("lineage")
        .select(
            pl.col("lineage")
            .str.split(";")
            .list.to_struct(n_field_strategy="max_width", fields=ranks)
            .alias("taxonomy")
        )
        .unnest("taxonomy")
        .collect()
    )

    # Update nodes table
    for parent_rank, rank in itertools.pairwise(["root", *ranks]):
        nodes_table_update = (  # noqa: F841
            tax_table.lazy()
            .with_columns(  # insert root column and map everything to root
                root=pl.lit("root")
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
                pl.col(rank).str.slice(3).alias("name"),  # remove rank prefix from name
            )
            .unique()
            .collect()
        )
        for row in nodes_table_update.iter_rows(named=True):
            builder.insert(
                identifier=row["taxid"],
                name=row["name"],
                rank=row["rank"],
                parent=row["parent_taxid"],
                overwrite=False,
            )

    print(f"Loaded taxonomic tree nodes from {len(gtdb_df)} {domain} entries")


@staticmethod
def _gtdb_add_assemblies(
    builder: _rust.TaxonomyTreeBuilder,
    domain: str,
    gtdb_df: pl.DataFrame,
    db_prefix: str,
):
    print(f"Add assemblies from GTDB for domain {domain}")

    assemblies_with_taxa = (
        gtdb_df.lazy()
        .select(  # noqa: F841
            # GTDB assembly accessions are prefixed (with RS_ or GB_ depending on source (RefSeq/GenBank)) NCBI assembly
            # accessions. By slicing the first three characters, we get the NCBI assembly accession.
            pl.col("assembly_id").str.slice(3),
            pl.concat_str(
                [
                    pl.lit(f"{db_prefix}:"),
                    pl.col("lineage").str.split(";").list.last(),
                ]
            ).alias("species"),
        )
        .unique("assembly_id")
        .collect()
    )

    # upsert assemblies in GTDB to the assemblies table
    # overwrite existing entries with the same assembly accession (from NCBI)
    for row in assemblies_with_taxa.iter_rows(named=True):
        builder.insert(
            identifier=row["assembly_id"],
            name=row["assembly_id"],
            rank="assembly",
            parent=row["species"],
            overwrite=True,
        )
    print(f"Loaded assemblies for {domain} with {len(assemblies_with_taxa)} entries")
