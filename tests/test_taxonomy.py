import pickle
from pathlib import Path

import pytest
from taxonomy_tree import Taxonomy, TaxonomyEntry

# https://goat.genomehubs.org
HUMAN_ASSEMBLY = "GCA_000001405.29"
HUMAN_NAME = "Homo sapiens"
HUMAN_TAXID = "9606"
HOMINIDAE_TAXID = "9604"


@pytest.fixture()
def taxonomy() -> Taxonomy:
    path = Path("taxonomy_ncbi.db")
    taxonomy = Taxonomy(path)
    if not path.exists():
        print("Creating taxonomy database")
        taxonomy.create_db()
    return taxonomy


@pytest.fixture()
def taxonomy_ncbi_and_gtdb() -> Taxonomy:
    path = Path("taxonomy_ncbi_and_gtdb.db")
    taxonomy = Taxonomy(path)
    if not path.exists():
        print("Creating taxonomy database")
        taxonomy.create_db(include_gtdb=True)
    return taxonomy


def test_pickle(taxonomy: Taxonomy) -> None:
    assert taxonomy.contains(HUMAN_ASSEMBLY)
    taxonomy2: Taxonomy = pickle.loads(pickle.dumps(taxonomy))
    assert taxonomy2.contains(HUMAN_ASSEMBLY)


def test_scientific_name(taxonomy: Taxonomy) -> None:
    scientific_name = taxonomy.find_scientific_name(HUMAN_TAXID)
    assert scientific_name == HUMAN_NAME


def test_scientific_name_gtdb(taxonomy_ncbi_and_gtdb: Taxonomy) -> None:
    scientific_name = taxonomy_ncbi_and_gtdb.find_scientific_name("gtdb:s__Escherichia coli")
    assert scientific_name == "Escherichia coli"
    scientific_name = taxonomy_ncbi_and_gtdb.find_scientific_name("gtdb:s__Palsa-295 sp036518755")
    assert scientific_name == "Palsa-295 sp036518755"
    scientific_name = taxonomy_ncbi_and_gtdb.find_scientific_name("gtdb:s__Methanocatella smithii")
    assert scientific_name == "Methanocatella smithii"
    scientific_name = taxonomy_ncbi_and_gtdb.find_scientific_name("gtdb:s__MGIIb-O2 sp002686525")
    assert scientific_name == "MGIIb-O2 sp002686525"


def test_search_scientific_name(taxonomy: Taxonomy) -> None:
    results = list(taxonomy.search_scientific_name("homo sapiens", lower_case=False, regex=False))
    assert len(results) == 0
    results = list(taxonomy.search_scientific_name("Homo sapiens", lower_case=False, regex=False))
    assert len(results) == 1
    assert results[0].identifier == HUMAN_TAXID
    assert results[0].name == HUMAN_NAME
    results = list(taxonomy.search_scientific_name("homo sapiens", lower_case=True, regex=False))
    assert len(results) == 1
    assert results[0].identifier == HUMAN_TAXID
    assert results[0].name == HUMAN_NAME
    results = list(taxonomy.search_scientific_name("homo sap", lower_case=True, regex=False))
    assert len(results) == 0
    results = list(
        taxonomy.search_scientific_name("[Homo]{4} sapiens", lower_case=False, regex=False)
    )
    assert len(results) == 0


def test_search_scientific_name_regex(taxonomy: Taxonomy) -> None:
    results = list(
        taxonomy.search_scientific_name("^homo [sapiens]$", lower_case=False, regex=True)
    )
    assert len(results) == 0
    results = list(
        taxonomy.search_scientific_name("^[Homo]{4} sapiens$", lower_case=False, regex=True)
    )
    assert len(results) == 1
    assert results[0].identifier == HUMAN_TAXID
    assert results[0].name == HUMAN_NAME
    results = list(
        taxonomy.search_scientific_name("^homo [SAPIENS]{7}$", lower_case=True, regex=True)
    )
    assert len(results) == 1
    assert results[0].identifier == HUMAN_TAXID
    assert results[0].name == HUMAN_NAME
    results = list(taxonomy.search_scientific_name("^homo sap$", lower_case=True, regex=True))
    assert len(results) == 0


def test_lineage(taxonomy: Taxonomy) -> None:
    lineage = list(taxonomy.find_lineage(HUMAN_ASSEMBLY, stop_rank="class"))
    assert all(isinstance(entry, TaxonomyEntry) for entry in lineage)
    identifiers = [entry.identifier for entry in lineage]
    assert "2759" not in identifiers, "Eukaryotes in lineage although beyond stop rank."
    expected = [HUMAN_ASSEMBLY, HOMINIDAE_TAXID, "9443", "40674"]
    index = -1
    for identifier in expected:
        find = identifiers.index(identifier)
        assert find > index, f"Identifier {identifier} not in correct order."
        index = find
    assert index == len(identifiers) - 1, "40674 not last in lineage despite being the stop rank."


def test_lineage_gtdb(taxonomy_ncbi_and_gtdb: Taxonomy) -> None:
    lineage = list(taxonomy_ncbi_and_gtdb.find_lineage("GCF_000744315.1", stop_rank="class"))
    assert all(isinstance(entry, TaxonomyEntry) for entry in lineage)
    identifiers = [entry.identifier for entry in lineage]
    assert "gtdb:d__Archaea" not in identifiers, (
        "gtdb:d__Archaea in lineage although beyond stop rank."
    )
    expected = [
        "GCF_000744315.1",
        "gtdb:s__Methanosarcina mazei",
        "gtdb:g__Methanosarcina",
        "gtdb:f__Methanosarcinaceae",
        "gtdb:o__Methanosarcinales",
        "gtdb:c__Methanosarcinia",
    ]
    index = -1
    for identifier in expected:
        find = identifiers.index(identifier)
        assert find > index, f"Identifier {identifier} not in correct order."
        index = find
    assert index == len(identifiers) - 1, (
        "gtdb:c__Methanosarcinia not last in lineage despite being the stop rank."
    )


def test_gtdb_root_linked_to_ncbi_root(taxonomy_ncbi_and_gtdb: Taxonomy) -> None:
    # Check if the GTDB root is linked to the NCBI root
    lineage = list(taxonomy_ncbi_and_gtdb.find_lineage("GCF_000744315.1"))
    identifiers = [entry.identifier for entry in lineage]
    assert "gtdb:root" in identifiers, "GTDB root not in lineage."
    assert "1" in identifiers, "NCBI root not in lineage."


def test_children(taxonomy: Taxonomy) -> None:
    children = list(taxonomy.find_children(HOMINIDAE_TAXID))
    identifiers = {child.identifier for child in children}
    assert HUMAN_ASSEMBLY in identifiers, "Human assembly not in children of Hominidae."
    assert HUMAN_TAXID in identifiers, "Human taxid not in children of Hominidae."
    assert HOMINIDAE_TAXID not in identifiers, "Hominidae taxid should not be in its own children."


def test_children_with_lineage(taxonomy: Taxonomy) -> None:
    children = list(taxonomy.find_children_with_lineage(HOMINIDAE_TAXID))
    identifiers = {lineage[0].identifier for lineage in children}
    assert HUMAN_ASSEMBLY in identifiers, "Human assembly not in children of Hominidae."
    assert HUMAN_TAXID in identifiers, "Human taxid not in children of Hominidae."
    assert HOMINIDAE_TAXID not in identifiers, "Hominidae taxid should not be in its own child."
    for lineage in children:
        reference_lineage = list(taxonomy.find_lineage(lineage[0].identifier, stop_rank="family"))
        assert all(e1.identifier == e2.identifier for e1, e2 in zip(lineage, reference_lineage)), (
            f"Child lineage {lineage} does not match reference lineage {reference_lineage}."
        )


def test_children_gtdb(taxonomy_ncbi_and_gtdb: Taxonomy) -> None:
    children = list(taxonomy_ncbi_and_gtdb.find_children("gtdb:o__Acidiferrobacterales"))
    identifiers = {child.identifier for child in children}
    assert "GCA_035292245.1" in identifiers, (
        "GCA_035292245.1 assembly not in children of gtdb:o__Acidiferrobacterales."
    )
    assert "gtdb:f__SPGG2" in identifiers, (
        "gtdb:f__SPGG2 taxid not in children of gtdb:o__Acidiferrobacterales."
    )


def test_nearest_neighbor(taxonomy: Taxonomy) -> None:
    neighbor = taxonomy.find_nearest_neighbor(HUMAN_TAXID, {"4897", "10090"})
    assert neighbor is not None, "No nearest neighbor found."
    assert neighbor == "10090", "Nearest neighbor not 10090."
