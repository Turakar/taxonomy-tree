import pickle
from pathlib import Path

import pytest

from taxonomy_tree import Taxonomy, TaxonomyEntry

# https://goat.genomehubs.org
HUMAN_ASSEMBLY = "GCA_000001405.29"
HUMAN_NAME = "Homo sapiens"
HUMAN_TAXID = "9606"


@pytest.fixture()
def taxonomy() -> Taxonomy:
    path = Path("taxonomy.db")
    taxonomy = Taxonomy(path)
    if not path.exists():
        print("Creating taxonomy database")
        taxonomy.create_db()
    return taxonomy


def test_pickle(taxonomy: Taxonomy) -> None:
    assert taxonomy.contains(HUMAN_ASSEMBLY)
    taxonomy2: Taxonomy = pickle.loads(pickle.dumps(taxonomy))
    assert taxonomy2.contains(HUMAN_ASSEMBLY)


def test_scientific_name(taxonomy: Taxonomy) -> None:
    scientific_name = taxonomy.find_scientific_name(HUMAN_TAXID)
    assert scientific_name == HUMAN_NAME


def test_lineage(taxonomy: Taxonomy) -> None:
    lineage = list(taxonomy.find_lineage(HUMAN_ASSEMBLY, stop_rank="class"))
    assert all(isinstance(entry, TaxonomyEntry) for entry in lineage)
    identifiers = [entry.identifier for entry in lineage]
    assert "2759" not in identifiers, "Eukaryotes in lineage although beyond stop rank."
    expected = [HUMAN_ASSEMBLY, "9604", "9443", "40674"]
    index = -1
    for identifier in expected:
        find = identifiers.index(identifier)
        assert find > index, f"Identifier {identifier} not in correct order."
        index = find
    assert index == len(identifiers) - 1, "40674 not last in lineage despite being the stop rank."


def test_children(taxonomy: Taxonomy) -> None:
    children = list(taxonomy.find_children("9604"))
    identifiers = {child.identifier for child in children}
    assert HUMAN_ASSEMBLY in identifiers, "Human assembly not in children of 9604."
    assert HUMAN_TAXID in identifiers, "Human taxid not in children of 9604."


def test_nearest_neighbor(taxonomy: Taxonomy) -> None:
    neighbor = taxonomy.find_nearest_neighbor(HUMAN_TAXID, {"4897", "10090"})
    assert neighbor is not None, "No nearest neighbor found."
    assert neighbor == "10090", "Nearest neighbor not 10090."
