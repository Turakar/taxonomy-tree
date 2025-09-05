from pathlib import Path
from typing import Any, Iterator

from .taxonomy_tree import TaxonomyEntry
from .taxonomy_tree import TaxonomyTree as _TaxonomyTree


class Taxonomy:
    def __init__(self, db_path: str | Path):
        self.db_path = str(db_path)
        self._db_: _TaxonomyTree | None = None

    def create_db(self, include_gtdb: bool = False) -> None:
        from .build import make_taxonomy_tree

        make_taxonomy_tree(self.db_path, include_gtdb=include_gtdb)

    @property
    def _db(self) -> _TaxonomyTree:
        if self._db_ is None:
            self._db_ = _TaxonomyTree(str(self.db_path))
        return self._db_

    def __getstate__(self) -> Any:
        return {"db_path": str(self.db_path)}

    def __setstate__(self, state: Any) -> None:
        self.db_path = Path(state["db_path"])
        self._db_ = None

    def contains(self, identifier: str) -> bool:
        return self._db.contains(identifier)

    def find_scientific_name(self, identifier: str) -> str:
        return self._db.find_scientific_name(identifier)

    def search_scientific_name(
        self, search: str, lower_case: bool = True, regex: bool = False
    ) -> list[TaxonomyEntry]:
        return self._db.search_scientific_name(
            search=search,
            lower_case=lower_case,
            regex=regex,
        )

    def find_lineage(
        self,
        identifier: str,
        stop_rank: str | None = None,
    ) -> Iterator[TaxonomyEntry]:
        return self._db.find_lineage(
            id=identifier,
            stop_rank=stop_rank,
        )

    def find_children(
        self,
        identifier: str,
        stop_rank: str | None = None,
        skip: set[str] | None = None,
    ) -> Iterator[TaxonomyEntry]:
        for lineage in self.find_children_with_lineage(
            identifier=identifier,
            stop_rank=stop_rank,
            skip=skip,
        ):
            if len(lineage) > 1:
                # Do not yield the lineage of the identifier itself, only the children
                yield lineage[0]

    def find_children_with_lineage(
        self,
        identifier: str,
        stop_rank: str | None = None,
        skip: set[str] | None = None,
    ) -> Iterator[list[TaxonomyEntry]]:
        for child in self._db.find_children(
            id=identifier,
            skip=skip,
            stop_rank=stop_rank,
        ):
            if child[0].identifier != identifier:
                yield child

    def depth_first_search(
        self,
        identifier: str,
        stop_rank_top: str | None = None,
        stop_rank_bottom: str | None = None,
    ) -> Iterator[TaxonomyEntry]:
        return self._db.depth_first_search(
            id=identifier,
            stop_rank_top=stop_rank_top,
            stop_rank_bottom=stop_rank_bottom,
        )

    def find_nearest_neighbor(self, identifier: str, targets: set[str]) -> str | None:
        for entry in self.depth_first_search(identifier):
            if entry.identifier in targets:
                return entry.identifier
        return None


if __name__ == "__main__":
    import pickle

    taxonomy = Taxonomy("taxonomy_ncbi.db")
    if not Path(taxonomy.db_path).exists():
        taxonomy.create_db()

    print(f"9606: {taxonomy.find_scientific_name('9606')}")

    taxonomy = pickle.loads(pickle.dumps(taxonomy))

    print(f"GCA_000001405.29: {taxonomy.find_scientific_name('GCA_000001405.29')}")
    print(f"and its lineage: {list(taxonomy.find_lineage('GCA_000001405.29', stop_rank='class'))}")
    print(f"Children of 4897: {list(taxonomy.find_children('4897'))}")
    children_314146 = list(taxonomy.find_children("314146", stop_rank="species"))
    print("10090 child of 314146:", "10090" in {c.identifier for c in children_314146})
    print(
        f"Nearest neighbor of 9606 in (4897, 10090): {taxonomy.find_nearest_neighbor('9606', {'4897', '10090'})}"
    )
