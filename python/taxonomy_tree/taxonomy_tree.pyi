from typing import Iterator

class TaxonomyEntry:
    identifier: str
    name: str
    rank: str | None

class TaxonomyTreeBuilder:
    def __init__(self) -> None: ...
    def insert(
        self,
        identifier: str,
        name: str,
        rank: str | None,
        parent: str | None,
        overwrite: bool,
    ) -> None: ...
    def write(self, path: str) -> None: ...

class TaxonomyTree:
    def __init__(self, path: str) -> None: ...
    def contains(self, identifier: str) -> bool: ...
    def find_scientific_name(self, identifier: str) -> str: ...
    def search_scientific_name(
        self, search: str, lower_case: bool, regex: bool
    ) -> list[TaxonomyEntry]: ...
    def find_lineage(self, id: str, stop_rank: str | None) -> Iterator[TaxonomyEntry]: ...
    def find_children(
        self, id: str, skip: set[str] | None, stop_rank: str | None
    ) -> Iterator[list[TaxonomyEntry]]: ...
    def depth_first_search(
        self,
        id: str,
        stop_rank_top: str | None,
        stop_rank_bottom: str | None,
    ) -> Iterator[TaxonomyEntry]: ...
