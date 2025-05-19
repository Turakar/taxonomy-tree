import pickle

from taxonomy_tree import Taxonomy


def main():
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


if __name__ == "__main__":
    main()
