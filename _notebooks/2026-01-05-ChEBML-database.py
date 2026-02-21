import marimo

__generated_with = "0.19.11"
app = marimo.App()


@app.cell
def _(mo):
    mo.md(r"""
    # ChEMBL Compounds, Targets, and Rule of 5
    """)
    return


@app.cell
def _(mo):
    mo.md(r"""
    When reviewing data to find pharma compounds for virtual screening, we might want to check what they target and rank candidates by how many [Lipinski's rule of five](https://en.wikipedia.org/wiki/Lipinski's_rule_of_five) violations they have--the fewer the better. This post uses the ChEMBL API and a SQLite database to do that.
    """)
    return


@app.cell
def _(mo):
    mo.md(r"""
    This post pulls data from ChEMBL using its chembl_webresource_client for Python. It's a helpful package which handles the API calls. It also provides caching so you won't accidentally run the same queries more than once. APIs often ask users to cache the results; I like that ChEMBL goes ahead and does that for you. (If it didn't, I would have used [DiskCache](https://pypi.org/project/diskcache/), which as the name implies caches results to disk so they persist across code runs, and which I've found works well for storing results from other API calls.)

    We write the results directly to a SQLite database. SQLite is file-based so its uptime is nearly 100%. That means we don't need to worry about its availability. Of course it being file-based is not ideal if users are distributed across the Internet, but that's not what we're doing here.
    """)
    return


@app.cell
def _(mo):
    mo.md(r"""
    ## Database schemas--conceptual
    """)
    return


@app.cell
def _(mo):
    mo.md(r"""
    ChEMBL has the connections as compounds ↔ activities ↔ targets. Let's plot that as an [entity-relationship diagram](https://en.wikipedia.org/wiki/Entity%E2%80%93relationship_model) (ERD) using [graphviz](https://graphviz.org/).
    """)
    return


@app.cell
def _():
    import logging
    import time
    from collections import defaultdict

    import pydot
    import sqlalchemy
    from chembl_webresource_client.new_client import new_client
    from graphviz import Digraph
    from IPython.display import SVG, display
    from sqlalchemy import (
        Column,
        Float,
        Integer,
        String,
        create_engine,
        func,
        select,
    )
    from sqlalchemy.exc import IntegrityError, SQLAlchemyError
    from sqlalchemy.orm import declarative_base, sessionmaker
    from sqlalchemy_schemadisplay import create_schema_graph

    return (
        Column,
        Digraph,
        Float,
        Integer,
        IntegrityError,
        SQLAlchemyError,
        SVG,
        String,
        create_engine,
        create_schema_graph,
        declarative_base,
        defaultdict,
        display,
        func,
        logging,
        new_client,
        pydot,
        select,
        sessionmaker,
        sqlalchemy,
        time,
    )


@app.cell
def _():
    # Scale down size of ERD relationship labels and arrowheads
    crow_fontsize = "8"  # smaller font for labels
    arrow_scale = "0.4"  # smaller arrowheads
    return arrow_scale, crow_fontsize


@app.cell
def _(Digraph, SVG, arrow_scale, crow_fontsize, display):
    # Create the diagram
    dot_chembl = Digraph(format="svg")
    dot_chembl.attr(rankdir="LR", splines="ortho")
    dot_chembl.attr(
        "node", shape="record", fontsize="10", style="filled", fillcolor="lightgrey"
    )

    # Nodes
    dot_chembl.node("Compound", "{Compound|compound_id (PK)}", fillcolor="#A3C1DA")
    dot_chembl.node("Activity", "{Activity|activity_id (PK)}", fillcolor="#A3C1DA")
    dot_chembl.node("Target", "{Target|target_id (PK)}", fillcolor="#A3C1DA")

    dot_chembl.edge(
        "Compound",
        "Activity",
        fontsize=crow_fontsize,
        arrowhead="crow",
        arrowtail="crow",
        dir="both",
        color="black",
        arrowsize=arrow_scale,
    )
    dot_chembl.edge(
        "Activity",
        "Target",
        fontsize=crow_fontsize,
        arrowhead="crow",
        arrowtail="crow",
        dir="both",
        color="black",
        arrowsize=arrow_scale,
    )

    # Render inline SVG
    svg_chembl = dot_chembl.pipe(format="svg").decode("utf-8")
    display(SVG(svg_chembl))
    return


@app.cell
def _(mo):
    mo.md(r"""
    That makes sense as a comprehensive schema; in this code, I wanted to simplify it by connecting compounds to targets directly. So my schema is simply compounds ↔ targets:
    """)
    return


@app.cell
def _(Digraph, SVG, arrow_scale, crow_fontsize, display):
    # Create the diagram
    dot_simple = Digraph(format="svg")
    dot_simple.attr(rankdir="LR", splines="ortho")
    dot_simple.attr(
        "node", shape="record", fontsize="10", style="filled", fillcolor="lightgrey"
    )

    # Define nodes
    dot_simple.node("Compound", "{Compound|compound_id (PK)}", fillcolor="#A3C1DA")
    dot_simple.node("Target", "{Target|target_id (PK)}", fillcolor="#A3C1DA")

    # Many-to-many edge (crow's foot at both ends)
    dot_simple.edge(
        "Compound",
        "Target",
        fontsize=crow_fontsize,
        arrowhead="crow",
        arrowtail="crow",
        dir="both",
        color="black",
        arrowsize=arrow_scale,
    )

    # Render inline SVG
    svg_simple = dot_simple.pipe(format="svg").decode("utf-8")
    display(SVG(svg_simple))
    return


@app.cell
def _(mo):
    mo.md(r"""
    ## Marimo pros and cons
    """)
    return


@app.cell
def _(mo):
    mo.md(r"""
    This is my first post using Marimo as the notebook: I had previously used just Jupyter. Initially Marimo wasn't great because

    - I discovered that publishing to a Markdown file directly from Marimo didn't lead to good formatting on my Jekyll, so I converted from Marimo to Jupyter and then to Markdown
    - I didn't have a reliable Internet connection when working on this and Marimo seemed to need a connection in VS Code
    - I sometimes am not allowed to rename a variable in VS Code (and automatically change the variable name wherever it's used)
    - The ruff VS Code extension doesn't work that well with Marimo notebooks for me
    - Marimo is very serious about not allowing you to re-use a variable name, even an iterator variable which is usually a throw-away. So for example if you're trying out two versions of a code block, you have to rename every variable, even the iterators. While I understand the need for this, perhaps there's a way to make it easier by specifying a suffix to append to each variable name when you clone a code block--seems like something an LLM could handle.

    However, it seemed worth it when, before committing via git, the diff was so much more readable than in Jupyter (which is a ton of TypeScript, metadata, etc.). With Marimo, the diff is just the actual code changes and a small amount of formatting in Python. With Jupyter Notebooks, in theory source control works, but in practice the diff is so large it's so difficult to tell what changes were made that I didn't find it useful for identifying or rolling back changes.
    """)
    return


@app.cell
def _(mo):
    mo.md(r"""
    ### Code setup
    """)
    return


@app.cell
def _(logging):
    # Set logger level to INFO
    logging.basicConfig(
        level=logging.INFO,
        format="[%(levelname)s] %(message)s",
        force=True,
    )
    logger = logging.getLogger(__name__)
    return (logger,)


@app.cell
def _(create_engine, declarative_base, sessionmaker):
    Base = declarative_base()

    engine = create_engine("sqlite:///compounds.db", echo=False)
    Session = sessionmaker(bind=engine)
    return Base, Session, engine


@app.cell
def _(mo):
    mo.md(r"""
    ## ChEMBL data fetching
    """)
    return


@app.cell
def _(mo):
    mo.md(r"""
    Now let's get the data from ChEMBL. The `chembl_webresource_client` `new_client` provides objects for molecule, activity, and target.

    We'll start with a set of molecules, then fetch associated data. Here we fetch a number of compounds `n_compounds` starting at ChEMBL ID `start_id`. Because ChEMBL assigns a ChEMBL ID to not only compounds but also targets, we'll count how many IDs correspond to compounds.

    Then we'll fetch associated activities (interactions between compounds and targets) for those compounds. Let's say we're interested in older (before 2010) studies on humans with IC50 standards and binding assays. We specify those as filters when we ask for the activities. This will serve to limit us to a reasonable number of activities and thus targets.

    Then we get target data including its ChEMBL ID, name, type, and organism.

    Lastly we associate targets with compounds by creating a list of targets for each molecule. This will make it easier to populate our database tables.
    """)
    return


@app.cell
def _(mo):
    mo.md(r"""
    #NoteToSelf: Instead of determining all_target_ids twice, output it from get_chembl_molecules and import it into save_compounds_to_db?
    """)
    return


@app.cell
def _(defaultdict, logger, logging, new_client):
    def get_chembl_molecules(
        n_compounds: int = 2,
        start_id: int = 1,
    ):
        chembl_ids = [f"CHEMBL{id}" for id in range(start_id, start_id + n_compounds)]

        molecule = new_client.molecule
        activity = new_client.activity
        target = new_client.target

        # Suppress INFO logs from the chembl package while keeping global INFO level
        logging.getLogger("chembl_webresource_client").setLevel(logging.WARNING)

        # --------------------
        # 1) Fetch molecules
        # --------------------
        mols = list(
            molecule.filter(molecule_chembl_id__in=chembl_ids).only(
                [
                    "molecule_chembl_id",
                    "molecule_structures",
                    "pref_name",
                    "molecule_properties",
                ]
            )
        )

        # From the molecules in mols, extract the digits after "CHEMBL" to check which IDs were found
        chembl_ids_found = set()
        for mol in mols:
            chembl_id = mol.get("molecule_chembl_id", "")
            if chembl_id.startswith("CHEMBL"):
                chembl_ids_found.add(int(chembl_id.replace("CHEMBL", "")))
        chembl_ids_found_str = ", ".join(map(str, sorted(chembl_ids_found)))

        logger.info(
            f"Of the {n_compounds} ChEMBL IDs "
            f"({start_id}-{start_id + n_compounds - 1}), {len(mols)} are compounds: "
            f"ChEMBL IDs {chembl_ids_found_str}"
        )

        if not mols:
            return []

        mol_ids_present = [m["molecule_chembl_id"] for m in mols]

        # ---------------------------------
        # 2) Bulk fetch activities → targets
        # ---------------------------------
        activities = activity.filter(
            molecule_chembl_id__in=mol_ids_present,
            target_organism="Homo sapiens",
            standard_type="IC50",
            assay_type="B",
            document_year__lt=2010,
        ).only(
            [
                "molecule_chembl_id",
                "target_chembl_id",
            ]
        )
        logger.info(f"Activities: {len(activities)}")

        mol_to_target_ids = defaultdict(set)

        for act in activities:
            try:
                mol_to_target_ids[act["molecule_chembl_id"]].add(act["target_chembl_id"])
            except Exception as e:
                logger.warning(e)

        # ---------------------------------
        # 3) Fetch target metadata (bulk)
        # ---------------------------------
        all_target_ids = sorted(
            {tar_id for tar_ids in mol_to_target_ids.values() for tar_id in tar_ids}
        )

        # Create dictionary of ChEMLB ID: target entries
        targets = {}
        if all_target_ids:
            for t in target.filter(target_chembl_id__in=all_target_ids).only(
                [
                    "target_chembl_id",
                    "pref_name",
                    "target_type",
                    "organism",
                ]
            ):
                targets[t["target_chembl_id"]] = t
        logger.info(f"Fetched metadata for {len(targets)} targets.")

        # ---------------------------------
        # 4) Attach targets to molecules
        # ---------------------------------

        for m in mols:
            t_ids = mol_to_target_ids.get(m["molecule_chembl_id"], [])
            m["targets"] = [targets[tar_id] for tar_id in t_ids if tar_id in targets]

        return mols

    return (get_chembl_molecules,)


@app.cell
def _(mo):
    mo.md(r"""
    Now let's define a function to save our compounds and targets to our SQLite database. To avoid duplication, we start by preloading all the targets into that table. Then we create a dictionary, which is an O(1) lookup, between target ChEMBL ID and its database entry so we can quickly link the compound to the target. That saves us from having to query the database each time we want to associate a target with a compound.
    existing_target=('CHEMBL1855', <__main__.Target object at 0x11623bb90>)
    """)
    return


@app.cell
def _(
    Compound,
    CompoundTarget,
    IntegrityError,
    SQLAlchemyError,
    Session,
    Target,
    logger,
):
    def save_compounds_to_db(molecules: list[dict]) -> tuple[int, int, int]:
        """Save multiple compounds and their targets to the database efficiently avoiding duplicate Targets."""
        # collect all target ids present in incoming molecules
        all_target_ids = {
            t["target_chembl_id"]
            for m in molecules
            for t in m.get("targets", [])
            if t.get("target_chembl_id")
        }

        # build a dictionary of metadata keyed by target Chembl ID to deduplicate
        target_map: dict[str, dict] = {
            t["target_chembl_id"]: {
                "organism": t.get("organism"),
                "pref_name": t.get("pref_name"),
                "target_chembl_id": t.get("target_chembl_id"),
                "target_type": t.get("target_type"),
            }
            for m in molecules
            for t in m.get("targets", [])
            if t.get("target_chembl_id")
        }
        # convert back to list of dicts for bulk insert
        all_targets = list(target_map.values())
        if all_targets:
            print(f"first unique target={all_targets[0]}")

        logger.info(f"{len(all_target_ids)=}: {all_target_ids=}")
        n_mols_saved = 0
        n_targets_saved = 0
        n_compounds_targets_saved = 0

        with Session() as db_session:
            # Bulk insert targets into Target table
            db_session.bulk_insert_mappings(Target, all_targets)
            db_session.commit()

            try:
                # preload existing targets into a mapping chembl_id -> Target instance
                existing_targets = {}
                if all_target_ids:
                    rows = (
                        db_session.query(Target)
                        .filter(Target.target_chembl_id.in_(list(all_target_ids)))
                        .all()
                    )
                    logger.info(f"{len(rows)=}: {rows=}")
                    existing_targets = {r.target_chembl_id: r for r in rows}
                    for existing_target in existing_targets.items():
                        print(f"{existing_target=}")
                        break

                for mol in molecules:
                    chembl_id = mol.get("molecule_chembl_id")
                    pref_name = mol.get("pref_name")
                    props = mol.get("molecule_properties", {}) or {}
                    compound = Compound(
                        chembl_id=chembl_id,
                        sml=mol.get("molecule_structures", {}).get("canonical_smiles"),
                        pref_name=pref_name,
                        molwt=props.get("full_molweight"),
                        tpsa=props.get("tpsa"),
                        num_h_acceptors=props.get("num_h_acceptors"),
                        num_h_donors=props.get("num_h_donors"),
                        num_ro5=props.get("num_ro5_violations"),
                        mol_logp=props.get("alogp"),
                    )

                    with db_session.begin_nested():
                        db_session.add(compound)
                        db_session.flush()  # get compound.id

                        for target_data in mol.get("targets", []):
                            target_id = target_data.get("target_chembl_id")
                            if not target_id:
                                continue

                            # reuse existing Target if present
                            target_obj = existing_targets.get(target_id)
                            logger.info(f" {target_id=}, {target_obj=}")

                            if target_obj is None:
                                logger.info(f"Had to create {target_id=}")
                                # create new Target, add and flush to get id, then cache it to existing_targets
                                target_obj = Target(
                                    organism=target_data.get("organism"),
                                    pref_name=target_data.get("pref_name"),
                                    target_chembl_id=target_id,
                                    target_type=target_data.get("target_type"),
                                )
                                db_session.add(target_obj)
                                db_session.flush()  # populates target_obj.id
                                existing_targets[target_id] = target_obj
                                n_targets_saved += 1

                            # create association
                            compound_target = CompoundTarget(
                                compound_id=compound.id,
                                target_id=target_obj.id,
                            )
                            db_session.add(compound_target)
                            n_compounds_targets_saved += 1

                    n_mols_saved += 1

                # outer commit happens when exiting the with Session() context
            except IntegrityError as e:
                # handle unexpected integrity issues gracefully
                logger.info(f"IntegrityError while saving: {e}")
                db_session.rollback()
            except SQLAlchemyError:
                logger.exception("Database error saving compounds")
                db_session.rollback()
                raise

        return n_mols_saved, n_targets_saved, n_compounds_targets_saved

    return (save_compounds_to_db,)


@app.cell
def _(Base, engine):
    def init_db():
        """Create database tables (call once at app startup or from scripts/tests)."""
        Base.metadata.create_all(engine)

    return (init_db,)


@app.cell
def _(Base, engine):
    def reset_db():
        """Drop all database tables (use with caution)."""
        Base.metadata.drop_all(engine)

    return (reset_db,)


@app.cell
def _(Base, Column, Float, Integer, String):
    class Compound(Base):
        __tablename__ = "compound"

        id = Column(Integer, primary_key=True)
        chembl_id = Column(Integer, unique=True)
        sml = Column(String)
        pref_name = Column(String)
        molwt = Column(Float)  # MolWt
        tpsa = Column(Float)  # TPSA
        num_h_acceptors = Column(Integer)  # NumHAcceptors
        num_h_donors = Column(Integer)  # NumHDonors
        num_ro5 = Column(Integer)  # NumRo5
        mol_logp = Column(Float)  # MolLogP

    return (Compound,)


@app.cell
def _(Base, Column, Integer, String):
    class Target(Base):
        __tablename__ = "target"

        id = Column(Integer, primary_key=True)
        organism = Column(String)
        pref_name = Column(String)
        target_chembl_id = Column(String, unique=True)
        target_type = Column(String)

    return (Target,)


@app.cell
def _(Base, Column, Integer, sqlalchemy):
    # Join table
    class CompoundTarget(Base):
        __tablename__ = "compound_target"

        id = Column(Integer, primary_key=True)
        compound_id = Column(Integer, sqlalchemy.ForeignKey("compound.id"))
        target_id = Column(Integer, sqlalchemy.ForeignKey("target.id"))

    return (CompoundTarget,)


@app.cell
def _(init_db, reset_db):
    # Reset database (uncomment to start fresh)
    reset_db()

    # Ensure tables exist
    init_db()
    return


@app.cell
def _(pydot):
    def add_ordering_edges(graph, Base, exclude_tables=None):
        """
        Add invisible edges based on foreign key relationships to enforce left-to-right ordering

        Args:
            graph: A pydot.Dot graph object (e.g. ERD) to add edges to.
            Base: SQLAlchemy declarative base containing table metadata.
            exclude_tables: Set of table names to exclude from processing

        Returns:
            The modified graph object with invisible ordering edges added.
        """
        if exclude_tables is None:
            exclude_tables = set()

        # Get all table names
        tables = Base.metadata.tables.keys()

        # For each table, check for foreign keys
        for table_name in tables:
            if table_name in exclude_tables:
                continue

            table = Base.metadata.tables[table_name]

            for fk in table.foreign_key_constraints:
                parent_table = fk.referred_table.name

                if parent_table not in exclude_tables:
                    # Add invisible edge: parent -> child
                    graph.add_edge(pydot.Edge(parent_table, table_name, style="invis"))

        return graph

    return (add_ordering_edges,)


@app.cell
def _(Base, SVG, add_ordering_edges, create_schema_graph, display, engine):
    # Create the ERD graph
    graph_full = create_schema_graph(
        engine=engine,
        metadata=Base.metadata,
        show_datatypes=True,
        show_indexes=False,
        rankdir="LR",
        concentrate=False,
    )

    # Force strict left-to-right ordering with invisible edges;
    # add them programmatically by inspecting the SQLAlchemy model
    add_ordering_edges(graph_full, Base)

    graph_full.set("splines", "ortho")

    # Move FK labels horizontally away from edges
    for edge_full in graph_full.get_edges():

        head_full = edge_full.get_headlabel()
        tail_full = edge_full.get_taillabel()

        if head_full:
            clean_head_full = head_full.replace("+ ", "").replace("+", "")
            edge_full.set_headlabel(clean_head_full)

        if tail_full:
            clean_tail_full = tail_full.replace("+ ", "").replace("+", "")
            edge_full.set_taillabel(clean_tail_full)

            edge_full.set_label("")  # critical fix
            edge_full.set_labeldistance("2.5")
        
    # Increase horizontal spacing between tables
    graph_full.set("ranksep", "1.0")
    svg_content_full = graph_full.create_svg()
    display(SVG(svg_content_full))
    return (graph_full,)


@app.cell
def _(
    Base,
    SVG,
    add_ordering_edges,
    create_schema_graph,
    display,
    engine,
    graph_full,
    pydot,
):
    def detect_join_tables(Base):
        """
        Detect join tables (tables with only an id and two foreign keys).

        Returns:
            dict: Mapping of join_table_name -> (parent_table1, parent_table2)
        """
        join_tables = {}

        for table_name, table in Base.metadata.tables.items():
            # Get foreign key constraints
            fk_constraints = list(table.foreign_key_constraints)

            # A join table typically has:
            # - Only id as primary key (or composite PK of the two FKs)
            # - Exactly 2 foreign key columns
            # - Minimal or no other columns

            foreign_key_columns = []
            for fk in fk_constraints:
                foreign_key_columns.extend(fk.column_keys)

            # Check if this looks like a join table:
            # Has exactly 2 FK constraints pointing to different tables
            if len(fk_constraints) == 2:
                referred_tables = [fk.referred_table.name for fk in fk_constraints]

                # Make sure they point to different tables
                if referred_tables[0] != referred_tables[1]:
                    join_tables[table_name] = tuple(referred_tables)

        return join_tables

    def get_primary_key_name(Base, table_name):
        """
        Get the primary key column name for a table.

        Args:
            Base: SQLAlchemy declarative base
            table_name: Name of the table

        Returns:
            str: Primary key column name (or comma-separated list if composite)
        """
        table = Base.metadata.tables[table_name]
        pk_columns = [col.name for col in table.columns if col.primary_key]

        if pk_columns:
            return ", ".join(pk_columns)
        else:
            # Fallback
            return "id"

    def remove_join_tables_from_graph(graph, Base):
        """
        Remove join tables from graph and replace with direct many-to-many edges.

        Args:
            graph: pydot.Dot graph object
            Base: SQLAlchemy declarative base

        Returns:
            set: Names of removed join tables
        """
        join_tables = detect_join_tables(Base)

        for join_table, (table1, table2) in join_tables.items():
            # Remove the join table node
            graph.del_node(join_table)

            # Remove all edges connected to the join table
            for edge in list(graph.get_edges()):
                if (
                    edge.get_source() == join_table
                    or edge.get_destination() == join_table
                ):
                    graph.del_edge(edge.get_source(), edge.get_destination())

            # Get primary key names for both tables
            table1_pk = get_primary_key_name(Base, table1)
            table2_pk = get_primary_key_name(Base, table2)

            # Add a direct many-to-many edge between the two tables
            graph.add_edge(
                pydot.Edge(
                    table1,
                    table2,
                    taillabel=table1_pk,
                    headlabel=table2_pk,
                    arrowhead="crow",
                    arrowtail="crow",
                    dir="both",
                )
            )

        return set(join_tables.keys())

    # Create the ERD graph
    graph = create_schema_graph(
        engine=engine,
        metadata=Base.metadata,
        show_datatypes=True,
        show_indexes=False,
        rankdir="LR",
        concentrate=False,
    )

    # Automatically detect and remove join tables
    excluded_tables = remove_join_tables_from_graph(graph, Base)

    # Add ordering edges (excluding detected join tables)
    add_ordering_edges(graph, Base, exclude_tables=excluded_tables)

    graph_full.set("splines", "ortho")

    # Move FK labels horizontally away from edges
    for edge in graph.get_edges():
        head = edge.get_headlabel()
        tail = edge.get_taillabel()

        if head:
            edge.set_headlabel(head)

        if tail:
            edge.set_taillabel(tail)
            edge.set_labeldistance("2.5")

    graph_full.set("ranksep", "1.0")

    svg_content = graph.create_svg()
    display(SVG(svg_content))
    return


@app.cell
def _(get_chembl_molecules, logger, save_compounds_to_db, time):
    # Measure how long it takes to fetch ChEMBL molecules
    start = time.time()
    mols = get_chembl_molecules(
        n_compounds=30,
        start_id=1000,  # Has targets
        # start_id=3430873, # Not a molecule
    )

    end = time.time()
    logger.info(f"Fetched {len(mols)} molecules in {end - start:.2f} seconds.")

    start = time.time()

    n_mols_saved, n_targets_saved, n_compounds_targets_saved = save_compounds_to_db(
        mols
    )
    logger.info(
        f"Saved {n_mols_saved} compounds, "
        f"{n_targets_saved} targets, "
        f"{n_compounds_targets_saved} compound-target associations to the database, "
        f"in {time.time() - start:.2f} seconds."
    )
    return


@app.cell
def _(Compound, CompoundTarget, Session, Target, func, logger, select):
    def run_queries():
        """Run the required queries against the ChEMBL database and print the results."""
        with Session() as db_session:
            print("run_queries")
            # 1. Find the count of all distinct counts of compound targets (ex: Voltage-gated inwardly rectifying potassium channel KCNH2:14, Neuronal acetylcholine receptor subunit alpha-3/Neuronal acetylcholine receptor subunit alpha-7: 5).
            # For each compound, concatenate its targets with a slash.
            # Then, count how many distinct such tuples exist in the database.

            # Build per-compound target combo subquery:
            # as correlated scalar subquery: group_concat over an ordered selection of types to ensure consistent ordering.

            # correlated inner select returning pref_name for the current Compound, ordered
            inner = (
                select(Target.pref_name)
                .select_from(
                    Target.__table__.join(
                        CompoundTarget.__table__, CompoundTarget.target_id == Target.id
                    )
                )
                .where(CompoundTarget.compound_id == Compound.id)
                .order_by(func.lower(Target.pref_name))
                .correlate(Compound)
            )

            # name the derived table so the outer group_concat can select from it
            ordered_targets = inner.subquery("ordered_targets")

            # aggregate the ordered names with a '\' separator
            target_combo_subq = select(
                func.group_concat(ordered_targets.c.pref_name, "\\")
            ).scalar_subquery()

            # Create a subquery that selects each compound's id and its target combination.
            target_combinations = db_session.query(
                Compound.id.label("compound_id"),
                target_combo_subq.label("target_combo"),
            ).subquery()

            # Create a list of distinct type combinations and their counts
            # where each is a tuple like (target_combo, num_compounds)
            compound_targets = (
                db_session.query(
                    target_combinations.c.target_combo,
                    func.count().label("num_compounds"),
                    func.group_concat(Compound.chembl_id, ", ").label("chembl_ids"),
                )
                .join(Compound, Compound.id == target_combinations.c.compound_id)
                .group_by(target_combinations.c.target_combo)
                .order_by(target_combinations.c.target_combo)
                .all()
            )

            n_compound_by_target = 0
            logger.info("1. Distinct compound target combinations and their counts:")
            for target_combo, count, chembl_ids in compound_targets:
                logger.info(f"    {target_combo}: {count} (Compounds: {chembl_ids})")
                n_compound_by_target += count
            logger.info(
                f"    Total compounds counted by target combinations: {n_compound_by_target}"
            )

            # Query compounds grouped by type and ordered by ascending number of Rule of 5 violations
            compounds_by_ro5 = (
                db_session.query(
                    target_combinations.c.target_combo,
                    Compound.chembl_id,
                    Compound.pref_name,
                    Compound.num_ro5,
                )
                .join(Compound, Compound.id == target_combinations.c.compound_id)
                .order_by(
                    target_combinations.c.target_combo,
                    Compound.num_ro5,
                )
                .all()
            )
            logger.info(
                "2. Compounds grouped by target combination and ordered by descending number of Rule of 5 violations:"
            )
            current_target_combo = ""
            for target_combo, chembl_id, pref_name, num_ro5 in compounds_by_ro5:
                if target_combo != current_target_combo:
                    current_target_combo = target_combo
                    logger.info(f"    Target combination: {current_target_combo}")
                    logger.info("        Rule of 5 violation count")

                logger.info(f"        {num_ro5} for {pref_name} ({chembl_id})")

    return (run_queries,)


@app.cell
def _(run_queries):
    run_queries()
    return


if __name__ == "__main__":
    app.run()
