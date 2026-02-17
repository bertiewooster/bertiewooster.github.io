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
    ChEMBL has the connections as compounds ↔ activities ↔ targets.
    """)
    return


@app.cell
def _():
    from graphviz import Digraph
    from IPython.display import SVG, display
    from sqlalchemy_schemadisplay import create_schema_graph
    # from sqlalchemy import Table
    import pydot

    return Digraph, SVG, create_schema_graph, display, pydot


@app.cell
def _():
    # Scale down size of labels and arrowheads
    crow_fontsize = "8"   # smaller font for "1" labels
    arrow_scale = "0.4"   # smaller arrowheads
    return arrow_scale, crow_fontsize


@app.cell
def _(Digraph, SVG, arrow_scale, crow_fontsize, display):
    # Create the diagram
    dot_no_activity = Digraph(format="svg")
    dot_no_activity.attr(rankdir="LR", splines="ortho")
    dot_no_activity.attr('node', shape='record', fontsize='10', style='filled', fillcolor='lightgrey')

    # Nodes
    dot_no_activity.node("Compound", "{Compound|compound_id (PK)}", fillcolor="#A3C1DA")
    dot_no_activity.node("Assay", "{Assay|assay_id (PK)}", fillcolor="#A3C1DA")
    dot_no_activity.node("Target", "{Target|target_id (PK)}", fillcolor="#A3C1DA")

    # Invisible edges for layout
    dot_no_activity.edge("Compound", "Assay", style="invis")
    dot_no_activity.edge("Assay", "Target", style="invis")

    dot_no_activity.edge("Compound", "Assay", fontsize=crow_fontsize,
             arrowhead="crow", arrowtail="crow", dir="both", color="black", arrowsize=arrow_scale)
    dot_no_activity.edge("Assay", "Target", fontsize=crow_fontsize,
             arrowhead="crow", arrowtail="crow", dir="both", color="black", arrowsize=arrow_scale)

    # Render inline SVG (scalable, no file saved)
    svg_content_no_activity = dot_no_activity.pipe(format="svg").decode("utf-8")
    display(SVG(svg_content_no_activity))
    return


@app.cell
def _(mo):
    mo.md(r"""
    That makes sense as a comprehensive schema; here, I wanted to simplify it by connecting molecules to targets more directly. So my schema is simply Compound ↔ Target.
    """)
    return


@app.cell
def _(Digraph, SVG, arrow_scale, crow_fontsize, display):
    # Create the diagram
    dot_simple_no_join = Digraph(format="svg")
    dot_simple_no_join.attr(rankdir="LR", splines="ortho")
    dot_simple_no_join.attr('node', shape='record', fontsize='10', style='filled', fillcolor='lightgrey')

    # Define nodes
    dot_simple_no_join.node("Compound", "{Compound|compound_id (PK)}", fillcolor="#A3C1DA")
    dot_simple_no_join.node("Target", "{Target|target_id (PK)}", fillcolor="#A3C1DA")

    # Invisible edge to force left-to-right ordering
    dot_simple_no_join.edge("Compound", "Target", style="invis")

    # Many-to-many edge (crow's foot at both ends)
    dot_simple_no_join.edge(
        "Compound",
        "Target",
        fontsize=crow_fontsize,
        arrowhead="crow",
        arrowtail="crow",
        dir="both",
        color="black",
        arrowsize=arrow_scale
    )

    # Render inline SVG (no file saved)
    svg_content_simple_no_join = dot_simple_no_join.pipe(format="svg").decode("utf-8")
    display(SVG(svg_content_simple_no_join))
    return


@app.cell
def _(mo):
    mo.md(r"""
    This is my first post using Marimo as the notebook: I had previously used just Jupyter. Initially Marimo wasn't great because

    - I discovered that publishing to a Markdown file directly from Marimo didn't lead to good formatting on my Jekyll, so I converted from Marimo to Jupyter and then to Markdown
    - I didn't have a reliable Internet connection when working on this and Marimo seemed to need a connection in VS Code
    - I sometimes am not allowed to rename a variable in VS Code (and automatically change the variable name wherever it's used)
    - Marimo is very serious about not allowing you to re-use a variable name, even an iterator variable which is usually a throw-away. So for example if you're trying out two versions of a code block, you have to rename every variable, even the iterators. While I understand the need for this, perhaps there's a way to make it easier by specifying a suffix to append to each variable name when you clone a code block--seems like something an LLM could handle.

    However, it seemed worth it when, before committing via git, the diff was so much more readable than in Jupyter (which is a ton of TypeScript, metadata, etc.). With Marimo, the diff is just the actual code changes and a small amount of formatting in Python. With Jupyter Notebooks, in theory source control works, but in practice it's so difficult to tell what changes were made that I didn't find it useful for identifying or rolling back changes.
    """)
    return


@app.cell
def _():
    import logging
    import time
    from collections import defaultdict

    import sqlalchemy
    from sqlalchemy.exc import IntegrityError, SQLAlchemyError
    from sqlalchemy.orm import declarative_base, sessionmaker
    from sqlalchemy import (
        Column,
        Float,
        Integer,
        String,
        create_engine,
        func,
        select,
    )
    from chembl_webresource_client.new_client import new_client

    return (
        Column,
        Float,
        Integer,
        IntegrityError,
        SQLAlchemyError,
        String,
        create_engine,
        declarative_base,
        defaultdict,
        func,
        logging,
        new_client,
        select,
        sessionmaker,
        sqlalchemy,
        time,
    )


@app.cell
def _(create_engine, declarative_base, logging, sessionmaker):
    # Set logging level to INFO
    logging.getLogger().setLevel(logging.INFO)

    Base = declarative_base()

    engine = create_engine("sqlite:///compounds.db", echo=False)
    Session = sessionmaker(bind=engine)
    return Base, Session, engine


@app.cell
def _(defaultdict, logging, new_client):
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

        # From the molecules in mols, extract the digits after "CHEMBL" to see which IDs were found
        chembl_ids_found = set()
        for mol in mols:
            chembl_id = mol.get("molecule_chembl_id", "")
            if chembl_id.startswith("CHEMBL"):
                chembl_ids_found.add(int(chembl_id.replace("CHEMBL", "")))
        chembl_ids_found_str = ", ".join(map(str, sorted(chembl_ids_found)))

        logging.info(
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
        logging.info(f"Activities: {len(activities)}")

        mol_to_target_ids = defaultdict(set)

        for a in activities:
            try:
                mol_to_target_ids[a["molecule_chembl_id"]].add(a["target_chembl_id"])
            except Exception as e:
                print(e)

        # ---------------------------------
        # 3) Fetch target metadata (bulk)
        # ---------------------------------
        all_target_ids = sorted(
            {tid for tids in mol_to_target_ids.values() for tid in tids}
        )

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
        logging.info(f"Fetched metadata for {len(targets)} targets.")

        # ---------------------------------
        # 4) Attach targets to molecules
        # ---------------------------------

        for m in mols:
            t_ids = mol_to_target_ids.get(m["molecule_chembl_id"], [])
            m["targets"] = [targets[tid] for tid in t_ids if tid in targets]

        return mols

    return (get_chembl_molecules,)


@app.cell
def _(
    Compound,
    CompoundTarget,
    IntegrityError,
    SQLAlchemyError,
    Session,
    Target,
    logging,
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

        n_mols_saved = 0
        n_targets_saved = 0
        n_compounds_targets_saved = 0

        with Session() as db_session:
            try:
                # preload existing targets into a mapping chembl_id -> Target instance
                existing_targets = {}
                if all_target_ids:
                    rows = (
                        db_session.query(Target)
                        .filter(Target.target_chembl_id.in_(list(all_target_ids)))
                        .all()
                    )
                    existing_targets = {r.target_chembl_id: r for r in rows}

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
                            if target_obj is None:
                                # create new Target, add and flush to get id, then cache it
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
                logging.info(f"IntegrityError while saving: {e}")
                db_session.rollback()
            except SQLAlchemyError:
                logging.exception("Database error saving compounds")
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
        rankdir='LR',
        concentrate=False
    )

    # Force strict left-to-right ordering with invisible edges;
    # add them programmatically by inspecting the SQLAlchemy model
    add_ordering_edges(graph_full, Base)

    graph_full.set_splines("ortho")

    # Move FK labels horizontally away from edges
    for edge_full in graph_full.get_edges():
        head_full = edge_full.get_headlabel()
        tail_full = edge_full.get_taillabel()

        if head_full:
            # Remove the "+ " prefix
            clean_head_full = head_full.replace("+ ", "").replace("+", "")
            edge_full.set_headlabel(clean_head_full)

        if tail_full:
            # Remove the "+ " prefix
            clean_tail_full = tail_full.replace("+ ", "").replace("+", "")
            edge_full.set_taillabel(clean_tail_full)
            edge_full.set_labeldistance("2.5")

    # Increase horizontal spacing between tables
    graph_full.set_ranksep("1.0")

    svg_content_full = graph_full.create_svg()
    display(SVG(svg_content_full))
    return


@app.cell
def _(
    Base,
    SVG,
    add_ordering_edges,
    create_schema_graph,
    display,
    engine,
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
    
        if len(pk_columns) == 1:
            return pk_columns[0]
        elif len(pk_columns) > 1:
            return ', '.join(pk_columns)  # Composite key
        else:
            return 'id'  # Fallback


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
                if edge.get_source() == join_table or edge.get_destination() == join_table:
                    graph.del_edge(edge.get_source(), edge.get_destination())
        
            # Get primary key names for both tables
            table1_pk = get_primary_key_name(Base, table1)
            table2_pk = get_primary_key_name(Base, table2)
        
            # Add a direct many-to-many edge between the two tables
            graph.add_edge(pydot.Edge(table1, table2, 
                                      taillabel=table1_pk, 
                                      headlabel=table2_pk,
                                      arrowhead='crow',
                                      arrowtail='crow',
                                      dir='both'))
    
        return set(join_tables.keys())

    # Create the ERD graph
    graph = create_schema_graph(
        engine=engine,
        metadata=Base.metadata,
        show_datatypes=True,
        show_indexes=False,
        rankdir='LR',
        concentrate=False
    )

    # Automatically detect and remove join tables
    excluded_tables = remove_join_tables_from_graph(graph, Base)

    # Add ordering edges (excluding detected join tables)
    add_ordering_edges(graph, Base, exclude_tables=excluded_tables)

    graph.set_splines("ortho")

    # Move FK labels horizontally away from edges
    for edge in graph.get_edges():
        head = edge.get_headlabel()
        tail = edge.get_taillabel()

        if head:
            edge.set_headlabel(head)

        if tail:
            edge.set_taillabel(tail)
            edge.set_labeldistance("2.5")

    graph.set_ranksep("1.0")

    svg_content = graph.create_svg()
    display(SVG(svg_content))
    return


@app.cell
def _(get_chembl_molecules, logging, save_compounds_to_db, time):
    # Measure how long it takes to fetch ChEMBL molecules
    start = time.time()
    mols = get_chembl_molecules(
        n_compounds=30,
        start_id=1000,  # Has targets
        # start_id=3430873, # Not a molecule
    )

    end = time.time()
    logging.info(f"Fetched {len(mols)} molecules in {end - start:.2f} seconds.")

    start = time.time()

    n_mols_saved, n_targets_saved, n_compounds_targets_saved = save_compounds_to_db(
        mols
    )
    logging.info(
        f"Saved {n_mols_saved} compounds, "
        f"{n_targets_saved} targets, "
        f"{n_compounds_targets_saved} compound-target associations to the database, "
        f"in {time.time() - start:.2f} seconds."
    )
    return


@app.cell
def _(Compound, CompoundTarget, Session, Target, func, logging, select):
    def run_queries():
        """Run the required queries against the ChEMBL database and print the results."""
        with Session() as db_session:
            # 1. Find the count of all distinct counts of compound targets (ex: Voltage-gated inwardly rectifying potassium channel KCNH2:14, Neuronal acetylcholine receptor subunit alpha-3/Neuronal acetylcholine receptor subunit alpha-7: 5).
            # For each compound, concatenate its targets with a slash.
            # Then, count how many distinct such tuples exist in the database.

            # Build per-compound target combo subquery:
            # as correlated scalar subquery: group_concat over an ordered selection of types to ensure consistent ordering.

            # correlated inner select returning pref_name for the current Compound, ordered
            inner = (
                select(Target.pref_name)
                .select_from(Target.__table__.join(CompoundTarget.__table__, CompoundTarget.target_id == Target.id))
                .where(CompoundTarget.compound_id == Compound.id)
                .order_by(func.lower(Target.pref_name))
                .correlate(Compound)
            )

            # name the derived table so the outer group_concat can select from it
            ordered_targets = inner.subquery("ordered_targets")

            # aggregate the ordered names with a '\' separator
            target_combo_subq = select(func.group_concat(ordered_targets.c.pref_name, "\\")).scalar_subquery()

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
            logging.info("1. Distinct compound target combinations and their counts:")
            for target_combo, count, chembl_ids in compound_targets:
                logging.info(f"    {target_combo}: {count} (Compounds: {chembl_ids})")
                n_compound_by_target += count
            logging.info(
                f"    Total compounds counted by target combinations: {n_compound_by_target}"
            )

            # Query compounds grouped by type and ordered by descending number of Rule of 5 violations
            compounds_by_ro5 = (
                db_session.query(
                    target_combinations.c.target_combo,
                    Compound.chembl_id,
                    Compound.pref_name,
                    Compound.num_ro5,
                    # Compound.molwt,
                    # Compound.tpsa,
                    # Compound.num_h_acceptors,
                    # Compound.num_h_donors,
                    # Compound.mol_logp,
                )
                .join(Compound, Compound.id == target_combinations.c.compound_id)
                .order_by(
                    target_combinations.c.target_combo,
                    Compound.num_ro5,
                )
                .all()
            )
            logging.info("2. Compounds grouped by target combination and ordered by descending number of Rule of 5 violations:")
            current_target_combo = ""
            logging.info("        Rule of 5 violation(s)")
            for target_combo, chembl_id, pref_name, num_ro5 in compounds_by_ro5:
                if target_combo != current_target_combo:
                    current_target_combo = target_combo
                    logging.info(f"    Target combination: {current_target_combo}")

                logging.info(f"        {num_ro5} for {pref_name} ({chembl_id})")

    return (run_queries,)


@app.cell
def _(run_queries):
    run_queries()
    return


if __name__ == "__main__":
    app.run()
