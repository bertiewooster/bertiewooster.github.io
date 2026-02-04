import marimo

__generated_with = "0.19.2"
app = marimo.App()


@app.cell
def _():
    import matplotlib.pyplot as plt
    return (plt,)


@app.cell
def _():
    a = 7
    return (a,)


@app.cell
def _(a):
    b = a + 3
    b
    return (b,)


@app.cell
def _(a, b, plt):
    # Create the plot
    fig, ax = plt.subplots()
    ax.scatter(a, b)
    ax.set_xlabel("a")
    ax.set_ylabel("b")
    ax.set_title("Plot of b against a")

    fig
    return


@app.cell
def _(create_engine, declarative_base, logging, sessionmaker):
    logging.getLogger().setLevel(logging.INFO)

    Base = declarative_base()

    engine = create_engine("sqlite:///compounds.db", echo=False)
    Session = sessionmaker(bind=engine)
    return (Session,)


@app.cell
def _(
    Compound,
    CompoundTarget,
    IntegrityError,
    SQLAlchemyError,
    Session,
    Target,
    defaultdict,
    logging,
    new_client,
):
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
        acts = activity.filter(
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
        logging.info(f"Activities: {len(acts)}")

        mol_to_target_ids = defaultdict(set)

        for a in acts:
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
    return


@app.cell
def _():
    return


if __name__ == "__main__":
    app.run()
