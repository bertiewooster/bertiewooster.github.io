import asyncio
import concurrent.futures
import json
import logging
import time
import traceback
from collections import defaultdict
from typing import Optional

import aiohttp
import sqlalchemy

# import database
from aiohttp import ClientSession
from chembl_webresource_client.new_client import new_client
from chembl_webresource_client.utils import utils
from sqlalchemy import Column, Float, Integer, String, create_engine

# from queries import run_queries
from sqlalchemy.exc import IntegrityError, SQLAlchemyError
from sqlalchemy.orm import declarative_base, sessionmaker

# Set logging level to INFO
logging.getLogger().setLevel(logging.INFO)

Base = declarative_base()

engine = create_engine("sqlite:///compounds.db", echo=False)
Session = sessionmaker(bind=engine)


def get_chembl_molecules(
    n_compounds: int = 2,
    start_id: int = 1,
):
    fn_start = time.time()
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

    logging.info(
        f"Of the {n_compounds} ChEMBL IDs "
        f"({start_id}-{start_id + n_compounds - 1}), {len(mols)} are compounds."
    )

    if not mols:
        return []

    mol_ids_present = [m["molecule_chembl_id"] for m in mols]

    # ---------------------------------
    # 2) Bulk fetch activities → targets
    # ---------------------------------
    start = time.time()
    acts = activity.filter(
        molecule_chembl_id__in=mol_ids_present,
        target_organism="Homo sapiens",
        standard_type="IC50",
        assay_type="B",
        # document_year__gt=2010,
        # document_year__lt=1990, # 16 activities
        # document_year__lt=2000, # 46 activities
        document_year__lt=2010,  # 284 activities
    ).only(
        [
            "molecule_chembl_id",
            "target_chembl_id",
        ]
    )
    end = time.time()
    logging.info(
        f"For {len(acts)} activities, bulk fetch activities → targets took {end - start} seconds."
    )

    mol_to_target_ids = defaultdict(set)

    start = time.time()
    for a in acts:
        try:
            mol_to_target_ids[a["molecule_chembl_id"]].add(a["target_chembl_id"])
        except Exception as e:
            print(e)

    end = time.time()
    logging.info(
        f"For {len(acts)} activities, setting mol_to_target_ids with 'try' took {end - start} seconds."
    )

    # ---------------------------------
    # 3) Fetch target metadata (bulk)
    # ---------------------------------
    all_target_ids = sorted(
        {tid for tids in mol_to_target_ids.values() for tid in tids}
    )

    targets = {}
    if all_target_ids:
        start = time.time()

        for t in target.filter(target_chembl_id__in=all_target_ids).only(
            [
                "target_chembl_id",
                "pref_name",
                "target_type",
                "organism",
            ]
        ):
            targets[t["target_chembl_id"]] = t
    end = time.time()
    logging.info(
        f"For {len(targets)} targets, fetch target metadata took {end - start} seconds."
    )

    # ---------------------------------
    # 4) Attach targets to molecules
    # ---------------------------------
    start = time.time()

    for m in mols:
        t_ids = mol_to_target_ids.get(m["molecule_chembl_id"], [])
        m["targets"] = [targets[tid] for tid in t_ids if tid in targets]
    end = time.time()
    logging.info(
        f"Attach targets to molecules target metadata took {end - start} seconds."
    )

    logging.info(
        f"Total time for get_chembl_molecules: {time.time() - fn_start} seconds."
    )

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




def init_db():
    """Create database tables (call once at app startup or from scripts/tests)."""
    Base.metadata.create_all(engine)


def reset_db():
    """Drop all database tables (use with caution)."""
    Base.metadata.drop_all(engine)


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


class Target(Base):
    __tablename__ = "target"

    id = Column(Integer, primary_key=True)
    organism = Column(String)
    pref_name = Column(String)
    target_chembl_id = Column(String, unique=True)
    target_type = Column(String)


# Join table
class CompoundTarget(Base):
    __tablename__ = "compound_target"

    id = Column(Integer, primary_key=True)
    compound_id = Column(Integer, sqlalchemy.ForeignKey("compound.id"))
    target_id = Column(Integer, sqlalchemy.ForeignKey("target.id"))


def run_queries():
    """Run the required queries against the Pokemon database and print the results."""
    with Session() as db_session:
        # Create a list of distinct type combinations and their counts
        # where each is a tuple like (type_combo, num_pokemons)
        poke_types = (
            db_session.query(
                type_combinations.c.type_combo, func.count().label("num_pokemons")
            )
            .group_by(type_combinations.c.type_combo)
            .order_by(type_combinations.c.type_combo)
            .all()
        )

        n_poke_by_type = 0
        logging.info("1. Distinct Pokemon type combinations and their counts:")
        for type_combo, count in poke_types:
            logging.info(f"    {type_combo}: {count}")
            n_poke_by_type += count
        logging.info(
            f"    Total Pokemon counted by type combinations: {n_poke_by_type}"
        )


if __name__ == "__main__":
    # Reset database (uncomment to start fresh)
    reset_db()

    # Ensure tables exist
    init_db()

    # Measure how long it takes to fetch ChEMBL molecules
    start = time.time()
    mols = get_chembl_molecules(
        n_compounds=10,
        start_id=1,  # Has targets
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
        f"{n_compounds_targets_saved} compound-target associations to the database,"
        f"in {time.time() - start:.2f} seconds."
    )

    pass
