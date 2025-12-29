import asyncio
import concurrent.futures
import json
import logging
import time
import traceback
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
    chembl_ids = [f"CHEMBL{id}" for id in range(start_id, start_id + n_compounds)]
    molecule = new_client.molecule
    
    # Suppress INFO logs from the chembl package while keeping global INFO level
    logging.getLogger("chembl_webresource_client").setLevel(logging.WARNING)

    mols = molecule.filter(molecule_chembl_id__in=chembl_ids).only(
        [
            "molecule_chembl_id",
            "molecule_structures",
            "pref_name",
            "molecule_properties",
        ]
    )

    logging.info(
        f"Of the {n_compounds} ChEMBL IDs ({start_id}-{start_id + n_compounds - 1}), {len(mols)} are compounds."
    )
    return mols


def save_compound_to_db(
    cid: int, name: str, type_names: list[str], stats_map: dict[str, int]
) -> bool:
    """Save a single compound and its descriptors to the database.
    Args:
        cid (int): The PubChem CID.
        name (str): The compound name.
        type_names (list of str): List of type names for the Pokemon. #TODO fix
        stats_map (dict): Mapping of stat names to their values. #TODO fix
    Returns:
        bool: True if saved successfully.
    """
    compound = database.Compound(
        cid=cid,
        name=name,
        hp=stats_map.get("hp"),
        attack=stats_map.get("attack"),
        defense=stats_map.get("defense"),
        speed=stats_map.get("speed"),
        special_attack=stats_map.get("special-attack"),
        special_defense=stats_map.get("special-defense"),
    )
    # Create and use a database session here so each async task has its own session, to prevent conflicts
    with database.Session() as db_session:
        # Use transaction context manager to ensure atomicity and auto-rollback on error
        try:
            with db_session.begin():
                db_session.add(compound)
                db_session.flush()
                for type_name in type_names:
                    pokemon_type = database.Type(
                        pokemon_id=compound.id, type_name=type_name
                    )
                    db_session.add(pokemon_type)
        except IntegrityError as e:
            logging.info(f"IntegrityError saving PubChem CID {cid}: {e}")
            return False
        except SQLAlchemyError:
            logging.exception(f"Database error saving PubChem CID {cid}")
            raise
    return True


def main():
    """Main entry point to initialize/reset DB, fetch compound data, compute descriptors, and run queries."""
    # configure logging in the application entrypoint
    logging.basicConfig(level=logging.INFO)
    # Ensure tables exist
    database.init_db()

    # Get Pokemon data from the PubChem API and store it in the database.
    # Create and own the ThreadPoolExecutor from the synchronous main function so it
    # is shut down deterministically before interpreter teardown.
    try:
        with concurrent.futures.ThreadPoolExecutor(max_workers=1) as executor:
            n_compounds_added, _ = asyncio.run(
                get_pubchem_api(
                    n_compounds=151,
                    executor=executor,
                )
            )
    except ExceptionGroup as eg:
        logging.exception(
            f"\nCompleted with ExceptionGroup: {len(eg.exceptions)} sub-exception(s); "
            f"{getattr(eg, 'successful_tasks', 'unknown')} succeeded"
        )
        for i, sub in enumerate(eg.exceptions, 1):
            logging.exception(
                f"\n--- Sub-exception #{i}: {type(sub).__name__}: {sub} ---"
            )
            traceback.print_exception(type(sub), sub, sub.__traceback__)
        return
    else:
        logging.info(
            f"\nSuccessfully added {n_compounds_added} compounds to the database."
        )

    # Run the queries from this script; can alternatively run from queries.py directly
    run_queries()
    pass


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
    target_chembl_id = Column(String, unique=True)
    target_type = Column(String)
    organism = Column(String)
    pref_name = Column(String)


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
    # main()
    # Measure how long it takes to fetch ChEMBL molecules
    start = time.time()
    result = get_chembl_molecules(
        n_compounds=1000,
        # start_id=3430873,
    )
    end = time.time()
    logging.info(f"Fetched {len(result)} molecules in {end - start:.2f} seconds.")
    pass
