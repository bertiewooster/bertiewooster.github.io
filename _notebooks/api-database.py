import asyncio
import concurrent.futures
import json
import logging
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

Base = declarative_base()

engine = create_engine("sqlite:///compounds.db", echo=False)
Session = sessionmaker(bind=engine)

def get_chembl_molecules(
    n_compounds: int = 2,
):
    chembl_ids = [f'CHEMBL{id}' for id in range(1, n_compounds + 1)]
    molecule = new_client.molecule
    mols = molecule.filter(molecule_chembl_id__in=chembl_ids).only(['molecule_chembl_id', 'molecule_structures', 'pref_name'])
    for index, mol in enumerate(mols):
        canon_sml = mol.get('molecule_structures').get('canonical_smiles')
        mol_ctab = utils.smiles2ctab(canon_sml)
        descs = json.loads(utils.chemblDescriptors(mol_ctab))[0]
        mol['descs'] = descs
    return mols

async def get_pubchem_api(
    n_compounds: int = 2,
    executor: Optional[concurrent.futures.ThreadPoolExecutor] = None,
) -> tuple[int, Optional[Exception]]:
    """Fetch compound data via the PubChem API and store it in the database.
    Args:
        n_compounds (int): Number of compounds to fetch (from poke_num 1 to n_compounds).
    Returns:
        int: The number of compounds fetched.
        Optional[Exception]: An exception if one occurred, else None.
    """
    # Create queue for single-writer executor to serialize all DB writes
    # safe for SQLite as a file-based database
    queue: asyncio.Queue = asyncio.Queue()
    # executor must be provided by the synchronous caller (see main). this ensures
    # the ThreadPoolExecutor is shutdown deterministically in the main thread
    # rather than during interpreter teardown (which triggered the __del__ error).
    if executor is None:
        raise RuntimeError(
            "get_pubchem_api requires a ThreadPoolExecutor instance (pass executor=...)"
        )
    try:

        async def _db_writer(
            q: asyncio.Queue, exec: concurrent.futures.ThreadPoolExecutor
        ):
            loop = asyncio.get_running_loop()
            while True:
                item = await q.get()
                if item is None:
                    q.task_done()
                    break
                cid, name, type_names, stats_map = item
                try:
                    # run the blocking DB write in the single-thread executor (serializes writes)
                    await loop.run_in_executor(
                        exec, save_compound_to_db, cid, name, type_names, stats_map
                    )
                except Exception as e:
                    # log and continue; writer keeps processing other items
                    logging.exception(f"DB writer error for PubChem ID {cid}: {e}")
                    traceback.print_exception(type(e), e, e.__traceback__)
                finally:
                    q.task_done()

        # Limit concurrent connections to avoid overloading the API/server
        connector = aiohttp.TCPConnector(limit=5)
        client = aiohttp.ClientSession(connector=connector)
        async with client as asyncio_session:
            writer_task = asyncio.create_task(_db_writer(queue, executor))
            try:
                async with asyncio.TaskGroup() as tg:
                    for id in range(1, n_compounds + 1):
                        tg.create_task(
                            fetch_compound_data(
                                asyncio_session=asyncio_session,
                                id=id,
                                write_queue=queue,
                            )
                        )
                # all fetch tasks finished — wait for queued DB writes to complete
                await queue.join()
                # signal writer to stop
                await queue.put(None)
                await writer_task
                return n_compounds, None
            except ExceptionGroup as eg:
                # Compute and attach successful task count to the ExceptionGroup,
                # print full tracebacks for debugging, ensure writer exits, then re-raise.
                successes = n_compounds - len(eg.exceptions)
                setattr(eg, "successful_tasks", successes)
                logging.exception(
                    f"TaskGroup raised ExceptionGroup with {len(eg.exceptions)} sub-exception(s); "
                    f"{successes} succeeded"
                )
                for i, sub in enumerate(eg.exceptions, 1):
                    logging.exception(
                        f"\n--- Sub-exception #{i}: {type(sub).__name__}: {sub} ---"
                    )
                    traceback.print_exception(type(sub), sub, sub.__traceback__)
                # ensure writer exits
                await queue.put(None)
                await writer_task
                # Re-raise so callers get the ExceptionGroup (inspect eg.successful_tasks)
                raise
    finally:
        # nothing to shutdown here; executor lifetime is owned by caller (main)
        pass


async def fetch_compound_data(
    asyncio_session,
    #: ClientSession, 
    cid: int, 
    write_queue: asyncio.Queue
) -> Optional[bool]:
    """Fetch data for a single compound by ID from the PubChem API and enqueue it for DB write.
    Args:
        asyncio_session: An aiohttp ClientSession to use for the request.
        id (int): The Pokemon ID to fetch.
        write_queue: An asyncio.Queue to place the prepared data for the DB writer.
    Returns:
        Optional[bool]: True if saved successfully, False if database integrity error, None if 404
    """
    # Get the URL for this Pokemon
    pubchem_URL = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound/{cid}/XML?heading=Chemical-Vendors"

    max_attempts = 5
    backoff = 1
    last_status = None

    try:
        for attempt in range(1, max_attempts + 1):
            async with asyncio_session.get(pubchem_URL) as response:
                last_status = response.status

                if response.status == 200:
                    get_poke_json = await response.json()
                    break

                # Treat 404 as non-fatal: log and skip this Pokemon
                if response.status == 404:
                    logging.warning(f"PubChem CID {cid} not found (404), skipping")
                    return None

                # Retry on rate limiting or server errors
                if response.status == 429 or 500 <= response.status < 600:
                    if attempt < max_attempts:
                        await asyncio.sleep(backoff)
                        backoff *= 2
                        continue

                # Non-retriable status or out of attempts
                raise ConnectionError(
                    f"Failed to fetch data for PubChem CID {cid}, status code: {response.status}"
                )
        else:
            # Exhausted attempts
            raise ConnectionError(
                f"Failed to fetch data for PubChem CID {cid} after {max_attempts} attempts, last status: {last_status}"
            )

    except asyncio.CancelledError:
        # Let cancellations propagate
        raise
    except aiohttp.ClientError as e:
        # Wrap network errors with context
        raise ConnectionError(f"Network error fetching PubChem CID {id}") from e

    # # Continue processing if we have JSON
    # name = get_poke_json.get("name")
    # types = get_poke_json.get("types", [])
    # type_names = [t.get("type", {}).get("name") for t in types if t.get("type")]

    # stats = get_poke_json.get("stats", [])
    # stats_map = {
    #     stat.get("stat", {}).get("name"): stat.get("base_stat")
    #     for stat in stats
    #     if stat.get("stat")
    # }

    #TODO convert SMILES into RDKit molecule, calculate descriptors for

    # enqueue the prepared data for the single DB writer to consume
    await write_queue.put((id, name, type_names, stats_map))
    return True


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
        logging.info(f"\nSuccessfully added {n_compounds_added} compounds to the database.")

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
    molwt = Column(Float) # MolWt
    tpsa = Column(Float) # TPSA
    num_h_acceptors = Column(Integer) # NumHAcceptors
    num_h_donors = Column(Integer) # NumHDonors
    num_ro5 = Column(Integer) # NumRo5
    mol_logp = Column(Float) # MolLogP

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
    result = get_chembl_molecules(1)
    pass
