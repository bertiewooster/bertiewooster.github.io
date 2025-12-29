import sqlalchemy
from sqlalchemy import Column, Integer, String, create_engine
from sqlalchemy.orm import declarative_base, sessionmaker

Base = declarative_base()

engine = create_engine("sqlite:///compounds.db", echo=False)
Session = sessionmaker(bind=engine)


def init_db():
    """Create database tables (call once at app startup or from scripts/tests)."""
    Base.metadata.create_all(engine)


def reset_db():
    """Drop all database tables (use with caution)."""
    Base.metadata.drop_all(engine)


class Pokemon(Base):
    __tablename__ = "pokemon"

    id = Column(Integer, primary_key=True)
    poke_num = Column(Integer, unique=True)
    name = Column(String)
    hp = Column(Integer)
    attack = Column(Integer)
    defense = Column(Integer)
    speed = Column(Integer)
    special_attack = Column(Integer)
    special_defense = Column(Integer)

    # Relationship to type table:
    types = sqlalchemy.orm.relationship("Type", back_populates="pokemon")


class Type(Base):
    """Type-to-Pokemon mapping table."""

    __tablename__ = "type"

    id = Column(Integer, primary_key=True)
    pokemon_id = Column(Integer, sqlalchemy.ForeignKey("pokemon.id"))
    type_name = Column(String)

    # Relationship back to Pokemon
    pokemon = sqlalchemy.orm.relationship("Pokemon", back_populates="types")
