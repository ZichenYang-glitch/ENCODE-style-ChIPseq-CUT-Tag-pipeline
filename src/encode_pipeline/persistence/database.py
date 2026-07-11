"""SQLAlchemy engine and session composition."""

from __future__ import annotations

from pathlib import Path

from sqlalchemy import Engine, create_engine, event
from sqlalchemy.orm import Session, sessionmaker


def create_database_engine(database_url: str, *, echo: bool = False) -> Engine:
    """Create a SQLAlchemy engine with safe SQLite connection defaults."""
    if not isinstance(database_url, str) or not database_url.strip():
        raise ValueError("database_url must be a non-empty string")

    normalized_url = database_url.strip()
    connect_args: dict[str, object] = {}
    if normalized_url.startswith("sqlite"):
        connect_args = {"check_same_thread": False, "timeout": 30}
        _ensure_sqlite_parent(normalized_url)

    engine = create_engine(
        normalized_url,
        echo=echo,
        future=True,
        connect_args=connect_args,
    )
    if engine.dialect.name == "sqlite":
        _configure_sqlite(engine)
    return engine


def create_session_factory(engine: Engine) -> sessionmaker[Session]:
    """Return a session factory whose domain conversions survive commits."""
    if not isinstance(engine, Engine):
        raise ValueError("engine must be a SQLAlchemy Engine")
    return sessionmaker(bind=engine, class_=Session, expire_on_commit=False)


def _configure_sqlite(engine: Engine) -> None:
    @event.listens_for(engine, "connect")
    def set_sqlite_pragmas(dbapi_connection, _connection_record) -> None:
        cursor = dbapi_connection.cursor()
        cursor.execute("PRAGMA foreign_keys=ON")
        cursor.execute("PRAGMA busy_timeout=30000")
        if _is_file_backed_sqlite(engine):
            cursor.execute("PRAGMA journal_mode=WAL")
        cursor.close()


def _is_file_backed_sqlite(engine: Engine) -> bool:
    database = engine.url.database
    return database not in {None, "", ":memory:"}


def _ensure_sqlite_parent(database_url: str) -> None:
    if ":memory:" in database_url:
        return
    prefix = "sqlite:///"
    if not database_url.startswith(prefix):
        return
    raw_path = database_url[len(prefix) :].split("?", maxsplit=1)[0]
    if not raw_path:
        return
    path = Path(raw_path).expanduser()
    path.parent.mkdir(parents=True, exist_ok=True)
