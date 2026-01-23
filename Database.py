import sqlite3
from contextlib import closing


class Database:
    def __init__(self, db_path: str):
        """
        Initialize database connection.

        Example:
            db = Database("app.db")
        """
        self.db_path = db_path
        self.conn = sqlite3.connect(self.db_path)
        self.conn.row_factory = sqlite3.Row  # dict-like rows

    def execute(self, query: str, params: tuple = ()):
        """
        Execute INSERT / UPDATE / DELETE queries.
        """
        with closing(self.conn.cursor()) as cursor:
            cursor.execute(query, params)
            self.conn.commit()
            return cursor.rowcount

    def fetchone(self, query: str, params: tuple = ()):
        """
        Fetch a single row.
        """
        with closing(self.conn.cursor()) as cursor:
            cursor.execute(query, params)
            row = cursor.fetchone()
            return dict(row) if row else None

    def fetchall(self, query: str, params: tuple = ()):
        """
        Fetch all rows.
        """
        with closing(self.conn.cursor()) as cursor:
            cursor.execute(query, params)
            rows = cursor.fetchall()
            return [dict(row) for row in rows]

    def executemany(self, query: str, data: list[tuple]):
        """
        Execute many inserts at once.
        """
        with closing(self.conn.cursor()) as cursor:
            cursor.executemany(query, data)
            self.conn.commit()

    def close(self):
        """
        Close DB connection.
        """
        self.conn.close()
