import sqlite3
import dataclasses
from datetime import datetime

import typing
from typing import Set, Type, Any, TypeVar, List, Literal, Union, Dict, NoReturn, Optional
from types import TracebackType

Q = TypeVar("Q")

StoreID = typing.NewType("StoreID", int)

Operator = Literal['<', '<=', '=', '!=', '>=', '>', 'and', 'or']
JoinHow = Literal['inner', 'left']
Order = Literal['asc', 'desc']


def assert_never(value: NoReturn) -> NoReturn:
    assert False, f'Unhandled value: {value} ({type(value).__name__})'


class Column:

    def __init__(self, name: str):
        self.name = name

    def __lt__(self, other: "Column") -> 'Condition':
        return Condition(self, other, "<")

    def __le__(self, other: "Column") -> 'Condition':
        return Condition(self, other, "<=")

    def __eq__(self, other: Union["Column", Any]) -> Union[bool, 'Condition']: #type: ignore
        return Condition(self, other, "=")

    def __ne__(self, other: Union["Column", Any]) -> Union[bool, 'Condition']: #type: ignore
        return Condition(self, other, "!=")

    def __ge__(self, other: "Column") -> 'Condition':
        return Condition(self, other, ">=")

    def __gt__(self, other: "Column") -> 'Condition':
        return Condition(self, other, ">")


class Condition:

    def __init__(self, c1: Column, c2: Union[Any, Column], operator: Operator):
        self.c1 = c1
        self.c2 = c2
        self.operator = operator

    def __and__(self, other: "Condition") -> 'AndCondition':
        return AndCondition([self, other])

    def __or__(self, other: "Condition") -> 'OrCondition':
        return OrCondition([self, other])

    def __invert__(self) -> 'Condition':
        return ICondition(self)

    def __str__(self) -> str:
        return f"Condition({self.c1} {self.operator} {self.c2})"


class AndCondition(Condition):

    def __init__(self, cs: List[Condition]):
        self.cs = cs

    def __and__(self, other: Condition) -> 'AndCondition':
        return AndCondition(self.cs+[other])


class OrCondition(Condition):

    def __init__(self, cs: List[Condition]):
        self.cs = cs

    def __or__(self, other: Condition) -> 'OrCondition':
        return OrCondition(self.cs+[other])


class ICondition(Condition):

    def __init__(self, c: Condition):
        self.c = c

    def __invert__(self) -> Condition:
        return self.c


class Table:

    def __init__(self, T: type[Q]):
        self.T = T
        self._name = T.__name__

    @property
    def name(self) -> str:
        return self._name

    def join(self, other: 'Table', on: Column, how: JoinHow="inner") -> 'Join':
        return Join(self, other, on, on, how)

    def where(self, condition: Condition) -> 'ConditionedTable':
        return ConditionedTable(self, condition)

    def order_by(self, col: Column, order: Order="asc") -> 'OrderedTable':
        return OrderedTable(self, col, order)


class Join(Table):

    def __init__(self, t1: Table, t2: Table, on1: Column, on2: Column, how: JoinHow):
        self.t1 = t1
        self.t2 = t2
        self.on1 = on1
        self.on2 = on2
        self.how = how

    @property
    def name(self) -> str:
        return self.t2.name

    @property
    def inner_name(self) -> str:
        return self.t1.name

class ConditionedTable(Table):

    def __init__(self, table: Table, condition: Condition):
        self.table = table
        self.condition = condition

    @property
    def name(self) -> str:
        return self.table.name

    def where(self, condition: Condition) -> 'ConditionedTable':
        return ConditionedTable(self.table, self.condition & condition)


class OrderedTable(Table):

    def __init__(self, table: Table, column: Column, order: Order):
        self.table = table
        self.column = column
        self.order = order

    @property
    def name(self) -> str:
        return self.table.name


col = Column
table = Table



class Store:

    def __init__(self, path: str):
        self._path = path
        self._tables: Set[str] = set()
        self._types: Set[type] = set()

        self._initialized = False
        self._in_transaction = 0

    def __enter__(self) -> "Store":
        self._begin()
        return self

    def __exit__(self, exc_type: type[BaseException], exc_value: BaseException, exc_traceback: TracebackType) -> None:
        if self._in_transaction:
            if exc_type is None:
                self._commit()
            else:
                self._rollback()

    def _register_dataclass(self, T: type) -> None:
        self._types.add(T)

    def _get_primary(self, T: type[Q]) -> str:
        for field in dataclasses.fields(T):
            if field.type == StoreID:
                return field.name
            if field.metadata.get("store.type", None) == "primary":
                return field.name
        return "rowid"

    def _create_table_query(self, T: type[Q]) -> str:
        create = "CREATE TABLE IF NOT EXISTS"
        tablename = T.__name__
        fields = []
        for field in dataclasses.fields(T):
            python_field_type = field.type

            optional = False
            if typing.get_origin(python_field_type) == typing.Union:
                a1, a2 = typing.get_args(python_field_type)
                if a1 is None:
                    a1, a2 = a2, a1
                python_field_type = a1
                optional = True

            field_name = field.name
            if python_field_type == str:
                field_type = "TEXT"
            elif python_field_type == bytes:
                field_type = "BLOB"
            elif python_field_type == datetime:
                field_type = "timestamp"
            elif python_field_type == int:
                field_type = "INTEGER"
            elif python_field_type == float:
                field_type = "REAL"
            elif python_field_type == bool:
                field_type = "BOOL"
            elif python_field_type == StoreID:
                field_type = "INTEGER PRIMARY KEY"
            else:
                raise Exception(f"Can not detect mapping of python type {python_field_type} to sqlite type")

            if optional and python_field_type != StoreID:
                field_type = field_type + " NOT NULL"

            if field.metadata.get("store.type", None) == "primary":
                if field.type != StoreID:
                    field_type = field_type + " PRIMARY KEY"
            fields.append(field_name + ' ' + field_type)
        field_statement = ", ".join(fields)
        return f'{create} {tablename}({field_statement});'

    def _begin(self) -> None:
        if self._in_transaction == 0:
            self._connection.execute("BEGIN TRANSACTION;")
        self._in_transaction += 1

    def _commit(self) -> None:
        self._in_transaction -= 1
        if self._in_transaction == 0:
            self._connection.execute("COMMIT;")

    def _rollback(self) -> None:
        self._connection.execute("ROLLBACK;")
        self._in_transaction = 0

    def create_db(self) -> None:
        sqlite3.register_converter("BOOL", lambda v: v == b"1")
        self._connection = sqlite3.connect(self._path, detect_types=sqlite3.PARSE_DECLTYPES, isolation_level=None)
        cur = self._connection.cursor()
        cur.execute("pragma journal_mode = WAL;")
        cur.execute("pragma synchronous = normal;")
        cur.execute("pragma temp_store = memory;")
        cur.execute("pragma mmap_size = 30000000000;")
        # cur.execute("pragma page_size = 32768;")

        with self:
            for dc in self._types:
                q = self._create_table_query(dc)
                cur.execute(q)
            for dc in self._types:
                self._tables.add(dc.__name__)

        self._initialized = True

    def create_table(self, dc: type[Q]) -> None:
        q = self._create_table_query(dc)
        cur = self._connection.cursor()
        with self:
            cur.execute(q)
        self._tables.add(dc.__name__)

    def insert(self, t: List[Q], update: bool=False) -> None:
        if not isinstance(t, list):
            t = [t]
        t0 = t[0]
        T = type(t0)
        tablename = T.__name__

        if not self._initialized:
            self.create_db()

        if tablename not in self._tables:
            self.create_table(T)

        aliases = [field.name for field in dataclasses.fields(type(t0))]
        placeholders = ["?" for field in dataclasses.fields(type(t0))]

        q = f"INSERT INTO {tablename} ({', '.join(aliases)}) VALUES ({', '.join(placeholders)})"

        if update:
            primary = self._get_primary(T)
            upsert_clause_l = [f" ON CONFLICT({primary}) DO UPDATE SET"]

            upsert_clause_l_fields = []
            for field in dataclasses.fields(T):
                upsert_clause_l_fields.append(f"{field.name}=?")
            upsert_clause_l.append(", ".join(upsert_clause_l_fields))
            upsert_clause = " ".join(upsert_clause_l)

            q += upsert_clause
            fillers = []
            for instance in t:
                tup = dataclasses.astuple(instance)
                fillers.append(tup+tup)
        else:
            fillers = [dataclasses.astuple(instance) for instance in t]

        q += ";"

        cur = self._connection.cursor()
        with self:
            cur.executemany(q, fillers)

    def get_single(self, T: Type[Q], **kwargs: Any) -> Optional[Q]:
        ts = self.get(T, **kwargs)
        if len(ts) == 1:
            return ts[0]
        elif len(ts) == 0:
            return None
        else:
            raise KeyError("Multiple values for query")

    def get(self, T: Type[Q], **kwargs: Any) -> List[Q]:
        fields = [field.name for field in dataclasses.fields(T)]
        tablename = T.__name__

        if tablename not in self._tables:
            return []

        if not self._initialized:
            self.create_db()

        where_clause = " and ".join([f"{key} = ?" for key in kwargs])
        if where_clause:
            where_clause = "WHERE " + where_clause
        else:
            where_clause = ''
        select_clause = f"SELECT {', '.join(fields)} from {tablename} "

        q = ' '.join([select_clause, where_clause, ';'])

        cur = self._connection.cursor()
        cur.execute(q, tuple([value for value in kwargs.values()]))
        rows = cur.fetchall()

        result = [T(*args) for args in rows]
        return result

    def _get_table_postfix(self, table_counter: Dict[str, int], name: str, increment: bool = True) -> str:
        table_count = table_counter.get(name, 0)
        if increment:
            table_counter[name] = table_count + 1
        if table_count == 0:
            table_postfix = ''
        else:
            table_postfix = str(table_count)
        return table_postfix

    def _get_value_as_str(self, tableid: str, value: Union[str, int, Column], param_container: List[Union[str, int]]) -> str:
        if isinstance(value, str) or isinstance(value, int) or isinstance(value, bytes) or isinstance(value, bool):
            param_container.append(value)
            return "?"
        elif isinstance(value, Column):
            return f"{tableid}.{value.name}"
        else:
            raise ValueError(f"Cannot deal with value of type {type(value)}")

    def _get_cond_str(self, tableid: str, cond: Condition, param_container: List[Union[str, int]]) -> str:
        if isinstance(cond, OrCondition):
            return "(" + " or ".join([self._get_cond_str(tableid, c, param_container) for c in cond.cs]) + ")"
        elif isinstance(cond, AndCondition):
            return "(" + " and ".join([self._get_cond_str(tableid, c, param_container) for c in cond.cs]) + ")"
        elif isinstance(cond, ICondition):
            return "NOT " + self._get_cond_str(tableid, cond.c, param_container)
        elif isinstance(cond, Condition):
            left = self._get_value_as_str(tableid, cond.c1, param_container)
            if cond.c2 is None:
                return f"{left} IS NULL"
            right = self._get_value_as_str(tableid, cond.c2, param_container)
            return f"{left} {cond.operator} {right}"
        else:
            assert_never(cond)

    def select(self, query: Union[Table, ConditionedTable, OrderedTable, Join], limit: Optional[int]=None) -> List[Q]:

        table_clauses: List[str] = []
        where_clauses: List[str] = []
        order_by = ""

        param_container: List[Union[str, int]] = []
        table_counter: Dict[str, int] = {}

        table: Optional[Table] = query

        while table is not None:
            if isinstance(table, OrderedTable):
                assert order_by == ""

                table_name = table.name
                table_postfix = self._get_table_postfix(table_counter, table_name, increment=False)
                table_id = table_name + table_postfix

                column = table.column.name if isinstance(table.column, Column) else table.column
                order_by = f'ORDER BY {table_id}.{column}' + (' ASC' if table.order == 'asc' else ' DESC')
                table = table.table
            elif isinstance(table, ConditionedTable):
                condition = table.condition
                table = table.table

                table_name = table.name
                table_postfix = self._get_table_postfix(table_counter, table_name, increment=False)
                table_id = table_name + table_postfix

                where_clauses.append(self._get_cond_str(table_id, condition, param_container))

            elif isinstance(table, Join):
                table_name = table.t2.name
                table_postfix = self._get_table_postfix(table_counter, table.t2.name)
                table_id = table_name + table_postfix

                inner_table_name = table.inner_name
                inner_table_postfix = self._get_table_postfix(table_counter, inner_table_name, increment=False)
                inner_table_id = inner_table_name + inner_table_postfix

                if table.how == 'inner':
                    join_op = 'INNER JOIN'
                elif table.how == 'left':
                    join_op = 'LEFT JOIN'
                else:
                    assert_never(table.how)

                table_clauses.insert(0, f'{join_op} {table_name} {table_id} ON {inner_table_id}.{table.on1} = {table_id}.{table.on2}')

                table = table.t1
            elif isinstance(table, Table):
                T = table.T

                tablename = T.__name__

                if tablename not in self._tables:
                    self.create_table(T)

                table_postfix = self._get_table_postfix(table_counter, tablename)
                table_id = tablename + table_postfix

                fields = [f'{table_id}.{field.name}' for field in dataclasses.fields(T)]

                table_clauses.insert(0, f"SELECT {', '.join(fields)} from {tablename}{table_postfix}")
                table = None
            else:
                assert_never(table)

        select_clause = " ".join(table_clauses)
        where_clause = ("WHERE " + " and ".join(where_clauses) ) if where_clauses else ''

        limit_clause = f'LIMIT {limit}' if limit else ''

        clauses = [select_clause, where_clause, order_by, limit_clause]

        q = ' '.join([clause for clause in clauses if clause])

        cur = self._connection.cursor()
        cur.execute(q, param_container)
        rows = cur.fetchall()

        result = [T(*args) for args in rows]
        return result




