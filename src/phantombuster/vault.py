from functools import wraps
from hashlib import sha3_256
import inspect
from datetime import datetime
from typing import Any
import uuid
import sqlite3

from dataclasses import dataclass, field
import dataclasses

from phantombuster import store
from phantombuster.store import Store, StoreID

import pickle
import time

from typing import Any, Optional, Union, Callable, Dict, List

@dataclass
class FunctionDefinition:
    name: str
    source: str
    function_hash: bytes = field(metadata={"store.type": "primary"})
    datetime: datetime


@dataclass
class FunctionCall:
    call_hash: bytes = field(metadata={"store.type": "primary"})
    function_hash: bytes
    kwargs: bytes
    datetime: datetime


@dataclass
class FunctionCallResult:
    call_hash: bytes
    result_type: bytes
    result: bytes
    computation_time: float
    time: datetime
    overrideable: bool = False



Chatter = Union[FunctionDefinition, FunctionCall, FunctionCallResult]


@dataclass(slots=True)
class FunctionCallProxy:
    call_hash: bytes


def hash_function_source(f):
    source = inspect.getsource(f)
    s = sha3_256(source.encode("utf-8"))
    hash = s.digest()
    return hash


class Vault:

    NO_RESULT = object()

    def __init__(self, storepath: Union[Store, str]):
        if isinstance(storepath, Store):
            self._store = storepath
        else:
            self._store = Store(storepath)
        self._store._register_dataclass(FunctionDefinition)
        self._store._register_dataclass(FunctionCall)
        self._store._register_dataclass(FunctionCallResult)
        self._store.create_db()
        self._function_store = {}
        self._store_function = False

    def consume(self, x: Chatter) -> None:
        if isinstance(x, FunctionDefinition):
            self._store.insert(x)
        elif isinstance(x, FunctionCall):
            self._store.insert(x)
        elif isinstance(x, FunctionCallResult):
            self._store.insert(x)
        else:
            raise ValueError(f'Can only handle values of type FunctionDefinition, FunctionCall or FunctionCallResult, not {type(x)}')

    def _get_store_def(self):
        with self._store:
            if not self._store_function:
                function_hash = b"vault.store"
                fdef = self.recall_definition_from_hash(function_hash)
                if fdef:
                    self._store_function = fdef
                else:
                    fdef = FunctionDefinition("vault.store", "vault internal hack", function_hash, datetime.now())
                    self._store.insert(fdef)
                    self._store_function = fdef
        return self._store_function

    def _calculate_store_variable_call(self, name):
        store_def = self._get_store_def()

        pickled_kwargs = pickle.dumps({"name": name}, protocol=5)
        call_hash = calculate_call_hash(store_def, pickled_kwargs)

        now = datetime.now()
        fcall = FunctionCall(call_hash=call_hash, function_hash=store_def.function_hash, kwargs=pickled_kwargs, datetime=now)
        return fcall

    def store_variable(self, name, value, final=False):
        with self._store:
            fcall = self._calculate_store_variable_call(name)
            print(f"Storing variable {name} with value {value} on callhash {fcall.call_hash}")
            if not self._store.get(FunctionCall, call_hash=fcall.call_hash):
                self._store.insert(fcall)
            self.save_result(fcall, value, 0.0, not final)

    def load_variable(self, name):
        fcall = self._calculate_store_variable_call(name)

        result = store.table(FunctionCallResult)
        correct_callhash = (store.col("call_hash") == fcall.call_hash)
        creation_time = store.col("time")

        results = self._store.select(result.where(correct_callhash).order_by(creation_time, "desc"), limit=1)
        if results:
            result = results[0]
            return self.load_result(result)
        else:
            return None

    def define(self, f: Callable) -> FunctionDefinition:
        source = inspect.getsource(f)
        function_hash = hash_function_source(f)
        fun_def = self._store.get(FunctionDefinition, function_hash=function_hash)
        if fun_def:
            fun_def = fun_def[0]
        else:
            name = f'{f.__module__}.{f.__name__}'
            fun_def = FunctionDefinition(name, source, function_hash, datetime.now())
            self._store.insert(fun_def)
        return fun_def

    def recall_definition_from_hash(self, function_hash: bytes) -> FunctionDefinition:
        return self._store.get_single(FunctionDefinition, function_hash=function_hash)

    def recall_definition_from_name(self, function_name: str) -> FunctionDefinition:
        return self._store.get(FunctionDefinition, name=function_name)

    def call(self, function_definition: FunctionDefinition, kwargs: Dict[str, Any]) -> FunctionCall:

        kwargs = proxy_function_calls(kwargs)
        pickled_kwargs = pickle.dumps(kwargs, protocol=5)

        call_hash = calculate_call_hash(function_definition, pickled_kwargs)

        fhash = function_definition.function_hash

        now = datetime.now()

        fcall = FunctionCall(call_hash=call_hash, function_hash=fhash, kwargs=pickled_kwargs, datetime=now)

        try:
            self._store.insert(fcall)
        except sqlite3.IntegrityError:
            pass
        return fcall

    def get_newest_unhandled_calls(self, limit: int=1) -> List[FunctionCall]:
        call = store.table(FunctionCall)
        result = store.table(FunctionCallResult)
        result_present = (store.col("result") == None)
        creation_time = store.col("time")

        calls = self._store.select(call.join(result, on="call_hash", how="left").where(result_present).order_by(creation_time, "desc"), limit=limit)
        return calls

    def get_result(self, function_call: Union[FunctionCall, FunctionCallProxy]):
        results = self._store.get(FunctionCallResult, call_hash=function_call.call_hash)
        results = [result for result in results if result.result_type != 'intent']
        if len(results) > 1:
            raise Exception("Got multiple results for the same function call, something seems wrong")
        elif len(results) == 1:
            return results[0]
        else:
            return Vault.NO_RESULT

    def save_result(self, function_call: FunctionCall, result: Any, time: float, overrideable: bool=False) -> FunctionCallResult:
        result_type = 'pickle'
        b = pickle.dumps(result, protocol=5)
        call_hash = function_call.call_hash
        fun_call_result = FunctionCallResult(call_hash, result_type, b, time, datetime.now(), overrideable)

        previous_results = self._store.get(FunctionCallResult, call_hash=call_hash)

        if len(previous_results) > 0:
            if previous_results[0].overrideable:
                self._store.insert(fun_call_result, update=True)
                return fun_call_result
            else:
                print("A result was already present, discarding the new one")
                return previous_results[0]
        else:
            self._store.insert(fun_call_result)
            return fun_call_result

    def load_result(self, result: FunctionCallResult) -> Any:
        if result.result_type == 'pickle':
            return pickle.loads(result.result)
        else:
            return result

    def lock_call(self, function_call: FunctionCall) -> FunctionCallResult:
        result_type = 'intent'
        b = b""
        call_hash = function_call.call_hash
        dt = 0
        fun_call_intent = FunctionCallResult(call_hash, result_type, b, dt, datetime.now(), True)
        self._store.insert(fun_call_intent)
        return fun_call_intent


def proxy_function_calls(value):
    if isinstance(value, list):
        return [proxy_function_calls(inner) for inner in value]
    elif isinstance(value, dict):
        return {key: proxy_function_calls(inner) for key, inner in value.items()}
    elif isinstance(value, FunctionCall):
        return FunctionCallProxy(value.hash)
    else:
        return value


def _create_kwargs(f, *args: Any, **kwargs: Any):
    params = inspect.signature(f).parameters
    kwargs.update({param: proxy_function_calls(arg) for param, arg in zip(params, args)})
    return kwargs



def calculate_call_hash(function_definition, pickled_kwargs: bytes):
    hasher = sha3_256()
    hasher.update(function_definition.function_hash)
    hasher.update(pickled_kwargs)
    hash = hasher.digest()
    return hash


