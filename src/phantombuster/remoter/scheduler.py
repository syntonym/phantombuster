from typing import Any, Callable, Optional
import os.path
import sys
import inspect
from datetime import datetime
import logging
import uuid
import time
import threading

from . import globaladdress
from .persisters import PersisterMapper
from .task import Task, hash_function
from .store import Store
from .messages import Message, StoreMsg, JoinMsg, WaitMsg, QuitMsg, parse_json_to_message, message_to_json
from .async_ import async_gather, async_sleep
from .socket_ import Socket, ReqSocket
from .worker import Worker
from .server import Server

class ServerManager:

    def __init__(self, address, storepath=None):
        self.address = address
        self.server = None
        self.storepath = storepath

    def start(self):
        self.thread = threading.Thread(target=self._start)
        self.thread.start()

    def _start(self):
        self.server = Server(self.address, storepath=self.storepath)
        self.server.start_async()

    def wait(self):
        self.thread.join()

    def stop(self):
        self.thread.stop()


class SingleWorker:

    def __init__(self, address, name):
        self.address = address
        self.name = name

    def start(self):
        self.thread = threading.Thread(target=self._start)
        self.thread.start()

    def _start(self):
        self.worker = Worker(self.address, self.name)
        self.worker.start_async()

    def stop(self):
        self.thread.stop()

    def wait(self):
        self.thread.join()


class WorkerManager:

    def __init__(self, address):
        self.address = address
        self.workers = []

    def start(self):
        worker = SingleWorker(self.address, name=f'worker{len(self.workers)+1}')
        self.workers.append(worker)
        worker.start()

    def wait(self):
        for worker in self.workers:
            worker.wait()

    def stop(self):
        for worker in self.workers:
            worker.stop()


class Scheduler:

    DEFAULT_WORKERS = 0

    def __init__(self,
                 address: Optional[str] = None,
                 persister_mapper=PersisterMapper,
                 name: Optional[str] = None,
                 start_server: bool = False,
                 start_workers: int = -1,
                 save_functions: bool = False,
                 server_store_path: Optional[str] = None,
                 ):
        # If no address is given and no explicit address was ever given, use default address
        if address is None and globaladdress.global_address is None:
            address = "ipc://.remoter.ipc"
        # if explicit address was given, reuse it for later
        elif address is not None and globaladdress.global_address is None:
            globaladdress.global_address = address
        # if no explicit address is given, reuse last explicit address
        elif address is None and globaladdress.global_address is not None:
            address = globaladdress.global_address

        if start_workers == -1:
            start_workers = Scheduler.DEFAULT_WORKERS

        self.worker_manager = WorkerManager(address)
        self.server_manager = ServerManager(address, server_store_path)

        self._address = address
        if name:
            self._name = name
        else:
            self._name = uuid.uuid4().hex
        self._running = False
        self._persister_mapper = persister_mapper
        self._current_file_module = os.path.split(sys.argv[0])[-1].replace(".py", "")
        self._store = Store(":memory:", name=self._name)
        self._logger = logging.getLogger(f"remoter.scheduler.{self._name}")
        self._logger.debug(f"Connecting to {address}")
        self._socket: Socket = self._create_socket()
        self._socket.connect(address)
        self._joined = False
        self._entered = 0
        self._failed = False
        self._tasks: list["Task"] = []
        self._start_server = start_server
        self._start_workers = start_workers
        self._save_functions = save_functions




    def _create_task(self, f: Callable, *args: Any, **kwargs: Any) -> Task:

        params = inspect.signature(f).parameters
        kwargs.update({param: arg for param, arg in zip(params, args)})
        module = f.__module__ if f.__module__ != "__main__" else self._current_file_module
        code = inspect.getsource(f)
        task = Task(None, self._store, datetime.now(), module, f.__name__, hash_function(f), kwargs, code=code)
        task.persister = self._persister_mapper(f)
        task.task_hash = task.calc_hash()
        return task

    def _fullfill(self, task, result):
        task.result = result
        task.done = True
        task._save_force()

    def _create_socket(self):
        return ReqSocket()

    def _join(self):
        msg = JoinMsg(self._name, "scheduler")
        self._send(msg)
        msg = self._recv()
        self._joined = True

    def run(self, f: Callable, *args: Any, **kwargs: Any) -> Task:
        if not self._joined:
            self._join()
        task = self.schedule(f, *args, **kwargs)
        self.gather(task)
        return task.result

    def request(self, f: Callable, *args: Any, **kwargs: Any) -> Any:
        if not self._joined:
            self._join()
        task = self._create_task(f, *args, **kwargs)
        self._logger.info("scheduling task: %s", task.function_name)
        msg = task.to_task_message(name=self._name, request=False)
        self._send(msg)
        return_msg = self._recv()
        if isinstance(return_msg, StoreMsg):
            self._store._income(return_msg)
            task.calc_hash()
            task.load()
        return task

    def schedule(self, f: Callable, *args: Any, **kwargs: Any) -> Any:
        if not self._joined:
            self._join()
        task = self._create_task(f, *args, **kwargs)
        if self._save_functions:
            msg = task.to_fsave_message()
            self._send(msg)
            return_msg = self._recv()
        self._logger.info("scheduling task: %s", task.function_name)
        msg = task.to_task_message(name=self._name)
        self._send(msg)
        return_msg = self._recv()
        if isinstance(return_msg, StoreMsg):
            self._store._income(return_msg)
            task.calc_hash()
        self._tasks.append(task)
        return task

    def wait_all(self, tasks, unpack=False):
        return [self.wait(task, unpack=unpack) for task in tasks]

    def wait(self, task, unpack=False):
        return self._wait_for(task, unpack=unpack)

    def _wait_for(self, task, unpack=False):
        self._logger.debug("Waiting for completion of task: %s", task)
        task.calc_hash()
        while not task.is_loadable():
            msg = task.to_task_message(self._name)
            self._send(msg)
            msg = self._recv()

            if isinstance(msg, StoreMsg):
                self._store._income(msg)
            elif isinstance(msg, WaitMsg):
                time.sleep(msg.time)
            elif isinstance(msg, QuitMsg):
                self._failed = True
                self._joined = False
                return None
            else:
                raise Exception(f"Scheduler received unexpected message: {msg}")
        self._logger.debug("Successfully waited for completion of task: %s", task)
        if unpack:
            task.calc_hash()
            task.load()
            return task.result
        else:
            return task

    def gather(self, *tasks):
        if self._joined:
            for task in tasks:
                self._wait_for(task)
            rs = [task.load() for task in tasks]
            for task in tasks:
                if task.failed:
                    err, tb = task.exception
                    raise err.with_traceback(tb)
            return rs
        else:
            if all(task.is_loadable() for task in tasks):
                rs = [task.load() for task in tasks]
                for task in tasks:
                    if task.failed:
                        err, tb = task.exception
                        raise err.with_traceback(tb)
                return rs
            else:
                raise Exception("Not all tasks are loadable and I'm not connected to a server")

    def _recv(self) -> Message:
        if self._socket is not None:
            self._logger.debug("Waiting for message sync")
            msg_json = self._socket.recv_json()
            msg = parse_json_to_message(msg_json, self._store)
            self._logger.debug(f"Sync Received message: {msg}")
            return msg
        else:
            raise Exception("Scheduler has to connect before being able to receive messages")

    def _send(self, msg: Message) -> None:
        if self._socket:
            j = message_to_json(msg)
            self._logger.debug(f"Sync Sending message: {msg}")
            self._socket.send_json(j)
        else:
            raise Exception("Scheduler has to connect before being able to send messages")

    def __enter__(self):
        if self._entered == 0:
            if self._start_server:
                self.server_manager.start()
            for i in range(self._start_workers):
                self.worker_manager.start()
            self._join()
        self._entered += 1
        return self

    def __exit__(self, exc_type, exc, tb):
        if self._entered == 1:
            self._logger.debug("%s, %s, %s", exc_type, exc, tb)
            self._logger.debug("Scheduler exited, sending exit msg to server sync")
            self.exit(force=True)
            self._logger.debug("Notified server of exit, exiting now")
            self.server_manager.wait()
            self.worker_manager.wait()
        self._entered -= 1
        if exc:
            raise exc
        return self

    def exit(self, force=False):
        if not force:
            for task in self._tasks:
                self._wait_for(task)
        self._send(QuitMsg(self._name, "scheduler"))
        self._recv()
        self._joined = False
        self._running = False

    async def _async_join(self):
        msg = JoinMsg(self._name, "scheduler")
        await self._async_send(msg)
        msg = await self._async_recv()
        self._joined = True

    async def async_run(self, f: Callable, *args: Any, **kwargs: Any) -> Any:
        if not self._joined:
            await self._async_join()
        task = await self.async_schedule(f, *args, **kwargs)
        await task.wait()
        return task.result

    async def async_schedule(self, f: Callable, *args: Any, **kwargs: Any) -> Any:
        if not self._joined:
            await self._async_join()
        task = self._create_task(f, *args, **kwargs)
        if self._save_functions:
            msg = task.to_fsave_message()
            await self._async_send(msg)
            return_msg = await self._async_recv()
        self._logger.info("scheduling task: %s", task.function_name)
        msg = task.to_task_message(self._name)
        await self._async_send(msg)
        return_msg = await self._async_recv()
        if isinstance(return_msg, StoreMsg):
            self._store._income(return_msg)
            task.calc_hash()
            task.load()
        self._tasks.append(task)
        return task

    async def _async_wait_for(self, task):
        self._logger.debug("Waiting for completion of task: %s", task)
        task.calc_hash()

        while not task.is_loadable():
            msg = task.to_task_message(self._name)
            j = message_to_json(msg)
            self._logger.debug(f"Sending message: {msg}")
            await self._socket.async_send_json(j)

            self._logger.debug("Waiting for message async")
            msg_json = await self._socket.async_recv_json()
            msg = parse_json_to_message(msg_json, self._store)
            self._logger.debug(f"Received message: {msg}")

            if isinstance(msg, StoreMsg):
                self._store._income(msg)
                task.calc_hash()
                task.load()
            elif isinstance(msg, WaitMsg):
                await async_sleep(msg.time)
            elif isinstance(msg, QuitMsg):
                self._failed = True
                self._joined = False
                return None
            else:
                raise Exception(f"Scheduler received unexpected message: {msg}")
        self._logger.debug("Successfully waited for completion of task: %s", task)
        return task.result

    async def async_gather(self, *tasks):
        rs = []
        for task in tasks:
            r = await self._async_wait_for(task)
            rs.append(r)
        for task in tasks:
            if task.failed:
                raise Exception(task.result)
        return rs

    async def _async_recv(self) -> Message:
        if self._socket is not None:
            self._logger.debug("Waiting for message async")
            msg_json = await self._socket.async_recv_json()
            msg = parse_json_to_message(msg_json, self._store)
            self._logger.debug(f"Received message: {msg}")
            return msg
        else:
            raise Exception("Scheduler has to connect before being able to receive messages")

    async def _async_send(self, msg: Message) -> None:
        if self._socket:
            j = message_to_json(msg)
            self._logger.debug(f"Sending message: {msg}")
            await self._socket.async_send_json(j)
        else:
            raise Exception("Scheduler has to connect before being able to send messages")

    async def __aenter__(self):
        if self._entered == 0:
            if self._start_server:
                start_server(self._address)
            for i in range(self._start_workers):
                start_worker(self._address, i)
            self._logger.debug("Scheduler entered async, sending join msg to server async")
            await self._async_join()
            self._logger.debug("Joined async")
        self._entered += 1
        return self

    async def __aexit__(self, exc_type, exc, tb):
        if self._entered == 1:
            self._logger.error("%s, %s, %s", exc_type, exc, tb)
            self._logger.debug("Scheduler exited async, sending exit msg to server async")
            await self.async_exit()
            self._logger.debug("Notified server of exit, exiting now")
        self._entered -= 1
        if exc:
            raise exc
        return self

    async def async_exit(self):
        rs = []
        for task in self._tasks:
            r = await self._async_wait_for(task)
            rs.append(r)
        await self._async_send(QuitMsg(self._name, "scheduler"))
        await self._async_recv()
        self._joined = False
        self._running = False
        return rs
