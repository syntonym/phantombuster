import json
import logging
import zmq
from zmq.utils.strtypes import bytes, unicode, cast_bytes

LOG_RECORDS_ATTRIBUTES = ["args", "asctime", "created", "exc_info", "filename", "funcName", "levelname", "levelno", "lineno", "message", "module", "msecs", "msg", "name", "pathname", "process", "processName", "relativeCreated", "stack_info", "thread", "threadName"]

class PushHandler(logging.Handler):
    """A basic logging handler that emits log messages through a Push socket.
    Takes a PUB socket already bound to interfaces or an interface to bind to.
    Example::
        sock = context.socket(zmq.PUSH)
        sock.bind('inproc://log')
        handler = PushHandler(sock)
    Or::
        handler = PushHandler('inproc://loc')
    These are equivalent.
    """

    socket = None

    def __init__(self, interface_or_socket, context=None):
        logging.Handler.__init__(self)
        self.formatters = {
            logging.DEBUG: logging.Formatter(
                "%(levelname)s %(filename)s:%(lineno)d - %(message)s\n"
            ),
            logging.INFO: logging.Formatter("%(message)s\n"),
            logging.WARN: logging.Formatter(
                "%(levelname)s %(filename)s:%(lineno)d - %(message)s\n"
            ),
            logging.ERROR: logging.Formatter(
                "%(levelname)s %(filename)s:%(lineno)d - %(message)s - %(exc_info)s\n"
            ),
            logging.CRITICAL: logging.Formatter(
                "%(levelname)s %(filename)s:%(lineno)d - %(message)s\n"
            ),
        }
        if isinstance(interface_or_socket, zmq.Socket):
            self.socket = interface_or_socket
            self.ctx = self.socket.context
        else:
            self.ctx = context or zmq.Context()
            self.socket = self.ctx.socket(zmq.PUSH)
            self.socket.connect(interface_or_socket)

    def setFormatter(self, fmt, level=logging.NOTSET):
        """Set the Formatter for this handler.
        If no level is provided, the same format is used for all levels. This
        will overwrite all selective formatters set in the object constructor.
        """
        if level == logging.NOTSET:
            for fmt_level in self.formatters.keys():
                self.formatters[fmt_level] = fmt
        else:
            self.formatters[level] = fmt

    def format(self, record):
        """Format a record."""
        return self.formatters[record.levelno].format(record)

    def emit(self, record):
        """Emit a log message on my socket."""
        ei = record.exc_info
        if ei:
            # just to get traceback text into record.exc_text ...
            dummy = self.format(record)
        # See issue #14436: If msg or args are objects, they may not be
        # available on the receiving end. So we convert the msg % args
        # to a string, save it as msg and zap the args.
        d = dict(record.__dict__)
        d['msg'] = record.getMessage()
        d['args'] = None
        d['exc_info'] = None
        # Issue #25685: delete 'message' if present: redundant with 'msg'
        d.pop('message', None)

        self.socket.send_json(d)
        #self.socket.send(bmsg)



