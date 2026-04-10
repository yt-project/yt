import contextlib
import io
import sys
import types
from unittest import mock

import pytest

import yt.utilities.rpdb as rpdb


def _sample_traceback():
    try:
        raise RuntimeError("boom")
    except RuntimeError as exc:
        return exc.__traceback__


def test_pdb_xmlrpc_server_signal_shutdown_and_serve_forever():
    server = rpdb.PdbXMLRPCServer.__new__(rpdb.PdbXMLRPCServer)
    server.finished = False

    with mock.patch.object(rpdb.signal, "signal") as signal_mock:
        server.register_signal(15)

    signum, handler = signal_mock.call_args[0]
    assert signum == 15
    assert handler.__self__ is server
    assert handler.__func__ is rpdb.PdbXMLRPCServer.signal_handler

    assert server.shutdown() == 1
    assert server.finished is True

    server.finished = False
    stdout = io.StringIO()
    with contextlib.redirect_stdout(stdout):
        server.signal_handler(9, None)
    assert server.finished is True
    assert stdout.getvalue().strip() == "Caught signal 9"

    served = rpdb.PdbXMLRPCServer.__new__(rpdb.PdbXMLRPCServer)
    served.finished = False
    calls = []

    def _handle_request():
        calls.append("handled")
        served.finished = True

    served.handle_request = _handle_request
    stdout = io.StringIO()
    with contextlib.redirect_stdout(stdout):
        served.serve_forever()

    assert calls == ["handled"]
    assert stdout.getvalue().strip() == "DONE SERVING"


def test_pdb_handler_initializes_debugger_and_executes_commands():
    tb = _sample_traceback()

    class FakePdb:
        def __init__(self, stdin, stdout):
            self.stdin = stdin
            self.stdout = stdout
            self.reset_calls = 0
            self.setup_args = None

        def reset(self):
            self.reset_calls += 1

        def setup(self, frame, traceback_obj):
            self.setup_args = (frame, traceback_obj)

        def onecmd(self, line):
            self.stdout.write(f"executed {line}\n")

    original_stdin, original_stdout, original_stderr = sys.stdin, sys.stdout, sys.stderr
    try:
        with mock.patch.dict(sys.modules, {"pdb": types.SimpleNamespace(Pdb=FakePdb)}):
            handler = rpdb.pdb_handler(tb)

        assert handler.debugger.reset_calls == 1
        assert handler.debugger.setup_args == (tb.tb_frame, tb)
        assert handler.execute("where") == "executed where\n"
    finally:
        sys.stdin = original_stdin
        sys.stdout = original_stdout
        sys.stderr = original_stderr


def test_rpdb_excepthook_starts_server_and_optionally_barriers():
    tb = _sample_traceback()
    fake_handler = object()

    class FakeServer:
        created = []

        def __init__(self, address):
            self.address = address
            self.introspection_registered = False
            self.instance = None
            self.registered_function = None
            self.served = False
            self.closed = False
            self.finished = False
            type(self).created.append(self)

        def register_introspection_functions(self):
            self.introspection_registered = True

        def register_instance(self, instance):
            self.instance = instance

        def register_function(self, func):
            self.registered_function = func

        def shutdown(self):
            self.finished = True
            return 1

        def serve_forever(self):
            self.served = True

        def server_close(self):
            self.closed = True

    barrier = mock.Mock()
    fake_mpi = types.SimpleNamespace(
        MPI=types.SimpleNamespace(COMM_WORLD=types.SimpleNamespace(Barrier=barrier))
    )

    stdout = io.StringIO()
    with (
        contextlib.redirect_stdout(stdout),
        mock.patch.object(rpdb.traceback, "print_exception") as print_exception,
        mock.patch.object(rpdb, "pdb_handler", return_value=fake_handler) as handler_ctor,
        mock.patch.object(rpdb, "PdbXMLRPCServer", FakeServer),
        mock.patch.object(rpdb.ytcfg, "get", side_effect=[2, 3]),
        mock.patch.dict(sys.modules, {"mpi4py": fake_mpi}),
    ):
        rpdb.rpdb_excepthook(RuntimeError, RuntimeError("boom"), tb)

    print_exception.assert_called_once()
    handler_ctor.assert_called_once_with(tb)
    assert "Starting RPDB server on task 2" in stdout.getvalue()
    server = FakeServer.created[0]
    assert server.address == ("localhost", 8012)
    assert server.introspection_registered is True
    assert server.instance is fake_handler
    assert server.registered_function() == 1
    assert server.served is True
    assert server.closed is True
    barrier.assert_called_once_with()

    FakeServer.created.clear()
    stdout = io.StringIO()
    with (
        contextlib.redirect_stdout(stdout),
        mock.patch.object(rpdb.traceback, "print_exception"),
        mock.patch.object(rpdb, "pdb_handler", return_value=fake_handler),
        mock.patch.object(rpdb, "PdbXMLRPCServer", FakeServer),
        mock.patch.object(rpdb.ytcfg, "get", side_effect=[0, 1]),
    ):
        rpdb.rpdb_excepthook(RuntimeError, RuntimeError("boom"), tb)

    assert "Starting RPDB server on task 0" in stdout.getvalue()


def test_rpdb_cmd_prints_proxy_results_and_swallows_shutdown_errors():
    class FakeProxy:
        def __init__(self):
            self.calls = []

        def execute(self, command):
            self.calls.append(("execute", command))
            return f"result:{command}"

        def shutdown(self):
            self.calls.append(("shutdown", None))
            return 1

    proxy = FakeProxy()
    stdout = io.StringIO()
    with contextlib.redirect_stdout(stdout):
        shell = rpdb.rpdb_cmd(proxy)
        shell.default("next")
        shell.do_help("where")
        assert shell.do_shutdown("") is True

    assert proxy.calls == [
        ("execute", "bt"),
        ("execute", "next"),
        ("execute", "help where"),
        ("shutdown", None),
    ]
    assert stdout.getvalue().splitlines() == [
        "result:bt",
        "result:next",
        "result:help where",
        "1",
    ]
    assert shell.postcmd("stop", "line") == "stop"

    failing_shell = rpdb.rpdb_cmd.__new__(rpdb.rpdb_cmd)
    failing_shell.proxy = mock.Mock()
    failing_shell.proxy.shutdown.side_effect = RuntimeError("ignore")
    failing_shell.postloop()


def test_run_rpdb_connects_to_requested_task_and_handles_errors():
    class FakeShell:
        created = []

        def __init__(self, proxy):
            self.proxy = proxy
            self.headers = []
            type(self).created.append(self)

        def cmdloop(self, header):
            self.headers.append(header)

    with (
        mock.patch.object(rpdb, "ServerProxy", return_value="proxy") as server_proxy,
        mock.patch.object(rpdb, "rpdb_cmd", FakeShell),
    ):
        rpdb.run_rpdb(task=3)

    server_proxy.assert_called_once_with("http://localhost:8013/")
    assert "task 3" in FakeShell.created[0].headers[0]

    stdout = io.StringIO()
    with (
        contextlib.redirect_stdout(stdout),
        mock.patch.object(rpdb, "ServerProxy", return_value="proxy"),
        mock.patch.object(rpdb, "rpdb_cmd", side_effect=OSError),
        mock.patch.object(rpdb.sys, "exit", side_effect=SystemExit(1)),
        pytest.raises(SystemExit, match="1"),
    ):
        rpdb.run_rpdb(task=0)

    assert stdout.getvalue().strip() == "Connection refused.  Is the server running?"

    with mock.patch.object(rpdb.sys, "argv", ["yt", "rpdb", "4"]):
        with pytest.raises(TypeError):
            rpdb.run_rpdb(task=None)
