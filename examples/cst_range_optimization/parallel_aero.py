"""Independent Scirocco processes evaluated concurrently, with analytic Jacobians."""
from concurrent.futures import ThreadPoolExecutor
import json
import os
from pathlib import Path
import socket
import subprocess
import time
import grpc
import numpy as np
import openmdao.api as om
from philote_mdo.general import ExplicitClient


class AeroPool:
    def __init__(self, binary, panels, targets, directory, workers=8):
        self.directory = Path(directory)
        self.targets = np.asarray(targets)
        self.panels = panels
        self.processes = []
        self.channels = []
        self.clients = []
        self.logs = []
        self.executor = ThreadPoolExecutor(max_workers=workers)
        self.last = None
        self.outputs = None
        self.jac = None
        self.calls = 0
        self.derivative_calls = 0
        self.history = (self.directory / "evaluations.jsonl").open("a", buffering=1)
        try:
            for i, target in enumerate(targets):
                with socket.socket() as sock:
                    sock.bind(("127.0.0.1", 0))
                    address = f"127.0.0.1:{sock.getsockname()[1]}"
                log = (self.directory / f"server_{i}.log").open("w")
                self.logs.append(log)
                process = subprocess.Popen(
                    [str(binary), address], stdout=log, stderr=log
                )
                self.processes.append(process)
                self.write_manifest()
                channel = grpc.insecure_channel(address)
                self.channels.append(channel)
                grpc.channel_ready_future(channel).result(timeout=20)
                client = ExplicitClient(channel)
                client.get_available_options()
                client.send_options({"naca": "0012", "npanel": panels})
                client.run_setup()
                client.get_variable_definitions()
                client.get_partials_definitions()
                self.clients.append(client)
        except BaseException:
            self.close()
            raise

    def write_manifest(self):
        (self.directory / "processes.json").write_text(
            json.dumps(
                {
                    "parent_pid": os.getpid(),
                    "server_pids": [p.pid for p in self.processes],
                },
                indent=2,
            )
            + "\n"
        )

    def matches(self, x, alpha):
        return (
            self.last is not None
            and np.array_equal(x, self.last[0])
            and np.array_equal(alpha, self.last[1])
        )

    def requests(self, x, alpha):
        return [
            {
                "x": np.asarray(x).copy(),
                "alpha": np.array([a]),
                "Re": np.array([1e6]),
                "Ma": np.array([0.0]),
            }
            for a in alpha
        ]

    def settled(self, method, requests):
        futures = [
            self.executor.submit(getattr(c, method), v)
            for c, v in zip(self.clients, requests)
        ]
        results = []
        errors = []
        # Drain every request before another batch may touch the server states.
        for i, future in enumerate(futures):
            try:
                results.append(future.result())
            except Exception as exc:
                results.append(None)
                errors.append(f"point {i}: {exc}")
        if errors:
            raise om.AnalysisError("; ".join(errors))
        return results

    def evaluate(self, x, alpha):
        if self.matches(x, alpha):
            return self.outputs
        self.last = None
        self.jac = None
        start = time.monotonic()
        self.calls += 1
        record = {
            "evaluation": self.calls,
            "success": False,
            "error": "Evaluation interrupted",
        }
        try:
            results = self.settled("run_compute", self.requests(x, alpha))
            outputs = {
                k: np.array([float(np.asarray(r[k]).reshape(-1)[0]) for r in results])
                for k in ["cl", "cd", "cm"]
            }
            if not all(np.all(np.isfinite(v)) for v in outputs.values()):
                raise om.AnalysisError("Nonfinite force coefficient")
            self.last = (np.array(x, copy=True), np.array(alpha, copy=True))
            self.outputs = outputs
            record = {
                "evaluation": self.calls,
                "success": True,
                "seconds": time.monotonic() - start,
                "x": np.asarray(x).tolist(),
                "alpha": np.asarray(alpha).tolist(),
                **{k: v.tolist() for k, v in outputs.items()},
            }
            (self.directory / "last_successful_evaluation.json").write_text(
                json.dumps(record) + "\n"
            )
            return outputs
        except Exception as exc:
            record = {
                "evaluation": self.calls,
                "success": False,
                "seconds": time.monotonic() - start,
                "error": str(exc),
            }
            raise
        finally:
            self.history.write(json.dumps(record) + "\n")

    def derivatives(self, x, alpha):
        self.evaluate(x, alpha)
        if self.jac is None:
            start = time.monotonic()
            results = self.settled("run_compute_partials", self.requests(x, alpha))
            self.jac = {}
            for key in ["cl", "cd", "cm"]:
                self.jac[key, "x"] = np.vstack(
                    [np.asarray(r[key, "x"]).reshape(1, -1) for r in results]
                )
                self.jac[key, "alpha"] = np.array(
                    [float(np.asarray(r[key, "alpha"]).reshape(-1)[0]) for r in results]
                )
            self.derivative_calls += 1
            self.history.write(
                json.dumps(
                    {
                        "derivative_evaluation": self.derivative_calls,
                        "seconds": time.monotonic() - start,
                    }
                )
                + "\n"
            )
        return self.jac

    def close(self):
        for process in self.processes:
            if process.poll() is None:
                process.terminate()
        for process in self.processes:
            try:
                process.wait(timeout=5)
            except subprocess.TimeoutExpired:
                process.kill()
                process.wait()
        for channel in self.channels:
            channel.close()
        self.executor.shutdown(wait=True, cancel_futures=True)
        for log in self.logs:
            log.close()
        self.history.close()
        if self.processes:
            (self.directory / "processes.json").write_text(
                json.dumps(
                    {
                        "parent_pid": os.getpid(),
                        "server_pids": [],
                        "servers_stopped": True,
                    },
                    indent=2,
                )
                + "\n"
            )


class ParallelAero(om.ExplicitComponent):
    def initialize(self):
        self.options.declare("pool", recordable=False)

    def setup(self):
        pool = self.options["pool"]
        n = len(pool.targets)
        self.add_input("x", shape=2 * (pool.panels + 1), units="m")
        self.add_input("alpha", shape=n, units="deg")
        for key in ["cl", "cd", "cm"]:
            self.add_output(key, shape=n)
            self.declare_partials(key, "x")
            self.declare_partials(key, "alpha", rows=np.arange(n), cols=np.arange(n))

    def compute(self, inputs, outputs):
        values = self.options["pool"].evaluate(inputs["x"], inputs["alpha"])
        for key, value in values.items():
            outputs[key] = value

    def compute_partials(self, inputs, partials):
        for key, value in (
            self.options["pool"].derivatives(inputs["x"], inputs["alpha"]).items()
        ):
            partials[key] = value
