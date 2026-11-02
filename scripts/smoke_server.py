"""Exercise a running Scirocco server using the Philote-MDO 0.8 protocol."""

import argparse
import math

import grpc
from google.protobuf.empty_pb2 import Empty
from google.protobuf.struct_pb2 import Struct
from philote_mdo.generated import data_pb2, disciplines_pb2_grpc


def smoke(address, kind):
    with grpc.insecure_channel(address) as channel:
        grpc.channel_ready_future(channel).result(timeout=20)
        service = disciplines_pb2_grpc.DisciplineServiceStub(channel)
        explicit = disciplines_pb2_grpc.ExplicitServiceStub(channel)
        options = Struct()
        options.update(
            {"npanel": 4.5, **({"naca": "0012"} if kind == "analysis" else {})}
        )
        try:
            service.SetOptions(data_pb2.DisciplineOptions(options=options), timeout=10)
        except grpc.RpcError:
            pass
        else:
            raise RuntimeError("Fractional panel count was accepted")
        options["npanel"] = 60
        service.SetOptions(data_pb2.DisciplineOptions(options=options), timeout=10)
        service.Setup(Empty(), timeout=10)
        if kind == "analysis":
            lower_x = [(1 + math.cos(math.pi * i / 30)) / 2 for i in range(31)]
            xs = lower_x + lower_x[-2::-1]
            ys = []
            for i, x in enumerate(xs):
                thickness = 0.6 * (
                    0.2969 * math.sqrt(x)
                    - 0.126 * x
                    - 0.3516 * x**2
                    + 0.2843 * x**3
                    - 0.1015 * x**4
                )
                ys.append(-thickness if i <= 30 else thickness)
            inputs = {"x": xs + ys, "alpha": [4.0], "Re": [0.0], "Ma": [0.0]}
        else:
            inputs = {"cst_upper": [0.1] * 6, "cst_lower": [0.1] * 6}
        messages = [
            data_pb2.VariableMessage(
                continuous=data_pb2.Array(
                    name=name,
                    start=0,
                    end=len(values) - 1,
                    type=data_pb2.kInput,
                    data=values,
                )
            )
            for name, values in inputs.items()
        ]
        outputs = {}
        for message in explicit.ComputeFunction(iter(messages), timeout=20):
            chunk = message.continuous
            values = outputs.setdefault(chunk.name, [])
            if chunk.start != len(values):
                raise RuntimeError("Unexpected output chunk ordering")
            values.extend(chunk.data)
        if not outputs or not all(
            math.isfinite(v) for values in outputs.values() for v in values
        ):
            raise RuntimeError("Missing or nonfinite outputs")
        if kind == "analysis":
            if not (0.3 < outputs["cl"][0] < 0.7 and abs(outputs["cd"][0]) < 1e-12):
                raise RuntimeError(f"Unexpected analysis outputs: {outputs}")
        elif len(outputs["x"]) != 122 or not (
            min(outputs["x"][61:]) < 0 < max(outputs["x"][61:])
        ):
            raise RuntimeError("Unexpected Kulfan geometry")
        print(f"{kind} server passed malformed-request recovery and compute smoke test")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--address", required=True)
    parser.add_argument("--kind", choices=["analysis", "geometry"], required=True)
    args = parser.parse_args()
    smoke(args.address, args.kind)
