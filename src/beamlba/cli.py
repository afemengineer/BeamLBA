"""Command line for verification fixtures and reviewed matrix bundles."""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

from . import __version__
from .assessment import assess
from .benchmarks import crowded_frame, fixture_metadata, vlasov_beam
from .bundle import load_bundle, pack_matrix_market, write_bundle
from .numerics import Config, InvalidModel
from .report import write_report


def parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="BeamLBA research alpha - no design approval")
    p.add_argument("--version", action="version", version=__version__)
    sub = p.add_subparsers(dest="command", required=True)
    demo = sub.add_parser("demo", help="Generate and assess a verification fixture")
    demo.add_argument("--kind", choices=("beam", "crowded"), default="beam")
    demo.add_argument("--out", required=True)
    demo.add_argument("--target", type=int, default=281)
    run = sub.add_parser("run", help="Assess a prepared hash-checked matrix bundle")
    run.add_argument("bundle")
    run.add_argument("--out", required=True)
    run.add_argument("--config", help="JSON Config fields plus optional member-name shifts")
    run.add_argument("--diagnostic", action="store_true", help="Allow unreviewed ANSYS inputs, explicitly flagged")
    verify = sub.add_parser("verify", help="Check bundle integrity and native algebraic parity")
    verify.add_argument("bundle")
    pack = sub.add_parser("pack", help="Package Matrix Market files and an equation-space specification")
    for key in ("k", "g", "spec", "out"):
        pack.add_argument("--" + key, required=True)
    pack.add_argument("--g-sign", type=int, choices=(-1, 1), default=1)
    pack.add_argument("--triangle", choices=("none", "lower", "upper"), default="none")
    return p


def main(argv: list[str] | None = None) -> int:
    args = parser().parse_args(argv)
    try:
        if args.command == "pack":
            print(pack_matrix_market(args.k, args.g, args.spec, args.out, args.g_sign, args.triangle))
            return 0
        if args.command == "demo":
            if not 1 <= args.target <= 300:
                raise InvalidModel("Target must be between 1 and 300")
            root = Path(args.out)
            root.mkdir(parents=True, exist_ok=False)
            if args.kind == "beam":
                pencil, analytical = vlasov_beam()
                values, xs = pencil.dense()
                check = {"analytical_mcr_Nm": analytical, "fe_mcr_Nm": float(values[0]*20000),
                         "relative_error": abs(float(values[0]*20000/analytical)-1)}
            else:
                pencil = crowded_frame()
                values, xs = pencil.dense()
                target = pencil.members[args.target-1]
                participation = np.sum(xs[target.dofs, :]**2, axis=0)
                rank = int(np.argmax(participation))
                check = {"dense_reference_rank": rank+1, "dense_reference_factor": float(values[rank]),
                         "first_20_target_participation_max": float(np.max(participation[:20]))}
                pencil.members = [target]
            metadata = fixture_metadata(args.kind)
            write_bundle(root / "bundle", pencil.k, pencil.g, pencil.members, metadata)
            report, vectors = assess(pencil)
            report["provenance"], report["verification_fixture"] = metadata, check
            write_report(root / "assessment", report, vectors)
            print(json.dumps(check, indent=2))
            print(root / "assessment" / "report.html")
            return 0
        settings = {}
        if args.command == "run" and args.config:
            settings = json.loads(Path(args.config).read_text(encoding="utf-8"))
        overrides = settings.pop("shifts", {})
        pencil, manifest, native = load_bundle(args.bundle, Config(**settings))
        metadata = manifest["metadata"]
        parity = pencil.native_parity(*native) if native is not None else {"algebraic_parity": False}
        reviewed = (native is not None and parity["algebraic_parity"]
                    and len(native[0]) >= min(3, pencil.n)
                    and metadata.get("mapping_reviewed") is True
                    and metadata.get("element_formulation_reviewed") is True)
        if args.command == "verify":
            print(json.dumps({"integrity_valid": True, "metadata": metadata, "parity": parity,
                              "input_review_gate_satisfied": reviewed}, indent=2))
            return 0 if metadata["source"] == "synthetic" or reviewed else 2
        if metadata["source"] == "ansys" and not reviewed and not args.diagnostic:
            raise InvalidModel("ANSYS gate: require >=3 native eigenpairs where available, reviewed equation mapping and element formulation. --diagnostic is for integration debugging only")
        report, vectors = assess(pencil, overrides, initial_modes=native)
        report["provenance"], report["native_parity"] = metadata, parity
        report["input_review_gate_satisfied"] = bool(reviewed)
        if metadata["source"] == "ansys" and not reviewed:
            report["warnings"].insert(0, "UNVERIFIED ANSYS INPUT: diagnostic use only.")
            for row in report["members"]:
                row["status"] = "DIAGNOSTIC_ONLY"
        write_report(args.out, report, vectors)
        print(Path(args.out) / "report.html")
        return 0
    except (InvalidModel, ValueError, TypeError, OSError, KeyError) as exc:
        print(f"BeamLBA error: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
