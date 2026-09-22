"""Human-readable and machine-readable reports with no implied design pass."""
from __future__ import annotations

import html
import json
from pathlib import Path

import numpy as np


def write_report(directory: str | Path, report: dict, vectors: np.ndarray) -> Path:
    root = Path(directory)
    root.mkdir(parents=True, exist_ok=False)
    (root / "report.json").write_text(json.dumps(report, indent=2, allow_nan=False), encoding="utf-8")
    np.savez_compressed(root / "modes.npz", vectors=vectors,
                        factors=np.array([m["factor"] for m in report["modes"]]))
    rows = []
    for member in report["members"]:
        factors = ", ".join(f'{m["factor"]:.7g}' for m in member["candidate_modes"]) or "--"
        cells = [member["member"], member["status"], factors, "NOT CERTIFIED",
                 member["design_check"], "; ".join(member["notes"])]
        rows.append("<tr>" + "".join("<td>" + html.escape(str(c)) + "</td>" for c in cells) + "</tr>")
    page = '''<!doctype html><html lang="en"><meta charset="utf-8">
<title>BeamLBA numerical assessment</title><style>
body{font:15px system-ui,sans-serif;margin:3rem;line-height:1.5;max-width:1500px}
table{border-collapse:collapse;width:100%}th,td{border:1px solid #bbb;padding:.6rem;text-align:left}
th{background:#eee}aside{padding:1rem;border:2px solid #555;margin:1rem 0}
</style><h1>BeamLBA - numerical candidates</h1>
<aside><strong>RESEARCH ALPHA - NOT A DESIGN APPROVAL.</strong><br>'''
    page += html.escape(" ".join(report["warnings"])) + "</aside><p>"
    page += html.escape("Load case: " + str(report.get("provenance", {}).get("load_case", "unspecified")))
    page += "</p><table><thead><tr>"
    page += "".join(f"<th>{t}</th>" for t in ("Member", "Status", "Candidate factors", "Coverage", "Design check", "Notes"))
    page += "</tr></thead><tbody>" + "".join(rows) + "</tbody></table></html>"
    (root / "report.html").write_text(page, encoding="utf-8")
    return root
