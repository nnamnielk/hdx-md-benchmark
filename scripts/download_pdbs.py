#!/usr/bin/env python3
"""Download and clean bundled PDB references for the active benchmark."""

from pathlib import Path
from urllib.request import urlopen


ROOT = Path(__file__).resolve().parents[1]
PDB_DIR = ROOT / "pdb"

PROTEINS = {
    # Trim the HHHAM tag, then convert NuG2 to the NuG2b sequence reported
    # by Skinner et al. Mutated sites retain backbone atoms only.
    "nug2b": ("1MI0", "A", 6, {37: "ALA", 46: "ASP", 47: "ALA"}),
    "ubiquitin": ("1UBQ", "A", 1, {}),
}


def select_chain(
    pdb_text: str,
    chain_id: str,
    first_residue: int,
    mutations: dict[int, str],
) -> str:
    """Keep first-model protein ATOM records for one chain and one altloc."""
    selected = []
    in_first_model = True
    residue_numbers = {}

    for line in pdb_text.splitlines():
        record = line[:6].strip()
        if record == "MODEL":
            in_first_model = line[10:14].strip() in {"", "1"}
            continue
        if record == "ENDMDL" and in_first_model:
            break
        if record != "ATOM" or not in_first_model:
            continue
        if len(line) < 22 or line[21].strip() != chain_id:
            continue
        if len(line) > 16 and line[16] not in {" ", "A"}:
            continue
        residue_number = int(line[22:26])
        if residue_number < first_residue:
            continue
        residue_key = (residue_number, line[26])
        if residue_key not in residue_numbers:
            residue_numbers[residue_key] = len(residue_numbers) + 1
        renumbered = residue_numbers[residue_key]
        if renumbered in mutations:
            if line[12:16].strip() not in {"N", "CA", "C", "O"}:
                continue
            line = f"{line[:17]}{mutations[renumbered]:>3}{line[20:]}"
        selected.append(f"{line[:22]}{renumbered:4d} {line[27:]}")

    if not selected:
        raise ValueError(f"No ATOM records found for chain {chain_id}")

    return "\n".join(selected + ["TER", "END", ""])


def main() -> None:
    PDB_DIR.mkdir(parents=True, exist_ok=True)

    for protein, (pdb_id, chain_id, first_residue, mutations) in PROTEINS.items():
        url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
        print(f"Downloading {protein}: {pdb_id} chain {chain_id}")
        with urlopen(url, timeout=30) as response:
            pdb_text = response.read().decode("ascii")

        output = PDB_DIR / f"{protein}.pdb"
        output.write_text(
            select_chain(pdb_text, chain_id, first_residue, mutations),
            encoding="ascii",
        )
        print(f"  wrote {output}")


if __name__ == "__main__":
    main()
