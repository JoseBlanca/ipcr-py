from pathlib import Path
from enum import Enum
from subprocess import run, SubprocessError
from shutil import which

MAKEBLAST_BIN = "makeblastdb"


class BlastDbType(Enum):
    NUCL = 1
    PROT = 2


def prepare_blast(
    seqs_fasta: Path,
    db_type: BlastDbType,
    out_db_dir: Path,
    skip_if_exists=True,
):
    seqs_fasta = Path(seqs_fasta)
    out_db_dir = Path(out_db_dir)

    out_base_path = out_db_dir / seqs_fasta.stem

    if db_type == BlastDbType.NUCL:
        out_path = Path(str(out_base_path) + ".nhr")
    elif db_type == BlastDbType.PROT:
        out_path = Path(str(out_base_path) + ".phr")

    if skip_if_exists and out_path.exists():
        return {"db_path": out_base_path}

    cmd = [
        MAKEBLAST_BIN,
        "-in",
        str(seqs_fasta),
        "-out",
        str(out_base_path),
        "-dbtype",
    ]

    if db_type == BlastDbType.NUCL:
        cmd.append("nucl")
    elif db_type == BlastDbType.PROT:
        cmd.append("prot")

    try:
        run(cmd, check=True, capture_output=True)
    except SubprocessError:
        if which(MAKEBLAST_BIN) is None:
            msg = f"{MAKEBLAST_BIN} is not installed in your system and it is required to create the blast database"
            raise RuntimeError(msg)
        else:
            raise

    return {"db_path": out_base_path}
