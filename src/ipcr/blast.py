from pathlib import Path
from enum import Enum
from subprocess import run, SubprocessError
from shutil import which
from typing import Iterable
from tempfile import NamedTemporaryFile
from collections import defaultdict
from dataclasses import dataclass

MAKEBLAST_BIN = "makeblastdb"


class BlastDbType(Enum):
    NUCL = 1
    PROT = 2


class BlastProgram(Enum):
    BLASTN = 1
    BLASTP = 2
    BLASTX = 3
    TBLASTN = 4


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


TABBLAST_OUTFMT = "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qframe sframe"


@dataclass
class Seq:
    name: str
    seq: str


@dataclass
class HSP:
    identity: float
    ali_len: int
    mismatch: int
    gap_opens: int
    query_start: int
    query_end: int
    subject_start: int
    subject_end: int
    evalue: float
    score: float
    query_strand: int
    subject_strand: int


def create_fasta_file(seqs, fasta_fhand):
    for seq in seqs:
        fasta_fhand.write(f">{seq.name}\n{seq.seq}\n")

    fasta_fhand.flush()


def blast_seqs(
    seqs: Iterable[Seq],
    db_path: Path,
    blast_program: BlastProgram,
    tmp_dir: Path | None = None,
    evalue_threshold=1e-5,
    short=False,
):
    if blast_program == BlastProgram.BLASTN:
        blast_cmd = "blastn"
    elif blast_program == BlastProgram.BLASTP:
        blast_cmd = "blastp"
    elif blast_program == BlastProgram.BLASTX:
        blast_cmd = "blastx"
    elif blast_program == BlastProgram.TBLASTN:
        blast_cmd = "tblastn"

    with NamedTemporaryFile(suffix=".fasta", dir=tmp_dir, mode="wt") as tmp_fasta:
        create_fasta_file(seqs, tmp_fasta)

        cmd = [
            blast_cmd,
            "-query",
            tmp_fasta.name,
            "-db",
            str(db_path),
            "-evalue",
            str(evalue_threshold),
            "-outfmt",
            TABBLAST_OUTFMT,
        ]

        if short:
            if blast_program == BlastProgram.BLASTN:
                cmd.extend(["-task", "blastn-short"])
            else:
                raise ValueError(
                    f"short option can only be used with blastn, but blast program is: {blast_cmd}"
                )

        try:
            process = run(cmd, check=False, capture_output=True)
        except SubprocessError:
            if which(blast_cmd) is None:
                msg = f"{blast_cmd} is not installed in your system and it is required to create the blast database"
                raise RuntimeError(msg)
            else:
                raise
    lines = process.stdout.decode().splitlines()
    result = defaultdict(dict)
    for line in lines:
        items = line.strip().split()
        if len(items) == 14:
            (
                query,
                subject,
                identity,
                ali_len,
                mis,
                gap_opens,
                query_start,
                query_end,
                subject_start,
                subject_end,
                expect,
                score,
                qstrand,
                sstrand,
            ) = items
        else:
            raise RuntimeError("Wrong blast output")

        hsp = HSP(
            identity=float(identity),
            ali_len=int(ali_len),
            mismatch=int(mis),
            gap_opens=int(gap_opens),
            query_start=int(query_start),
            query_end=int(query_end),
            subject_start=int(subject_start),
            subject_end=int(subject_end),
            evalue=float(expect),
            score=float(score),
            query_strand=int(qstrand),
            subject_strand=int(sstrand),
        )

        try:
            hsps = result[query][subject]
        except KeyError:
            hsps = []
            result[query][subject] = hsps
        hsps.append(hsp)
    return result


def filter_hsps_by_align_len(hsps, min_len=None):
    if min_len is None:
        return hsps

    filtered_hsps = []
    for hsp in hsps:
        if hsp.ali_len >= min_len:
            filtered_hsps.append(hsp)
    return filtered_hsps


def filter_hsps_by_identity(hsps, min_identity=None):
    if min_identity is None:
        return hsps

    filtered_hsps = []
    for hsp in hsps:
        if hsp.identity >= min_identity:
            filtered_hsps.append(hsp)
    return filtered_hsps


def filter_hsps_by_mismatch(hsps, max_mismatch=None):
    if max_mismatch is None:
        return hsps

    filtered_hsps = []
    for hsp in hsps:
        if hsp.mismatch <= max_mismatch:
            filtered_hsps.append(hsp)
    return filtered_hsps
