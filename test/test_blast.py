from pathlib import Path
from tempfile import TemporaryDirectory
from collections import namedtuple

from ipcr.blast import prepare_blast, BlastDbType, blast_seqs, BlastProgram

TEST_DIR = Path(__file__).parent
TEST_DATA_DIR = TEST_DIR / "data"

Seq = namedtuple("Seq", ["name", "seq"])


def test_prepare_blast_db():
    seq_fasta = TEST_DATA_DIR / "fxn_human.fasta"
    seq = Seq("primer1", "TAGTGCTGTTTCTCCCACATATTC")
    with TemporaryDirectory() as temp_dir:
        db_path = prepare_blast(
            seq_fasta, db_type=BlastDbType.NUCL, out_db_dir=Path(temp_dir)
        )["db_path"]
        res = blast_seqs(
            [seq], db_path=db_path, blast_program=BlastProgram.BLASTN, short=True
        )
        hsps = res["primer1"]["NM_000144.5"]
        assert len(hsps) == 1
        hsp = hsps[0]
        assert hsp.mismatch == 0
