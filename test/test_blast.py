from pathlib import Path
from tempfile import TemporaryDirectory
from collections import namedtuple

from ipcr.blast import (
    prepare_blast,
    BlastDbType,
    blast_seqs,
    BlastProgram,
    HSP,
    filter_hsps_by_align_len,
    filter_hsps_by_identity,
    filter_hsps_by_mismatch,
)

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


def test_hsp_filtering():
    hsp1 = HSP(
        identity=100.0,
        ali_len=24,
        mismatch=1,
        gap_opens=0,
        query_start=1,
        query_end=24,
        subject_start=1268,
        subject_end=1291,
        evalue=3.78e-10,
        score=48.1,
        query_strand=1,
        subject_strand=1,
    )
    hsps = filter_hsps_by_align_len([hsp1], min_len=24)
    assert len(hsps) == 1
    hsps = filter_hsps_by_align_len([hsp1], min_len=25)
    assert not hsps

    hsps = filter_hsps_by_identity([hsp1], min_identity=100.0)
    assert len(hsps) == 1
    hsps = filter_hsps_by_align_len([hsp1], min_len=110.0)
    assert not hsps

    hsps = filter_hsps_by_mismatch([hsp1], max_mismatch=1)
    assert len(hsps) == 1
    hsps = filter_hsps_by_mismatch([hsp1], max_mismatch=0)
    assert not hsps
