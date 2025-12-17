from pathlib import Path
from tempfile import TemporaryDirectory

from ipcr.blast import prepare_blast, BlastDbType

TEST_DIR = Path(__file__).parent
TEST_DATA_DIR = TEST_DIR / "data"


def test_prepare_blast_db():
    seq_fasta = TEST_DATA_DIR / "fxn_human.fasta"
    with TemporaryDirectory() as temp_dir:
        prepare_blast(seq_fasta, db_type=BlastDbType.NUCL, out_db_dir=Path(temp_dir))
