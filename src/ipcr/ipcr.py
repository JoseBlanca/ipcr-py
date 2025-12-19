from pathlib import Path
from tempfile import TemporaryDirectory
import itertools
from dataclasses import dataclass

from ipcr.blast import Seq, blast_seqs, prepare_blast, BlastDbType, BlastProgram

MAX_PRODUCT_SIZE = 2000
MAX_PRIMER_MISMATCHES = 0


class _GenomeLocation:
    def __init__(self, chrom, start, end, fwd_strand, primer_name):
        if start > end:
            start, end = end, start
        self.start = start
        self.end = end
        self.chrom = chrom
        self.fwd_strand = fwd_strand
        self.primer_name = primer_name

    @classmethod
    def from_hsp(
        cls,
        chrom,
        subject_start,
        subject_end,
        query_strand,
        subject_strand,
        primer_name,
    ):
        if query_strand == 1 and subject_strand == 1:
            fwd_strand = True
        elif query_strand == -1 and subject_strand == 1:
            fwd_strand = False
        elif query_strand == 1 and subject_strand == -1:
            fwd_strand = False
        elif query_strand == -1 and subject_strand == -1:
            fwd_strand = True
        else:
            raise RuntimeError("We should not be here")

        return cls(chrom, subject_start, subject_end, fwd_strand, primer_name)

    def __str__(self):
        return f"{self.chrom}:{self.start}-{self.end}({self.fwd_strand})"

    def __repr__(self):
        return str(self)


@dataclass
class PCRProduct:
    chrom: str
    start: int
    end: int


class OverlappingPrimerProduct(RuntimeError): ...


class OnePrimerProduct(RuntimeError): ...


def _is_valid_product(
    loc1: _GenomeLocation,
    loc2: _GenomeLocation,
    max_product_size=MAX_PRODUCT_SIZE,
    allow_overlapping_primer_products=False,
    allow_one_primer_products=False,
):
    if loc1.chrom != loc2.chrom:
        return False
    if loc1.fwd_strand == loc2.fwd_strand:
        return False

    fwd_primer, rev_primer = (loc1, loc2) if loc1.fwd_strand else (loc2, loc1)

    product_size = rev_primer.end - fwd_primer.start
    if product_size > max_product_size:
        return False
    if product_size < 0:
        return False

    if fwd_primer.end > rev_primer.end:
        return False

    if (
        not allow_overlapping_primer_products
        and (fwd_primer.end >= rev_primer.start)
        and (fwd_primer.start <= rev_primer.start)
    ):
        raise OverlappingPrimerProduct(
            f"Two primers overlap and would create a product: {fwd_primer.chrom}:{fwd_primer.start}-{rev_primer.end}"
        )

    if not allow_one_primer_products and (
        fwd_primer.primer_name == rev_primer.primer_name
    ):
        raise OnePrimerProduct(
            f"Primer {fwd_primer.primer_name} generates products by itself"
        )

    return True


def filter_hsps_by_mismatch(hsps, max_mismatch, query_len):
    filtered_hsps = []
    for hsp in hsps:
        mismatches = hsp.mismatch + (query_len - hsp.ali_len)
        if mismatches <= max_mismatch:
            filtered_hsps.append(hsp)
    return filtered_hsps


def _get_valid_locations_for_primer(primer, db_path, max_mismatches):
    hsps = blast_seqs(
        [primer], db_path=db_path, blast_program=BlastProgram.BLASTN, short=True
    )
    if primer.name not in hsps:
        return []

    valid_locations = []
    for chrom, chrom_hsps in hsps[primer.name].items():
        chrom_hsps = filter_hsps_by_mismatch(
            chrom_hsps, max_mismatch=max_mismatches, query_len=len(primer.seq)
        )
        for hsp in chrom_hsps:
            valid_locations.append(
                _GenomeLocation.from_hsp(
                    chrom,
                    hsp.subject_start,
                    hsp.subject_end,
                    hsp.query_strand,
                    hsp.subject_strand,
                    primer.name,
                )
            )
    return valid_locations


def _create_pcr_products(
    primer1_locations,
    primer2_locations,
    max_product_size,
    allow_overlapping_primer_products,
    allow_one_primer_products,
):
    locations = [
        (loc1, loc2)
        for loc1, loc2 in itertools.combinations_with_replacement(
            [*primer1_locations, *primer2_locations], 2
        )
        if _is_valid_product(
            loc1,
            loc2,
            max_product_size,
            allow_overlapping_primer_products,
            allow_one_primer_products,
        )
    ]

    products = []
    for location in locations:
        chrom = location[0].chrom
        poss = [location[0].start, location[0].end, location[1].start, location[1].end]
        start = min(poss)
        end = max(poss)
        products.append(PCRProduct(chrom, start, end))
    return products


def run_ipcr(
    primer1: str,
    primer2: str,
    genome_fasta: Path | None,
    genome_blast_db: Path | None = None,
    max_mismatches=MAX_PRIMER_MISMATCHES,
    max_product_size=MAX_PRODUCT_SIZE,
    allow_one_primer_products=False,
    allow_overlapping_primer_products=False,
):
    if genome_fasta is not None and genome_blast_db is not None:
        raise ValueError(
            "You should provide either genome_fasta or genome_blast_db, but not both"
        )
    elif genome_fasta is None and genome_blast_db is None:
        raise ValueError("You should provide either genome_fasta or genome_blast_db")

    p1 = "primer1"
    p2 = "primer2"
    primer1 = Seq(p1, primer1)
    primer2 = Seq(p2, primer2)
    with TemporaryDirectory() as temp_dir:
        if genome_blast_db:
            db_path = genome_blast_db
        else:
            db_path = prepare_blast(
                genome_fasta, db_type=BlastDbType.NUCL, out_db_dir=Path(temp_dir)
            )["db_path"]

        primer1_valid_locations = _get_valid_locations_for_primer(
            primer1,
            db_path,
            max_mismatches,
        )
        primer2_valid_locations = _get_valid_locations_for_primer(
            primer2,
            db_path,
            max_mismatches,
        )
    products = _create_pcr_products(
        primer1_valid_locations,
        primer2_valid_locations,
        max_product_size,
        allow_overlapping_primer_products=allow_overlapping_primer_products,
        allow_one_primer_products=allow_one_primer_products,
    )
    return products
