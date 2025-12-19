from test_blast import TEST_DATA_DIR
from ipcr import run_ipcr, OverlappingPrimerProduct, OnePrimerProduct
from ipcr.ipcr import _GenomeLocation, _is_valid_product

import pytest

COMPLEMENTARY_NUCL = {"A": "T", "C": "G", "G": "C", "T": "A"}


def rev_comp(seq):
    return "".join(COMPLEMENTARY_NUCL[nucl] for nucl in reversed(seq))


def test_genome_location():
    loc = _GenomeLocation.from_hsp("chrom1", 10, 20, 1, 1, "primer")
    assert loc.fwd_strand
    loc = _GenomeLocation.from_hsp("chrom1", 10, 20, -1, 1, "primer")
    assert not loc.fwd_strand
    loc = _GenomeLocation.from_hsp("chrom1", 10, 20, 1, -1, "primer")
    assert not loc.fwd_strand
    loc = _GenomeLocation.from_hsp("chrom1", 10, 20, -1, -1, "primer")
    assert loc.fwd_strand


def test_is_valid_product():
    # ---->
    #        <------
    loc1 = _GenomeLocation("chrom1", 10, 20, True, "primer1")
    loc2 = _GenomeLocation("chrom1", 30, 40, False, "primer2")
    assert _is_valid_product(loc1, loc2)
    assert _is_valid_product(loc2, loc1)

    # ---->
    #                                  <------
    loc1 = _GenomeLocation("chrom1", 10, 20, True, "primer1")
    loc2 = _GenomeLocation("chrom1", 11000, 11010, False, "primer2")
    assert not _is_valid_product(loc1, loc2)

    # <----
    #        ------>
    loc1 = _GenomeLocation("chrom1", 10, 20, False, "primer1")
    loc2 = _GenomeLocation("chrom1", 30, 40, True, "primer2")
    assert not _is_valid_product(loc1, loc2)

    #              ------->
    # <------
    loc2 = _GenomeLocation("chrom1", 10, 20, False, "primer1")
    loc1 = _GenomeLocation("chrom1", 30, 40, True, "primer2")
    assert not _is_valid_product(loc1, loc2)

    # ------>
    #                  --------->
    loc1 = _GenomeLocation("chrom1", 10, 20, True, "primer1")
    loc2 = _GenomeLocation("chrom1", 30, 40, True, "primer2")
    assert not _is_valid_product(loc1, loc2)

    # <-------
    #                   <--------
    loc1 = _GenomeLocation("chrom1", 10, 20, False, "primer1")
    loc2 = _GenomeLocation("chrom1", 30, 40, False, "primer2")
    assert not _is_valid_product(loc1, loc2)

    # ------->
    #      <--------
    loc1 = _GenomeLocation("chrom1", 10, 20, True, "primer1")
    loc2 = _GenomeLocation("chrom1", 15, 30, False, "primer2")
    with pytest.raises(OverlappingPrimerProduct):
        _is_valid_product(loc1, loc2)

    #    ------->
    #  <-------
    loc1 = _GenomeLocation("chrom1", 10, 20, True, "primer1")
    loc2 = _GenomeLocation("chrom1", 5, 15, False, "primer2")
    assert not _is_valid_product(loc1, loc2)


def test_prepare_blast_db():
    seq_fasta = TEST_DATA_DIR / "fxn_human.fasta"
    primer1 = "TAGTGCTGTTTCTCCCACATATTC"
    primer2 = "CTGCACAGGAGGCTGAGACAGGAGGA"
    primer2_rev = rev_comp(primer2)

    products = run_ipcr(primer1, primer2_rev, seq_fasta)
    assert len(products) == 1


def test_one_primer_product():
    seq_fasta = TEST_DATA_DIR / "seq_with_duplication.fasta"
    primer1 = "TAGTGCTGTTTCTCCCACATATTC"
    primer2 = "CTGCACAGGAGGCTGAGACAGGAGGA"

    with pytest.raises(OnePrimerProduct):
        run_ipcr(primer1, primer2, seq_fasta)
