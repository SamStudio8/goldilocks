import pytest

import numpy as np

from goldilocks import Goldilocks
from goldilocks.strategies import BaseStrategy, MotifCounterStrategy

DATA_FAI = {
    "HOOT": {
        "file": "tests/dat/hoot.fa.fai"
    },
}

@pytest.fixture(params=[48, 96, 24])
def g(request):
    length=request.param
    if length == 48:
        MOTIFS = [
            "CGTGGGACATTTACGTGAGAGTTGTGGTTCATGCAAACCTAAGAAAGC",
            "CCCTTGGGACCGAGAAAGCGCATCAGTGCACTGCTCTTTGCGGACCGA",
            "TAAGATCGTCACATGTCCCACTGTCGCGAGTAGAAATAGCTTAGGAAA",
            "GAACGCGCAAATTACTTATAGGGAGTTAGCATGGTATCCTAAAAATCC",
            "CTGCTTTCTACAGCTCATCACGAGTACCGTATATCAAGGTACCTTTGC",
            "ACCCTTTCGCAAGGAACTCGCTCGCTCTCCGCTGGATATTACTGCAAG",
            "GGTGCCGATGGCGGTACGGCTCTAGCAAGGTACGCCTAGACGCAGGAC",
            "GGGTTCGTCTCGACATCCAGGTGTATAGATCCCCTGGATCTCAACGTT",
        ]
    elif length == 96:
        MOTIFS = [
            "CGTGGGACATTTACGTGAGAGTTGTGGTTCATGCAAACCTAAGAAAGCCCCTTGGGACCGAGAAAGCGCATCAGTGCACTGCTCTTTGCGGACCGA",
            "TAAGATCGTCACATGTCCCACTGTCGCGAGTAGAAATAGCTTAGGAAAGAACGCGCAAATTACTTATAGGGAGTTAGCATGGTATCCTAAAAATCC",
            "CTGCTTTCTACAGCTCATCACGAGTACCGTATATCAAGGTACCTTTGCACCCTTTCGCAAGGAACTCGCTCGCTCTCCGCTGGATATTACTGCAAG",
            "GGTGCCGATGGCGGTACGGCTCTAGCAAGGTACGCCTAGACGCAGGACGGGTTCGTCTCGACATCCAGGTGTATAGATCCCCTGGATCTCAACGTT",
        ]
    elif length == 24:
        MOTIFS = [
            "CGTGGGACATTTACGTGAGAGTTG",
            "TGGTTCATGCAAACCTAAGAAAGC",
            "CCCTTGGGACCGAGAAAGCGCATC",
            "AGTGCACTGCTCTTTGCGGACCGA",
            "TAAGATCGTCACATGTCCCACTGT",
            "CGCGAGTAGAAATAGCTTAGGAAA",
            "GAACGCGCAAATTACTTATAGGGA",
            "GTTAGCATGGTATCCTAAAAATCC",
            "CTGCTTTCTACAGCTCATCACGAG",
            "TACCGTATATCAAGGTACCTTTGC",
            "ACCCTTTCGCAAGGAACTCGCTCG",
            "CTCTCCGCTGGATATTACTGCAAG",
            "GGTGCCGATGGCGGTACGGCTCTA",
            "GCAAGGTACGCCTAGACGCAGGAC",
            "GGGTTCGTCTCGACATCCAGGTGT",
            "ATAGATCCCCTGGATCTCAACGTT",
        ]
    return Goldilocks(MotifCounterStrategy(MOTIFS), DATA_FAI, length=length, stride=length, is_faidx=True)

def test_group_track_counts_contents(g):
    """Ensure group_counts held by region metadata for each group-track
    combination (including the total-default group-track) match the
    expected number of bases.
    This test is somewhat cheating in that it fetches region metadata from
    the regions dictionary."""

    number_comparisons = 0
    ggroups = g.groups.keys()# + ["total"]
    ttracks = g.strategy.TRACKS# + ["default"]
    for group in ggroups:
        for t_id, track in enumerate(ttracks):
            for region_i, value in enumerate(g.counter_matrix[g._get_group_id(group), g._get_track_id(track),]):
                if region_i == t_id:
                    assert 1 == value
                else:
                    assert 0 == value

                number_comparisons += 1
    assert (len(g.groups.keys()) * len(g.strategy.TRACKS) * len(g.strategy.TRACKS)) == number_comparisons

def test_by_export_meta(g):
    g.query("max").export_meta()

def test_by_export_fasta(g, capsys):
    g.export_fasta()
    expected_lines = open("tests/dat/expected/hoot"+str(g.LENGTH)+".fa").readlines()+["\n"]
    assert len(expected_lines) > 0
    out, err = capsys.readouterr()
    num_lines = 0
    for i, line in enumerate(out.split("\n")):
        assert line+"\n" == expected_lines[i]
        num_lines += 1
    assert len(expected_lines) == num_lines

