from __future__ import absolute_import

from .goldilocks import Goldilocks
from .util import parse_si_bp
from . import strategies

import argparse
import re
import sys

def main():
    STRATEGIES = {
        "nuc": strategies.NucleotideCounterStrategy,
        "motif": strategies.MotifCounterStrategy,
        "gc": strategies.GCRatioStrategy,
        "ref": strategies.ReferenceConsensusStrategy,
    }
    FORMATS = ["bed", "circos", "melt", "table"] #TODO Gross duplication
    SORTS = ["min", "max", "mean", "median", "none"]

    if len(sys.argv) == 2:
        if sys.argv[1].lower() == "list":
            print("Available Strategies")
            for s in STRATEGIES:
                print("  * %s" % s)
            print
            print("Available Output Formats")
            for f in FORMATS:
                print("  * %s" % f)
            sys.exit(0)

    parser = argparse.ArgumentParser(description='Wrapper script for Goldilocks library.')
    parser.add_argument('strategy', choices=STRATEGIES.keys())
    parser.add_argument('sort', choices=SORTS)
    parser.add_argument('faidx', nargs='+')
    parser.add_argument('-t', '--tracks', nargs='+')
    parser.add_argument('-f', '--format', default="table", choices=FORMATS)
    parser.add_argument('-l', '--length', required=True, type=parse_si_bp)
    parser.add_argument('-s', '--stride', required=True, type=parse_si_bp)
    parser.add_argument('-@', '--processes', default=1, type=int)
    args = parser.parse_args()

    tracks = args.tracks
    if ',' in args.tracks:
        tracks = tracks.split(',')

    sequence_faidx = {}
    for i, idx in enumerate(args.faidx):
        sequence_faidx[i] = { "idx": idx }

    g = Goldilocks(STRATEGIES[args.strategy](tracks), sequence_faidx, length=args.length, stride=args.stride, processes=args.processes, is_faidx=True)

    if args.sort != "none":
        g.query(args.sort).export_meta(sep="\t")
    else:
        g.export_meta(sep="\t", fmt=args.format)
