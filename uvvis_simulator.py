# -*- coding: utf-8 -*-
"""Command-line interface for the UV-Vis simulator."""
import argparse
import glob
import sys
from uvvis_api import simulate_uvvis
sys.dont_write_bytecode = True


def build_parser():
    parser = argparse.ArgumentParser(description="UV-Vis Spectrum Simulator")
    parser.add_argument("files", nargs="*", help="QM output files")
    parser.add_argument("-s", "--sigma", type=float, default=0.33)
    parser.add_argument("-n", "--states", type=int, default=3)
    parser.add_argument("-o", "--out", default="UVVis_Simulation")
    parser.add_argument("-r", "--range", type=float, nargs=2, default=[200, 800],
                        metavar=("START_NM", "END_NM"))
    parser.add_argument("--ncomp", type=int, default=0)
    parser.add_argument("--mode",  choices=["individual", "average", "combined"], default="individual",
        help=("Plot mode: "
        "'individual' Plot separate files together (default), "
        "'average' (mean only), "
        "'combined' (mean + standard deviation)"
        ))
    profile = parser.add_mutually_exclusive_group()
    profile.add_argument("-G", "--GAU", action="store_true")
    profile.add_argument("-L", "--LOREN", action="store_true")
    profile.add_argument("-P", "--PSVoigt", type=float)
    return parser


def main(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)
    if args.GAU:
        weight = 1.0
    elif args.LOREN:
        weight = 0.0
    elif args.PSVoigt is not None:
        if not 0 <= args.PSVoigt <= 1:
            parser.error("-P/--PSVoigt must be between 0 and 1")
        weight = args.PSVoigt
    else:
        weight = 0.85
    files = args.files or glob.glob("*.log") + glob.glob("*.out")
    if not files:
        parser.error("No input files supplied or found")
    try:
        simulate_uvvis(files, sigma=args.sigma, n_states=args.states, weight=weight,
                       wavelength_range=args.range, n_components=args.ncomp,
                       plot_mode=args.mode, excel_out=f"{args.out}.xlsx",
                       png_out=f"{args.out}.png")
    except ValueError as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
