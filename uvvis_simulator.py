# -*- coding: utf-8 -*-
"""
UV-Vis Spectrum Simulator
Created on Wed Apr 29 2026
@author: ralosant
"""

import argparse
from datetime import datetime
import glob
import logging
from pathlib import Path
import sys
import cclib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


class LogReader:
    """
    Reads and extracts data from files 
    """
    @staticmethod
    def detect_software(file):
        """
        Identifies the software by reading the file header.

        Args:
            file (str): Path to the file.

        Returns:
            str: 'orca', 'molcas', or 'cclib_target'.
        """
        with open(file, 'r', errors='ignore', encoding="utf-8") as f:
            h = "".join([f.readline() for _ in range(50)])
        if "* O   R   C   A *" in h:
            return "orca"
        if "MOLCAS" in h or "OpenMolcas" in h:
            return "molcas"
        return "cclib_target"

    def read(self, file, n_states):
        """
        Calls the proper reader for each file.

        Args:
            file (str): Path to the file.

        Returns:
            run qm.parser.
        """
        sw = self.detect_software(file)
        if sw == "orca":
            return self._read_orca(file, n_states)
        if sw == "molcas":
            return self._read_molcas(file, n_states)
        return self._read_cclib(file)

    def _read_orca(self, file, n_states):
        """
        ORCA log reader with column detection logic (nm vs eV).

        Args:
            file (str): Path to the ORCA .log file.
            n_states (int): Number of excited states to read.

        Returns:
            tuple: (energies, oscillator_strengths) lists, or (None, None) if failed.
        """
        try:
            with open(file, encoding="utf-8") as data:
                lines = data.read().split("\n")

            target = (
                "ABSORPTION SPECTRUM VIA TRANSITION ELECTRIC DIPOLE MOMENTS"
            )
            n1_line = next(
                (i for i, line in enumerate(lines) if target in line), None
            )
            if n1_line is None:
                return None, None

            # Initial attempt: Columns 5 and 6
            energy = [
                float(lines[n1_line + 5 + i].split()[5])
                for i in range(n_states)
            ]
            f = [
                float(lines[n1_line + 5 + i].split()[6])
                for i in range(n_states)
            ]

            # Check < 100: If the first value is too low (old version)
            # Switch to alternative columns (2 and 3).
            if energy[0] < 100:
                energy = [
                    float(lines[n1_line + 5 + i].split()[2])
                    for i in range(n_states)
                ]
                f = [
                    float(lines[n1_line + 5 + i].split()[3])
                    for i in range(n_states)
                ]

            return energy, f
        except (pd.errors.ParserError, IOError, IndexError, ValueError) as e:
            logging.error("Error parsing ORCA in %s: %s", file, e)
            return None, None

    def _read_molcas(self, file, n_states):
        """
        Manual reader for Molcas/OpenMolcas output files.

        Args:
            file (str): Path to the Molcas output file.
            n_states (int): Number of excited states to read.

        Returns:
            tuple: (energies, oscillator_strengths) or (None, None) if failed.
        """
        try:
            with open(file, encoding="utf-8") as data:
                lines = list(data.read().split())

            n0_line = 0
            # 1. Search for energies (nm)
            # Uses the occurrence after 'D:o,' for nm
            n1_line = lines.index("D:o,", n0_line + 1)
            energy = [
                (1 / float(lines[n1_line + 9 + i * 4])) * 10**7
                for i in range(n_states)
            ]

            # 2. Search for oscillator strengths (f)
            # Skips '(sec-1)' headers according to original logic
            n1_line = lines.index("(sec-1)", n1_line + 1)
            n1_line = lines.index("(sec-1)", n1_line + 1)

            exc_num = [
                lines[n1_line + 2 + i * 7] for i in range(n_states)
            ]
            _n_states_rel = n_states
            for i in range(n_states):
                if exc_num[i] != "1":
                    _n_states_rel -= 1

            # Extract valid oscillator strengths
            f_vals = [
                float(lines[n1_line + 4 + i * 7]) for i in range(_n_states_rel)
            ]

            # Pad with zeros if there are fewer relative states than n_states
            f = f_vals + [0.0] * (n_states - _n_states_rel)

            # Trim or adjust to match exact length
            return energy[: n_states], f[: n_states]

        except (pd.errors.ParserError, IOError, IndexError, ValueError) as e:
            print(f"Error parsing Molcas in ({file}): {e}")
            return None, None

    def _read_cclib(self, file):
        """
        General reader using cclib for multiple software formats.

        Args:
            file (str): Path to the log file.

        Returns:
            tuple: (energies, oscillator_strengths) or (None, None).
        """
        try:
            data = cclib.io.ccread(file)
            if hasattr(data, 'etoscs'):
                return (1e7 / data.etenergies).tolist(), data.etoscs.tolist()
        except (pd.errors.ParserError, IOError, IndexError, ValueError):
            pass
        return None, None


class UVVisSimulator:
    """
    # =============================================================================
    # UV-Vis simulator to plot QM data parsed with cclib or manually
    # =============================================================================
    """

    def __init__(self, sigma, n_states, weight, x_range):
        """
        Initializes the UV-Vis simulator.

        Args:
            sigma (float): Bandwidth (Half-Width at Half-Maximum) in eV.
            weight_g (float): Gaussian weight (0 to 1) for the Pseudo-Voigt profile.
            x_range (tuple): Range for the X-axis (start_nm, end_nm, num_points).
        """
        self.config = {
            'sigma': sigma, 
            'weight': weight,
            'n_states': n_states
        }
        start, end = x_range
        num_points = int(abs(end - start) + 1)
        self.x = np.linspace(start, end, num_points)
        self.data = {
            'results': [],
            'total_list': [],
            'avg': None,
            'std': None
        }
        self.reader = LogReader()
        logging.basicConfig(level=logging.INFO, format="%(message)s")

    def _pseudo_voigt(self, energy_nm, f_osc):
        """
        Vectorized convolution engine for a single peak.

        Args:
            energy_nm (float): Peak position in nm.
            f_osc (float): Oscillator strength.

        Returns:
            np.array: Intensity values across the X range.
        """
        dev = self.config['sigma'] * 8065.54
        const = 1.3062974e8
        diff = (1e7 / self.x) - (1e7 / energy_nm)
        gau = (const * f_osc / dev) * np.exp(-((diff / dev) ** 2))
        loren = (const * (2 / np.pi) * f_osc * dev) / (4 * diff**2 + dev**2)
        return (self.config['weight'] * gau) + ((1 - self.config['weight']) * loren)

    def process_files(self, file_list):
        """
        Processes a complete list of files.

        Args:
            file_list (list): List of file paths.
        """

        for file_path in file_list:
            en, fosc = self.reader.read(file_path, self.config['n_states'])

            if en and fosc:
                # Ensure fosc is not longer than en
                fosc = fosc[: len(en)]
                peaks = [self._pseudo_voigt(e, f) for e, f in zip(en, fosc)]
                total = np.sum(peaks, axis=0)
                self.data['total_list'].append(total)
                self.data['results'].append(
                    {
                        "name": file_path,
                        "total": total,
                        "peaks": peaks,
                        "en_raw": en,
                        "f_raw": fosc,
                    }
                )

        if not self.data['total_list']:
            return

        print("Generating final reports...")
        matrix = np.array(self.data['total_list'])
        if matrix.ndim == 1:
            matrix = matrix[np.newaxis, :]
        self.data['avg'] = np.mean(matrix, axis=0)
        self.data['std'] = np.std(matrix, axis=0)

    def plotting(self, png_out, n_comp, avg_only):
        """
        Plotting results (PNG plot).

        Args:
            png_out (str): Output PNG filename.
            n_comp (int): Number of individual components to plot per series.
            avg_only (bool): If True, only the global average is plotted.
        """
        plt.figure(figsize=(10, 6))
        if not avg_only:
            if len(self.data['results']) > 5:
                alpha, lw = 0.1, 0.3
                n_comp = False
                show_labels = False
            else:
                alpha, lw = 0.4, 1.5
                show_labels = True

            for res in self.data['results']:
                label = (
                    f"Total {Path(res['name']).stem}" if show_labels else None
                )
                (l,) = plt.plot(
                    self.x, res["total"], alpha=alpha, lw=lw, label=label
                )
                color = l.get_color()

                # Plot individual components
                if n_comp:
                    for p in res["peaks"][:n_comp]:
                        if p.max() >= 200:
                            plt.plot(self.x, p, ls=":", lw=0.8,
                                     color=color, alpha=0.3)

        # Plot Average and Deviation
        plt.plot(
            self.x,
            self.data['avg'],
            color="black",
            lw=3,
            label="GLOBAL AVERAGE",
            zorder=10,
        )
        plt.fill_between(
            self.x,
            self.data['avg'] - self.data['std'],
            self.data['avg'] + self.data['std'],
            color="black",
            alpha=0.15,
            label="Standard Deviation (σ)",
            zorder=5,
        )

        # plt.title("UV-Vis Simulation: Average, Dispersion, and Contributions")
        plt.xlabel("Wavelength (nm)")
        plt.ylabel(r"Extinction Coeff. ($\epsilon$)")
        plt.legend(loc="upper right", ncol=2, fontsize="small")
        plt.grid(False)
        plt.tight_layout()

        plt.savefig(png_out, transparent=True, dpi=300, bbox_inches="tight")

    def _prepare_excel_data(self, n_comp):
        """Helper method to reduce local variable count in export."""
        dict_details = {"Wavelength (nm)": self.x}
        dict_raw = {}

        for i, res in enumerate(self.data['results']):
            tag = f"S{i + 1}"
            dict_details[f"{tag}_TOTAL"] = res["total"]
            for j, p in enumerate(res["peaks"][:n_comp]):
                dict_details[f"{tag}_E{j + 1}"] = p

            dict_raw[f"{tag}_E_nm"] = pd.Series(res["en_raw"])
            dict_raw[f"{tag}_f_osc"] = pd.Series(res["f_raw"])

        return pd.DataFrame(dict_details), pd.DataFrame(dict_raw)

    def export(self, excel_out="Results.xlsx", png_out="Spectra.png", n_comp=3, avg_only=False):
        """Generates final reports (Excel and PNG plot)."""
        if not self.data['total_list']:
            logging.warning("No data found to export.")
            return

        # 1. Statistical Calculations
        matrix = np.atleast_2d(self.data['total_list'])
        self.data['avg'] = np.mean(matrix, axis=0)
        self.data['std'] = np.std(matrix, axis=0)

        # 2. Get DataFrames from helper
        df_details, df_raw = self._prepare_excel_data(n_comp)

        # 3. Determine Model String
        w = self.config['weight']
        model_str = "Pure Gaussian" if w == 1.0 else "Pure Lorentzian" if w == 0.0 \
                    else f"Pseudo-Voigt ({w * 100:.1f}% G)"

        # 4. Write to Excel
        with pd.ExcelWriter(excel_out, engine="xlsxwriter") as writer:
            # Configuration Sheet
            pd.DataFrame({
                "Parameter": ["Date", "Sigma (eV)", "Model", "Files"],
                "Value": [datetime.now().strftime("%Y-%m-%d %H:%M"), 
                          self.config['sigma'], model_str, len(self.data['results'])]
            }).to_excel(writer, sheet_name="Configuration", index=False)

            # Summary and Data Sheets
            pd.DataFrame({
                "Wavelength (nm)": self.x, 
                "Average": self.data['avg'], 
                "SD": self.data['std']}
            ).to_excel(writer, sheet_name="Summary", index=False)

            df_details.to_excel(writer, sheet_name="Series_and_Components", index=False)
            df_raw.to_excel(writer, sheet_name="Raw_Data", index=False)

        # 5. Generate Plot
        self.plotting(png_out, n_comp, avg_only)
        print(f"Simulation saved to {excel_out} and {png_out}")



def main():
    """
    Main function to run the UV-Vis Simulator from the command line.
    """
    parser = argparse.ArgumentParser(
        description="UV-Vis Spectrum Simulator from ORCA, Gaussian, and Molcas outputs."
    )

    # Input files
    parser.add_argument(
        "files",
        nargs="*",
        help="Input log files (e.g., *.log, *.out). If empty, searches current folder.",
    )

    # Simulation parameters
    parser.add_argument(
        "-s",
        "--sigma",
        type=float,
        default=0.33,
        help="Bandwidth (Half-bandwidth) in eV. Default: 0.33",
    )
    parser.add_argument(
        "-n",
        "--states",
        type=int,
        default=3,
        help="Number of excited states to extract per file. Default: 3",
    )

    # Output options
    parser.add_argument(
        "-o",
        "--out",
        default="UVVis_Simulation",
        help="Base name for output files (Excel and PNG). Default: UVVis_Simulation",
    )
    parser.add_argument(
        "-r",
        "--range",
        type=float,
        nargs=2,
        default=[200, 800],
        help="X-axis range: start and end points (e.g., 200 800). Resolution is 1 nm by default.",
    )

    # Plotting options
    parser.add_argument(
        "--ncomp",
        type=int,
        default=3,
        help="Number of individual components to plot per series. Default: 3",
    )
    parser.add_argument(
        "--avg-only",
        action="store_true",
        help="Plot only the global average, hiding individual series.",
    )

    # Broadening model (Mutually exclusive group)
    group = parser.add_mutually_exclusive_group()
    group.add_argument(
        "-G", "--GAU", action="store_true", help="Use pure Gaussian functions."
    )
    group.add_argument(
        "-L",
        "--LOREN",
        action="store_true",
        help="Use pure Lorentzian functions.",
    )
    group.add_argument(
        "-P",
        "--PSVoigt",
        type=float,
        help="Use Pseudo-Voigt model with specified Gaussian weight (0.0 to 1.0).",
    )

    args = parser.parse_args()

    # Determine Gaussian weight logic
    if args.GAU:
        weight = 1.0
    elif args.LOREN:
        weight = 0.0
    elif args.PSVoigt is not None:
        if not 0.0 <= args.PSVoigt <= 1.0:
            parser.error(
                "Pseudo-Voigt weight (-P) must be between 0.0 and 1.0"
            )
        weight = args.PSVoigt
    else:
        weight = 0.85  # Default weight

    # File selection logic
    file_list = (
        args.files if args.files else (glob.glob("*.log") + glob.glob("*.out"))
    )

    if not file_list:
        print(
            "Error: No files provided and no .log or .out files found in the current directory."
        )
        sys.exit(1)

    # 1. Initialize simulator
    sim = UVVisSimulator(sigma=args.sigma, n_states=args.states,
                         weight=weight, x_range=args.range)

    # 2. Process all files
    sim.process_files(file_list)

    # 3. Export results (Excel and Plot)
    sim.export(
        excel_out=f"{args.out}.xlsx",
        png_out=f"{args.out}.png",
        n_comp=args.ncomp,
        avg_only=args.avg_only,
    )


if __name__ == "__main__":
    main()
