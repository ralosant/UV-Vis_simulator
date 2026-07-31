# -*- coding: utf-8 -*-
"""Public Python API for simulating UV-Vis spectra from QM output files."""

from datetime import datetime
import logging
from pathlib import Path

try:
    import cclib
except ImportError:  # ORCA/Molcas manual parsers remain available
    cclib = None
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

__all__ = ["LogReader", "UVVisSimulator", "simulate_uvvis"]


class LogReader:
    """Extract excitation wavelengths and oscillator strengths."""

    @staticmethod
    def detect_software(file):
        with open(file, "r", errors="ignore", encoding="utf-8") as stream:
            header = "".join(stream.readline() for _ in range(50))
        if "* O   R   C   A *" in header:
            return "orca"
        if "MOLCAS" in header or "OpenMolcas" in header:
            return "molcas"
        return "cclib_target"

    def read(self, file, n_states):
        software = self.detect_software(file)
        if software == "orca":
            return self._read_orca(file, n_states)
        if software == "molcas":
            return self._read_molcas(file, n_states)
        return self._read_cclib(file, n_states)

    @staticmethod
    def _read_orca(file, n_states):
        try:
            lines = Path(file).read_text(encoding="utf-8", errors="ignore").splitlines()
            target = "ABSORPTION SPECTRUM VIA TRANSITION ELECTRIC DIPOLE MOMENTS"
            start = next((i for i, line in enumerate(lines) if target in line), None)
            if start is None:
                return None, None
            rows = [lines[start + 5 + i].split() for i in range(n_states)]
            energy = [float(row[5]) for row in rows]
            oscillator = [float(row[6]) for row in rows]
            if energy and energy[0] < 100:
                energy = [float(row[2]) for row in rows]
                oscillator = [float(row[3]) for row in rows]
            return energy, oscillator
        except (OSError, IndexError, StopIteration, ValueError) as exc:
            logging.error("Error parsing ORCA file %s: %s", file, exc)
            return None, None

    @staticmethod
    def _read_molcas(file, n_states):
        try:
            tokens = Path(file).read_text(encoding="utf-8", errors="ignore").split()
            position = tokens.index("D:o,")
            energy = [1e7 / float(tokens[position + 9 + i * 4]) for i in range(n_states)]
            position = tokens.index("(sec-1)", position + 1)
            position = tokens.index("(sec-1)", position + 1)
            excitation_numbers = [tokens[position + 2 + i * 7] for i in range(n_states)]
            relevant = sum(value == "1" for value in excitation_numbers)
            oscillator = [float(tokens[position + 4 + i * 7]) for i in range(relevant)]
            oscillator.extend([0.0] * (n_states - relevant))
            return energy[:n_states], oscillator[:n_states]
        except (OSError, IndexError, ValueError) as exc:
            logging.error("Error parsing Molcas file %s: %s", file, exc)
            return None, None

    @staticmethod
    def _read_cclib(file, n_states):
        if cclib is None:
            logging.error("cclib is required to parse %s", file)
            return None, None
        try:
            data = cclib.io.ccread(file)
            if data is not None and hasattr(data, "etoscs") and hasattr(data, "etenergies"):
                energy = (1e7 / data.etenergies).tolist()[:n_states]
                oscillator = data.etoscs.tolist()[:n_states]
                return energy, oscillator
        except (OSError, AttributeError, TypeError, ValueError) as exc:
            logging.error("Error parsing %s with cclib: %s", file, exc)
        return None, None


class UVVisSimulator:
    """Low-level, reusable simulator for programmatic workflows."""

    def __init__(self, sigma=0.33, n_states=3, weight=0.85, x_range=(200, 800)):
        if sigma <= 0:
            raise ValueError("sigma must be greater than zero")
        if n_states <= 0:
            raise ValueError("n_states must be greater than zero")
        if not 0 <= weight <= 1:
            raise ValueError("weight must be between 0 and 1")
        start, end = x_range
        if start <= 0 or end <= 0 or start == end:
            raise ValueError("x_range must contain two distinct positive wavelengths")
        self.config = {"sigma": sigma, "weight": weight, "n_states": n_states}
        self.x = np.linspace(start, end, int(abs(end - start) + 1))
        self.data = {"results": [], "total_list": [], "avg": None, "std": None}
        self.reader = LogReader()

    def _pseudo_voigt(self, energy_nm, oscillator_strength):
        deviation = self.config["sigma"] * 8065.54
        constant = 1.3062974e8
        difference = (1e7 / self.x) - (1e7 / energy_nm)
        gaussian = (constant * oscillator_strength / deviation) * np.exp(
            -((difference / deviation) ** 2)
        )
        lorentzian = (
            constant * (2 / np.pi) * oscillator_strength * deviation
        ) / (4 * difference**2 + deviation**2)
        weight = self.config["weight"]
        return weight * gaussian + (1 - weight) * lorentzian

    def process_files(self, files):
        """Parse files and populate ``data``; return ``self`` for chaining."""
        for file in files:
            wavelengths, strengths = self.reader.read(str(file), self.config["n_states"])
            if wavelengths and strengths:
                strengths = strengths[:len(wavelengths)]
                peaks = [self._pseudo_voigt(e, f) for e, f in zip(wavelengths, strengths)]
                total = np.sum(peaks, axis=0)
                result = {
                    "name": str(file), "total": total, "peaks": peaks,
                    "en_raw": wavelengths, "f_raw": strengths,
                }
                self.data["results"].append(result)
                self.data["total_list"].append(total)
        return self

    def calculate_statistics(self):
        """Calculate and return global average and standard deviation."""
        if not self.data["total_list"]:
            raise ValueError("No spectra have been processed")
        matrix = np.atleast_2d(self.data["total_list"])
        self.data["avg"] = np.mean(matrix, axis=0)
        self.data["std"] = np.std(matrix, axis=0)
        return self.data["avg"], self.data["std"]

    def plot(self, output, n_components, mode="individual"):
        """Plot in ``combined``, ``average``, or ``individual`` mode."""
        if mode not in {"combined", "average", "individual"}:
            raise ValueError("mode must be 'combined', 'average', or 'individual'")
    
        if not self.data["results"]:
            raise ValueError("No spectra have been processed")
    
        if mode != "individual" and self.data["avg"] is None:
            self.calculate_statistics()
    
        plt.figure(figsize=(10, 6))
    
        if mode != "average":
            for result in self.data["results"]:
                label = Path(result["name"]).stem
                line, = plt.plot(
                    self.x,
                    result["total"],
                    lw=1.5,
                    alpha=0.9,
                    label=label
                )
    
                for peak in result["peaks"][:n_components]:
                    if peak.max() >= 200:
                        plt.plot(
                            self.x,
                            peak,
                            ":",
                            lw=0.8,
                            color=line.get_color(),
                            alpha=0.3
                        )
    
        if mode in {"average", "combined"}:
    
            plt.plot(
                self.x,
                self.data["avg"],
                color="black",
                lw=3,
                label="GLOBAL AVERAGE",
                zorder=10,
            )
    
            if mode == "combined":
                plt.fill_between(
                    self.x,
                    self.data["avg"] - self.data["std"],
                    self.data["avg"] + self.data["std"],
                    color="black",
                    alpha=0.15,
                    label="Standard Deviation (sigma)",
                    zorder=5,
                )
    
        plt.xlabel("Wavelength (nm)")
        plt.ylabel(r"Extinction Coeff. ($\epsilon$)")
        plt.legend(loc="upper right", ncol=2, fontsize="small")
        plt.tight_layout()
        plt.savefig(output, transparent=True, dpi=300, bbox_inches="tight")
        plt.close()
    
        return Path(output)
    
    
    def export_excel(self, output, n_components, mode="individual"):
        """Export configuration, spectra, components, and raw transitions."""
    
        if not self.data["results"]:
            raise ValueError("No spectra have been processed")
    
        if mode != "individual" and self.data["avg"] is None:
            self.calculate_statistics()
    
        details = {"Wavelength (nm)": self.x}
        raw = {}
    
        for result in self.data["results"]:
    
            tag = Path(result["name"]).stem
    
            details[f"{tag}_TOTAL"] = result["total"]
    
            for peak_index, peak in enumerate(
                result["peaks"][:n_components],
                start=1
            ):
                details[f"{tag}_E{peak_index}"] = peak
    
            raw[f"{tag}_E_nm"] = pd.Series(result["en_raw"])
            raw[f"{tag}_f_osc"] = pd.Series(result["f_raw"])
    
        weight = self.config["weight"]
    
        model = (
            "Pure Gaussian"
            if weight == 1
            else "Pure Lorentzian"
            if weight == 0
            else f"Pseudo-Voigt ({weight * 100:.1f}% G)"
        )
    
        with pd.ExcelWriter(output, engine="xlsxwriter") as writer:
    
            pd.DataFrame({
                "Parameter": [
                    "Date",
                    "Sigma (eV)",
                    "Model",
                    "Files",
                    "Plot mode",
                ],
                "Value": [
                    datetime.now().strftime("%Y-%m-%d %H:%M"),
                    self.config["sigma"],
                    model,
                    len(self.data["results"]),
                    mode,
                ],
            }).to_excel(
                writer,
                sheet_name="Configuration",
                index=False,
            )
    
            if mode == "average":
    
                pd.DataFrame({
                    "Wavelength (nm)": self.x,
                    "Average": self.data["avg"],
                }).to_excel(
                    writer,
                    sheet_name="Summary",
                    index=False,
                )
    
            elif mode == "combined":
    
                pd.DataFrame({
                    "Wavelength (nm)": self.x,
                    "Average": self.data["avg"],
                    "SD": self.data["std"],
                }).to_excel(
                    writer,
                    sheet_name="Summary",
                    index=False,
                )
    
            pd.DataFrame(details).to_excel(
                writer,
                sheet_name="Series_and_Components",
                index=False,
            )
    
            pd.DataFrame(raw).to_excel(
                writer,
                sheet_name="Raw_Data",
                index=False,
            )
    
        return Path(output)


def simulate_uvvis(files, *, sigma=0.33, n_states=3, weight=0.85,
                   wavelength_range=(200, 800), n_components=0,
                   plot_mode="individual", excel_out=None, png_out=None):
    """High-level Python API; return the populated :class:`UVVisSimulator`."""

    files = [str(file) for file in files]

    if not files:
        raise ValueError("At least one input file is required")

    simulator = UVVisSimulator(
        sigma,
        n_states,
        weight,
        wavelength_range,
    )

    simulator.process_files(files)

    if not simulator.data["results"]:
        raise ValueError("No spectra could be parsed from the supplied files")

    if plot_mode not in {"combined", "average", "individual"}:
        raise ValueError(
            "plot_mode must be 'combined', 'average', or 'individual'"
        )

    if plot_mode != "individual":
        simulator.calculate_statistics()

    if excel_out is not None:
        simulator.export_excel(
            excel_out,
            n_components,
            plot_mode,
        )

    if png_out is not None:
        simulator.plot(
            png_out,
            n_components,
            plot_mode,
        )

    return simulator
