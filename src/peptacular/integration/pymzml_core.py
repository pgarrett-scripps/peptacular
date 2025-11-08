from typing import Generator

from ..spectra.spectrum import Ms1Spectrum, Ms2Spectrum, Spectrum, SpectrumPeak


def parse_pymzml_spectrum(spectrum) -> Spectrum:
    """Parse pymzML spectrum object into Spectrum."""

    mzs = []
    intensities = []
    for mz, intensity in spectrum.peaks("centroided"):
        mzs.append(float(mz))
        intensities.append(float(intensity))

    if spectrum.ms_level == 1:
        return Ms1Spectrum(
            retention_time=float(spectrum.scan_time_in_minutes() or 0),
            ms1_scan_number=spectrum.ID,
            ion_mobility=float(spectrum.get("ion mobility drift time (ms)", 0) or 0)
            if spectrum.get("ion mobility drift time (ms)")
            else None,
            peaks=[
                SpectrumPeak(mz, intensity) for mz, intensity in zip(mzs, intensities)
            ],
        )
    elif spectrum.ms_level == 2:
        return Ms2Spectrum(
            retention_time=float(spectrum.scan_time_in_minutes() or 0),
            ms1_scan_number=int(spectrum.get("ms1 scan number", 0) or 0)
            if spectrum.get("ms1 scan number")
            else None,
            ion_mobility=float(spectrum.get("ion mobility drift time (ms)", 0) or 0)
            if spectrum.get("ion mobility drift time (ms)")
            else None,
            precursor_mz=float(spectrum.selected_precursors[0]["mz"])
            if spectrum.selected_precursors and "mz" in spectrum.selected_precursors[0]
            else None,
            precursor_charge=int(spectrum.selected_precursors[0]["charge"])
            if spectrum.selected_precursors
            and "charge" in spectrum.selected_precursors[0]
            and isinstance(spectrum.selected_precursors[0]["charge"], int)
            else None,
            ms2_scan_number=spectrum.ID,
            activation_method=spectrum.get("activation method"),
            isolation_window=(
                float(spectrum.get("isolation window target m/z", 0) or 0)
                - float(spectrum.get("isolation window lower offset", 0) or 0),
                float(spectrum.get("isolation window target m/z", 0) or 0)
                + float(spectrum.get("isolation window upper offset", 0) or 0),
            )
            if all(
                k in spectrum
                for k in [
                    "isolation window target m/z",
                    "isolation window lower offset",
                    "isolation window upper offset",
                ]
            )
            else None,
            peaks=[
                SpectrumPeak(mz, intensity) for mz, intensity in zip(mzs, intensities)
            ],
        )
    else:
        return Spectrum((mzs, intensities))


def parse_mzml(
    mzml_file: str, ms_level: int | None = None
) -> Generator[Spectrum, None, None]:
    """Parse mzML file and yield Spectrum objects."""
    import pymzml

    for spectrum in pymzml.run.Reader(mzml_file):
        if spectrum.ms_level != ms_level and ms_level is not None:
            continue

        spec = parse_pymzml_spectrum(spectrum)
        yield spec
