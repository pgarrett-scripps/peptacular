import streamlit as st
import peptacular as pt
from peptacular.plotting.tmp import spectrum
from peptacular.plotting.spectrum import generate_sprectrum_table, plot_spectrum_table

spectrum = [(float(peak.split(' ')[0]), float(peak.split(' ')[1])) for peak in spectrum.split('\n')]

fragments = pt.Fragmenter("[164.0700]-FDSFGDLSSASAIM[16]GNPK", True).fragment(ion_types=list("abcxyz"),
                                                    charges=[1,2,3],
                                                    isotopes=[0, 1, 2],
                                                    water_loss=False,
                                                    ammonia_loss=False,
                                                    losses=None,
                                                    max_losses=1,
                                                    return_type="fragment",
                                                    precision=None)

spectrum_df = generate_sprectrum_table(spectrum, fragments, 'ppm', 500.0)

print(spectrum_df[spectrum_df['charge'] == 1])

fig = plot_spectrum_table(spectrum_df, title="Spectrum for ACDEFGHIKLMNPQRSTVWY (b)", xaxis_title="M/Z", yaxis_title="Intensity")

st.plotly_chart(fig, use_container_width=True)