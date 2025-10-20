"""
This Streamlit application allows you to deconvolute mass spectra using the msdecon package.
It takes input spectra in (mz intensity) format, filters by intensity, and then uses the
deconvolute function from msdecon to identify monoisotopic peaks and assign charge states.
"""

import numpy as np
import streamlit as st
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import plotly.graph_objects as go
import pandas as pd


import peptacular as pt
st.set_page_config(
    page_title="MsDecon",
    page_icon=":bar_chart:",
    menu_items={
        'Get Help': 'https://github.com/pgarrett-scripps/MsDecon',
        'Report a bug': "https://github.com/pgarrett-scripps/MsDecon/issues",
        'About': "# This is a Streamlit app for deconvoluting mass spectra."
    }
)

with st.sidebar:
    st.title('Deconvolute Mass Spectra :bar_chart:')

    st.caption("""
    **MsDecon** deconvolutes mass spectra by identifying the base peak of an isotopic distribution and assigning charge states. 
    Paste your centroided spectra, adjust the settings, and visualize the results. Source code: https://github.com/pgarrett-scripps/MsDecon
    """)

    # read default spectra 'default_spectra.txt'
    default_spectra = open('default_spectra.txt', 'r').read()

    # mz intensity\n
    c1, c2 = st.columns(2)

    spectra = c1.text_area('Paste spectra here', default_spectra, height=160,
                           help='Paste the spectra in MS2 format: mz intensity\n')


    c2.caption("Tolerance")
    tolerance_type = c2.selectbox('Type', ['ppm', 'da'], index=0,
                                  help='The type of tolerance to use for deconvolution')

    if tolerance_type == 'ppm':
        deconvolute_tolerance = c2.number_input('Value', value=25, step=1, key='ppm',
                                                help='The tolerance to apply when assigning isotopic peaks')
    else:
        deconvolute_tolerance = c2.number_input('Value', value=0.005, step=0.001, key='da',
                                                help='The tolerance to apply when assigning isotopic peaks')

    min_intensity = st.number_input('Min intensity', value=0)


    min_charge, max_charge = st.slider('Charge range', 1, 10, (1, 2),
                                       help='The min and max charge to consider for deconvolution.')

    with st.expander('Advanced Options'):
        st.caption(
            "Advanced deconvolution settings. The 'Min Intensity' filter is applied before deconvolution to exclude"
            " low-intensity peaks. The 'Min Left/Right Intensity Decrease' controls how the algorithm traces "
            "isotopic envelopes, starting from the most intense peak. Generally, peak intensities decrease "
            "gradually to the left (toward lower m/z), while larger drops are expected to the right ("
            "toward higher m/z).")
        c1, c2 = st.columns(2)
        max_left_intensity_decrease = c1.number_input('Max intensity decrease (Left)', value=0.65,
                                                      help='Maximum intensity decrease when navigating left (lower m/z)')
        max_right_intensity_decrease = c2.number_input('Max intensity decrease (Right)', value=0.95,
                                                       help='Maximum intensity decrease when navigating right (higher m/z)')
        charge_carrier = c1.number_input('Charge carrier mass', value=pt.PROTON_MASS, step=None, format="%0.6f",
                                         help='The mass of the charge carrier (e.g., proton for ESI-MS)')
        isotope_mass_offset = c2.number_input('Isotope mass', value=pt.C13_NEUTRON_MASS, step=None, format="%0.6f",
                                                help='The expected mass difference between isotopes '
                                                     '(Less than neutron due to differences in nuclear binding '
                                                     'energies of different elements).')

        show_gap_statistic = st.checkbox('Show Isotope Gap Statistics', value=False,
                                         help='Display a histogram of the isotope gaps between peaks in the deconvoluted spectrum')

    st.caption("""    
    #### Interactive Plot:  
    - You can **Zoom**, **Pan**, and **Hover** over peaks for more information. **Double-click** on the plot to reset the view.    
    - The legend is also interactive: **Double-click** on a category to isolate it. **Single-Click** on a category to hide/show it.  
    """)
    #deconvolute_min_intensity = st.number_input('Deconvolution min intensity', value=0.1)

# Validate and parse input
lines = [line.strip() for line in spectra.split('\n') if line.strip()]
try:
    mz_array = [float(line.split()[0]) for line in lines]
    intensity_array = [float(line.split()[1]) for line in lines]
except (ValueError, IndexError):
    st.error("Invalid input format. Please ensure each line contains 'mz intensity' separated by a space.")
    st.stop()

peaks = [(mz, intensity) for mz, intensity in zip(mz_array, intensity_array)]

# Filter by minimum intensity and sort
peaks = list(filter(lambda x: x[1] >= min_intensity, peaks))
if len(peaks) == 0:
    st.stop()

if len(peaks) > 10_000:
    st.warning('Too many peaks. Please reduce the number of peaks')
    st.stop()

# sort peaks by mz
peaks = sorted(peaks, key=lambda x: x[0])


def run_deconvolution(peaks, charge_range, tolerance, tolerance_type, max_left_decrease, max_right_decrease):
    return pt.deconvolute(
        peaks,
        charge_range=charge_range,
        tolerance=tolerance,
        tolerance_type=tolerance_type,
        max_left_decrease=max_left_decrease,
        max_right_decrease=max_right_decrease,
        isotope_mass=isotope_mass_offset
    )


dpeaks = run_deconvolution(
    peaks,
    (min_charge, max_charge),
    deconvolute_tolerance,
    tolerance_type,
    max_left_intensity_decrease,
    max_right_intensity_decrease
)

peaks_df = pd.DataFrame()
peaks_df['base_mz'] = [dpeak.base_peak.mz for dpeak in dpeaks]
peaks_df['total_intensity'] = [dpeak.total_intensity for dpeak in dpeaks]
peaks_df['base_neutral_mass'] = [peak.base_peak_neutral_mass for peak in dpeaks]
peaks_df['charge'] = [peak.charge for peak in dpeaks]
peaks_df['num_peaks'] = [peak.num_peaks for peak in dpeaks]
peaks_df['peak_mzs'] = [';'.join([f'{p.mz:.4f}' for p in peak.peaks]) for peak in dpeaks]
peaks_df['peak_intensities'] = [';'.join([f'{p.intensity:.1f}' for p in peak.peaks]) for peak in dpeaks]
peaks_df['gap_ppm_errors'] = [';'.join([f'{gap:.2f}' for gap in peak.isotope_gaps_ppm_error]) for peak in dpeaks]

# Handle missing charges by filling with 0 or a placeholder value
peaks_df['charge'] = peaks_df['charge'].fillna(0)  # Use 0 or a designated value for None

fig = go.Figure()

fig.add_trace(go.Scatter(
    x=sum([[mz, mz, None] for mz, intensity in peaks], []),
    y=sum([[0, intensity, None] for mz, intensity in peaks], []),
    mode='lines',
    name='Raw Spectrum',
    line=dict(color='gray'),
    opacity=0.3,
))

# Normalize charge states to map to a colormap
unique_charges = np.unique(peaks_df['charge'])
norm = mcolors.Normalize(vmin=min(unique_charges), vmax=max(unique_charges))
colormap = plt.get_cmap('viridis')  # Choose a colormap ('viridis', 'plasma', etc.)

# Assign a distinct color for None or 0 (unassigned charges)
special_color = '#efc7a0'  # Light gray for None/NaN charges

# Generate colors for each charge state
charge_colors = {charge: mcolors.rgb2hex(colormap(norm(charge))) for charge in unique_charges}
charge_colors[0] = special_color  # Explicitly assign color for charge = 0 (was None)

# Map colors to charge states
colors = [charge_colors.get(charge, special_color) for charge in peaks_df['charge']]

# Plot each charge state as a separate trace
for charge in unique_charges:
    mask = peaks_df['charge'] == charge
    mz_values = peaks_df['base_mz'][mask]
    intensity_values = peaks_df['total_intensity'][mask]

    fig.add_trace(go.Scatter(
        x=sum([[mz, mz, None] for mz in mz_values], []),
        y=sum([[0, intensity, None] for intensity in intensity_values], []),
        mode='lines',
        name=f'+{int(charge)} Peaks' if charge > 0 else '+? Peaks',
        line=dict(color=charge_colors[charge], width=2),
        hovertext=sum([
            [
                f'Charge: {charge}<br>Base m/z: {mz}<br>Total Intensity: {intensity}',
                f'Charge: {charge}<br>Base m/z: {mz}<br>Total Intensity: {intensity}',
                None
            ]
            for mz, intensity in zip(
                mz_values,
                intensity_values,
            )
        ], [])
    ))

fig.update_layout(
    title='Deconvoluted Spectrum',
    width=800,
    height=450,
    legend=dict(
        x=0.02,  # Horizontal position (0 = left, 1 = right)
        y=0.98,  # Vertical position (0 = bottom, 1 = top)
        xanchor='left',
        yanchor='top',
        bgcolor='rgba(255,255,255,0.5)',  # Semi-transparent background
    )

)

# show peaks in original spectrum (count) and then in the deconvoluted spectrum

st.subheader('Deconvolution Results', divider=True)

c1, c2, c3 = st.columns(3)
c1.metric('Original Peaks', len(peaks),
          help='Number of peaks in the original spectrum (After filtering by min intensity)')
c2.metric('Deconvoluted Peaks', len(peaks_df), help='Number of peaks in the deconvoluted spectrum')

# peaks with valid charge
valid_peaks = peaks_df[peaks_df['charge'] != 0]
c3.metric('Valid Peaks', len(valid_peaks), help='Number of peaks with a valid charge state')


@st.fragment
def display():
    min_peak_mz = min([peak[0] for peak in peaks])
    max_peak_mz = max([peak[0] for peak in peaks])

    filter_min_mz, filter_max_mz = st.slider('Select a range of m/z', min_peak_mz, max_peak_mz,
                                             (min_peak_mz, max_peak_mz))

    peaks_df_filtered = peaks_df[
        (peaks_df['base_mz'] >= filter_min_mz) & (peaks_df['base_mz'] <= filter_max_mz)]

    if len(peaks_df_filtered) == 0:
        st.error('No peaks found in the selected m/z range')
        st.stop()

    max_intensity = max(peaks_df_filtered['total_intensity'])
    # zoom into plotly plot and filter df

    fig.update_layout(xaxis=dict(range=[filter_min_mz, filter_max_mz]), yaxis=dict(range=[0, max_intensity]))
    st.plotly_chart(fig, use_container_width=True)

    # set 0 chareg state to None:
    peaks_df_filtered['charge'] = peaks_df_filtered['charge'].replace(0, None)

    all_gap_errors = []
    for errors in peaks_df_filtered['gap_ppm_errors']:
        if errors:
            all_gap_errors.extend([float(error) for error in errors.split(';')])

    st.dataframe(peaks_df_filtered, use_container_width=True, hide_index=True)

    # download table
    st.download_button(
        label="Download Table",
        data=peaks_df_filtered.to_csv(index=False),
        file_name="deconvoluted_peaks.csv",
        mime="text/csv",
        use_container_width=True,
        type='primary'
    )

    # plot histogram of isotope gaps
    if show_gap_statistic:
        st.subheader('Isotope Gap Statistics', divider=True)
        # write a metric for eman, median and std of the gap errors
        c1, c2, c3 = st.columns(3)
        c1.metric('Mean Gap Error', f"{np.mean(all_gap_errors):.2f} ppm")
        c2.metric('Median Gap Error', f"{np.median(all_gap_errors):.2f} ppm")
        c3.metric('Std. Dev. Gap Error', f"{np.std(all_gap_errors):.2f} ppm")

        fig_hist = go.Figure()
        fig_hist.add_trace(go.Histogram(
            x=all_gap_errors,
            nbinsx=50,
            marker_color='blue',
            opacity=0.75
        ))

        st.plotly_chart(fig_hist)


display()