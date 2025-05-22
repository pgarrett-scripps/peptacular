# try to import plotly and pandas
from typing import List, Tuple, Union
import peptacular as pt

try:
    import plotly.express as px
except ImportError:
    raise ImportError(
        "Plotly is not installed. Please install it with `pip install plotly`."
    )

try:
    import pandas as pd
except ImportError:
    raise ImportError(
        "Pandas is not installed. Please install it with `pip install pandas`."
    )

SPECTRUM_INPUT = Union[Tuple[List[float], List[float]], List[Tuple[float, float]]]

# FRAGMENT COLUMNS
MZ_COLUMN = "mz"
INTENSITY_COLUMN = "intensity"
CHARGE_COLUMN = "charge"
ISOTOPE_COLUMN = "isotope"
LOSS_COLUMN = "loss"
ERROR_COLUMN = "error"

# SCATTER COLUMNS
NAME_COLUMN = "scatter_name"
VISIBLE_COLUMN = "scatter_visible"
LEGENDGROUP_COLUMN = "scatter_legendgroup"
OPACITY_COLUMN = "scatter_opacity"
ZORDER_COLUMN = "scatter_zorder"
TEXT_COLUMN = "scatter_text"
TEXTPOSITION_COLUMN = "scatter_textposition"

# LINE COLUMNS
LINE_BACKOFF_COLUMN = "scatter_line_backoff"
LINE_COLOR_COLUMN = "scatter_line_color"
LINE_DASH_COLUMN = "scatter_line_dash"
LINE_WIDTH_COLUMN = "scatter_line_width"

# TEXTFONTS COLUMNS
TEXTFONT_COLOR_COLUMN = "scatter_textfont_color"
TEXTFONT_FAMILY_COLUMN = "scatter_textfont_family"
TEXTFONT_LINEPOSITION_COLUMN = "scatter_textfont_lineposition"
TEXTFONT_SHADOW_COLUMN = "scatter_textfont_shadow"
TEXTFONT_SIZE_COLUMN = "scatter_textfont_size"
TEXTFONT_STYLE_COLUMN = "scatter_textfont_style"
TEXTFONT_TEXTCASE_COLUMN = "scatter_textfont_textcase"
TEXTFONT_VARIANT_COLUMN = "scatter_textfont_variant"
TEXTFONT_WEIGHT_COLUMN = "scatter_textfont_weight"

### DEFAULT VALUES

# DEFAULT FRAGMENT VALUES
DEFAULT_MZ = 0.0
DEFAULT_INTENSITY = 0.0
DEFAULT_CHARGE = 0
DEFAULT_IOSTOPE = 0
DEFAULT_LOSS = 0.0
DEFAULT_ERROR = 0.0

# DEFAULT SCATTER VALUES
DEFAULT_NAME = ""
DEFAULT_VISIBLE = True
DEFAULT_LEGENDGROUP = "unassigned"
DEFAULT_OPACITY = 1.0
DEFAULT_ZORDER = 0
DEFAULT_TEXT = ""
DEFAULT_TEXTPOSITION = "top center"

# DEFAULT LINE VALUES
DEFAULT_LINE_BACKOFF = 0.0
DEFAULT_LINE_COLOR = 'gray'
DEFAULT_LINE_DASH = 'solid'
DEFAULT_LINE_WIDTH = 2.0

# DEFAULT VALUES
DEFAULT_TEXTFONT_COLOR_COLUMN = "gray"
DEFAULT_TEXTFONT_FAMILY_COLUMN = "Arial"
DEFAULT_TEXTFONT_LINEPOSITION_COLUMN = "none"
DEFAULT_TEXTFONT_SHADOW_COLUMN = "none"
DEFAULT_TEXTFONT_SIZE_COLUMN = 20
DEFAULT_TEXTFONT_STYLE_COLUMN = "normal"
DEFAULT_TEXTFONT_TEXTCASE_COLUMN = "normal"
DEFAULT_TEXTFONT_VARIANT_COLUMN = "normal"
DEFAULT_TEXTFONT_WEIGHT_COLUMN = 20

# FRAGMENT COLORS
DEFAULT_COLORS = {
    'a': 'red',
    'b': 'blue',
    'c': 'green',
    'x': 'purple',
    'y': 'orange',
    'z': 'pink',
}


def generate_sprectrum_table(
        spectrum: SPECTRUM_INPUT,
                  fragments: List[pt.Fragment],
                  tolerance_type: str,
                  tolerance_value: float) -> pd.DataFrame:
    
    mz_array, intensity_array = return_spectrum_values(spectrum)

    if len(mz_array) == 0 or len(intensity_array) == 0:
        raise ValueError("Spectrum must contain at least one peak.")
    
    if len(mz_array) != len(intensity_array):
        raise ValueError("Spectrum must contain the same number of mz and intensity values.")

    # Create a DataFrame from the spectrum data
    spectrum_table = pd.DataFrame({
        MZ_COLUMN: mz_array,
        INTENSITY_COLUMN: intensity_array,
    })

    # sum the intensities if there are multiple peaks with the same mz
    spectrum_table = spectrum_table.groupby(MZ_COLUMN).sum().reset_index()

    # set default values for the columns
    spectrum_table[CHARGE_COLUMN] = DEFAULT_CHARGE
    spectrum_table[ISOTOPE_COLUMN] = DEFAULT_IOSTOPE
    spectrum_table[LOSS_COLUMN] = DEFAULT_LOSS
    spectrum_table[ERROR_COLUMN] = DEFAULT_ERROR

    spectrum_table[NAME_COLUMN] = DEFAULT_NAME
    spectrum_table[VISIBLE_COLUMN] = DEFAULT_VISIBLE
    spectrum_table[LEGENDGROUP_COLUMN] = DEFAULT_LEGENDGROUP
    spectrum_table[OPACITY_COLUMN] = DEFAULT_OPACITY
    spectrum_table[ZORDER_COLUMN] = DEFAULT_ZORDER
    spectrum_table[TEXT_COLUMN] = DEFAULT_TEXT
    spectrum_table[TEXTPOSITION_COLUMN] = DEFAULT_TEXTPOSITION

    spectrum_table[LINE_BACKOFF_COLUMN] = DEFAULT_LINE_BACKOFF
    spectrum_table[LINE_COLOR_COLUMN] = DEFAULT_LINE_COLOR
    spectrum_table[LINE_DASH_COLUMN] = DEFAULT_LINE_DASH
    spectrum_table[LINE_WIDTH_COLUMN] = DEFAULT_LINE_WIDTH

    spectrum_table[TEXTFONT_COLOR_COLUMN] = DEFAULT_TEXTFONT_COLOR_COLUMN
    spectrum_table[TEXTFONT_FAMILY_COLUMN] = DEFAULT_TEXTFONT_FAMILY_COLUMN
    spectrum_table[TEXTFONT_LINEPOSITION_COLUMN] = DEFAULT_TEXTFONT_LINEPOSITION_COLUMN
    spectrum_table[TEXTFONT_SHADOW_COLUMN] = DEFAULT_TEXTFONT_SHADOW_COLUMN
    spectrum_table[TEXTFONT_SIZE_COLUMN] = DEFAULT_TEXTFONT_SIZE_COLUMN
    spectrum_table[TEXTFONT_STYLE_COLUMN] = DEFAULT_TEXTFONT_STYLE_COLUMN
    spectrum_table[TEXTFONT_TEXTCASE_COLUMN] = DEFAULT_TEXTFONT_TEXTCASE_COLUMN
    spectrum_table[TEXTFONT_VARIANT_COLUMN] = DEFAULT_TEXTFONT_VARIANT_COLUMN
    spectrum_table[TEXTFONT_WEIGHT_COLUMN] = DEFAULT_TEXTFONT_WEIGHT_COLUMN

    
    matches = pt.get_fragment_matches(fragments=fragments, 
                                      mz_spectra=spectrum_table[MZ_COLUMN].tolist(),
                                      intensity_spectra=spectrum_table[INTENSITY_COLUMN].tolist(),
                                      tolerance_value=tolerance_value,
                                      tolerance_type=tolerance_type,
                                      mode='closest')

    # update the spectrum table with the matches
    for match in matches:
        # get the mz value
        mz = match.mz

        index = spectrum_table[MZ_COLUMN] == mz

        # update the color and group columns
        spectrum_table.loc[index, LINE_COLOR_COLUMN] = DEFAULT_COLORS.get(match.ion_type, DEFAULT_LINE_COLOR)
        spectrum_table.loc[index, LEGENDGROUP_COLUMN] = match.ion_type
        spectrum_table.loc[index, TEXT_COLUMN] = match.label
        spectrum_table.loc[index, CHARGE_COLUMN] = match.charge
        spectrum_table.loc[index, ISOTOPE_COLUMN] = match.isotope
        spectrum_table.loc[index, LOSS_COLUMN] = match.loss
        spectrum_table.loc[index, ERROR_COLUMN] = match.error
        spectrum_table.loc[index, TEXTFONT_COLOR_COLUMN] = DEFAULT_COLORS.get(match.ion_type, DEFAULT_LINE_COLOR)

    return spectrum_table

def spectrum_to_line(mz_array: List[float], intensity_array: List[float]) -> Tuple[List[float], List[float]]:
    """
    Convert a spectrum to a line representation for plotting.
    Each peak is represented by two points: (mz, 0) and (mz, intensity).
    """
    x_values = []
    y_values = []

    for mz, intensity in zip(mz_array, intensity_array):
        x_values.extend([mz, mz, None])
        y_values.extend([0, intensity, None])

    return x_values, y_values


def get_group_value(group_df, column: str) -> Union[str, float, int]:
    """
    Get the value of a column for a group in the DataFrame.
    This is used to assert that all values in the group are the same.
    """
    assert len(group_df[column].unique()) == 1, f"All values in the group must be the same for {column}"
    return group_df[column].unique()[0]

def plot_spectrum_table(
    spectrum_table: pd.DataFrame,
    title: str = "Spectrum",
    xaxis_title: str = "M/Z",
    yaxis_title: str = "Intensity",
):
    """
    Plot a spectrum table using Plotly with vertical lines from y=0 to intensity values.
    Uses an efficient approach for large datasets by using scatter with connecting lines.
    """
    # Get the data from the DataFrame
    
    # Create a figure using go.Scatter for more control
    import plotly.graph_objects as go
    fig = go.Figure()

    # group by group and plot each group with a different color
    for group, group_df in spectrum_table.groupby(LEGENDGROUP_COLUMN):

        mz_values, intensity_values = spectrum_to_line(group_df[MZ_COLUMN].tolist(), 
                                                       group_df[INTENSITY_COLUMN].tolist())
        
        # assert that the colors are the same
        name = get_group_value(group_df, NAME_COLUMN)
        visible = get_group_value(group_df, VISIBLE_COLUMN)
        legendgroup = get_group_value(group_df, LEGENDGROUP_COLUMN)
        opacity = get_group_value(group_df, OPACITY_COLUMN)
        zorder = get_group_value(group_df, ZORDER_COLUMN)

        line_color = get_group_value(group_df, LINE_COLOR_COLUMN)
        line_dash = get_group_value(group_df, LINE_DASH_COLUMN)
        line_width = get_group_value(group_df, LINE_WIDTH_COLUMN)


        fig.add_trace(
            go.Scattergl(
                x=mz_values,
                y=intensity_values,
                name=name,
                visible=visible,
                legendgroup=legendgroup,
                mode='lines',
                line=dict(
                          color=line_color, 
                          dash=line_dash,
                          width=line_width),
                hoverinfo='none',  # Disable hover info on connecting lines
                opacity=opacity,
            )
        )


        # Add annotations for each peak label
        for _, row in group_df.iterrows():
            if row[TEXT_COLUMN]:  # Only add annotation if there's text
                fig.add_annotation(
                    x=row[MZ_COLUMN],
                    y=row[INTENSITY_COLUMN],
                    text=row[TEXT_COLUMN],
                    showarrow=False,
                    textangle=-90,
                    yshift=10,  # Shift the text up a bit from the peak
                    font=dict(
                    color=row[TEXTFONT_COLOR_COLUMN],
                    family=row[TEXTFONT_FAMILY_COLUMN],
                    size=row[TEXTFONT_SIZE_COLUMN],
                    ),
                    xanchor='center',
                    yanchor='bottom'
                )

        # Add hover data as separate scatter trace with no visible markers
        fig.add_trace(
            go.Scatter(
            x=group_df[MZ_COLUMN],
            y=group_df[INTENSITY_COLUMN],
            mode='text',
            text="",
            name=name,
            hovertemplate=
            'M/Z: %{x:.4f}<br>' +
            'Intensity: %{y:.2f}<br>' +
            'Ion: ' + legendgroup + '<br>',
            legendgroup=legendgroup,
            showlegend=False
            )
        )

    # Update layout
    fig.update_layout(
        title=title,
        xaxis_title=xaxis_title,
        yaxis_title=yaxis_title,
        showlegend=False,
        hovermode='closest'
    )
    
    return fig


def return_spectrum_values(spectrum: SPECTRUM_INPUT) -> Tuple[List[float], List[float]]:
    mz_array = []
    intensity_array = []

    if isinstance(spectrum, tuple):
        # If the spectrum is a tuple, unpack it
        mz_array, intensity_array = spectrum

    elif isinstance(spectrum, list):
        for mz, intensity in spectrum:
            mz_array.append(mz)
            intensity_array.append(intensity)

    else:
        raise ValueError("Spectrum must be a tuple of (mz_array, intensity_array) or a list of (mz, intensity) tuples.")

    return mz_array, intensity_array


def plot_spectrum(spectrum: SPECTRUM_INPUT,
                  sequence: str, 
                  ion_types: str, 
                  charge_range: Tuple[int, int],
                  error_type: str,
                  error_tolerance: float,
                  monoisotopic: bool):
    
    spectrum_df = generate_sprectrum_table(spectrum, sequence, ion_types, charge_range, error_type, error_tolerance, monoisotopic)
    fig = plot_spectrum_table(spectrum_df, title=f"{sequence}")

    return fig
    
    

if __name__  == "__main__":
    from peptacular.plotting.tmp import spectrum

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
    fig.show()
