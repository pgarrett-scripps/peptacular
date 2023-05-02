from peptacular.constants import ENZYMES_OPTIONS


def search_components():
    import streamlit as st

    tolerance = st.number_input(label='Tolerance')
    tolerance_type = st.radio(label='Tolerance Type', options=('PPM', 'TH', 'Dalton'))
    return {'tolerance': tolerance, 'tolerance_type': tolerance_type}


def digestion_components():
    import streamlit as st

    non_enzymatic = st.checkbox(label='Non Enzymatic', value=False)

    col1, col2 = st.columns(2)
    min_peptide_length = col1.number_input(label='Min Peptide Length', value=7)
    max_peptide_length = col2.number_input(label='Max Peptide Length', value=30)

    enzyme_regex_location_pairs, missed_cleavages, semi_enzymatic = None, None, None
    if non_enzymatic is False:

        missed_cleavages = st.number_input(label='Number of Missed Cleavages', value=0)
        semi_enzymatic = st.checkbox(label='Semi Enzymatic', value=False)

        # Enzyme Parser
        enzyme_input = st.radio(label='Enzyme Input', options=('Regex', 'Enzyme'))
        enzyme_regex_location_pairs = None
        if enzyme_input == 'Regex':
            col1, col2 = st.columns(2)
            enzyme_regex = col1.text_input(label='Enzyme Regex', value='([KR])([^P])')
            enzyme_regex_cleavage = col2.number_input(label='Cleavage Site', value=1)
            enzyme_regex_location_pairs = ([(enzyme_regex, enzyme_regex_cleavage)], [])
        elif enzyme_input == 'Enzyme':
            enzyme_regex_location_pairs = ENZYMES_OPTIONS[st.selectbox(label='Enzyme', options=ENZYMES_OPTIONS.keys())]

    return {'non_enzymatic': non_enzymatic,
            'min_peptide_length': min_peptide_length,
            'max_peptide_length': max_peptide_length,
            'enzyme_regexes': enzyme_regex_location_pairs,
            'missed_cleavages': missed_cleavages,
            'semi_enzymatic': semi_enzymatic}
