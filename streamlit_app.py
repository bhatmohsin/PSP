import streamlit as st
import requests
import re
import py3Dmol
from stmol import showmol
import biotite.structure.io as bsio

# Advanced protein analysis (optional)
try:
    from Bio.SeqUtils.ProtParam import ProteinAnalysis
    BIOPYTHON_AVAILABLE = True
except ImportError:
    BIOPYTHON_AVAILABLE = False

# Page configuration
st.set_page_config(
    page_title="SPCFold Protein Structure Predictor", 
    page_icon="🧬", 
    layout="wide"
)

# Sidebar introduction
st.sidebar.title('🧬 SPCFold Protein Structure Predictor')
st.sidebar.markdown("""
    SPCFold is an advanced protein structure prediction tool 
    based on the ESM-2 language model. 
    [Learn More](https://esmatlas.com/about)
""")

# Default sequence
DEFAULT_SEQ = "MGSSHHHHHHSSGLVPRGSHMRGPNPTAASLEASAGPFTVRSFTVSRPSGYGAGTVYYPTNAGGTVGAIAIVPGYTARQSSIKWWGPRLASHGFVVITIDTNSTLDQPSSRSSQQMAALRQVASLNGTSSSPIYGKVDTARMGVMGWSMGGGGSLISAANNPSLKAAAPQAPWDSSTNFSSVTVPTLIFACENDSIAPVNSSALPIYDSMSRNAKQFLEINGGSHSCANSGNSNQALIGKKGVAWMKRFMDNDTRYSTFACENPNSTRVSDFRTANCSLEDPAANKARKEAELAAATAEQ"

# Utility Functions
def validate_sequence(sequence):
    """
    Validate protein sequence input
    """
    # Remove whitespace and convert to uppercase
    cleaned_sequence = sequence.replace(' ', '').upper()
    
    # Validate sequence
    if not re.match(r'^[ACDEFGHIKLMNPQRSTVWY]+$', cleaned_sequence):
        st.error('Invalid protein sequence. Use standard amino acid letters.')
        return None
    
    # Check sequence length
    if len(cleaned_sequence) < 10:
        st.error('Sequence too short. Minimum 10 amino acids required.')
        return None
    
    if len(cleaned_sequence) > 1000:
        st.error('Sequence too long. Maximum 1000 amino acids supported.')
        return None
    
    return cleaned_sequence

@st.cache_data(ttl=3600)  # Cache for 1 hour
def predict_structure(sequence):
    """
    Predict protein structure using ESMFold API with error handling
    """
    try:
        headers = {
            'Content-Type': 'application/x-www-form-urlencoded',
        }
        response = requests.post(
            'https://api.esmatlas.com/foldSequence/v1/pdb/', 
            headers=headers, 
            data=sequence,
            timeout=30  # Add timeout
        )
        response.raise_for_status()  # Raise exception for bad responses
        return response.content.decode('utf-8')
    
    except requests.RequestException as e:
        st.error(f"API Request Failed: {e}")
        return None

def render_mol(pdb, style='cartoon'):
    """
    Render 3D molecular visualization with multiple style options
    """
    pdbview = py3Dmol.view()
    pdbview.addModel(pdb, 'pdb')
    
    # Configurable visualization styles
    style_options = {
        'Cartoon (Spectrum)': {'cartoon': {'color': 'spectrum'}},
        'Stick': {'stick': {}},
        'Sphere': {'sphere': {'scale': 0.5}},
        'Cross': {'cross': {'scale': 1}},
        'Ribbon': {'ribbon': {}}
    }
    
    pdbview.setStyle(style_options.get(style, {'cartoon': {'color': 'spectrum'}}))
    pdbview.setBackgroundColor('white')
    pdbview.zoomTo()
    pdbview.zoom(2, 800)
    pdbview.spin(True)
    showmol(pdbview, height=500, width=800)

def analyze_protein_sequence(sequence):
    """
    Provide basic protein sequence analysis
    """
    if BIOPYTHON_AVAILABLE:
        try:
            protein_analysis = ProteinAnalysis(sequence)
            
            col1, col2 = st.columns(2)
            with col1:
                st.metric('Sequence Length', len(sequence))
                st.metric('Molecular Weight', 
                          f"{protein_analysis.molecular_weight():.2f} Da")
                st.metric('Isoelectric Point', 
                          f"{protein_analysis.isoelectric_point():.2f}")
            
            with col2:
                # Amino acid composition
                composition = protein_analysis.count_amino_acids()
                most_common = max(composition, key=composition.get)
                st.metric('Most Common Amino Acid', 
                          f"{most_common} (x{composition[most_common]})")
                
                # Hydrophobicity
                st.metric('Hydrophobicity', 
                          f"{protein_analysis.gravy():.2f}")
        
        except Exception as e:
            st.warning(f"Protein analysis error: {e}")
    else:
        st.info("Install Biopython for advanced protein analysis")

def main():
    # Sequence input
    txt = st.sidebar.text_area(
        'Input Protein Sequence', 
        DEFAULT_SEQ, 
        height=275, 
        help="Enter a valid protein sequence using standard amino acid letters"
    )

    # Visualization style selector
    visualization_style = st.sidebar.selectbox(
        'Visualization Style', 
        ['Cartoon (Spectrum)', 'Stick', 'Sphere', 'Cross', 'Ribbon'],
        help="Choose how the protein structure will be displayed"
    )

    # Prediction button
    if st.sidebar.button('Predict Structure'):
        # Validate sequence
        validated_seq = validate_sequence(txt)
        
        if validated_seq:
            # Predict structure
            pdb_string = predict_structure(validated_seq)
            
            if pdb_string:
                # Structure prediction successful
                st.success("Structure prediction completed!")
                
                # Render molecular structure
                st.subheader('Predicted Protein Structure')
                render_mol(pdb_string, visualization_style)
                
                # Load structure for analysis
                struct = bsio.load_structure('predicted.pdb', extra_fields=["b_factor"])
                b_value = round(struct.b_factor.mean(), 4)
                
                # Confidence and download
                st.subheader('Prediction Confidence')
                st.info(f'plDDT Score: {b_value} (0-100 scale)')
                st.caption('plDDT: Per-residue confidence of structure prediction')
                
                # Download PDB
                st.download_button(
                    label="Download PDB File",
                    data=pdb_string,
                    file_name='predicted_structure.pdb',
                    mime='text/plain',
                )
                
                # Sequence analysis
                st.subheader('Sequence Analysis')
                analyze_protein_sequence(validated_seq)
            
    else:
        st.warning('👈 Enter a protein sequence and click Predict!')