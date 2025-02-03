import streamlit as st
import requests
import re
import py3Dmol
from stmol import showmol

AA_MAP = {
   'a': 'A', 'b': 'C', 'c': 'C', 'd': 'D', 'e': 'E', 'f': 'F', 
   'g': 'G', 'h': 'H', 'i': 'I', 'j': 'L', 'k': 'K', 'l': 'L', 
   'm': 'M', 'n': 'N', 'o': 'O', 'p': 'P', 'q': 'Q', 'r': 'R', 
   's': 'S', 't': 'T', 'u': 'U', 'v': 'V', 'w': 'W', 'x': 'X', 
   'y': 'Y', 'z': 'F'
}

def name_to_amino_acid_sequence(name):
   name = name.replace(' ', '').lower()
   sequence = ''.join(AA_MAP.get(char, 'G') for char in name)
   repeat_count = max(10, len(name))
   return (sequence * repeat_count)[:1000]

def validate_sequence(sequence):
   if not re.match(r'^[ACDEFGHIKLMNPQRSTVWY]+$', sequence):
       st.error('Invalid protein sequence')
       return None
   return sequence

def predict_structure(sequence):
   try:
       headers = {'Content-Type': 'application/x-www-form-urlencoded'}
       response = requests.post(
           'https://api.esmatlas.com/foldSequence/v1/pdb/', 
           headers=headers, 
           data=sequence,
           timeout=30
       )
       response.raise_for_status()
       return response.content.decode('utf-8')
   except requests.RequestException as e:
       st.error(f"API Request Failed: {e}")
       return None

def render_mol(pdb, style='cartoon'):
   pdbview = py3Dmol.view()
   pdbview.addModel(pdb, 'pdb')
   
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

def main():
   st.set_page_config(page_title="Name to Protein Structure", page_icon="🧬", layout="wide")
   st.sidebar.title('🧬 Name to Protein Structure')
   
   name = st.sidebar.text_input('Enter Your Name')
   visualization_style = st.sidebar.selectbox(
       'Visualization Style', 
       ['Cartoon (Spectrum)', 'Stick', 'Sphere', 'Cross', 'Ribbon']
   )

   if st.sidebar.button('Predict Structure'):
       if name:
           amino_acid_sequence = name_to_amino_acid_sequence(name)
           validated_seq = validate_sequence(amino_acid_sequence)
           
           if validated_seq:
               pdb_string = predict_structure(validated_seq)
               
               if pdb_string:
                   st.success(f"Structure prediction for {name} completed!")
                   st.subheader('Predicted Protein Structure')
                   render_mol(pdb_string, visualization_style)
                   
                   st.download_button(
                       label="Download PDB File",
                       data=pdb_string,
                       file_name=f'{name.replace(" ", "_")}_structure.pdb',
                       mime='text/plain'
                   )
       else:
           st.warning('Please enter a name')

if __name__ == "__main__":
   main()
