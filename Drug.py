import streamlit as st
import pandas as pd
import requests
from rdkit import Chem
from rdkit.Chem import AllChem

def fetch_molecule_json(url):
    """Fetch molecule data from ChEMBL JSON URL."""
    try:
        response = requests.get(url)
        response.raise_for_status()
        return response.json()
    except requests.exceptions.RequestException as e:
        st.error(f"Failed to fetch molecule data from {url}: {e}")
        return None

def generate_fingerprint(smiles):
    """Generate Morgan fingerprint from SMILES."""
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            return AllChem.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=2048)
        else:
            st.error(f"Invalid SMILES: {smiles}")
            return None
    except Exception as e:
        st.error(f"Failed to generate fingerprint: {e}")
        return None

def main():
    st.title("Drug Discovery: Fingerprint Generation")

    # File uploader
    uploaded_file = st.file_uploader("Upload your dataset (CSV format)", type="csv")

    if uploaded_file is not None:
        # Load the dataset
        try:
            data = pd.read_csv(uploaded_file, header=None, names=["ChEMBL_ID", "URL", "Target"])
            st.write("### Uploaded Data:", data)
        except Exception as e:
            st.error(f"Error reading the file: {e}")
            return

        fingerprints = []

        for index, row in data.iterrows():
            chembl_id = row["ChEMBL_ID"]
            url = row["URL"]

            # Fetch molecule data
            molecule_data = fetch_molecule_json(url)
            if molecule_data and "molecule_structures" in molecule_data:
                smiles = molecule_data["molecule_structures"].get("canonical_smiles", None)
                if smiles:
                    fingerprint = generate_fingerprint(smiles)
                    fingerprints.append((chembl_id, smiles, fingerprint))
                else:
                    st.warning(f"No SMILES available for {chembl_id}")
            else:
                st.warning(f"No molecule data for {chembl_id}")

        # Display the results
        if fingerprints:
            st.write("### Generated Fingerprints:")
            fingerprint_df = pd.DataFrame(
                [(chembl_id, smiles, list(fp)) for chembl_id, smiles, fp in fingerprints if fp is not None],
                columns=["ChEMBL_ID", "SMILES", "Fingerprint"]
            )
            st.write(fingerprint_df)

            # Option to download results
            csv = fingerprint_df.to_csv(index=False)
            st.download_button(
                label="Download Fingerprints as CSV",
                data=csv,
                file_name="fingerprints.csv",
                mime="text/csv",
            )

if __name__ == "__main__":
    main()
