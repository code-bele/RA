import urllib.request
import xml.etree.ElementTree as ET
import pandas as pd

# Sample DataFrame with a column of Gene IDs
df = pd.read_csv("output_genbank_id.csv")
df = df.dropna(subset=["gene_id"])
df['gene_id'] = df['gene_id'].astype(int)



def fetch_gene_symbol(gene_id):
    """Fetch the gene symbol for a given gene ID from NCBI Entrez."""
    query = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=gene&id={gene_id}&retmode=xml"

    try:
        with urllib.request.urlopen(query) as response:
            content = response.read().decode()
            root = ET.fromstring(content)
            gene_symbol = root.find(".//Gene-ref_locus")
            return gene_symbol.text if gene_symbol is not None else None
    except Exception as e:
        print(f"Error fetching gene {gene_id}: {e}")
        return None


# Apply function to DataFrame
df['gene_symbol'] = df['gene_id'].apply(fetch_gene_symbol)
df.to_csv("genbank_with_symbol.csv", index=False)

print(df)
