import xml.dom.minidom
import urllib.request
import pandas as pd


def _getText(node):
    if node.firstChild and node.firstChild.nodeType == node.TEXT_NODE:
        return node.firstChild.data
    return ""


def fetch_gene_id(genbank_id):
    query = f"http://www.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gene&term={genbank_id}[Nucleotide+Accession]"
    with urllib.request.urlopen(query) as response:
        doc = xml.dom.minidom.parse(response)
    idlist = list(map(_getText, doc.getElementsByTagName('Id')))
    return idlist[0] if idlist else None


# Load your DataFrame (replace with actual file or data source)
df = pd.read_csv("male_early_ra_genbank.csv")  # Ensure this file has a 'genbank_id' column

# Apply the function to fetch gene IDs
df["gene_id"] = df["gene"].apply(fetch_gene_id)

# Save or display the result
df.to_csv("output_genbank_id.csv", index=False)
print(df.head())
