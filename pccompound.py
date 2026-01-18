from Bio import Entrez

Entrez.email = "furkanakcay@proton.me"


def pccompound(query):
    # TODO her hit için istek atmak mı, yoksa ilaç listesini çekip depolayıp onu kullanmak mı?

    # relevance – default sort order, (“Best Match”) on web PubMed. https://www.ncbi.nlm.nih.gov/books/NBK25499/#chapter4.ESearch:~:text=sort
    handle = Entrez.esearch(db="pccompound", term=f"{query}[Synonym]", retmax=10)
    result = Entrez.read(handle)
    handle.close()

    meshlist = dict()

    for i in result.get("IdList", []):
        meshlist["CID"] = i
        meshlist["mesh_list"] = get_mesh_term_directly(query)

    return meshlist


def get_mesh_term_directly(chemical_name):
    # Step 1: Search the MeSH database for the chemical name
    with Entrez.esearch(db="mesh", term=chemical_name) as search_handle:
        search_results = Entrez.read(search_handle)

    ids = search_results.get("IdList", [])

    if not ids:
        return f"No MeSH terms found for {chemical_name}."

    # Step 2: Retrieve the specific MeSH descriptor names
    with Entrez.esummary(db="mesh", id=",".join(ids)) as summary_handle:
        summaries = Entrez.read(summary_handle)

    return [s["DS_MeshTerms"] for s in summaries]
