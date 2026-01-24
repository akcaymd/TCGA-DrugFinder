from Bio import Entrez

Entrez.email = "furkanakcay@proton.me"


def pccompound(drug_title):
    # relevance – default sort order, (“Best Match”) on web PubMed. https://www.ncbi.nlm.nih.gov/books/NBK25499/#chapter4.ESearch:~:text=sort
    handle = Entrez.esearch(db="pccompound", term=f"{drug_title}[Synonym]", retmax=10)
    result = Entrez.read(handle)
    handle.close()
    id_list = result.get("IdList", [])
    meshlist = []
    for i in get_mesh_ids(id_list):
        handle = Entrez.esummary(db="mesh", id=i["mesh_id"])

        mesh = Entrez.read(handle)
        handle.close()
        meshlist.append(
            {
                "CID": i["CID"],
                "mesh_list": mesh[0]["DS_MeshTerms"],
            }
        )

    return meshlist


"""
def get_mesh_term(chemical_name):
    # Step 1: Search the MeSH database for the chemical name
    search_handle = Entrez.esearch(db="mesh", term=f"{chemical_name}[MeSH Terms]")
    search_results = Entrez.read(search_handle)
    ids = search_results.get("IdList", [])

    if not ids:
        return None

    # Step 2: Retrieve the specific MeSH descriptor names
    summary_handle = Entrez.esummary(db="mesh", id=",".join(ids))
    summaries = Entrez.read(summary_handle)

    return [s["DS_MeshTerms"] for s in summaries][0]
"""


def get_mesh_ids(ids):
    handle = Entrez.elink(dbfrom="pccompound", db="mesh", id=ids)
    records = Entrez.read(handle)
    handle.close()

    mesh_ids = []

    for linkset in records:
        cid = linkset.get("IdList")[0]
        for dbset in linkset.get("LinkSetDb", []):
            for link in dbset.get("Link", []):
                mesh_ids.append({"CID": cid, "mesh_id": link["Id"]})
    return mesh_ids
