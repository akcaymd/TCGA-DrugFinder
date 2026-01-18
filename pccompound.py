from Bio import Entrez, SeqIO


def pccompound(query):
    # her hit için istek atmak mı, yoksa ilaç listesini çekip depolayıp onu kullanmak mı?
    Entrez.email = "furkanakcay@proton.me"

    # relevance – default sort order, (“Best Match”) on web PubMed. https://www.ncbi.nlm.nih.gov/books/NBK25499/#chapter4.ESearch:~:text=sort
    handle = Entrez.esearch(db="pccompound", term=f"{query}[Synonym]", retmax=10)
    result = Entrez.read(handle)
    handle.close()

    id_list = ",".join(result.get("IdList", []))

    if id_list == "":
        return []

    meshlist = []
    drugs_handle = Entrez.esummary(db="pccompound", id=id_list)
    records = Entrez.read(drugs_handle)
    drugs_handle.close()

    for i in records:
        if i["MeSHHeadingList"]:
            meshlist.append(
                {
                    "CID": i["CID"],
                    "MeSHHeadingList": i["MeSHHeadingList"],
                    "MeSHTermList": i["MeSHTermList"],
                }
            )

    return meshlist
