from Bio import Entrez, SeqIO

Entrez.email = "furkanakcay@proton.me"

# relevance – default sort order, (“Best Match”) on web PubMed. https://www.ncbi.nlm.nih.gov/books/NBK25499/#chapter4.ESearch:~:text=sort
handle = Entrez.esearch(db="pccompound", term="Cisplatin[Synonym]", retmax=10)
result = Entrez.read(handle)
handle.close()

id_list = ",".join(result.get("IdList", []))

drugs_handle = Entrez.esummary(db="pccompound", id=id_list)
records = Entrez.read(drugs_handle)
drugs_handle.close()

for i in records:
    if i["MeSHTermList"]:
        print(i["CID"])
        print(i["MeSHTermList"])
        print("*" * 20)
