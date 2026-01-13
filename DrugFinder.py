import requests
import json
from Bio import Entrez, SeqIO


class DrugFinder:

    GDC_URL = "https://api.gdc.cancer.gov/cases"

    def __init__(self):
        self.raw_data = None
        self.hash = None

    def case_fetch(self, query):
        """
        Fetches clinical data for a specific project (e.g., 'LUAD').
        """
        fields = [
            "submitter_id",
            "case_id",
            "primary_site",
            "diagnoses.treatments.therapeutic_agents",
            "diagnoses.treatments.treatment_type",
            "exposures.tobacco_smoking_status",
            "demographic.days_to_death",
        ]

        filters = {
            "op": "and",
            "content": [
                {
                    "op": "in",
                    "content": {
                        "field": "project.project_id",
                        "value": [query],
                    },
                },
                {
                    "op": "in",
                    "content": {
                        "field": "diagnoses.treatments.treatment_type",
                        "value": [
                            "Chemotherapy",
                            "Immunotherapy (Including Vaccines)",
                            "Pharmaceutical Therapy, NOS",
                        ],
                    },
                },
                {
                    "op": "not",
                    "content": {
                        "field": "diagnoses.treatments.therapeutic_agents",
                    },
                },
                {
                    "op": "not",
                    "content": {
                        "field": "demographic.days_to_death",
                    },
                },
                {
                    "op": "not",
                    "content": {
                        "field": "exposures.tobacco_smoking_status",
                    },
                },
                {
                    "op": "!=",
                    "content": {
                        "field": "exposures.tobacco_smoking_status",
                        "value": "Not Reported",
                    },
                },
                {
                    "op": "!=",
                    "content": {
                        "field": "diagnoses.treatments.therapeutic_agents",
                        "value": "None",
                    },
                },
            ],
        }
        params = {
            "filters": json.dumps(filters),
            "fields": ",".join(fields),
            "format": "JSON",
            "size": "100",
        }

        try:
            response = requests.get(self.GDC_URL, params=params)
            response.raise_for_status()
            self.raw_data = response.json()
        except requests.exceptions.RequestException as e:
            print(f"Error fetching data: {e}")
            self.raw_data = None

    def get_cleaned_data(self):
        if not self.raw_data or "data" not in self.raw_data:
            return []

        patients = []
        _hash = []

        for hit in self.raw_data["data"]["hits"]:
            patient_id = hit.get("submitter_id")

            has_valid_drug = False
            for diagnosis in hit.get("diagnoses", []):
                for treatment in diagnosis.get("treatments", []):
                    agent = treatment.get("therapeutic_agents")

                    if (
                        agent
                        and agent.lower()
                        not in [
                            "unknown",
                            "none",
                        ]
                        and hit.get("demographic")
                    ):
                        has_valid_drug = True
                        smoker = hit.get("exposures")[0].get("tobacco_smoking_status")
                        _hash.append(
                            {
                                "submitter_id": hit["submitter_id"],
                                "treatment_type": treatment.get("treatment_type"),
                                "therapeutic_agents": treatment.get(
                                    "therapeutic_agents"
                                ),
                                "tobacco_smoking_status": smoker,
                                "is_smoker": (
                                    1
                                    if smoker == "Current Smoker"
                                    else 0 if smoker == "Lifelong Non-Smoker" else 2
                                ),
                                "days_to_death": hit.get("demographic").get(
                                    "days_to_death"
                                ),
                                "cid": "",
                                "mesh_term": "",
                            }
                        )

            if has_valid_drug:
                patients.append(hit)
        self.hash = _hash
        return patients

    def pccompound(self, query):
        # her hit için istek atmak mı, yoksa ilaç listesini çekip depolayıp onu kullanmak mı?
        Entrez.email = "furkanakcay@proton.me"

        # relevance – default sort order, (“Best Match”) on web PubMed. https://www.ncbi.nlm.nih.gov/books/NBK25499/#chapter4.ESearch:~:text=sort
        handle = Entrez.esearch(db="pccompound", term=f"{query}[Synonym]", retmax=10)
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


data = DrugFinder()
data.case_fetch("TCGA-LUAD")
patients = data.get_cleaned_data()

for i in data.hash:
    print(i)
    print("*" * 20)
print(len(data.hash))
