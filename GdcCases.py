import requests
import json
from Bio import Entrez, SeqIO
from pccompound import pccompound
import time


class GdcCases:

    GDC_URL = "https://api.gdc.cancer.gov/cases"

    def __init__(self):
        self.raw_data = None
        self.cases = None
        self.query = ""

    def fetch(self, query):
        """
        Fetches clinical data for a specific project (e.g., 'LUAD').
        """
        self.query = query
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
                        "value": [self.query],
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

    def get_data(self):
        if not self.raw_data or "data" not in self.raw_data:
            return []

        cases = []

        # Loop in cases
        for hit in self.raw_data["data"]["hits"]:
            patient_id = hit.get("submitter_id")

            for diagnosis in hit.get("diagnoses", []):
                for treatment in diagnosis.get("treatments", []):
                    agent = treatment.get("therapeutic_agents")
                    demographic = hit.get("demographic")

                    if (
                        agent
                        and agent.lower() not in ["unknown", "none", "clinical trial"]
                        and demographic
                    ):
                        cases.append(
                            {
                                "submitter_id": hit["submitter_id"],
                                "treatment_type": treatment.get("treatment_type"),
                                "therapeutic_agents": agent,
                                "is_smoker": self.is_smoker(hit),
                                "days_to_death": demographic.get("days_to_death"),
                            }
                        )

        self.cases = cases

    def is_smoker(self, hit):
        status = hit.get("exposures")[0].get("tobacco_smoking_status")
        return (
            "smoker"
            if status == "Current Smoker"
            else ("nonsmoker" if status == "Lifelong Non-Smoker" else "reformed smoker")
        )

    def get_mesh_list(self):
        hash_list = []
        # for i in range(len(self.cases)):
        for i in range(10):
            hit = self.cases[i]

            agent = hit.get("therapeutic_agents")

            time.sleep(3)
            mesh_list = pccompound(agent)
            for mesh in mesh_list:
                hash_list.append(
                    {
                        **hit,
                        **{
                            "CID": mesh.get("CID"),
                            "mesh_list": mesh.get("mesh_list"),
                        },
                    }
                )
        return hash_list
