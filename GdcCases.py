import requests
import json
from Bio import Entrez, SeqIO
from pccompound import pccompound
import time


class GdcCases:

    GDC_URL = "https://api.gdc.cancer.gov/cases"

    def __init__(self):
        self.raw_data = None
        self.hash = None

    def fetch(self, query):
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

        # Loop in cases
        for hit in self.raw_data["data"]["hits"]:
            patient_id = hit.get("submitter_id")

            has_valid_drug = False

            # Loop in case's treatments
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
                        tobacco_smoking_status = hit.get("exposures")[0].get(
                            "tobacco_smoking_status"
                        )
                        _hash.append(
                            {
                                "submitter_id": hit["submitter_id"],
                                "treatment_type": treatment.get("treatment_type"),
                                "therapeutic_agents": agent,
                                "tobacco_smoking_status": tobacco_smoking_status,
                                "is_smoker": (
                                    "smoker"
                                    if tobacco_smoking_status == "Current Smoker"
                                    else (
                                        "nonsmoker"
                                        if tobacco_smoking_status
                                        == "Lifelong Non-Smoker"
                                        else "reformed smoker"
                                    )
                                ),
                                "days_to_death": hit.get("demographic").get(
                                    "days_to_death"
                                ),
                            }
                        )

            if has_valid_drug:
                patients.append(hit)
        self.hash = _hash
        self.patients = patients

    def pretty(self):
        print(json.dumps(self.hash, indent=4, sort_keys=True))

    def get_mesh_list(self):
        for i in self.hash:
            meshlist = pccompound(i["therapeutic_agents"])
            print(meshlist)
