import requests
import json
from Bio import Entrez, SeqIO
from pccompound import pccompound
import time
from Database import Database


class GdcCases:

    GDC_URL = "https://api.gdc.cancer.gov/cases"

    def __init__(self):
        self.db = Database("./db.sqlite")
        self.raw_data = None
        self.hash = None
        self.query = ""

        self.db.execute(
            """CREATE TABLE IF NOT EXISTS drugs ( id INTEGER PRIMARY KEY AUTOINCREMENT, title TEXT NOT NULL, CID TEXT, mesh_list TEXT )"""
        )

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

        patients = []
        _hash = []

        # Loop in cases
        for hit in self.raw_data["data"]["hits"]:
            patient_id = hit.get("submitter_id")
            has_valid_drug = False

            for diagnosis in hit.get("diagnoses", []):
                for treatment in diagnosis.get("treatments", []):
                    agent = treatment.get("therapeutic_agents")
                    demographic = hit.get("demographic")

                    if (
                        agent
                        and agent.lower() not in ["unknown", "none", "clinical trial"]
                        and demographic
                    ):
                        has_valid_drug = True
                        _hash.append(
                            {
                                "submitter_id": hit["submitter_id"],
                                "treatment_type": treatment.get("treatment_type"),
                                "therapeutic_agents": agent,
                                "is_smoker": self.is_smoker(hit),
                                "days_to_death": demographic.get("days_to_death"),
                            }
                        )
            if has_valid_drug:
                patients.append(hit)
        self.hash = _hash
        self.patients = patients

    def is_smoker(self, hit):
        status = hit.get("exposures")[0].get("tobacco_smoking_status")
        return (
            "smoker"
            if status == "Current Smoker"
            else ("nonsmoker" if status == "Lifelong Non-Smoker" else "reformed smoker")
        )

    def get_mesh_list(self):
        for i in range(len(self.hash)):
            hit = self.hash[i]

            agent = hit.get("therapeutic_agents")
            mesh_list = self.db.fetchone("SELECT * FROM drugs WHERE title=?", (agent,))

            if mesh_list == None:
                time.sleep(1)
                mesh_list = pccompound(agent)
                self.db.execute(
                    "INSERT INTO drugs (title, CID, mesh_list) VALUES (?, ?, ?)",
                    (
                        agent,
                        mesh_list.get("CID"),
                        json.dumps(mesh_list.get("mesh_list")),
                    ),
                )
            self.hash[i]["CID"] = mesh_list.get("CID")
            self.hash[i]["mesh_list"] = mesh_list.get("mesh_list")
