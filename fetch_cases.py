import requests
import json

cases_endpt = "https://api.gdc.cancer.gov/cases"

fields = [
    "submitter_id",
    "case_id",
    "primary_site",
    "diagnoses.treatments.therapeutic_agents",
    "diagnoses.treatments.treatment_type",
    "diagnoses.days_to_last_follow_up",
    "exposures.tobacco_smoking_status",
    "demographic.days_to_death",
]

fields = ",".join(fields)

filters = {
    "op": "and",
    "content": [
        {
            "op": "in",
            "content": {
                "field": "project.project_id",
                "value": ["TCGA-LUAD"],
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
    ],
}

params = {
    "filters": json.dumps(filters),
    "fields": fields,
    "format": "JSON",
    "size": "100",
}

response = requests.get(cases_endpt, params=params)
data = response.json()

# clean data and create a list
# patients containing drug information
patients = []

for i in data["data"]["hits"]:
    for treatments in i["diagnoses"]:
        if not treatments.get("treatments"):
            continue

        for k in treatments["treatments"]:
            if not k.get("therapeutic_agents"):
                continue

            if (
                k.get("therapeutic_agents") == "Unknown"
                or k.get("therapeutic_agents") == None
            ):
                continue

    patients.append(i)


with open("file.json", "w", encoding="utf-8") as f:
    json.dump(data, f)

print(len(patients))
for i in patients:
    print(i["submitter_id"])
    print("*" * 20)
