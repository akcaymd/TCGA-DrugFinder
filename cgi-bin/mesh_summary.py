from collections import defaultdict
import numpy as np


def build_mesh_summary(mesh_hash_list):
    summary = {}

    for row in mesh_hash_list:

        patient = row["submitter_id"]
        mesh_terms = row["mesh_list"]
        survival = row["days_to_death"]
        smoker = row["is_smoker"]
        treatment = row["treatment_type"]

        for mesh in mesh_terms:

            if mesh not in summary:
                summary[mesh] = {
                    "patients": set(),
                    "survival_days": [],
                    "smokers": 0,
                    "nonsmokers": 0,
                    "therapy": defaultdict(int),
                }

            summary[mesh]["patients"].add(patient)
            summary[mesh]["survival_days"].append(survival)

            # smoking status
            if smoker == "smoker":
                summary[mesh]["smokers"] += 1
            else:
                summary[mesh]["nonsmokers"] += 1

            # therapy type
            if treatment:
                summary[mesh]["therapy"][treatment] += 1

    return summary


def mesh_summary_table(summary):
    table = []

    for mesh, data in summary.items():

        mean_survival = int(sum(data["survival_days"]) / len(data["survival_days"]))

        table.append(
            {
                "mesh_term": mesh,
                "patients": len(data["patients"]),
                "smokers": data["smokers"],
                "nonsmokers": data["nonsmokers"],
                "therapy": dict(data["therapy"]),
                "mean_survival": mean_survival,
            }
        )

    # 🔥 sort by survival descending
    table.sort(key=lambda x: x["mean_survival"], reverse=True)

    return table


def print_mesh_summary(table, limit=20):
    print("\n=== Drug MeSH Survival Summary ===\n")

    for row in table[:limit]:
        print(f"MeSH Term: {row['mesh_term']}")
        print(f" Patients: {row['patients']}")
        print(f" Smokers: {row['smokers']}")
        print(f" Nonsmokers: {row['nonsmokers']}")
        print(f" Mean survival: {row['mean_survival']} days")

        print(" Therapy breakdown:")
        for k, v in row["therapy"].items():
            print(f"   {k}: {v}")

        print("-" * 60)
