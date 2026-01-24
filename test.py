from collections import defaultdict

mesh_stats = defaultdict(
    lambda: {"patients": set(), "survival": [], "therapy": defaultdict(int)}
)

for _, row in df_exploded.iterrows():

    mesh = row["mesh_list"]
    patient = row["submitter_id"]
    survival = row["days_to_death"]
    therapy = row["treatment_type"]

    mesh_stats[mesh]["patients"].add(patient)

    if pd.notna(survival):
        mesh_stats[mesh]["survival"].append(survival)

    mesh_stats[mesh]["therapy"][therapy] += 1

summary = []

for mesh, data in mesh_stats.items():

    mean_survival = (
        sum(data["survival"]) / len(data["survival"]) if data["survival"] else 0
    )

    summary.append(
        {
            "mesh_term": mesh,
            "patients": len(data["patients"]),
            "mean_survival": round(mean_survival, 2),
            "chemotherapy": data["therapy"].get("Chemotherapy", 0),
            "immunotherapy": data["therapy"].get("Immunotherapy", 0),
            "targeted_therapy": data["therapy"].get("Targeted Molecular Therapy", 0),
        }
    )

summary_df = pd.DataFrame(summary)
summary_df = summary_df.sort_values(by="mean_survival", ascending=False)
print(summary_df.head(10))
