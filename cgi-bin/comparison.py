from pccompound import pccompound
import statistics
import pandas as pd


def summarize_mesh_terms(mesh_terms, mesh_summary):
    """
    Combine statistics of multiple MeSH terms
    into a single compound-level summary.
    """

    patients = set()
    survival = []
    smokers = 0
    nonsmokers = 0
    therapy = {}

    for mesh in mesh_terms:
        if mesh not in mesh_summary:
            continue

        data = mesh_summary[mesh]

        patients |= data["patients"]
        survival += data["survival_days"]
        smokers += data["smokers"]
        nonsmokers += data["nonsmokers"]

        for k, v in data["therapy"].items():
            therapy[k] = therapy.get(k, 0) + v

    if survival:
        mean_survival = int(sum(survival) / len(survival))
        median_survival = int(statistics.median(survival))
    else:
        mean_survival = 0
        median_survival = 0

    return {
        "mesh_count": len(mesh_terms),
        "patients": len(patients),
        "smokers": smokers,
        "nonsmokers": nonsmokers,
        "mean_survival": mean_survival,
        "median_survival": median_survival,
        "therapy": therapy,
    }


def compare_compounds(compound1, compound2, mesh_summary) -> pd.DataFrame:
    """
    Returns side-by-side compound comparison as pandas DataFrame
    """

    # ── get mesh terms from PubChem
    c1_mesh = pccompound(compound1)
    c2_mesh = pccompound(compound2)

    mesh1 = set()
    mesh2 = set()

    for m in c1_mesh:
        mesh1.update(m["mesh_list"])

    for m in c2_mesh:
        mesh2.update(m["mesh_list"])

    # ── summarize
    stat1 = summarize_mesh_terms(mesh1, mesh_summary)
    stat2 = summarize_mesh_terms(mesh2, mesh_summary)

    # ── build dataframe
    metrics = sorted(set(stat1.keys()) | set(stat2.keys()))

    data = []

    for metric in metrics:
        data.append(
            {
                "Metric": metric,
                compound1: stat1.get(metric, 0),
                compound2: stat2.get(metric, 0),
            }
        )

    df = pd.DataFrame(data).set_index("Metric")

    return df
