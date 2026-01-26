from pccompound import pccompound
import statistics


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


def compare_compounds(compound1, compound2, mesh_summary):
    """
    Prints side-by-side TCGA statistics for two compounds.
    """

    c1_mesh = pccompound(compound1)
    c2_mesh = pccompound(compound2)

    mesh1 = set()
    mesh2 = set()

    for m in c1_mesh:
        mesh1.update(m["mesh_list"])

    for m in c2_mesh:
        mesh2.update(m["mesh_list"])

    stat1 = summarize_mesh_terms(mesh1, mesh_summary)
    stat2 = summarize_mesh_terms(mesh2, mesh_summary)

    print("\n" + "=" * 70)
    print("                SIDE BY SIDE COMPOUND COMPARISON")
    print("=" * 70)

    header = f"{'Metric':30s}{compound1:20s}{compound2:20s}"
    print(header)
    print("-" * 70)

    def row(label, v1, v2):
        print(f"{label:30s}{str(v1):20s}{str(v2):20s}")

    row("MeSH terms", stat1["mesh_count"], stat2["mesh_count"])
    row("Patients", stat1["patients"], stat2["patients"])
    row("Smokers", stat1["smokers"], stat2["smokers"])
    row("Non-smokers", stat1["nonsmokers"], stat2["nonsmokers"])
    row("Mean survival (days)", stat1["mean_survival"], stat2["mean_survival"])
    row("Median survival (days)", stat1["median_survival"], stat2["median_survival"])

    print("-" * 70)

    all_therapies = set(stat1["therapy"]) | set(stat2["therapy"])

    for therapy in sorted(all_therapies):
        row(
            therapy,
            stat1["therapy"].get(therapy, 0),
            stat2["therapy"].get(therapy, 0),
        )

    print("=" * 70)
