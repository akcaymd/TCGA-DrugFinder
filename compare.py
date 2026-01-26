from GdcCases import GdcCases
from mesh_summary import build_mesh_summary
from comparison import compare_compounds


gdc = GdcCases()

project = input("Enter cancer project (e.g. TCGA-LUAD): ")

gdc.fetch(project)
gdc.get_data()

mesh_hash = gdc.get_mesh_list()

mesh_summary = build_mesh_summary(mesh_hash)

compound1 = input("\nEnter first compound name: ")
compound2 = input("Enter second compound name: ")

compare_compounds(compound1, compound2, mesh_summary)
