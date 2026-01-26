#!/usr/bin/env python3

import os
from urllib.parse import parse_qs
from GdcCases import GdcCases
from mesh_summary import build_mesh_summary
from comparison import compare_compounds

query = os.environ.get("QUERY_STRING", "")
params = parse_qs(query)

project = params.get("query", [""])[0]
compound1 = params.get("comp1", [""])[0]
compound2 = params.get("comp2", [""])[0]

gdc = GdcCases()
gdc.fetch(project)
gdc.get_data()
mesh_hash = gdc.get_mesh_list()

mesh_summary = build_mesh_summary(mesh_hash)


print("Content-Type: text/html\n")

print(compare_compounds(compound1, compound2, mesh_summary))
