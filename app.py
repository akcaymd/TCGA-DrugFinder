from GdcCases import GdcCases


data = GdcCases()
data.fetch("TCGA-LUAD")
data.get_data()
data.get_mesh_list()

print(data.hash)
print(len(data.hash))
