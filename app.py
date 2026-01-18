from GdcCases import GdcCases


data = GdcCases()
data.fetch("TCGA-LUAD")
data.get_cleaned_data()


print(data.get_mesh_list())
print(len(data.patients))
print(len(data.hash))
