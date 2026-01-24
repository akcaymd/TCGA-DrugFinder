from GdcCases import GdcCases
import pandas as pd


data = GdcCases()
data.fetch("TCGA-LUAD")
data.get_data()
hash_list = data.get_mesh_list()
df = pd.DataFrame(hash_list)
df_exploded = df.explode("mesh_list")
