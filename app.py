from GdcCases import GdcCases
from pccompound import pccompound
import time


data = GdcCases()
data.fetch("TCGA-LUAD")
patients = data.get_cleaned_data()

for i in data.hash:
    time.sleep(1)
    meshlist = pccompound(i["therapeutic_agents"])
    print(i["therapeutic_agents"], meshlist)
    print("*" * 20)
print(len(data.hash))
