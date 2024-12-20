import matplotlib.pyplot as plt
import pandas as pd

from upsetplot import from_contents
from upsetplot import UpSet
from upsetplot import plot

df = pd.read_csv('data/Rarefied_SNV_Table_for_Upset_plot.csv',index_col=0)

location = {}
location['Surface'] =  df[df.index.str.contains('BdS')].index.to_list() 
location['Grimsel'] =  df[df.index.str.contains('Gts')].index.to_list()
location['TM-0444'] = df[df.index.str.contains('BdTfrac444-')].index.to_list()
location['TM-1306'] = df[df.index.str.contains('BdTfrac1306C-')].index.to_list()
location['TM-1494'] = df[df.index.str.contains('BdTfrac1494-')].index.to_list()
location['TM-2848'] = df[df.index.str.contains('BdTfrac2848A-')].index.to_list() 
location['TM-4652'] = df[df.index.str.contains('BdTfrac4652-')].index.to_list()

print(location)

overlap_metadata = open('data/delos-surface-grimsel-sites-labels.csv','w')
overlap_metadata.write('full_name,location_id\n')

d = {}
for key in location:
    tmp = df.loc[location[key]]
    for id in location[key]:
        overlap_metadata.write(id+','+key+'\n')
    d[key] = tmp.columns[tmp.sum() > 1].tolist()
    print(key, len(d[key]))

overlap_metadata.close()

md = pd.read_csv('data/delos-surface-grimsel-sites-labels.csv')




dlocation = from_contents(d)
upset = UpSet(dlocation, subset_size="count",sort_by='cardinality')

upset.plot()

fig = plt.gcf()  

fig.set_size_inches(8.27*2,2.5*2.8) 

plt.yscale('log')
plt.xlim(-1,36)
# plt.xlim(6,36)
# plt.ylim(0,250)
plt.tight_layout()
plt.show()

