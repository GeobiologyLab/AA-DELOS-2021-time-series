import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import pandas as pd
import matplotlib.patches as mpatches
from wqchartpy import triangle_piper



df = pd.read_csv('data/Non-Bedretto-wqchartpy-piper.csv')
df = df[~(df['location'] == 'stripa')]

bd = pd.read_csv('data/Bedretto_data_wqchartpy_piper.csv')

ions = ['Ca', 'Na', 'Cl', 'HCO3','SO4','Mg','K','CO3']

bd = bd.dropna(subset=ions)


frac = bd[bd['Label']=='DELOS Gneiss Fracture']
frac = frac.groupby('Tunnel_Location_tm')[ions].mean().reset_index()
frac['Alpha'] = 0.6
frac['Color'] = 'tab:orange'
frac['Size']   = 50
frac['Marker'] = '^'
frac['Label'] = 'DELOS Gneiss Fracture'

gran = bd[bd['Label']=='DELOS Granite Fracture']
gran = gran.groupby('Tunnel_Location_tm')[ions].mean().reset_index()
gran['Alpha'] = 0.6
gran['Color'] = 'tab:red'
gran['Size']   = 50
gran['Marker'] = 's'
gran['Label'] = 'DELOS Granite Fracture'




combined_df = pd.concat([frac, gran], ignore_index=True)
combined_df = pd.concat([combined_df, df], ignore_index=True)

print(combined_df)

triangle_piper.plot(combined_df, unit='mg/L', figname='piper-diagram-chalk-aspo-test', figformat='svg')
 

