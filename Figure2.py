### Original Code by Moira Arnet

import pandas as pd
import numpy as np 
import matplotlib.pyplot as plt

## load datasets
df = pd.read_csv('data/chemistry_and_qpcr_from_combine_monitoring_data_v1.1.csv')
mf = pd.read_csv('data/tunnel-outline.csv')


## subselect the dna samples
sequenced_samples = df[df['qPCR_mean'].notna()]
dna_fractures = sequenced_samples[sequenced_samples['Excel_Category']=='subsurface'].Excel_TM.unique()
dna_fractures_colors = ['tab:gray'] * len(dna_fractures)

## define the coloring for samples in the time series
for i,t in enumerate(dna_fractures):
    if t == 444:
        dna_fractures_colors[i] = 'tab:blue'
    elif t == 1306:
        dna_fractures_colors[i] = 'tab:orange'
    elif t == 1494:
        dna_fractures_colors[i] = 'tab:green'
    elif t == 2848:
        dna_fractures_colors[i] = 'tab:red'
    elif t == 4652:
        dna_fractures_colors[i] = 'tab:purple'
    

special_d = {444.0: 'tab:blue', 1306.0: 'tab:orange', 1494.0: 'tab:green', 2848.0: 'tab:red', 4652.0: 'tab:purple'}

## calculation of flow by M.A.
massetts = mf['Massets_cum_dis']
lutz = mf['lutz_total_spacing']
x = mf['X']
tunnel_height = mf['Y']


## Define the figure parameters
fig = plt.figure()
gs = fig.add_gridspec(5,1, hspace=0)

ax1 = fig.add_subplot(gs[0, 0])
ax2 = fig.add_subplot(gs[1, 0], sharex= ax1)
ax3 = fig.add_subplot(gs[2, 0], sharex= ax1)
ax4 = fig.add_subplot(gs[3, 0], sharex= ax1)
ax5 = fig.add_subplot(gs[4, 0], sharex= ax1)
ax6 =ax5.twinx()




### temperature plot
temperature = df.groupby('Excel_TM')['Temperature'].agg(['mean', 'std'])
temperature.columns = ['T_mean','T_std']

filtered_temperature = temperature[~temperature.index.isin(set(special_d.keys()))]
ax1.errorbar(filtered_temperature.index,filtered_temperature.T_mean,filtered_temperature.T_std, linestyle='None', marker='o', c='#FFADAD',mec='k', ecolor='k', capsize=5, elinewidth=1, markersize=8, alpha = 0.9)

for key in special_d:
    ax1.errorbar(key, temperature[temperature.index==key].T_mean, temperature[temperature.index==key].T_std, marker='*', markersize=18, c=special_d[key],mec='k',ecolor='k', elinewidth=1,linestyle='None', capsize=5)

ax1.set_ylabel('T (°C)', fontsize=13)
ax1.tick_params(which='both',labelsize=13)


## EC plot
ec = df.groupby('Excel_TM')['EC'].agg(['mean', 'std'])
ec.columns = ['EC_mean','EC_std']
filtered_ec = ec[~ec.index.isin(set(special_d.keys()))]

ax2.errorbar(filtered_ec.index,filtered_ec.EC_mean,filtered_ec.EC_std, linestyle='None', marker='o', c='#FDFFB6',mec='k', ecolor='k', capsize=5, elinewidth=1, markersize=8, alpha = 0.9)

for key in special_d:
    ax2.errorbar(key, ec[ec.index==key].EC_mean, ec[ec.index==key].EC_std, marker='*', markersize=18, c=special_d[key],mec='k',ecolor='k', elinewidth=1,linestyle='None', capsize=5)


ax2.set_ylabel('EC ($\mu$S cm$^{-1}$)',fontsize=13)
ax2.tick_params(which='both',labelsize=13)
ax2.set_yscale('log')
ax2.set_ylim(10,2000)

## oxygen isotope plot
ox = df.groupby('Excel_TM')['d18O_ppt'].agg(['mean', 'std'])
ox.columns = ['ox_mean','ox_std']
filtered_ox = ox[~ox.index.isin(set(special_d.keys()))]

ax3.errorbar(filtered_ox.index,filtered_ox.ox_mean,filtered_ox.ox_std, linestyle='None', marker='o', c='#CEFFC4',mec='k', ecolor='k', capsize=5, elinewidth=1, markersize=8, alpha = 0.9)

for key in special_d:
    ax3.errorbar(key, ox[ox.index==key].ox_mean, ox[ox.index==key].ox_std, marker='*', markersize=18, c=special_d[key],mec='k',ecolor='k', elinewidth=1,linestyle='None', capsize=5)


ax3.set_ylabel(u'$\delta^{18}$O (‰)', fontsize=13)
ax3.tick_params(which='both',labelsize=13)


## nitrate plot
nit = df[df['Nitrate']>0].groupby('Excel_TM')['Nitrate'].agg(['mean', 'std'])
nit.columns = ['nit_mean','nit_std']

filtered_nit = nit[~nit.index.isin(set(special_d.keys()))]

ax4.errorbar(filtered_nit.index,filtered_nit.nit_mean,filtered_nit.nit_std, linestyle='None', marker='o', c='#FFD6A5',mec='k', ecolor='k', capsize=5, elinewidth=1, markersize=8, alpha = 0.9)

for key in special_d:
    try: 
        nit[nit.index==key].nit_mean
        ax4.errorbar(key, nit[nit.index==key].nit_mean, nit[nit.index==key].nit_std, marker='*', markersize=18, c=special_d[key],mec='k',ecolor='k', elinewidth=1,linestyle='None', capsize=5)
    except:
        print('nitrate below detection at TM', key)

ax4.set_ylabel('Nitrate (ppm)', fontsize=13)
ax4.tick_params(which='both',labelsize=13)


## flow rate data 
index = np.arange(len(massetts))
bar_width = 100
bar_height = massetts
cd=ax6.bar(x,bar_height,bar_width,color="#79B8F4",alpha=0.5)

fr = df.groupby('Excel_TM')['Q'].agg(['mean', 'std'])
fr.columns = ['Q_mean','Q_std']
filtered_fr = fr[~fr.index.isin(set(special_d.keys()))]

ax5.errorbar(filtered_fr.index,filtered_fr.Q_mean,filtered_fr.Q_std, linestyle='None', marker='o', c='k',mec='k', ecolor='k', capsize=5, elinewidth=1, markersize=8, alpha = 0.9)

for key in special_d.keys():
    ax5.errorbar(key, fr[fr.index==key].Q_mean, fr[fr.index==key].Q_std, marker='*', markersize=18, c=special_d[key],mec='k',ecolor='k', elinewidth=1,linestyle='None', capsize=5)

ax5.set_xlabel('Tunnel meter (m)', fontsize=13)
ax5.set_ylabel('Discharge (L s$^{-1}$)', c='k',fontsize=13)
ax5.tick_params(colors = 'k',which='both',labelsize=13, color='k')
ax5.set_yscale('log')
ax5.set_ylim(0.0005,8000)


ax6.tick_params(labelsize=12, colors= '#79B8F4')
ax6.yaxis.label.set_color('#79B8F4')
ax6.set_ylabel('Cum. Inflow (L s$^{-1}$)', fontsize=13, color='#79B8F4')





plt.xlabel('Tunnelmeter',fontsize=12)
plt.xlim(0, 5220)

fig.align_labels()
fig.set_size_inches(10, 15)
plt.show()
