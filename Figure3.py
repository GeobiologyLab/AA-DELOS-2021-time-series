import pandas as pd
import numpy as np 
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib as mpl
from matplotlib.colors import ListedColormap


## prepare DELOS Data
df = pd.read_csv('data/chemistry_and_qpcr_from_combine_monitoring_data_v1.1.csv')
df = df[~df['Excel_Category'].isin(['experimental'])]
precip = df[df['type']=='precipitation']
precip = precip.dropna(subset=['d18O_ppt','dD_ppt'])
df = df[~df['type'].isin(['precipitation'])]

subsurface = df[df['Excel_Category']=='subsurface']
snow = df[df['type']=='snow']


# Get meteoric water line

swissprecip = pd.read_csv('data/regional_precipdata.csv',comment='#') #data from Ofterdinger (2001)
localprecip = pd.read_csv("data/wateriso.uta.edu_estimate.csv",comment='#') #estimates from waterisotopes.org
localprecip = localprecip.drop('Elevation', axis=1)


xmin = df[df['Excel_Category']=='surface']['d18O_ppt'].min()
xmax = df[df['Excel_Category']=='surface']['d18O_ppt'].max()
# xmax = precip['d18O_ppt'].max()

x = np.linspace(xmin-0.2, xmax+0.2, 100)

# mean fit
slope, intercept = np.polyfit(localprecip['d18O_ppt'], localprecip['dD_ppt'], 1)
y = slope * x + intercept
plt.plot(x, y, 'k--', label='LMWL')


y = slope*localprecip['d18O_ppt'] + intercept
y_upper = y + localprecip['dD_CI']
y_lower = y - localprecip['dD_CI']


# plot shaded region
plt.fill_between(localprecip['d18O_ppt'], y_lower, y_upper, color='gray', alpha=0.3)

# plot DELOS
n = 256
new_colors = cm.BuPu(np.linspace(0.3, 1, n)) #GnBu is green blue
truncated_GnBu = ListedColormap(new_colors)
norm = mpl.colors.Normalize(vmin=subsurface['Excel_TM'].min(), vmax=subsurface['Excel_TM'].max())
subsurface_colors = truncated_GnBu(norm(subsurface['Excel_TM']))
sc = plt.scatter(subsurface['d18O_ppt'], subsurface['dD_ppt'], alpha=0.25, label='DELOS Fractures', c=subsurface_colors)
cbar = plt.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=truncated_GnBu), alpha=.5, ax=plt.gca())
cbar.set_label('Tunnel Meter (TM)')



surface_df = df[df['Excel_Category']=='surface']
plt.plot(surface_df['d18O_ppt'], surface_df['dD_ppt'],marker = 's',fillstyle='none',markersize=6,label='DELOS Surface',c='k',linestyle='None')


## add mean and var of time series points
dfmv = subsurface.groupby('Excel_TM').agg({'d18O_ppt': ['mean', 'std'], 'dD_ppt': ['mean', 'std']})
dfmv.columns = ['d18O_ppt_mean', 'd18O_ppt_std', 'dD_ppt_mean', 'dD_ppt_std']
print(dfmv.dropna())

d = {444.0:'tab:blue',1306.0:'tab:orange',1494.0:'tab:green',2848.0:'tab:red',4652.0:'tab:purple'}

for key in d:
    print(dfmv.loc[key])
    if key == 1306:
        # keytm = '1306C'
        keytm = '1306'
    elif key == 2848:
        # keytm = '2848C'
        keytm='2848'
    else:
        keytm = str(key)
        keytm = keytm.split('.')[0]
    plt.errorbar(dfmv.loc[key, 'd18O_ppt_mean'], dfmv.loc[key, 'dD_ppt_mean'], xerr=dfmv.loc[key, 'd18O_ppt_std'], yerr=dfmv.loc[key, 'dD_ppt_std'], fmt='o',markersize=8, color=d[key], label=f'TM-{keytm}',markeredgecolor='k')


### Load data from Schneeberger, Raphael, Urs K. Mäder, and H. Niklaus Waber. "Hydrochemical and isotopic (δ2H, δ18O, 3H) characterization of fracture water in crystalline rock (Grimsel, Switzerland)." Procedia Earth and Planetary Science 17 (2017): 738-741.
gts = pd.read_csv('data/Schneeberger_waterIsotopes_Gts2017.csv',comment='#')

plt.scatter(gts['d18O'],gts['d2H'],marker = 'x',c='k',label='Grimsel Test Site')

plt.legend(loc='lower right')
plt.xlabel(u'$\delta^{18}$O (‰)', fontsize=18)
plt.ylabel(u'$\delta^{2}$H (‰)', fontsize=18)
plt.xlim(xmin-.2,xmax+.2)
plt.show()
