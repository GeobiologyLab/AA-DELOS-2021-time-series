import pandas as pd
import matplotlib.pyplot as plt

# Load the total environmental dataset
df = pd.read_csv('data/supplementary_data_2.csv')

df['date'] = pd.to_datetime(df['Date_ddmmmyy'],  format='%d.%m.%Y')

# Load the DNA sequence metadata
perSampleMDfiles = {"TM-444":"data/TM444_time_series_metadata.csv","TM-1306":"data/TM1306C_time_series_metadata.csv", "TM-1494":"data/TM1494_time_series_metadata.csv", "TM-2848":"data/TM2848A_time_series_metadata.csv", "TM-4652":"data/TM4652_time_series_metadata.csv"}

fracName = {"TM-444":"BdTfrac444","TM-1306":"BdTfrac1306C","TM-1494":"BdTfrac1494","TM-2848":"BdTfrac2848A","TM-4652":"BdTfrac4652"}
plt.rcParams['font.size'] = 14


def plot_environmental_factor(measurement,yaxislabel):
    if measurement == 'EC_field_ySpcm':
        plt.scatter(x= pd.to_datetime(['2021-05-04','2021-05-23'],format='%Y-%m-%d'),y=[55,55],label='Stimulation',color='tab:gray',marker='s',alpha=0.5,s=1)
    if measurement == 'T_degC':
        plt.scatter(x= pd.to_datetime(['2021-05-04','2021-05-23'],format='%Y-%m-%d'),y=[12,12],label='Stimulation',color='tab:gray',marker='s',alpha=0.5,s=1)
    for key in perSampleMDfiles:
        tmp = df[df['location_name']==fracName[key]]
        m = pd.read_csv(perSampleMDfiles[key])
        m['date'] = pd.to_datetime(m['date'],format='%Y%m%d')
        cmb = pd.concat([m,tmp],axis=0)
        print(cmb)
        cmb = cmb[['date', measurement]].dropna()
        cmb = cmb.drop_duplicates()
        cmb = cmb.sort_values('date')
        cmb = cmb[cmb['date'] >= '2020-11-01']
        if measurement == 'pH_field':
            cmb = cmb[cmb[measurement] < 12]

        plt.scatter(cmb['date'], cmb[measurement],label=key,alpha = 0.5)
        plt.plot(cmb['date'], cmb[measurement], '-', alpha=0.5)
    
    if measurement == 'EC_field_ySpcm':
        plt.yscale('log')
        plt.ylim(50,1000)

    plt.xlabel('Date')
    plt.ylabel(yaxislabel)
    plt.legend(loc='upper center', bbox_to_anchor=(0.5, 1.15), ncol=3)


# First subplot: Temperature
plt.subplot(2, 1, 1)
plot_environmental_factor('T_degC','Temperature (°C)')

# Second subplot: EC
plt.subplot(2, 1, 2)
plot_environmental_factor('EC_field_ySpcm','EC (µS/cm)')
plt.tight_layout()
plt.show()
