import os
import statistics

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import matplotlib as mpl

from scipy import optimize

from Cyanopackage.ActivityAnalysisClass import MetalActivityData, michaelis_menten, transform_substratedf_productdf

mpl.rcParams['font.size'] = 25
mpl.rcParams['font.family'] = 'Arial'
fig, ax = plt.subplots(2, 2, figsize=(22, 18), dpi=450)
# %% Metal figure

os.chdir('/Users/u2176312/OneDrive - University of Warwick/Thesis/Paper Draft/')
TableMetals = pd.read_excel('/Users/u2176312/OneDrive - University of Warwick/Thesis/Paper '
                            'Draft/Andrew_PsipMetals.xlsx', sheet_name='Sheet3', header=0, index_col=0)
x_value = list(TableMetals.index)
ax[0, 0].spines['top'].set_visible(False)
ax[0, 0].spines['right'].set_visible(False)

shapedict = {
    'Calcium': 's',
    'Iron': 'o',
    'Calcium + Iron': '*'
}
for column in TableMetals.columns:
    if '_std' in str(column):
        continue

    if 'Calcium +' in str(column):
        y = list(TableMetals[column])
        error = list(TableMetals[column + '_std'])
        ax[0, 0].errorbar(x_value, y, yerr=error, capsize=3, capthick=1,
                          fmt='--' + str(shapedict[str(column)]),
                          markerfacecolor='black',
                          color='black',
                          markeredgecolor='black',
                          label=str(column),
                          markersize=16)

    elif 'Calcium' in str(column):
        y = list(TableMetals[column])
        error = list(TableMetals[column + '_std'])
        ax[0, 0].errorbar(x_value, y, yerr=error, capsize=3, capthick=1,
                          fmt='--' + str(shapedict[str(column)]),
                          markerfacecolor='black',
                          color='black',
                          markeredgecolor='black',
                          label=str(column),
                          markersize=12)

    else:
        y = list(TableMetals[column])
        error = list(TableMetals[column + '_std'])
        ax[0, 0].errorbar(x_value, y, yerr=error, capsize=3, capthick=1,
                          fmt='-' + str(shapedict[str(column)]),
                          markerfacecolor='none',
                          color='black',
                          markeredgecolor='black',
                          label=str(column),
                          markersize=14)
ax[0, 0].set_xlabel('Time (mins)')
ax[0, 0].set_ylabel('Abs(405nm)')
ax[0, 0].set_title('A', loc='left', weight='bold', y=1.1, fontsize=28)
ax[0, 0].legend(loc='upper left', fontsize=20)

# %% pH figure
os.chdir('/Users/u2176312/Downloads')

trial = pd.read_excel('Psip1_experiment_pag53_table.xlsx', sheet_name='Study_noOut', index_col=0)
experiment = MetalActivityData(trial)
curves, rvalue = experiment.fit_curve()

pH = list([list(curves['pH6.8'].values()), list(curves['pH7.5'].values()), list(curves['pH8.8'].values()),
           list(curves['pH9.4'].values()), list(curves['pH9.8'].values()), list(curves['pH10.4'].values()),
           list(curves['pH11.2'].values())])

r_score = list([list(rvalue['pH6.8'].values()), list(rvalue['pH7.5'].values()), list(rvalue['pH8.8'].values()),
                list(rvalue['pH9.4'].values()), list(rvalue['pH9.8'].values()), list(rvalue['pH10.4'].values()),
                list(rvalue['pH11.2'].values())])
means = [statistics.mean(x) for x in pH]
error = [statistics.stdev(x) for x in pH]
x = [6.8, 7.5, 8.8, 9.4, 9.8, 10.4, 11.2]

ax[0, 1].errorbar(x, means, yerr=error, fmt='-o', capsize=3, capthick=1, color='black', markersize=10)
ax[0, 1].set_xlabel('pH values')
ax[0, 1].set_ylabel(r'$\Delta$Abs(405nm) min$^{-1}$ mg protein$^{-1}$')
ax[0, 1].spines['top'].set_visible(False)
ax[0, 1].spines['right'].set_visible(False)
ax[0, 1].spines['bottom'].set_color('black')
ax[0, 1].spines['left'].set_color('black')
ax[0, 1].set_title('B', loc='left', weight='bold', y=1.1, fontsize=28)

# %% Km MUF-P
os.chdir('/Users/u2176312/OneDrive - University of Warwick/Thesis/Paper Draft/')
Data_MUFP = pd.read_excel('MUFPkinetics.xlsx', sheet_name='MUFP', index_col=0)


collection = {}

for index, row in Data_MUFP.iterrows():
    if index not in collection.keys():
        collection.update({index: [row['Rate (nmoles/min/mg_protein)']]})
    else:
        collection[index].append(row['Rate (nmoles/min/mg_protein)'])

Y_values = [statistics.mean(collection[x]) for x in collection.keys()]
Y_values_2 = Data_MUFP['Rate (nmoles/min/mg_protein)']
error_values = [statistics.stdev(collection[x]) for x in collection.keys()]

amounts = list(collection.keys())
amounts_2 = Data_MUFP.index.tolist()

# --- Fitting of the
p0 = [0.00000001, 0.00000001]
params, cv = optimize.curve_fit(michaelis_menten,amounts_2, Y_values_2)
vmax, km = params
print("Vmax = ", vmax)
print("Km = ", km)
print('Estimated variance (Vm, Km) = ' + str(cv[0, 0]) + ', ' + str(cv[1, 1]))
print('Estimated standard devitation (Vm, Km) = ', np.sqrt(np.diag(cv)))

residuals = Y_values - michaelis_menten(amounts, *params)
ss_res = np.sum(residuals ** 2)
ss_tot = np.sum((Y_values - np.mean(Y_values)) ** 2)
r_squared = 1 - (ss_res / ss_tot)

ax[1, 0].errorbar(amounts, Y_values, yerr=error_values, fmt='o', capsize=3, capthick=1, color='black', markersize=10)
ax[1, 0].plot(np.linspace(0, 10, 1000), michaelis_menten(np.linspace(0, 15, 1000), vmax, km),
              "k--", label='Fitted Michaelis-Menten equation')
ax[1, 0].text(5, 1, u"R\u00b2= {:0.2f}".format(r_squared), style='italic', fontsize=22)
ax[1, 0].spines['top'].set_visible(False)
ax[1, 0].spines['right'].set_visible(False)
ax[1, 0].spines['bottom'].set_color('black')
ax[1, 0].spines['left'].set_color('black')
ax[1, 0].set_xlabel('MUFP (μM)')
ax[1, 0].set_ylabel(r'Velocity (nmoles min$^{-1}$ mg protein$^{-1}$)')
ax[1, 0].set_title('C', loc='left', weight='bold', y=1.1, fontsize=28)

ax[1, 0].legend(loc='lower right', fontsize=22)


#%% Km Bis-MUFP

os.chdir('/Users/u2176312/OneDrive - University of Warwick/Thesis/Paper Draft/')
Data_MUFP = pd.read_excel('MUFPkinetics.xlsx', sheet_name='BisMUFP', index_col=0)


collection = {}

for index, row in Data_MUFP.iterrows():
    if index not in collection.keys():
        collection.update({index: [row['Rate (nmoles/min/mg_protein)']]})
    else:
        collection[index].append(row['Rate (nmoles/min/mg_protein)'])

Y_values = [statistics.mean(collection[x]) for x in collection.keys()]
Y_values_2 = Data_MUFP['Rate (nmoles/min/mg_protein)']
error_values = [statistics.stdev(collection[x]) for x in collection.keys()]

amounts = list(collection.keys())
amounts_2 = Data_MUFP.index.tolist()

# --- Fitting of the
p0 = [1, 1]
params, cv = optimize.curve_fit(michaelis_menten,amounts_2, Y_values_2, p0=p0)
vmax, km = params
print("Vmax = ", vmax)
print("Km = ", km)
print('Estimated variance (Vm, Km) = ' + str(cv[0, 0]) + ', ' + str(cv[1, 1]))
print('Estimated standard devitation (Vm, Km) = ', np.sqrt(np.diag(cv)))

residuals = Y_values - michaelis_menten(amounts, *params)
ss_res = np.sum(residuals ** 2)
ss_tot = np.sum((Y_values - np.mean(Y_values)) ** 2)
r_squared = 1 - (ss_res / ss_tot)

ax[1, 1].errorbar(amounts, Y_values, yerr=error_values, fmt='o', capsize=3, capthick=1, color='black', markersize=10)
ax[1, 1].plot(np.linspace(0, 200, 1000), michaelis_menten(np.linspace(0, 200, 1000), vmax, km),
              "k--", label='Fitted Michaelis-Menten equation')
ax[1, 1].text(6, 1, u"R\u00b2= {:0.2f}".format(r_squared), style='italic', fontsize=22)
ax[1, 1].spines['top'].set_visible(False)
ax[1, 1].spines['right'].set_visible(False)
ax[1, 1].spines['bottom'].set_color('black')
ax[1, 1].spines['left'].set_color('black')
ax[1, 1].set_xlabel('Bis-MUFP (μM)')
ax[1, 1].set_ylabel(r'Velocity (nmoles min$^{-1}$ mg protein$^{-1}$)')
ax[1, 1].set_title('D', loc='left', weight='bold', y=1.1, fontsize=28)

ax[1, 1].legend(loc='lower right', fontsize=22)


plt.tight_layout()
plt.savefig('/Users/u2176312/OneDrive - University of Warwick/'
            'Thesis/Paper Draft/FinalDocumentsPaper/ReviewManuscript/To_submit_reviewApril/Fig_1_Final_Mar24.pdf', dpi=1000)
