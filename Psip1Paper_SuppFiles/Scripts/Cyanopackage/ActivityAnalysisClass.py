import statistics

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import optimize
from sklearn.metrics import r2_score


class MetalActivityData:
    """ Provided a dataframe with the different columns for the different metals, it will fit the curves of the
    different abs/time, and it will calculate the R coefficient of the element"""

    def __init__(self, dataframe):
        self.df = dataframe

    def fit_curve(self):
        index = np.array(self.df.index)  # Generate the X values (x-axis)

        slopes = {}
        r_squares = {}
        for row in self.df:  # This will fit a curve using the x-axis for each column
            y = np.array(self.df[row])
            idx = np.isfinite(index) & np.isfinite(y)  # Test to prevent NaNs
            ab = np.polyfit(index[idx], y[idx], 1)  # Combine those elements with no NaNs
            predict = np.poly1d(ab)
            r_score = r2_score(y[idx], predict(index[idx]))
            if str(row).rsplit("_")[0] in slopes:
                slopes[str(row).rsplit("_")[0]].update({str(row).rsplit("_")[1]: ab[0]})
                r_squares[str(row).rsplit("_")[0]].update({str(row).rsplit("_")[1]: r_score})
            else:
                slopes.update({str(row).rsplit("_")[0]: {str(row).rsplit("_")[1]: ab[0]}})
                r_squares.update({str(row).rsplit("_")[0]: {str(row).rsplit("_")[1]: r_score}})

        return slopes, r_squares

    def plot_metal(self, metal, filename1, filename2):

        collector, rscore = self.fit_curve()
        metal_x = np.array(list(collector[metal].keys()), dtype='float')
        metal_y = np.array(list(collector[metal].values()), dtype='float')

        df_metal = self.df[self.df.filter(like=metal).columns]
        df_metal.plot()

        plt.savefig(filename1)
        plt.show()

        plt.plot(metal_x, metal_y, 'b-o')
        plt.xlabel('Concentration (mM)')
        plt.ylabel('R coefficient')

        plt.savefig(filename2)
        plt.show()

    def plot_curve_errorbar(self, metal):
        df = self.df[self.df.filter(like=str(metal)).columns]
        x = np.array(df.index)
        y = []
        z = []
        for index, row in df.iterrows():
            y.append(statistics.mean(row))
            z.append(statistics.stdev(row))

        plt.errorbar(x, y, yerr=z, fmt='-', capsize=2, capthick=1)

    def plot_raw_data_rcoefficient(self, list_to_plot, outfile1, outfile2):
        plt.figure(figsize=(8, 7))

        for element in list_to_plot:
            self.plot_curve_errorbar(element)

        plt.xlabel('Time(mins)')
        plt.ylabel('Abs 405nm')
        plt.legend(list_to_plot, loc='upper left', prop={'size': 18})
        plt.savefig(outfile1)
        plt.show()

        experiment = MetalActivityData(self.df)
        curves, rscore = experiment.fit_curve()

        metal = np.array([list(curves[x].values()) for x in list_to_plot])

        means = [statistics.mean(x) for x in metal]
        error = [statistics.stdev(x) for x in metal]

        x_pos = np.arange(len(list_to_plot))

        figure, ax = plt.subplots()
        ax.bar(x_pos, means, yerr=error, align='center', alpha=0.5, ecolor='black', capsize=10)
        ax.set_ylabel('Ratio of transformation per minute')
        ax.set_xticks(x_pos)
        ax.set_xticklabels(list_to_plot)
        ax.yaxis.grid(True)

        # Save the figure and show
        plt.tight_layout()
        plt.savefig(outfile2)
        plt.show()


def michaelis_menten(x, Vm, Km):
    return Vm * (x / (Km + x))


def plot_curve_errorbar(dataframe: pd.DataFrame, metal):
    df = dataframe[dataframe.filter(like=str(metal)).columns]
    x = np.array(df.index)
    y = []
    z = []
    for index, row in df.iterrows():
        y.append(statistics.mean(row))
        z.append(statistics.stdev(row))

    plt.errorbar(x, y, yerr=z, fmt='-', capsize=2, capthick=1)


def plot_abs_time(metal, dataframe):
    average = dataframe[dataframe.filter(like='_mean').columns]
    std = dataframe[dataframe.filter(like='_std').columns]
    values = average[average.filter(like=str(metal)).columns]
    for column in values:
        xvalues = np.array(values.index)
        y = np.array(values[column])
        errorvalues = np.array(std[str(column).replace('mean', 'std')])
        plt.errorbar(xvalues, y, yerr=errorvalues, fmt='-o', capsize=2, capthick=1)
        plt.xlabel('Time (mins)')
        plt.ylabel('Abs 405nm')
        plt.legend(list(values.columns), loc='upper left')


def r_coefficient_values(array1, array2, average_list, variance_list):
    import statistics
    for element in range(len(array1)):
        average_list.append(statistics.mean([array1[element], array2[element]]))
        variance_list.append(statistics.stdev([array1[element], array2[element]]))


def transform_substratedf_productdf(dfsubstrate, standard_df, proteinMg):
    pnp = np.array(standard_df.index)
    abs_y = np.array(standard_df['Mean_abs'])
    error = np.array(standard_df['std'])

    def curve1degree(x, a):
        return x * a

    p0 = [0.000000001]

    slope, cv = optimize.curve_fit(curve1degree, pnp, abs_y, absolute_sigma=True, p0=p0, sigma=error)

    product_df = pd.DataFrame(index=dfsubstrate.index, columns=list(dfsubstrate))
    for row in dfsubstrate:
        product_df[row] = (dfsubstrate[row]) / slope

    product_df = product_df / proteinMg
    return product_df
