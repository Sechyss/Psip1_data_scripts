import matplotlib.pyplot as plt
import numpy as np
from scipy import stats


def plot_lineweaverburk(s, v):
    # take reciprocal of Michaelis-Menten to plot in Lineweaver-Burke
    recip_S = 1 / np.array(s)
    recip_V = 1 / np.array(v)

    # Calculate linear regression
    slope, intercept, r_value, p_value, std_err = stats.linregress(recip_S, recip_V)

    # Draw linear regression
    x = np.linspace(-.8, .5, 1000)
    y = []
    for i in x:
        y = (slope * x) + intercept

    # Plot 1/Vmax
    plt.scatter(0, intercept, color='red')
    print("\n1/Vmax = ", intercept)
    print("Vmax = ", 1 / intercept)

    # Plot -1/Km
    xintercept = ((0 - intercept) / slope)
    plt.scatter(xintercept, 0, color='red')
    print("\n-1/Km = ", xintercept)
    Km = (-1 / xintercept)
    print("Km = ", Km)
    print("\nKm/Vmax (slope): ", slope)

    # Draw x & y origins
    plt.axhline(0, color='black')
    plt.axvline(0, color='black')

    # Graph scatter points
    plt.scatter(recip_S, recip_V)

    # Graph linear regresison
    plt.plot(x, y)

    # Titles and labels
    plt.xlabel('1/[S] ($\mu$M)')
    plt.ylabel('1/v (nmoles/min/mg of protein)')
    plt.title('Lineweaver-Burk')

