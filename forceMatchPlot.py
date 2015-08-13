__author__ = 'Greg Poisson'

import numpy
import matplotlib.pyplot as plt

def makePlots(plots):
    print "FMPlot - Make plots"
    print len(plots)
    print len(plots[1])
    print len(plots[1][0])

# Draw subplots for each data set
def plotAllData(vars, color='r'):
    plt.close()

    for pair in range(0, len(vars.plots)):
        if pair % 2 == 0:
            xAxis = numpy.arange(len(vars.plots[pair+1][0])) * vars.binSize
            axes = []       # List of the requested plot order

            for ax in range(1, len(vars.plotOrder)+1):
                # go through plotOrder and arrange these plots.
                if ax in vars.plotOrder:
                    axes.append(vars.plotOrder.index(ax))


            # Create plots and labels
            if len(axes) == 0:        # If zero plots are requested, end the program
                break
            elif len(axes) > 1:       # Multiple plots
                f, yAxes = plt.subplots(len(axes), sharex=True)

                # Draw all requested plots
                for i in range(0, len(axes)):
                    yAxes[i].scatter(xAxis, vars.plots[pair+1][axes[i]], c=color, s=15)
                    yAxes[i].plot(xAxis, vars.plots[pair+1][axes[i]])
                    yAxes[i].grid(True)
                    yAxes[i].set_title(vars.titles[axes[i]])
                    yAxes[i].set_ylabel(vars.yLabels[axes[i]], fontsize=10)
                    #coef0 = numpy.polyfit(xAxis, plots[pair+1][0], 3)
                    #poly0 = numpy.poly1d(coef0)
                    #ax0fit = poly0(xAxis)
                    #ax0.plot(ax0fit)

                plt.suptitle("{} Interactions".format(vars.plots[pair]))

            elif len(axes) == 1:           # Single plot
                plt.scatter(xAxis, vars.plots[pair+1][axes[0]], c=color, s=15)
                plt.plot(xAxis, vars.plots[pair+1][axes[0]])
                plt.grid(True)
                plt.suptitle("{} -- {} Interactions".format(vars.titles[axes[0]], vars.plots[pair]))
                plt.ylabel(vars.yLabels[axes[0]], fontsize=10)

            plt.xlabel("Intermolecular Distance (Angstroms)", fontsize=10)

    plt.show()