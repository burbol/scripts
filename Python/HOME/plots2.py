#!/usr/bin/python

for k in range(0, 20, 2):
    number % 2 == 0
    
    def plotGraph(X,Y):
      fig = plt.figure()
      ### Plotting arrangements ###
      return fig
------ plotting module ------

----- mainModule ----

from matplotlib.backends.backend_pdf import PdfPages

plot1 = plotGraph(tempDLstats, tempDLlabels)
plot2 = plotGraph(tempDLstats_1, tempDLlabels_1)
plot3 = plotGraph(tempDLstats_2, tempDLlabels_2)

pp = PdfPages('foo.pdf')
pp.savefig(plot1)
pp.savefig(plot2)
pp.savefig(plot3)
pp.close()

def plot_circ(z, R, xc, yc, R_2):
    theta_fit = linspace(-pi, pi, 180)
    x_fit2 = xc + R_2*cos(theta_fit)
    y_fit2 = yc + R_2*sin(theta_fit)
    pl.plot(x_fit2, y_fit2, color="blue", label=method, lw=2)
    xlim(0,4.0)
    ylim(0,3)
    pl.plot(z, R, 'k.')
