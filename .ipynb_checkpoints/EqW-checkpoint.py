#!/usr/bin/env python
# encoding: utf-8

import pyspeckit as p
import matplotlib.pyplot as plt, matplotlib.mlab as mlab
import numpy as np
from scipy.stats import norm
from astropy import log

"""Calculate the equivalent width of an absorption or emission line for a given spectrum using PySpecKit. By: Munazza Alam
=======
Args:
=======
xmin,xmax - the specified interval of the spectrum to plot
excludemin, excludemax - the specified interval (in wavelength space) of the absorption feature
n - the number of Monte Carlo iterations

=======
Returns:
=======
- the mean and standard deviation of the equivalent width measured n times
- the spectrum plotted with the Voigt profile line fit (blue), the pseudo-continuum (yellow),
  and the approximated rectangle (green)
- a histogram of the EqW distribution
"""

def measure_equivalent_width(filename,xmin,xmax,exclude_min,exclude_max,n,
                             xunit='micron',xtype='wavelength',to_plot=True,
                             filebase="EQW.pdf"):
    sp = p.Spectrum(filename)
    sp.xarr.units = xunit
    sp.xarr.xtype = xtype

    if np.any(sp.error!=0):
        errs = sp.error
    else:
        w = sp.xarr.value
        continuum = np.append(sp.data[(w>=6558) & (w<=6562)],
                              sp.data[(w>=6566) & (w<=6569)])
        errs = np.std(continuum)
    # print("errs",errs)

    if to_plot is True:
        sp.plotter(xmin=xmin, xmax=xmax, ymin=0, errstyle='bars',
                   color='grey')
        ax1 = plt.gca()
    sp.baseline(xmin=xmin, xmax=xmax, exclude=[exclude_min,exclude_max],
                subtract=False, highlight_fitregion=False,
                selectregion=True, order=0)
    sp.specfit(plot=True, fittype='voigt', color='blue',
               guesses='moments', vheight=True)
    sp.specfit.EQW(plot=True, plotcolor='g', fitted=True, components=False,
                   annotate=True, loc='lower left', xmin=None, xmax=None)
    plt.savefig(filebase.replace(".","_initialfit."),bbox_inches="tight")

    sp2 = sp.copy()
    EQWs = np.zeros(n)

    if to_plot is True:
        fig, ax = plt.subplots(1,1)
        ax.step(sp.xarr,sp.data,color="C3",zorder=100)
        ax.set_xlim(xmin,xmax)
        ax.set_ylim(ax1.get_ylim())
    ds = sp.data.size
    for w in range(n):
        sp2.data = sp.data + np.random.randn(ds)*errs
        if to_plot is True:
            ax.plot(sp2.xarr,sp2.data,'.',alpha=0.01,color="Grey",zorder=-1)
        # print("data")
        sp2.baseline(xmin=xmin, xmax=xmax, exclude=[exclude_min,exclude_max],
                     subtract=False, highlight_fitregion=False,
                     selectregion=True, order=0)
        # print("base")
        sp2.specfit(fittype='voigt', guesses=sp.specfit.parinfo.values,
                    vheight=True,plot=False)
        # print("voight")
        dist = sp2.specfit.EQW(plotcolor='g', fitted=True, components=False,
                               annotate=True, loc='lower left', xmin=None, xmax=None)
        # print("fit")
        # print(dist)
        EQWs[w] = dist
    plt.savefig(filebase.replace(".","_mc."),bbox_inches="tight")

    mu,sigma = norm.fit(EQWs)
    print( mu, sigma)
    # Plot the median and 68th percentile values for comparison
    perc = np.percentile(EQWs, [16, 50, 84])
    print(perc)
    
    target_name1 = filename.split('/')[9]
    target_name = filename.split('.')[1:2]
    date = filename.split('/')[7]
    
    #save perc data to .csv file
    import csv
  
    fields=[str(target_name), str(date), str(perc)]
    with open('/Users/amandaash/Desktop/Research/EW.csv', 'a') as f:
        writer = csv.writer(f)
        writer.writerow(fields)

    if to_plot is True:
        plt.figure()
        num_bins = max(10,np.int(np.log10(n)) * 10)
        nn,bins,patches = plt.hist(EQWs, num_bins, normed=True, facecolor='green', histtype='stepfilled')
        y = mlab.normpdf(bins,mu,sigma)
        plt.plot(bins,y,'r--',linewidth=2)
        plt.grid(True)
        plt.ylabel('Probability')
        plt.xlabel('EQW')
        # plt.show()

        # Light blue histogram for contrast and for R/G colorblind folks
        n,bins,patches = plt.hist(EQWs, num_bins, normed=True,
                                  facecolor='lightblue',
                                  histtype='stepfilled')
        y = mlab.normpdf(bins,mu,sigma)
        plt.plot(bins,y,'r--',linewidth=2)
        ax = plt.gca()
        for p_value in perc:
            ax.axvline(p_value,linestyle=":",color="k",lw=1.5)

        plt.ylabel('Probability')
        plt.xlabel('EQW')
        # plt.show()
    plt.savefig(filebase.replace(".","_eqwhist."),bbox_inches="tight")

    plt.close("all")

    # Return the percentiles, since these aren't necessarily Gaussian distributions
    return perc


if __name__=="__main__":
    # xmin,xmax,exclude_min,exclude_max,n
    equivalent_width('2m1821_61_08jun11.txt',1.250, 1.255, 1.2514,
                     1.2538, 1000)

def final_plots(filename,xmin,xmax,exclude_min,exclude_max):
    vf = p.spectrum.models.inherited_voigtfitter.voigt_fitter()

    sp = p.Spectrum(filename)
    sp.xarr.units = 'micron'
    sp.xarr.xtype = 'wavelength'
    sp.plotter(xmin=xmin, xmax=xmax, ymin=0, errstyle='fill',color='grey')
    sp.baseline(xmin=xmin, xmax=xmax,exclude=[exclude_min,exclude_max],subtract=False,
                reset_selection=False,hightlight_fitregion=False,order=0)
    sp.specfit(plot=True, fittype='voigt', color='magenta', guesses='moments',
               vheight=True)
    fwhm = sp.specfit.measure_approximate_fwhm(threshold='error', emission=False,
                                               interpolate_factor=1024, plot=True,
                                               grow_threshold=1, color='magenta')
    sp.plotter.refresh()
    sp.specfit.EQW(plot=True, plotcolor='g', fitted=False, continuum=0.5, components=False, annotate=True, loc='lower left', xmin=None,
    xmax=None)
    sp.plotter.refresh()
    xarr_fit_units = 'microns'
    plt.ylabel('Normalized Flux')
    plt.xlabel('Wavelength ($\mu m$)')
