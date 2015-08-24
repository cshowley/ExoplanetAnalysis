import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import csv
import os
from collections import defaultdict


columns = defaultdict(list)
with open("planetData.csv","rU") as data:
    reader = csv.DictReader(data,delimiter=',')
    for row in reader:
        for (h,v) in row.items():
            columns[h].append(v)

topFluxes = defaultdict(list)
with open("topFluxes.csv","rU") as data:
    reader = csv.DictReader(data,delimiter=',')
    for row in reader:
        for (h,v) in row.items():
            topFluxes[h].append(v)

if not os.path.exists('Correlation Plots'):
    os.makedirs('Correlation Plots')

def histograms():
    e = []
    for i in columns['Eccentricity']:
        e.append(float(i))
    plt.figure()
    plt.hist(e,bins=100)
    plt.xlabel('Eccentricity')
    plt.ylabel('Counts')
    plt.title('Eccentricity')
    plt.yscale('log',nonposy='clip')
    plt.savefig(os.path.join("Correlation Plots/","Eccentricity1"))

    e = []
    for i in columns['Eccentricity']:
        if i == '0.0':
            continue
        else:
            e.append(float(i))
    plt.figure()
    plt.hist(e,bins=100)
    plt.xlabel('Eccentricity')
    plt.ylabel('Counts')
    plt.title('Eccentricity')
    plt.yscale('log',nonposy='clip')
    plt.savefig(os.path.join("Correlation Plots/","Eccentricity2"))

    e = []
    for i in columns['Argument of Periastron']:
        if i == '90.0':
            continue
        else:
            e.append(float(i))
    plt.figure()
    plt.hist(e,bins=100)
    plt.xlabel('Argument of Periastron')
    plt.ylabel('Counts')
    plt.title('Argument of periastron')
    plt.xlim([0, 360])
    plt.yscale('log',nonposy='clip')
    plt.savefig(os.path.join("Correlation Plots/","Arg_periastron"))

    e = []
    for i in columns['Orbital Inclination']:
        if i == '90.0':
            continue
        else:
            e.append(float(i))
    plt.figure()
    plt.hist(e,bins=100)
    plt.xlabel('Orbital Inclination')
    plt.ylabel('Counts')
    plt.title('Inclination')
    plt.yscale('log',nonposy='clip')
    plt.savefig(os.path.join("Correlation Plots/","Inclination"))

    e = []
    for i in columns['Planetary Radius']:
        if i == '1.0':
            continue
        else:
            e.append(float(i))
    plt.figure()
    plt.hist(e,bins=100)
    plt.xlabel('Radius (Jupiter-radii)')
    plt.ylabel('Counts')
    plt.title('Planetary Radius')
    plt.yscale('log',nonposy='clip')
    plt.savefig(os.path.join("Correlation Plots/","Radius"))
    
    e = []
    for i in columns['Period']:
        e.append(float(i))
    plt.figure()
    plt.hist(e,bins=100)
    plt.xlabel('Period')
    plt.title('Period')
    plt.yscale('log',nonposy='clip')
    plt.savefig(os.path.join("Correlation Plots/","Period"))

    e = []
    for i in columns['Semi-major Axis']:
        e.append(float(i))
    plt.figure()
    plt.hist(e,bins=100)
    plt.xlabel('Semi-major Axis')
    plt.ylabel('Counts')
    plt.title('Semimajor axis')
    plt.yscale('log',nonposy='clip')
    plt.savefig(os.path.join("Correlation Plots/","semimajor_axis"))

    e = []
    for i in columns['Albedo']:
        e.append(float(i))
    plt.figure()
    plt.hist(e,bins=100)
    plt.xlabel('Albedo')
    plt.ylabel('Counts')
    plt.title('Albedo')
    plt.yscale('log',nonposy='clip')
    plt.savefig(os.path.join("Correlation Plots/","albedo"))

    e = []
    for i in columns['Albedo at Semi-Major Axis']:
        e.append(float(i))
    plt.figure()
    plt.hist(e,bins=100)
    plt.xlabel('Albedo at Semi-Major Axis')
    plt.ylabel('Counts')
    plt.title('Albedo at semi-major axis')
    plt.yscale('log',nonposy='clip')
    plt.savefig(os.path.join("Correlation Plots/","smaxis_albedo"))


    # For the fluxes, the vast majority of points are very close to zero (but =/=0)

    e = []
    for i in columns['Flux']:
        e.append(float(i))
    plt.figure()
    plt.hist(e,bins=100)
    plt.xlabel('Flux')
    plt.ylabel('Counts')
    plt.title('Mean Flux')
    plt.yscale('log',nonposy='clip')
    plt.savefig(os.path.join("Correlation Plots/","flux"))

    e = []
    for i in columns['Peak Flux']:
        e.append(float(i))
    plt.figure()
    plt.hist(e,bins=100)
    plt.xlabel('Peak Flux')
    plt.ylabel('Counts')
    plt.title('Peak flux')
    plt.yscale('log',nonposy='clip')
    plt.savefig(os.path.join("Correlation Plots/","peak_flux"))



# Period, Semi-major Axis, Planetary Radius, Eccentricity, Argument of Periastron, Orbital Inclination, Albedo, Albedo at Semi-major Axis, Hilton Phase Function, Flux, Peak Flux


def scatterPlots():
    
    plt.figure()
    a = []
    b = []
    for i in range(len(columns['V-Magnitude'])):
        if columns['V-Magnitude'][i] == '30.0':
            continue
        else:
            a.append(columns['Peak Flux'][i])
            b.append(columns['V-Magnitude'][i])
    c = []
    d = []
    for i in range(len(topFluxes['V-Magnitude'])):
        if topFluxes['V-Magnitude'][i] == '30.0':
            continue
        else:
            c.append(topFluxes['Peak Flux'][i])
            d.append(topFluxes['V-Magnitude'][i])
    plt.semilogx(a,b,'.')
    plt.semilogx(c,d,'r.')
    plt.title('V-Magnitude v. Peak Flux')
    plt.xlabel('Peak Flux')
    plt.ylabel('V-Magnitude')
    plt.savefig(os.path.join("Correlation Plots/","V_PeakFlux"))


    plt.figure()
    plt.loglog(columns['Period'],columns['Albedo at Semi-Major Axis'],'.')
    plt.loglog(topFluxes['Period'],topFluxes['Albedo at Semi-Major Axis'],'r.')
    plt.title('Albedo at Semi-major axis v. Period')
    plt.xlabel('Period (days)')
    plt.ylabel('Albedo')
    plt.savefig(os.path.join("Correlation Plots/","Period_Albedo"))
    
    # Eliminates pre-assigned radius values
    colors = []
    flux = []
    sm = []
    for i in range(len(columns['Planetary Radius'])):
        colors.append(float(columns['Planetary Radius'][i]))
        flux.append(float(columns['Peak Flux'][i]))
        sm.append(float(columns['Semi-major Axis'][i]))
    fig = plt.figure()
    ax = fig.add_subplot(111)
    sc = ax.scatter(sm,flux,c=colors)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_ylim([1e-10,1e-3])
    ax.set_xlabel('Semi-major Axis (AU)')
    ax.set_ylabel('Peak Flux Ratio')
    ax.set_title('Peak Flux Ratio v. Semi-major Axis')
    plt.colorbar(sc)
    plt.savefig(os.path.join("Correlation Plots/","SMAxis_PeakFlux_heatmap"))

    plt.figure()
    plt.loglog(columns['Semi-major Axis'],columns['Peak Flux'],'.')
    plt.title('Peak Flux Ratio v. Semi-major Axis')
    plt.xlabel('Semi-major Axis (AU)')
    plt.ylabel('Peak Flux Ratio')
    plt.savefig(os.path.join("Correlation Plots/","SMAxis_PeakFlux"))
    
    plt.figure()
    a = []
    b = []
    for i in range(len(columns['Planetary Radius'])):
        if columns['Planetary Radius'][i] == '1.0':
            continue
        else:
            a.append(columns['Planetary Radius'][i])
            b.append(columns['Semi-major Axis'][i])
    c = []
    d = []
    for i in range(len(topFluxes['Planetary Radius'])):
        if columns['Planetary Radius'][i] == '1.0':
            continue
        else:
            c.append(topFluxes['Planetary Radius'][i])
            d.append(topFluxes['Semi-major Axis'][i])
    plt.loglog(a,b,'.')
    plt.loglog(c,d,'r.')
    plt.title('Semi-major Axis v. Planetary Radius')
    plt.xlabel('Radius (Jupiter-radii)')
    plt.ylabel('Semi-Major Axis (m)')
    plt.savefig(os.path.join("Correlation Plots/","Radius_SMAxis"))
    
    plt.figure()
    plt.plot(columns['Eccentricity'],columns['Planetary Radius'],'.')
    plt.title('Planetary Radius v. Eccentricity')
    plt.xlabel('Eccentricity')
    plt.ylabel('Radius (Jupiter-radii)')
    plt.savefig(os.path.join("Correlation Plots/","Eccentricity_Radius1"))

    plt.figure()
    a = []
    b = []
    for i in range(len(columns['Planetary Radius'])):
        if columns['Planetary Radius'][i] == '1.0':
            continue
        elif columns['Eccentricity'][i] == '0.0':
            continue
        else:
            a.append(columns['Eccentricity'][i])
            b.append(columns['Planetary Radius'][i])
    plt.plot(a,b,'.')
    plt.title('Planetary Radius v. Eccentricity')
    plt.xlabel('Eccentricity')
    plt.ylabel('Radius (Jupiter-radii)')
    plt.savefig(os.path.join("Correlation Plots/","Eccentricity_Radius2"))

    plt.figure()
    plt.loglog(columns['Period'],columns['Albedo'],'.')
    plt.title('Mean Albedo v. Period')
    plt.xlabel('Period (days)')
    plt.ylabel('Albedo')
    plt.savefig(os.path.join("Correlation Plots/","Period_Albedo2"))

    plt.figure()
    plt.loglog(columns['Semi-major Axis'],columns['Albedo at Semi-Major Axis'],'.')
    plt.loglog(topFluxes['Semi-major Axis'],topFluxes['Albedo at Semi-Major Axis'],'r.')
    plt.title('Albedo at Semi-Major Axis v. Semi-major Axis')
    plt.xlabel('Semi-major Axis (AU)')
    plt.ylabel('Albedo')
    plt.savefig(os.path.join("Correlation Plots/","SMAxis_Albedo"))

    plt.figure()
    a = []
    b = []
    for i in range(len(columns['Eccentricity'])):
        if columns['Eccentricity'][i] == '0.0':
            continue
        else:
            a.append(columns['Eccentricity'][i])
            b.append(columns['Albedo'][i])
    plt.loglog(a,b,'.')
    plt.title('Albedo v. Eccentricity')
    plt.xlabel('Eccentricity')
    plt.ylabel('Albedo')
    plt.savefig(os.path.join("Correlation Plots/","Eccentricity_Albedo"))

    plt.figure()
    a = []
    b = []
    for i in range(len(columns['Eccentricity'])):
        if columns['Eccentricity'][i] == '0.0':
            continue
        else:
            a.append(columns['Eccentricity'][i])
            b.append(columns['Peak Flux'][i])
    plt.loglog(a,b,'.')
    plt.title('Peak Flux v. Eccentricity')
    plt.xlabel('Eccentricity')
    plt.ylabel('Flux')
    plt.savefig(os.path.join("Correlation Plots/","Eccentricity_PeakFlux"))

    plt.figure()
    a = []
    b = []
    for i in range(len(columns['Eccentricity'])):
        if columns['Eccentricity'][i] == '0.0':
            continue
        else:
            a.append(columns['Eccentricity'][i])
            b.append(columns['Albedo at Semi-Major Axis'][i])
    plt.loglog(columns['Eccentricity'],columns['Albedo at Semi-Major Axis'],'.')
    plt.title('Albedo at Semi-Major Axis v. Eccentricity')
    plt.xlabel('Eccentricity')
    plt.ylabel('Albedo')
    plt.savefig(os.path.join("Correlation Plots/","Eccentricity_SMAxis"))
    
    plt.figure()
    a = []
    b = []
    for i in range(len(columns['V-Magnitude'])):
        if columns['V-Magnitude'][i] == '30.0':
            continue
        else:
            a.append(columns['Peak Flux'][i])
            b.append(columns['V-Magnitude'][i])
    c = []
    d = []
    for i in range(len(topFluxes['V-Magnitude'])):
        if topFluxes['V-Magnitude'][i] == '30.0':
            continue
        else:
            c.append(topFluxes['Peak Flux'][i])
            d.append(topFluxes['V-Magnitude'][i])
    plt.semilogx(a,b,'.')
    plt.semilogx(c,d,'r.')
    plt.title('V-Magnitude v. Peak Flux')
    plt.xlabel('Peak Flux')
    plt.ylabel('V-Magnitude')
    plt.savefig(os.path.join("Correlation Plots/","V_PeakFlux"))


def scatter3d():
    # 3-D scatter plot
    xs = [float(i) for i in columns['Semi-major Axis']]
    ys = [float(i) for i in columns['Albedo']]
    zs = [float(i) for i in columns['Period']]

    fig = plt.figure()
    ax = fig.add_subplot(111,projection='3d')
    ax.scatter(xs,ys,zs)
    ax.set_xlabel('Semi-major axis')
    ax.set_ylabel('Albedo')
    ax.set_zlabel('Period')




histograms()
scatterPlots()
scatter3d()

