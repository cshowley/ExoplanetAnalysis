# Plot phase variations for all confirmed exoplanets' systems
# Works by pulling a spreadsheet of all the raw data from the internet then manipulates and reorganizes it, outputting a streamlined spreadsheet and a plot for each system (947 plots)

import numpy as np
import matplotlib.pyplot as plt
import math
import csv
from collections import defaultdict
from itertools import chain
import os
import re
from sortedcontainers import SortedDict
import wget

# Define the functions required for flux ratio calculation

def calculateStartTime(T,e,w,stepSize):     # Find t_o (starting point) for Kepler equation solution
    x = int(1)
    if w > 270.0:
        w = w - 360.0
    w = math.radians(w)
    findStartDate = []
    if T < 1:
        T = 1.0
    for i in range(10**x * int(T)):
        M = (2 * np.pi / T) * (stepSize + 10**(-x) * i)
        M = M%(2*np.pi) - 2*np.pi
        Ei = np.pi
        for j in range(0,100): # Arbitrary number of steps chosen; it should find a solution for E after only a few iterations
            E = Ei - (Ei - e*np.sin(Ei) - M)/(1 - e*np.cos(Ei))
            if np.abs((E-Ei)/Ei) < 1e-4:
                break
            else:
                Ei = E
        cos_f = (np.cos(E) - e) / (1.0 - e * np.cos(E))
        sin_f = (np.sqrt(1.0 - e**2) * (np.sin(E) / (1.0 - e * np.cos(E))))
        if sin_f < 0:
            f = 2*np.pi - np.arccos(cos_f)
            findStartDate.append(abs(math.radians(270.0)-(w+f)))
        else:
            f = np.arccos(cos_f)
            findStartDate.append(abs(math.radians(270.0)-(w+f)))
    for i in range(len(findStartDate)):
        if findStartDate[i] == min(findStartDate):
            t_o = float(i) / len(findStartDate) * T
            break
    return t_o

def kepler(T,t_o,e,n,w,stepSize):      # Solve for true anomaly
    true_anomaly = []
    for i in range(0,n):
        M = (2 * np.pi / T) * ((i) * stepSize + t_o)
        M = M%(2*np.pi) - 2*np.pi
        E = 1.5  # Initial guess of pi doesn't work but 1.5 has empirically been determined to work for all but three systems; reason unknown
        for j in range(0,100):
            Ei = E
            E = Ei - (Ei - e*np.sin(Ei) - M)/(1.0 - e*np.cos(Ei))
            if np.abs((E-Ei)/Ei) < 1e-4:
                break
            else:
                continue
        cos_f = (np.cos(E) - e) / (1.0 - e * np.cos(E))
        sin_f = (np.sqrt(1.0 - e**2) * (np.sin(E) / (1.0 - e * np.cos(E))))
        if sin_f < 0:
            f = 2*np.pi - np.arccos(cos_f)
            true_anomaly.append(f)
        else:
            f = np.arccos(cos_f)
            true_anomaly.append(f)
    return true_anomaly

def phaseVariation(days,n,T,t_o,R_E,R,e,w,a,I,stepSize,name):
    T = np.array(T)
    R = R_E * np.array(R) # Convert R to km
    e = np.array(e)
    w = np.array(w)
    a = np.array(a)
    I = np.array(I)
    fluxRatio = np.zeros(shape=(len(R),n)) # Create matrix for planets' fluxes
    A = []      # Albedo
    smA = []    # Albedo at semi-major axis
    flux = []   # Flux (of individual planets)
    sps = []    # Star-planet separation
    timeRange = []
    for i in range(n):
        s = stepSize*i
        timeRange.append(s)
    # Calculate flux ratio
    for i in range(len(R)):
        f = np.array(kepler(T[i],t_o,e[i],n,w[i],stepSize)) # True anomaly
        f = f.tolist()
        alpha = []
        for j in range(len(f)):
            angle = np.arccos(np.sin(np.radians(w[i]) + f[j] + np.pi) * np.sin(np.radians(I[i])))
            alpha.append(angle)
        r = a[i] * ( 1 - (e[i]) ** 2) / (1 + e[i] * np.cos(f)) # Star-planet separation
        sps.append(r)
        albedo = (np.exp(r-1) - np.exp(1-r))/(5 * (np.exp(r-1) + np.exp(1-r))) + 3.0/10 # Calculate geometric albedo
        smA.append((np.exp(a[i]-1) - np.exp(1-a[i]))/(5 * (np.exp(a[i]-1) + np.exp(1-a[i]))) + 3.0/10)
        A.append(albedo)
        ma = .09 * (alpha/np.radians(100)) + 2.39 * (alpha/np.radians(100)) ** 2 - .65 * (alpha/np.radians(100)) ** 3 # Hilton phase function
        g = 10.0 ** (-.4*ma)
        fluxRatio[i] = albedo * g * (R[i] / (r*1.496e8))**2 # Flux ratio
        flux.append(fluxRatio[i])
    totalFlux = np.zeros(n)
    for i in range(len(R)):
        totalFlux += fluxRatio[i]
    plt.figure()
    plt.plot(timeRange/max(T),totalFlux,'k',markersize=2)
    plt.ticklabel_format(style = 'sci',axis='y',scilimits=(0,0))
    plt.xlim([0,timeRange[-1]/max(T)])
    plt.xlabel('Orbital phase of outermost planet')
    plt.ylabel('Flux ratio')
    plt.title(name+', '+str(len(R))+' planet(s)')
    return A,smA,flux,totalFlux,sps

#  Pull constants from Kepler data

# Read in the EOD column: if EOD == 1, add it to dictionary. If EOD == 0 (the candidates), skip it
# If no value for inclination, use I == 90
# Read in a column for planet mass
# If no value for planetary radius, use equation sent in Kane's email: if planet mass > .3, use R = 1 jupiter radius. If planet mass < .3, use the power law relation
# If no eccentricity, use e == 0
# If no argument of periastron, use w == 90
# If no period AND semi-major axis, ignore

exoDict = {}
columns = defaultdict(list)
if not os.path.exists('exoplanets.csv'):
    url = 'http://exoplanets.org/csv-files/exoplanets.csv'
    filename = wget.download(url)
with open("exoplanets.csv","rU") as data:
    reader = csv.DictReader(data,delimiter=',')
    for row in reader:
        for (h,v) in row.items():
            columns[h].append(v)
i = 0
for name in columns['NAME']:
    if columns['EOD'][columns['NAME'].index(name)] == '1':
        #if not re.search('KOI',name):
        planetName = name[-1]
        name = name[0:-2]
        if exoDict.has_key(name):
            if columns['I'][i] == '':
                columns['I'][i] = 90.0
            if columns['R'][i] == '':
                if columns['MASS'][i] == '':
                    columns['R'][i] = 0.001
                elif columns['MASS'][i] >= '0.3':
                    columns['R'][i] = 1.0
                else:
                    columns['R'][i] = 1.91286 * float(columns['MASS'][i]) ** 0.513388
            if columns['ECC'][i] == '':
                columns['ECC'][i] = 0.0
            if columns['OM'][i] == '':
                columns['OM'][i] = 90.0
            if columns['MSTAR'][i] == '':
                columns['MSTAR'][i] = 0.0
            if columns['V'][i] == '':
                columns['V'][i] = 30.0
            if columns['A'][i] == '' and columns['PER'][i] == '':
                i = i + 1
                continue
            exoDict[name]['Planet'].append(planetName)
            exoDict[name]['Period'].append(columns['PER'][i])
            exoDict[name]['Semi-major Axis'].append(columns['A'][i])
            exoDict[name]['Planetary Radius'].append(columns['R'][i])
            exoDict[name]['Eccentricity'].append(columns['ECC'][i])
            exoDict[name]['Argument of Periastron'].append(columns['OM'][i])
            exoDict[name]['Orbital Inclination'].append(columns['I'][i])
            exoDict[name]['Mass'].append(columns['MASS'][i])
            exoDict[name]['Star Mass'].append(columns['MSTAR'][i])
            exoDict[name]['V-Mag'].append(columns['V'][i])
            i = i + 1
        else:
            if columns['I'][i] == '':
                columns['I'][i] = 90.0
            if columns['R'][i] == '':
                if columns['MASS'][i] == '':
                    columns['R'][i] = 0.001
                elif columns['MASS'][i] >= '0.3':
                    columns['R'][i] = 1.0
                else:
                    columns['R'][i] = 1.91286 * float(columns['MASS'][i]) ** 0.513388
            if columns['ECC'][i] == '':
                columns['ECC'][i] = 0.0
            if columns['OM'][i] == '':
                columns['OM'][i] = 90.0
            if columns['MSTAR'][i] == '':
                columns['MSTAR'][i] = 0.0
            if columns['V'][i] == '':
                columns['V'][i] = 30.0
            if columns['A'][i] == '' and columns['PER'][i] == '':
                i = i + 1
                continue
            exoDict[name] = {'Planet':[planetName],'Period':[columns['PER'][i]],'Semi-major Axis':[columns['A'][i]],'Star-Planet Separation':[],'Planetary Radius':[columns['R'][i]],'Eccentricity':[columns['ECC'][i]],'Argument of Periastron':[columns['OM'][i]],'Orbital Inclination':[columns['I'][i]],'Mass':[columns['MASS'][i]],'Star Mass':[columns['MSTAR'][i]],'V-Mag':[columns['V'][i]],'Mean Albedo':[],'Semi-major Albedo':[],'Mean Flux':[],'Peak Flux':[],'Total Flux':[]}
            i = i + 1
    else:
        i = i + 1
    """ # This section contains the Kepler candidate planets; for now they can be ignored
        else:
        if exoDict.has_key(name):
        exoDict[name]['Planet'].append('N/A')
        exoDict[name]['Period'].append(columns['PER'][i])
        exoDict[name]['Semi-major Axis'].append(columns['A'][i])
        exoDict[name]['Planetary Radius'].append(columns['R'][i])
        exoDict[name]['Eccentricity'].append(columns['ECC'][i])
        exoDict[name]['Argument of Periastron'].append(columns['OM'][i])
        exoDict[name]['Orbital Inclination'].append(columns['I'][i])
        i = i + 1
        else:
        exoDict[name] = {'Planet':['N/A'],'Period':[columns['PER'][i]],'Semi-major Axis':[columns['A'][i]],'Star-Planet Separation':[],'Planetary Radius':[columns['R'][i]],'Eccentricity':[columns['ECC'][i]],'Argument of Periastron':[columns['OM'][i]],'Orbital Inclination':[columns['I'][i]],'Mass':[columns['MASS'][i]],'Star Mass':[columns['MSTAR'][i]],'V-Mag':[columns['V'][i]],'Mean Albedo':[],'Semi-major Albedo':[],'Mean Flux':[],'Peak Flux':[],'Total Flux':[]}
        i = i + 1
        """

# Run code to create phase curves for Kepler systems

if not os.path.exists('Plots'):
    os.makedirs('Plots')
for name in exoDict:
    print name
    T = []
    R = []
    e = []
    w = []
    I = []
    a = []
    for j in range(len(exoDict[name]['Planet'])):
        T.append(float(exoDict[name]['Period'][j]))
        R.append(float(exoDict[name]['Planetary Radius'][j]))
        e.append(float(exoDict[name]['Eccentricity'][j]))
        w.append(float(exoDict[name]['Argument of Periastron'][j]))
        I.append(float(exoDict[name]['Orbital Inclination'][j]))
        a.append(float(exoDict[name]['Semi-major Axis'][j]))
    n = 10000
    R_E = 71492.0                   # Jupiter-radius (km)
    maxPeriod = T.index(max(T))     # The number of elements into the array the maximum period occurs; use to calculate the start date
    days = max(T)
    stepSize = days/n               # Dividing outermost planet's period into n equal parts
    t_o = calculateStartTime(T[maxPeriod],e[maxPeriod],w[maxPeriod],stepSize)
    A,smA,flux,totalFlux,sps = phaseVariation(days,n,T,t_o,R_E,R,e,w,a,I,stepSize,name)
    for i in range(len(exoDict[name]['Planet'])):
        exoDict[name]['Mean Albedo'].append(np.mean(A[i]))
        exoDict[name]['Semi-major Albedo'].append(smA[i])
        exoDict[name]['Mean Flux'].append(np.mean(flux[i]))
        exoDict[name]['Peak Flux'].append(max(flux[i]))
        exoDict[name]['Total Flux'].append(max(totalFlux))
        exoDict[name]['Star-Planet Separation'].append(min(sps[i]))
    plt.savefig(os.path.join("Plots/",name))
    plt.close()
exoDict = SortedDict(exoDict)
with open("planetData.csv","w") as data:
    data.write('%s' % 'System' +","+ '%s' % 'Planet' +","+ '%s' % 'Period' +","+ '%s' % 'Semi-major Axis' +","+ '%s' % 'Minimum Star-Planet Separation' +","+ '%s' % 'Planetary Radius' +","+ '%s' % 'Eccentricity' +","+ '%s' % 'Argument of Periastron' +","+ '%s' % 'Orbital Inclination' +","+ '%s' % 'Albedo' +","+ '%s' % 'Albedo at Semi-Major Axis' +","+ '%s' % 'Flux' +","+ '%s' % 'Peak Flux' +","+ '%s' % 'Total Flux' +","+ '%s' % 'Star Mass' +","+ '%s\n' % 'V-Magnitude')
for names in exoDict:
    for i in range(len(exoDict[names]['Planet'])):
        with open("planetData.csv","a") as data:
            data.write('%s' % names +","+ '%s' % exoDict[names]['Planet'][i] +","+ '%s' % exoDict[names]['Period'][i] +","+ '%s' % exoDict[names]['Semi-major Axis'][i] +","+ '%s' % exoDict[names]['Star-Planet Separation'][i] +","+ '%s' % exoDict[names]['Planetary Radius'][i] +","+ '%s' % exoDict[names]['Eccentricity'][i] +","+ '%s' % exoDict[names]['Argument of Periastron'][i] +","+ '%s' % exoDict[names]['Orbital Inclination'][i] +","+ '%s' % exoDict[names]['Mean Albedo'][i] +","+ '%s' % exoDict[names]['Semi-major Albedo'][i] +","+ '%s' % exoDict[names]['Mean Flux'][i] +","+ '%s' % exoDict[names]['Peak Flux'][i] +","+ '%s' % exoDict[names]['Total Flux'][i] +","+ '%s' % exoDict[names]['Star Mass'][i] +","+ '%s\n' % exoDict[names]['V-Mag'][i])