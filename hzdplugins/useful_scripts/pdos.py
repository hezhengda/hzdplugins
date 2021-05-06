import re
import os
import sys
import numpy as np
import json

list_file = os.listdir('./')

if sys.argv[1] == 'spin':
    print("Your simulation is spin-polarized.")
elif sys.argv[1] == 'nspin':
    print("Your simulation is non-spin-polarized.")

p = re.compile('.+#([0-9]+)\(([a-zA-Z]+)\)_wfc#[0-9]\(([a-z])\)')

results = {}

json_file = open('pdos_results.json', 'w')

for file in list_file:
    m = p.match(file)
    if m is None:
        pass
    else:
        
        m1 = m.group(1)
        m2 = m.group(2)
        m3 = m.group(3)

        if m2 not in results:
            results[m2] = {}
        if m1 not in results[m2]:
            results[m2][m1] = {}
        if m3 not in results[m2][m1]:
            results[m2][m1][m3]={}

        if sys.argv[1] == 'spin': 
            if m3 == 's':
                energy = np.loadtxt(file, usecols=0)
                results[m2][m1][m3]['energy'] = np.ndarray.tolist(energy)

                ldosup = np.loadtxt(file, usecols=1)
                results[m2][m1][m3]['ldosup'] = np.ndarray.tolist(ldosup)

                ldosdw = np.loadtxt(file, usecols=2)
                results[m2][m1][m3]['ldosdw'] = np.ndarray.tolist(ldosdw)

                pdosup = np.loadtxt(file, usecols=3)
                results[m2][m1][m3]['pdosup'] = np.ndarray.tolist(pdosup)

                pdosdw = np.loadtxt(file, usecols=4)
                results[m2][m1][m3]['pdosdw'] = np.ndarray.tolist(pdosdw)

            if m3 == 'p':
                energy = np.loadtxt(file, usecols=0)
                results[m2][m1][m3]['energy'] = np.ndarray.tolist(energy)

                ldosup = np.loadtxt(file, usecols=1)
                results[m2][m1][m3]['ldosup'] = np.ndarray.tolist(ldosup)

                ldosdw = np.loadtxt(file, usecols=2)
                results[m2][m1][m3]['ldosdw'] = np.ndarray.tolist(ldosdw)

                p1up = np.loadtxt(file, usecols=3)
                results[m2][m1][m3]['p1up'] = np.ndarray.tolist(p1up)

                p1dw = np.loadtxt(file, usecols=4)
                results[m2][m1][m3]['p1dw'] = np.ndarray.tolist(p1dw)

                p2up = np.loadtxt(file, usecols=5)
                results[m2][m1][m3]['p2up'] = np.ndarray.tolist(p2up)

                p2dw = np.loadtxt(file, usecols=6)
                results[m2][m1][m3]['p2dw'] = np.ndarray.tolist(p2dw)

                p3up = np.loadtxt(file, usecols=7)
                results[m2][m1][m3]['p3up'] = np.ndarray.tolist(p3up)

                p3dw = np.loadtxt(file, usecols=8)
                results[m2][m1][m3]['p3dw'] = np.ndarray.tolist(p3dw)

                sumpup = p1up + p2up + p3up
                sumpdw = p1dw + p2dw + p3dw

                results[m2][m1][m3]['ptotup'] = np.ndarray.tolist(sumpup)
                results[m2][m1][m3]['ptotdw'] = np.ndarray.tolist(sumpdw)

            if m3 == 'd':
                energy = np.loadtxt(file, usecols=0)
                results[m2][m1][m3]['energy'] = np.ndarray.tolist(energy)

                ldosup = np.loadtxt(file, usecols=1)
                results[m2][m1][m3]['ldosup'] = np.ndarray.tolist(ldosup)

                ldosdw = np.loadtxt(file, usecols=2)
                results[m2][m1][m3]['ldosdw'] = np.ndarray.tolist(ldosdw)

                d1up = np.loadtxt(file, usecols=3)
                results[m2][m1][m3]['d1up'] = np.ndarray.tolist(d1up)

                d1dw = np.loadtxt(file, usecols=4)
                results[m2][m1][m3]['d1dw'] = np.ndarray.tolist(d1dw)

                d2up = np.loadtxt(file, usecols=5)
                results[m2][m1][m3]['d2up'] = np.ndarray.tolist(d2up)

                d2dw = np.loadtxt(file, usecols=6)
                results[m2][m1][m3]['d2dw'] = np.ndarray.tolist(d2dw)

                d3up = np.loadtxt(file, usecols=7)
                results[m2][m1][m3]['d3up'] = np.ndarray.tolist(d3up)

                d3dw = np.loadtxt(file, usecols=8)
                results[m2][m1][m3]['d3dw'] = np.ndarray.tolist(d3dw)

                d4up = np.loadtxt(file, usecols=5)
                results[m2][m1][m3]['d4up'] = np.ndarray.tolist(d4up)

                d4dw = np.loadtxt(file, usecols=6)
                results[m2][m1][m3]['d4dw'] = np.ndarray.tolist(d4dw)

                d5up = np.loadtxt(file, usecols=7)
                results[m2][m1][m3]['d5up'] = np.ndarray.tolist(d5up)

                d5dw = np.loadtxt(file, usecols=8)
                results[m2][m1][m3]['d5dw'] = np.ndarray.tolist(d5dw)

                sumdup = d1up + d2up + d3up + d4up + d5up
                sumddw = d1dw + d2dw + d3dw + d4dw + d5dw

                results[m2][m1][m3]['dtotup'] = np.ndarray.tolist(sumdup)
                results[m2][m1][m3]['dtotdw'] = np.ndarray.tolist(sumddw)

json.dump(results, json_file) 

