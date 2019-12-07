import subprocess
commands = ['rm *.dat','rm fort.*',
'ifort dvode_f90_m.f90 lyap.f90 -o lyap',
'./lyap',
'cp fort.1 period.dat',
'cp fort.2 lag.dat',
'cp fort.11 max_y1.dat',
'cp fort.21 min_y1.dat',
'cp fort.12 max_y2.dat',
'cp fort.22 min_y2.dat',
'cp fort.13 max_y3.dat',
'cp fort.23 min_y3.dat',
'cp fort.14 max_y4.dat',
'cp fort.24 min_y4.dat',
'cp fort.15 max_y5.dat',
'cp fort.25 min_y5.dat',
'cp fort.16 max_y6.dat',
'cp fort.26 min_y6.dat',
'cp fort.17 max_y7.dat',
'cp fort.27 min_y7.dat',
'cp fort.32 lyap_exp.dat',
'rm lyap',
'rm *.mod',
'rm fort.*']
for command in commands:
    subprocess.call((command), shell=True)
#subprocess.call((commands[2]), shell=True)
