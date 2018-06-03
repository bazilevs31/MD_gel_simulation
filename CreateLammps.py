#!/usr/bin/env python

# Program: CreateLammps.py
# Purpose: creates Lammps input files
# Author:  Triandafilidi Vasiliy , MSc student at CHBE UBC, Vancouver
# e-mail:  vtriandafilidi(at)chbe(dot)ubc(dot)ca
# Syntax:  python CreateLammps.py for help,
# example: python CreateLammps.py
# Requires:

# Copyright (c) 2014 Vasiliy Triandafilidi
# Released under the GNU Public Licence, v2 or any higher version

import numpy as np
import os
import read_parameters

def CreateTable():
    """creates cgpva.table file"""
    text = """\
#Angle potential for CG-PVA coarse grain model
#Col-1: index
#Col-2: angle  (degree)
#Col-3: energy (LAMMPS unit)
#Col-4: force  (LAMMPS unit/degree)

CG_PVA
N 181 FP 0 0

1    0.00000  200.000000000000000    0.748411898253894
2    1.00000  198.000000000000000    0.748411898253894
3    2.00000  195.500000000000000    0.748411898253894
4    3.00000  193.000000000000000    0.748411898253894
5    4.00000  190.500000000000000    0.748411898253894
6    5.00000  188.000000000000000    0.748411898253894
7    6.00000  185.500000000000000    0.748411898253894
8    7.00000  183.000000000000000    0.748411898253894
9    8.00000  180.500000000000000    0.748411898253894
10    9.00000  178.000000000000000    0.748411898253894
11   10.00000  175.500000000000000    0.748411898253894
12   11.00000  173.000000000000000    0.748411898253894
13   12.00000  170.500000000000000    0.748411898253894
14   13.00000  168.000000000000000    0.748411898253894
15   14.00000  165.500000000000000    0.748411898253894
16   15.00000  163.000000000000000    0.748411898253894
17   16.00000  160.500000000000000    0.748411898253894
18   17.00000  158.000000000000000    0.748411898253894
19   18.00000  155.500000000000000    0.748411898253894
20   19.00000  153.000000000000000    0.748411898253894
21   20.00000  150.500000000000000    0.748411898253894
22   21.00000  148.000000000000000    0.748411898253894
23   22.00000  145.500000000000000    0.748411898253894
24   23.00000  143.000000000000000    0.748411898253894
25   24.00000  140.500000000000000    0.748411898253894
26   25.00000  138.000000000000000    0.748411898253894
27   26.00000  135.500000000000000    0.748411898253894
28   27.00000  133.000000000000000    0.748411898253894
29   28.00000  130.500000000000000    0.748411898253894
30   29.00000  128.000000000000000    0.748411898253894
31   30.00000  125.500000000000000    0.748411898253894
32   31.00000  123.000000000000000    0.748411898253894
33   32.00000  120.500000000000000    0.748411898253894
34   33.00000  118.000000000000000    0.748411898253894
35   34.00000  115.500000000000000    0.748411898253894
36   35.00000  113.000000000000000    0.748411898253894
37   36.00000  110.500000000000000    0.748411898253894
38   37.00000  108.000000000000000    0.748411898253894
39   38.00000  105.500000000000000    0.748411898253894
40   39.00000  103.000000000000000    0.748411898253894
41   40.00000  100.500000000000000    0.748411898253894
42   41.00000   98.000000000000000    0.748411898253894
43   42.00000   95.500000000000000    0.748411898253894
44   43.00000   93.000000000000000    0.748411898253894
45   44.00000   90.500000000000000    0.748411898253894
46   45.00000   88.000000000000000    0.748411898253894
47   46.00000   85.500000000000000    0.748411898253894
48   47.00000   83.000000000000000    0.748411898253894
49   48.00000   80.500000000000000    0.748411898253894
50   49.00000   78.000000000000000    0.748411898253894
51   50.00000   75.500000000000000    0.748411898253894
52   51.00000   73.000000000000000    0.748411898253894
53   52.00000   70.500000000000000    0.748411898253894
54   53.00000   68.000000000000000    0.748411898253894
55   54.00000   65.500000000000000    0.748411898253894
56   55.00000   63.000000000000000    0.748411898253894
57   56.00000   60.500000000000000    0.748411898253894
58   57.00000   58.000000000000000    0.748411898253894
59   58.00000   55.500000000000000    0.748411898253894
60   59.00000   53.000000000000000    0.748411898253894
61   60.00000   50.500000000000000    0.748411898253894
62   61.00000   48.000000000000000    0.748411898253894
63   62.00000   45.500000000000000    0.748411898253894
64   63.00000   43.000000000000000    0.748411898253894
65   64.00000   40.500000000000000    0.748411898253894
66   65.00000   38.000000000000000    0.748411898253894
67   66.00000   35.500000000000000    0.748411898253894
68   67.00000   33.000000000000000    0.748411898253894
69   68.00000   30.500000000000000    0.748411898253894
70   69.00000   28.000000000000000    0.748411898253894
71   70.00000   25.500000000000000    0.748411898253894
72   71.00000   23.000000000000000    0.748411898253894
73   72.00000   20.500000000000000    0.748411898253894
74   73.00000   18.000000000000000    0.748411898253894
75   74.00000   15.500000000000000    0.748411898253894
76   75.00000   13.000000000000000    0.748411898253894
77   76.00000   10.944100000000001    0.748411898253894
78   77.00000    9.914809999999999    0.748411898253894
79   78.00000    9.146540000000000    0.748411898253894
80   79.00000    8.418080000000000    0.689934941528673
81   80.00000    7.766780000000000    0.593418872078838
82   81.00000    7.231370000000000    0.490805949169356
83   82.00000    6.785300000000000    0.420907037640622
84   83.00000    6.389700000000000    0.374501661702760
85   84.00000    6.036450000000000    0.323318447123692
86   85.00000    5.743220000000000    0.267716370218545
87   86.00000    5.501190000000000    0.236263339520293
88   87.00000    5.270870000000000    0.231830323209402
89   88.00000    5.037700000000000    0.225814912355779
90   89.00000    4.819410000000000    0.205821209851153
91   90.00000    4.626230000000000    0.180259921899849
92   91.00000    4.459070000000000    0.155173799322168
93   92.00000    4.316070000000000    0.132824498373964
94   93.00000    4.193600000000000    0.113196309253767
95   94.00000    4.089860000000000    0.093035712222347
96   95.00000    4.007720000000000    0.069904727165061
97   96.00000    3.950240000000000    0.044143041589219
98   97.00000    3.919630000000000    0.016145492264417
99   98.00000    3.918150000000000   -0.012976646218291
100   99.00000    3.945780000000000   -0.041948532741234
101  100.00000    4.002250000000000   -0.068037425912495
102  101.00000    4.082070000000000   -0.088050176311657
103  102.00000    4.178570000000000   -0.103343240931075
104  103.00000    4.288980000000000   -0.115848258138394
105  104.00000    4.410480000000000   -0.124016544251576
106  105.00000    4.537230000000000   -0.125316276719848
107  106.00000    4.661340000000000   -0.121711846157289
108  107.00000    4.780870000000000   -0.116428137366422
109  108.00000    4.894410000000000   -0.105063418139290
110  109.00000    4.991220000000000   -0.081636493853637
111  110.00000    5.057900000000000   -0.050740145697482
112  111.00000    5.092910000000000   -0.018679050114887
113  112.00000    5.095460000000000    0.023956384642703
114  113.00000    5.045190000000000    0.084211145942055
115  114.00000    4.927220000000000    0.148633362365101
116  115.00000    4.748090000000000    0.204093192752710
117  116.00000    4.519190000000000    0.254943964619320
118  117.00000    4.238360000000000    0.297942318225201
119  118.00000    3.923460000000000    0.315123795860061
120  119.00000    3.608250000000000    0.305968741909619
121  120.00000    3.311660000000000    0.282461058579403
122  121.00000    3.043480000000000    0.251911191863801
123  122.00000    2.808000000000000    0.215012064506227
124  123.00000    2.613630000000000    0.172503674404277
125  124.00000    2.463180000000000    0.129236250833180
126  125.00000    2.355350000000000    0.087074073563153
127  126.00000    2.289240000000000    0.046587371312115
128  127.00000    2.262400000000000    0.007687957559803
129  128.00000    2.274100000000000   -0.029440205547854
130  129.00000    2.321520000000000   -0.063995992631555
131  130.00000    2.402340000000000   -0.094567404700382
132  131.00000    2.510920000000000   -0.118691644526227
133  132.00000    2.640000000000000   -0.137810754984722
134  133.00000    2.786830000000000   -0.154510354982402
135  134.00000    2.949320000000000   -0.170331369351336
136  135.00000    3.127800000000000   -0.186334314372135
137  136.00000    3.322300000000000   -0.199965216165372
138  137.00000    3.528060000000000   -0.206310944322056
139  138.00000    3.735260000000000   -0.205801233394001
140  139.00000    3.940010000000000   -0.203592059361818
141  140.00000    4.142800000000000   -0.204218429040293
142  141.00000    4.348810000000000   -0.212950103105970
143  142.00000    4.569090000000000   -0.230450095976196
144  143.00000    4.810110000000000   -0.252695646776890
145  144.00000    5.074900000000000   -0.268238420032804
146  145.00000    5.347030000000000   -0.255631534779026
147  146.00000    5.586610000000000   -0.213358733491144
148  147.00000    5.774210000000000   -0.161911155927266
149  148.00000    5.910910000000000   -0.102560511282532
150  149.00000    5.979820000000000   -0.026484033298657
151  150.00000    5.964370000000000    0.059520944001749
152  151.00000    5.861280000000000    0.138934403972244
153  152.00000    5.687030000000000    0.195542067033090
154  153.00000    5.470740000000000    0.227554436835173
155  154.00000    5.232500000000000    0.249985577932544
156  155.00000    4.971390000000000    0.270308838222814
157  156.00000    4.692540000000000    0.280302593626241
158  157.00000    4.411500000000000    0.280183278965039
159  158.00000    4.132960000000000    0.284160304839920
160  159.00000    3.844040000000000    0.295264853759676
161  160.00000    3.543370000000000    0.300803694103935
162  161.00000    3.243460000000000    0.296905950849982
163  162.00000    2.950700000000000    0.290910006469021
164  163.00000    2.662930000000000    0.286648223731388
165  164.00000    2.378850000000000    0.282645724198492
166  165.00000    2.099270000000000    0.277745075536802
167  166.00000    1.825230000000000    0.271179061346769
168  167.00000    1.559080000000000    0.262424585842835
169  168.00000    1.302910000000000    0.249936838503170
170  169.00000    1.062180000000000    0.232935119979385
171  170.00000    0.840613000000000    0.211275428744556
172  171.00000    0.644005000000000    0.184970629587489
173  172.00000    0.476140000000000    0.156906211931453
174  173.00000    0.337224000000000    0.128283193039742
175  174.00000    0.228953000000000    0.098249593811098
176  175.00000    0.153872000000000    0.068082128021545
177  176.00000    0.112585000000000    0.064465125573232
178  177.00000    0.058306100000000    0.034906585039887
179  178.00000    0.020000000000000    0.017453292519943
180  179.00000    0.010000000000000    0.000000000000000
181  180.00000    0.000000000000000    0.000000000000000 """
    with open('cgpva.table', 'w') as table:
        table.write(text+'\n')

def CreateKremer():
    """creates inkremer.txt"""

    NameOfSim = "kremer"
    dcdfile = NameOfSim
    lammpsscript = "in"+dcdfile+".txt"
    filedumpskip = "dumpskip.txt"
    tmpdata = "tmp." + NameOfSim
    indata = "init"
    outdata = "equil"
    runpart = 1000000
    dumptraj = runpart/200
    dumpthermo = runpart/5000
    timestep = 0.01


    with open(filedumpskip, 'w') as f:
        f.write("log"+NameOfSim+"skip "+str(dumpthermo)+'\n')
        f.write("dump"+NameOfSim+"skip "+str(dumptraj)+'\n')
        f.write("time"+NameOfSim+"step "+str(timestep)+'\n')

    parameters = {'dcdfile':dcdfile,
                'tmpdata':tmpdata,
                'indata':indata,
                'outdata':outdata,
                'runpart':int(runpart),
                'dumptraj':int(dumptraj),
                'dumpthermo':int(dumpthermo),
                "runnvehot": runpart*1,
                'runnvehotmed' : 2*runpart,
                'runnvemed': 2*runpart,
                'runnvemednorm': 2*runpart,
                'runnvenorm': 2*runpart,
                'runnptrest': 5*runpart,
                'runquenchsim': 20*runpart,
                'kremerequil': 2*runpart,
                'timestep': timestep
    }


    script = """\

    # Kremer-Grest model.

    variable        Th equal 1.0   #-> hot

    units lj
    atom_style bond

    special_bonds lj/coul 0 1 1

    read_data {0[indata]}.data

    neighbor 0.4 bin
    neigh_modify every 1 delay 1
    comm_modify vel yes


    bond_style fene
    bond_coeff * 30.0 1.5 1.0 1.0

    dump            mydump all dcd {0[dumptraj]} {0[dcdfile]}.dcd
    timestep {0[timestep]}
    thermo {0[dumpthermo]}
    thermo_modify norm no




    pair_style dpd 1.0 1.0 122347
    pair_coeff * * 25 4.5 1.0


    velocity all create ${{Th}} 17786140

    fix 1 all nve/limit 0.001
    run 500
    fix 1 all nve/limit 0.01
    run 500
    fix 1 all nve/limit 0.05
    run 500
    fix 1 all nve/limit 0.1
    run 500
    unfix 1
    fix 1 all nve
    velocity all create ${{Th}} 17786140
    run 50000


    write_data {0[tmpdata]}.data

    pair_coeff * * 50.0 4.5 1.0
    velocity all create ${{Th}} 15086120
    run 50
    pair_coeff * * 100.0 4.5 1.0
    velocity all create ${{Th}} 15786120
    run 50
    pair_coeff * * 150.0 4.5 1.0
    velocity all create ${{Th}} 15486120
    run 50
    pair_coeff * * 200.0 4.5 1.0
    velocity all create ${{Th}} 17986120
    run 100
    pair_coeff * * 250.0 4.5 1.0
    velocity all create ${{Th}} 15006120
    run 100
    pair_coeff * * 500.0 4.5 1.0
    velocity all create ${{Th}} 15087720
    run 100
    pair_coeff * * 1000.0 4.5 1.0
    velocity all create ${{Th}} 15086189
    run 100
    write_data {0[tmpdata]}1.data

    pair_style hybrid/overlay lj/cut 1.122462 dpd/tstat 1.0 1.0 1.122462 122347
    pair_modify shift yes
    pair_coeff * * lj/cut 1.0 1.0 1.122462
    pair_coeff * * dpd/tstat 4.5 1.122462
    velocity all create ${{Th}} 1508612013
    run 50
    velocity all create ${{Th}} 15021
    run 50
    velocity all create ${{Th}} 2086111
    run 50
    velocity all create ${{Th}} 126111
    run 50
    velocity all create ${{Th}} 1286111
    run 50
    velocity all create ${{Th}} 15021
    run 50
    velocity all create ${{Th}} 26111
    run 50
    velocity all create ${{Th}} 126111
    run 50
    velocity all create ${{Th}} 1286111
    run 50
    write_data {0[tmpdata]}_push.data
    velocity all create ${{Th}} 15086125

    reset_timestep 0

    run {0[kremerequil]}

    write_data {0[outdata]}.data
    """.format(parameters)

    with open(lammpsscript, 'w') as lmp:
        lmp.write(script+'\n')

def CreateCoolPrep():
    """creates initial cooling preparation lammpsscript"""

    NameOfSim = "coolprep"
    dcdfile = NameOfSim
    lammpsscript = "in"+dcdfile+".txt"
    filedumpskip = "dumpskip.txt"
    tmpdata = "tmp." + NameOfSim
    indata = "tmp"
    outdata = "polymer_rest"
    runpart = 1000000
    dumptraj = runpart/10
    dumpthermo = runpart/200
    timestep = 0.001


    with open(filedumpskip, 'a') as f:
        f.write("log"+NameOfSim+"skip "+str(dumpthermo)+'\n')
        f.write("dump"+NameOfSim+"skip "+str(dumptraj)+'\n')
        f.write("time"+NameOfSim+"step "+str(timestep)+'\n')

    parameters = {'dcdfile':dcdfile,
                'tmpdata':tmpdata,
                'indata':indata,
                'outdata':outdata,
                'runpart':int(runpart),
                'dumptraj':int(dumptraj),
                'dumpthermo':int(dumpthermo),
                "runnvehot": runpart*1,
                'runnvehotmed' : 2*runpart,
                'runnvemed': 2*runpart,
                'runnvemednorm': 2*runpart,
                'runnvenorm': 2*runpart,
                'runnptrest': 5*runpart,
                'runquenchsim': 20*runpart,
                'timestep': timestep
    }


    # Run_times = {"Run_nvehot": runpart*1, "Run_nvecool": 56,"Run_nverest": 6*runpart, "Run_nptcool" : 4*runpart}

    script = """\



    variable        Th equal 1.2   #-> hot
    variable        Tm equal 1.0   #-> med
    variable        Tn equal 0.82   #-> norm
    variable        Ps equal 8.0
    variable        dump_traj equal {0[dumptraj]}

    units           lj
    boundary        p p p
    atom_style      angle


    pair_style  lj96/cut 1.0188
    pair_modify shift yes
    bond_style      harmonic
    angle_style   table spline 181
    read_data   {0[indata]}.data
    pair_coeff * * 0.37775 0.89
    bond_coeff * 1352.0 0.5
    angle_coeff * cgpva.table CG_PVA

    special_bonds   lj 0.0 0.0 1.0

    neighbor        0.4 bin
    neigh_modify    every 1 once no cluster yes
    timestep        {0[timestep]}
    thermo          {0[dumpthermo]}


    # this part is used for long polymer chains and their equilibration
    # velocity all create ${{Th}} 17786140
    # fix bswap all bond/swap 50 0.5 1.3 598934
    # run 100000
    # unfix bswap

    # fix 1 all nve/limit 0.001
    # run 500
    # fix 1 all nve/limit 0.01
    # run 500
    # fix 1 all nve/limit 0.05
    # run 500
    # fix 1 all nve/limit 0.1
    # run 500
    # unfix 1
    # fix 1 all nve
    # run 50000


    # here comes the initial equilibration at T = 1.2, with nve/limit


    velocity        all create   ${{Th}}  1231231

    fix             CENTER_of_mass all recenter INIT INIT INIT

    fix             1 all nve/limit 0.001
    fix             2 all langevin   ${{Th}}    ${{Th}}  10.0 904297
    thermo_style    custom step temp press vol epair ebond eangle etotal
    run             1000
    unfix           1
    unfix           2

    fix             1 all nve/limit 0.05
    fix             2 all langevin   ${{Th}}    ${{Th}}  10.0 904297
    thermo_style    custom step temp press vol epair ebond eangle etotal
    run             1000
    unfix           1
    unfix           2

    fix             1 all nve/limit 0.1
    fix             2 all langevin   ${{Th}}    ${{Th}}  10.0 904297
    thermo_style    custom step temp press vol epair ebond eangle etotal
    run             1000
    unfix           1
    unfix           2

    velocity        all create   ${{Th}}  2342322

    dump                    dump_traj all dcd {0[dumptraj]} {0[dcdfile]}.dcd
    dump_modify     dump_traj sort id unwrap yes

    # nve at T = 1.2

    fix             1 all nve
    fix             2 all langevin   ${{Th}}    ${{Th}}  10.0 904297
    thermo_style    custom step temp press vol epair ebond eangle etotal
    run             {0[runnvehot]}
    unfix           1
    unfix           2
    write_data      {0[tmpdata]}nvehot.data

    # nve cooling at T = 1.2 -> 1.0


    fix             1 all nve
    fix             2 all langevin   ${{Th}}    ${{Tm}}  10.0 904297
    thermo_style    custom step temp press vol epair ebond eangle etotal
    run             {0[runnvehotmed]}
    unfix           1
    unfix           2
    write_data      {0[tmpdata]}nvecooling.data


    # nve at T = 1.0

    velocity        all create   ${{Tm}}  2867876

    fix             1 all nve
    fix             2 all langevin   ${{Tm}}    ${{Tm}}  10.0 904297
    thermo_style    custom step temp press vol epair ebond eangle etotal
    run             {0[runnvemed]}
    unfix           1
    unfix           2
    write_data      {0[tmpdata]}nverest.data



    # nve cooling at T = 1.0 -> 0.82
    velocity        all create   ${{Tm}}  2867876

    fix             1 all nve
    fix             2 all langevin   ${{Tm}}    ${{Tn}}  10.0 904297
    thermo_style    custom step temp press vol epair ebond eangle etotal
    run             {0[runnvemednorm]}
    unfix           1
    unfix           2
    write_data      {0[tmpdata]}nvecool2.data


    # nve at T = 0.82
    velocity        all create   ${{Tn}}  2867876

    fix             1 all nve
    fix             2 all langevin   ${{Tn}}    ${{Tn}}  10.0 904297
    thermo_style    custom step temp press vol epair ebond eangle etotal
    run             {0[runnvenorm]}
    unfix           1
    unfix           2
    write_data      {0[outdata]}.data


    """.format(parameters)

    with open(lammpsscript, 'w') as lmp:
        lmp.write(script+'\n')
    return None



def CreateStretchPrep():
    """creates initial cooling preparation lammpsscript"""

    NameOfSim = "stretch"
    dcdfile = NameOfSim
    lammpsscript = "in"+dcdfile+".txt"
    filedumpskip = "dumpskip.txt"

    tmpdata = "tmp." + NameOfSim
    indata = "polymer_0.8"
    outdata = "polymer_stretch"
    runpart = 1000000
    dumptraj = runpart/50
    dumpthermo = runpart/200
    timestep = 0.005
    my_trate = np.log(2.)/float(runpart)
    drag = 2.0

    with open(filedumpskip, 'a') as f:
        f.write("log"+NameOfSim+"skip "+str(dumpthermo)+'\n')
        f.write("dump"+NameOfSim+"skip "+str(dumptraj)+'\n')
        f.write("time"+NameOfSim+"step "+str(timestep)+'\n')

    parameters = {'dcdfile': dcdfile,
                  'tmpdata': tmpdata,
                  'indata': indata,
                  'outdata': outdata,
                  'runpart': int(runpart),
                  'dumptraj': int(dumptraj),
                  'dumpthermo': int(dumpthermo),
                  "runnvehot": runpart*1,
                  'runnvehotmed': 2*runpart,
                  'runnvemed': 2*runpart,
                  'runnvemednorm': 2*runpart,
                  'runnvenorm': 2*runpart,
                  'runnptrest': 5*runpart,
                  'runnptstretch': 1*runpart,
                  'my_trate': my_trate,
                  'drag': drag,
                  'runquenchsim': 20*runpart,
                  'timestep': timestep}


    # Run_times = {"Run_nvehot": runpart*1, "Run_nvecool": 56,"Run_nverest": 6*runpart, "Run_nptcool" : 4*runpart}

    script = """\



    variable        Th equal 1.2
    variable        Tm equal 1.0
    variable        Tn equal 0.82
    variable        Ps equal 8.0
    variable        dump_traj equal {0[dumptraj]}

    units           lj
    boundary        p p p
    atom_style      angle


    pair_style  lj96/cut 1.0188
    pair_modify shift yes
    bond_style      harmonic
    angle_style   table spline 181
    read_data   {0[indata]}.data
    pair_coeff * * 0.37775 0.89
    bond_coeff * 1352.0 0.5
    angle_coeff * cgpva.table CG_PVA

    special_bonds   lj 0.0 0.0 1.0

    neighbor        0.4 bin
    neigh_modify    every 1 once no cluster yes
    timestep        {0[timestep]}
    thermo          {0[dumpthermo]}


    # here comes the initial equilibration at T = 1.2, with nve/limit
    dump                    dump_traj all dcd {0[dumptraj]} {0[dcdfile]}.dcd
    dump_modify     dump_traj sort id unwrap yes



    velocity all create ${{Tn}} 4928459 dist gaussian

    fix    1 all npt temp ${{Tn}} ${{Tn}} 10 y 0 0 1000 z 0 0 1000 drag {0[drag]}
    fix         2 all deform 1 x trate {0[my_trate]}  units box remap x

    thermo_style    custom step temp press density vol pe ke epair ebond eangle etotal
    reset_timestep  0
    run             {0[runnptstretch]}
    unfix           1
    unfix           2
    write_data      {0[outdata]}.data

    """.format(parameters)

    with open(lammpsscript, 'w') as lmp:
        lmp.write(script+'\n')


def CreateRestPrep():
    """creates initial cooling preparation lammpsscript"""

    NameOfSim = "restprep"
    dcdfile = NameOfSim
    lammpsscript = "in"+dcdfile+".txt"
    filedumpskip = "dumpskip.txt"

    tmpdata = "tmp." + NameOfSim
    indata = "polymer_rest"
    outdata = "polymer_0.8"
    runpart = 1000000
    dumptraj = runpart/50
    dumpthermo = runpart/200
    timestep = 0.005


    with open(filedumpskip, 'a') as f:
        f.write("log"+NameOfSim+"skip "+str(dumpthermo)+'\n')
        f.write("dump"+NameOfSim+"skip "+str(dumptraj)+'\n')
        f.write("time"+NameOfSim+"step "+str(timestep)+'\n')

    parameters = {'dcdfile': dcdfile,
                  'tmpdata': tmpdata,
                  'indata': indata,
                  'outdata': outdata,
                  'runpart': int(runpart),
                  'dumptraj': int(dumptraj),
                  'dumpthermo': int(dumpthermo),
                  "runnvehot": runpart*1,
                  'runnvehotmed': 2*runpart,
                  'runnvemed': 2*runpart,
                  'runnvemednorm': 2*runpart,
                  'runnvenorm': 2*runpart,
                  'runnptrest': 5*runpart,
                  'runquenchsim': 20*runpart,
                  'timestep': timestep}


    # Run_times = {"Run_nvehot": runpart*1, "Run_nvecool": 56,"Run_nverest": 6*runpart, "Run_nptcool" : 4*runpart}

    script = """\



    variable        Th equal 1.2
    variable        Tm equal 1.0
    variable        Tn equal 0.82
    variable        Ps equal 8.0
    variable        dump_traj equal {0[dumptraj]}

    units           lj
    boundary        p p p
    atom_style      angle


    pair_style  lj96/cut 1.0188
    pair_modify shift yes
    bond_style      harmonic
    angle_style   table spline 181
    read_data   {0[indata]}.data
    pair_coeff * * 0.37775 0.89
    bond_coeff * 1352.0 0.5
    angle_coeff * cgpva.table CG_PVA

    special_bonds   lj 0.0 0.0 1.0

    neighbor        0.4 bin
    neigh_modify    every 1 once no cluster yes
    timestep        {0[timestep]}
    thermo          {0[dumpthermo]}


    # here comes the initial equilibration at T = 1.2, with nve/limit
    dump                    dump_traj all dcd {0[dumptraj]} {0[dcdfile]}.dcd
    dump_modify     dump_traj sort id unwrap yes



	velocity all create ${{Tn}} 4928459 dist gaussian
	velocity all zero linear
	velocity all zero angular

# initial relaxation which will make the pressure from 0.0 -> 8.0 lj units
# and will relax the veloctiy at Temp = Tn = 0.82

	fix 1 all nvt temp ${{Tn}} ${{Tn}} 10
	fix 2 all press/berendsen iso 0.0 ${{Ps}} 10.0
	thermo_style    custom step temp press density vol pe ke epair ebond eangle etotal
	run 50000
	unfix 1
	unfix 2
	write_data {0[tmpdata]}.berendsen.data

# Now lets adapt that for the system that is simulated under NPT ensemble

	velocity        all create   ${{Tn}}  324231
	fix             1 all npt temp   ${{Tn}}  ${{Tn}}  10 iso   ${{Ps}} ${{Ps}} 100 drag 1.0
	fix             2 all momentum 100 linear 1 1 1 angular
	thermo_style    custom step temp press density vol pe ke epair ebond eangle etotal
	run             500000
	unfix           1
	unfix           2
	write_data      {0[tmpdata]}.drag.data


# now the actual simulation preparation - after that only quench
	velocity        all create   ${{Tn}}  324231
    fix             1 all npt temp   ${{Tn}}  ${{Tn}}  10 iso   ${{Ps}}    ${{Ps}}  1000
    fix             2 all momentum 10 linear 1 1 1 angular
    thermo_style    custom step temp press density vol pe ke epair ebond eangle etotal
    reset_timestep  0
    run             {0[runnptrest]}
    unfix           1
    unfix           2
    write_data      {0[outdata]}.data


    """.format(parameters)

    with open(lammpsscript, 'w') as lmp:
        lmp.write(script+'\n')
    return None


def CreateQuenchSim(NameOfSim='quench', NumOfRuns=20):
    """creates initial cooling preparation lammpsscript
    input: NameOfSim
    which can be one of ('quench', 'cooling', 'underquench')
    """

    dcdfile = NameOfSim
    lammpsscript = "in"+NameOfSim+".txt"
    filedumpskip = "dumpskip.txt"
    tmpdata = "tmp." + NameOfSim
    indata = "polymer_0.8"
    outdata = "polymer_" + NameOfSim

    runpart = 3000000
    # dumptraj = runpart/200
    dumptraj = runpart/40  # here I modified it to dump less steps
    dumpthermo = runpart/1000  # here as well  less steps
    timestep = 0.005

    with open(filedumpskip, 'a') as f:
        f.write("log"+NameOfSim+"skip "+str(dumpthermo)+'\n')
        f.write("dump"+NameOfSim+"skip "+str(dumptraj)+'\n')
        f.write("time"+NameOfSim+"step "+str(timestep)+'\n')

    parameters = {'dcdfile': dcdfile,
                  'tmpdata': tmpdata,
                  'indata': indata,
                  'outdata': outdata,
                  'runpart': int(runpart),
                  'dumptraj': int(dumptraj),
                  'dumpthermo': int(dumpthermo),
                  "runnvehot": runpart*1,
                  'runnvehotmed': 2*runpart,
                  'runnvemed': 2*runpart,
                  'runnvemednorm': 2*runpart,
                  'runnvenorm': 2*runpart,
                  'runnptrest': 5*runpart,
                  'runquenchsim': NumOfRuns*runpart,
                  'timestep': timestep}

    script = """\


    variable        Th equal 1.2
    variable        Tm equal 1.0
    variable        Tn equal 0.82
    variable        Tlunder equal 0.78
    variable        Tl equal 0.72
    variable        Ps equal 8.0
    variable        dump_traj equal {0[dumptraj]}

    units           lj
    boundary        p p p
    atom_style      angle


    pair_style  lj96/cut 1.0188
    pair_modify shift yes
    bond_style      harmonic
    angle_style   table spline 181
    read_data   ./figures/{0[indata]}.data
    pair_coeff * * 0.37775 0.89
    bond_coeff * 1352.0 0.5
    angle_coeff * cgpva.table CG_PVA

    special_bonds   lj 0.0 0.0 1.0

    neighbor        0.4 bin
    neigh_modify    every 1 once no cluster yes
    timestep        {0[timestep]}
    thermo          {0[dumpthermo]}

    dump                    dump_traj all dcd {0[dumptraj]} {0[dcdfile]}.dcd
    dump_modify     dump_traj sort id unwrap yes

    # velocity        all create   ${{Tl}}  1231231


    """.format(parameters)

    if NameOfSim is 'quench':
        script += """
        # npt quenching at T = 0.82 ---|--|--|-> 0.72, P = 8.0

        fix 1 all npt temp ${{Tl}} ${{Tl}}  10 iso ${{Ps}} ${{Ps}} 1000
        fix 2 all momentum 10 linear 1 1 1 angular
        thermo_style custom step temp press density vol pe ke epair ebond eangle etotal
        reset_timestep  0
    """.format(parameters)
    elif NameOfSim is 'cooling':
        script += """
        # npt cooling at T = 0.82 ---|--|--|-> 0.72, P = 8.0

        fix 1 all npt temp   ${{Tn}}  ${{Tl}}  10 iso   ${{Ps}}    ${{Ps}}  1000
        fix             2 all momentum 10 linear 1 1 1 angular
        thermo_style    custom step temp press density vol pe ke epair ebond eangle etotal
        reset_timestep  0
    """.format(parameters)
    elif NameOfSim is 'underquench':
        script += """
        # npt quenching at T = 0.82 ---|--|--|-> 0.78, P = 8.0
        fix  1 all npt temp ${{Tlunder}} ${{Tlunder}}  10 iso ${{Ps}} ${{Ps}}  1000
        fix             2 all momentum 10 linear 1 1 1 angular
        thermo_style    custom step temp press density vol pe ke epair ebond eangle etotal
        reset_timestep  0
    """.format(parameters)
    else:
        raise ValueError("you need to specify what\
         type of sim you want to run")

    string = """\
    run             {0[runpart]} start 0 stop {0[runquenchsim]}
    write_data  {0[tmpdata]}.*.data
    """.format(parameters)

    finish = """\
    write_data      {0[outdata]}.data
    """.format(parameters)


    for i in range(NumOfRuns):
        script += string
    script += finish

    with open(lammpsscript, 'w') as lmp:
        lmp.write(script+'\n')


def main():
    """creates lammps files"""
    args = read_parameters.read_createlammps()
    if args.quench:
        CreateKremer()
        CreateTable()
        CreateCoolPrep()
        CreateRestPrep()
        CreateQuenchSim(NameOfSim='quench', NumOfRuns=args.NumOfRuns)
    elif args.cooling:
        CreateTable()
        CreateQuenchSim(NameOfSim='cooling', NumOfRuns=args.NumOfRuns)
    elif args.underquench:
        CreateTable()
        CreateQuenchSim(NameOfSim='underquench', NumOfRuns=args.NumOfRuns)
    elif args.stretch:
        CreateTable()
        CreateStretchPrep()
    else:
        raise ValueError("you need to specify what\
         type of sim you want to run")
if __name__ == '__main__':
    main()
