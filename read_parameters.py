#!/usr/bin/env python

import os
import argparse
from termcolor import colored
import matplotlib

# modules for readparameters

# modify read_pbs to read, and CreatePbs to use the dumpskip.txt file

def read_traj_vmd():
    """
    read parameters from the commandline
    input
    output datafile,traj, trajskip, startframe, endframe, psffile
    uses vmd pluging to create a psf topology file
    """
    # parser = argparse.ArgumentParser()

    parser = argparse.ArgumentParser(description=__doc__,
                            formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-f", "--psf", dest="psffile",
                     help="Name of the future files, all other files will start with FILE",
                     metavar="FILE")
    parser.add_argument("-d", "--data", dest="datafile",
                    default="./figures/polymer_0.8.data",
                    # type=lambda x: is_valid_file(parser, x),
                    help="read datafile and if exists then convert it to psf file by invoking a vmd script",
                    metavar="FILE")

    parser.add_argument("-t", "--trajectroy", dest="traj",
                        default="quench.dcd",
                        type=lambda x: is_valid_file(parser, x),
                        help="Input trajectory file)", metavar="FILE")

    parser.add_argument("-e", "--endframe", dest="endframe",
                    default=-1,
                    type=int,
                    help="End frame of the trajectory file type (default: %(default)s)")
    parser.add_argument("-st",
                    "--startframe",
                    dest="startframe",
                    default=0,
                    type=int,
                    help="Start frame of the trajectory file type (default: %(default)s)")


    ######## skipping info

    parser.add_argument("-autos",
                    "--auto_traj_skip",
                    action="store_true",
                    dest="auto_traj_skip",
                    default=False,
                    help="---Skipping---: do you want skipping automaticallty configured as 120 frames of whole trajectory. If provided --trajskip becomes a dummy variable(default: %(default)s)")
    parser.add_argument("-s",  "--trajskip", dest="trajskip",
                    default=40,
                    type=int,
                    help="How many steps are to be skipped when trajectory \
                                        file is being red\
                                        (needs to be > 1, < number of frames) \
                                        type (default: %(default)s)")
    parser.add_argument("-fd",
                    "--filedumpskip",
                    action="store_true",
                    dest="filedumpskip",
                    default=False,
                    help="---Skipping---: Use dumpskip info from dumpskip.txt(default: %(default)s)")


    parser.add_argument("-kp",
                    "--keyparameter",
                    dest="keyparameter",
                    default="quench",
                    help="keyparameter, can be: kremer,coolprep,restprep,quench")


    parser.add_argument("-ds",
                    "--dumptrajskip",
                    dest="dumpskip",
                    default=200000,
                    type=float,
                    help="---Skipping---: This is a factor for transforming trajectory frames into lj units, numframes*dumpskip = simtime [lj] (default: %(default)s)")

    parser.add_argument("-ts",
                    "--timestep",
                    dest="timestep",
                    default=0.001,
                    type=float,
                    help="in case filelogskip=False, This is a timestep for transforming trajectory frames into lj units, numframes*logskip = simtime [lj] (default: %(default)s)")



    ## ----------Crystallization ----------------

    parser.add_argument("-ne",
                    "--neighbor",
                    dest="neigh",
                    default=10,
                    type=int,
                    help="ICC Crystalinity: Number of neighbours to consider for AnalyzeChain crystalinity parameter program (default: %(default)s)")
    parser.add_argument("-th",
                    "--threshold",
                    dest="threshold",
                    default=0.95,
                    type=float,
                    help="ICC Crystalinity: threshold for ICC crystalinity parameter, chains with p2 higher than threshold are \
                    considered to be crystaline(default: %(default)s)")

    parser.add_argument("-na",
                    "--nAtomsPerBox",
                    dest="nAtomsPerBox",
                    default=3,
                    type=int,
                    help="---YAMAMOTO---: Mesh size for Yamamoto crystalinity parameter program, How many atoms per little mesh box do you want? (default: %(default)s)")

    parser.add_argument("-wr",
                    "--wrap",
                    action="store_true",
                    dest="wrap",
                    default=False,
                    help="---YAMAMOTO---: wraping trajectory or no? Yamamoto crystalinity parameter program (default: %(default)s)")

    ## --------- for other programs like AnalyzeRDF, AnalyzeEnd2End

    parser.add_argument("-noff",
                    "--noffset",
                    dest="Noff",
                    default=1,
                    type=int,
                    help="End2End: Number of last points that we don't want to see (default: %(default)s)")


    parser.add_argument("-calc",
                    "--calculateRDF",
                    action="store_true",
                    dest="calc",
                    default=False,
                    help="(---RDF---: Do you want to write .vmd script that can analyze RDF (default: %(default)s)")
    parser.add_argument("-p",
                    "--plot",
                    action="store_true",
                    dest="plot",
                    default=False,
                    help="---RDF---: Do you want to plot your results? (default: %(default)s)")
    parser.add_argument("-lp",
                    "--logplot",
                    action="store_true",
                    dest="logplot",
                    default=False,
                    help="Do you want time(or N) to be plotted with log(t) scale? (default: %(default)s)")




    args = parser.parse_args()
    args = create_psf(args)
    if args.wrap==True:
        args = create_wrapdcd(args)
    else:
        print " wrap is set to False, using non-wrapped coordinates"

    print args
    print colored('parameters have been red', 'green')

    return args

def read_log():
    """
    provide information for logfile
    output: args , containing logfile, name of the output file
    """

    parser = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("-l", "--log", dest="logfile",
                        default="./figures/log.kremer",
                        help="read logfile to analyze")

    parser.add_argument("-s",
                    "--stride",
                    dest="initoffset",
                    default=1,
                    type=int,
                    help="How many first steps are to be skipped, because they are ususally too big")
    parser.add_argument("-nv",
                    "--Nevery",
                    dest="Nevery",
                    default=10,
                    type=int,
                    help="Plot every this step")

    parser.add_argument("-o",
                    "--output",
                    dest="outfile",
                    default="sim",
                    help="results of the logfile thermo info ploting will be saved here")

    parser.add_argument("-al",
                    "--all",
                    action="store_true",
                    dest="all",
                    default=False,
                    help="---Plot all---: do you want to plot all log files in ./figures(default: %(default)s)")

    ######## skipping info

    parser.add_argument("-fd",
                    "--filelogskip",
                    action="store_true",
                    dest="filelogskip",
                    default=False,
                    help="---Skipping---: do you want to get logskip,timestep from a dumpskip.txt file(default: %(default)s)")

    parser.add_argument("-kp",
                    "--keyparameter",
                    dest="keyparameter",
                    default="kremer",
                    help="keyparameter, can be: kremer,coolprep,restprep,quench")


    parser.add_argument("-ls",
                    "--logskip",
                    dest="logskip",
                    default=200000,
                    type=float,
                    help="---Skipping---: (in case filelogskip=False, This is a factor for transforming trajectory frames into lj units, numframes*logskip = simtime [lj] (default: %(default)s)")

    parser.add_argument("-ts",
                    "--timestep",
                    dest="timestep",
                    default=0.001,
                    type=float,
                    help="in case filelogskip=False, This is a timestep for transforming trajectory frames into lj units, numframes*logskip = simtime [lj] (default: %(default)s)")

    ######### log plotting

    parser.add_argument("-lp",
                    "--logplot",
                    action="store_true",
                    dest="logplot",
                    default=False,
                    help="Do you want time(or N) to be plotted with log(t) scale? (default: %(default)s)")

    args = parser.parse_args()

    return args

def read_pbs():
    """
    read information about the pbs file you want to create
    this read pbs can read what you want to do
    -about creating polymers
    -running the simulation including the equilibration
    -analazying the results, using certain program

    TODO: for analyzing part choose what kind of analysis you do want to do

    """
    parser = argparse.ArgumentParser(description=__doc__,
                            formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # group = parser.add_mutually_exclusive_group(help='What do you want to do create melt or run sim')
    group = parser.add_mutually_exclusive_group()
    group.add_argument('--create', action='store_true',help='Create and equilibrate the run')
    group.add_argument('--run', action='store_true', help='Run the simulation')
    group.add_argument('--analyze', action='store_true', help='Analyze the alignment using AnalyzeChain.py')
    group.add_argument('--cont', action='store_true', help='continue the simulation using in.continue')
    group.add_argument('--stretch', action='store_true', help='do stretch')

    parser.add_argument("-p",
                        "--procs",
                        dest="procs",
                        default=12,
                        type=int,
                        help="Number of processors")
    parser.add_argument("-n",
                        "--nodes",
                        dest="nodes",
                        default=5,
                        type=int,
                        help="Number of nodes")
    parser.add_argument("-m",
                        "--memmory",
                        dest="memmory",
                        default=600,
                        type=int,
                        help="Memmory to use")
    parser.add_argument("-wm",
                        "--wallminutes",
                        dest="wallminutes",
                        default=0,
                        type=int,
                        help="Walltime minutes")
    parser.add_argument("-wh",
                        "--wallhours",
                        dest="wallhours",
                        default=8,
                        type=int,
                        help="Number of nodes, now for Simulation 1500 atom/\
                         proc one needs  Equil = 7 hours, min, for run it will\
                          be 30 hours")
    parser.add_argument("-name",
                        "--pbsname",
                        dest="pbsname",
                        default='poly',
                        type=str,
                        help="Name of pbs file")
    parser.add_argument("-sim",
                        "--simtype",
                        dest="simtype",
                        default='quench',
                        type=str,
                        help="Name of the lammps file")



    parser.add_argument("-bg",
                    "--bugaboo",
                    action="store_true",
                    dest="bugaboo",
                    default=False,
                    help="Do you want to make a lammps file on BUGABOO? (default: %(default)s)")
    parser.add_argument("-orc",
                    "--orcinus",
                    action="store_true",
                    dest="orcinus",
                    default=False,
                    help="Do you want to make a lammps file on ORCINUS? (default: %(default)s)")
    parser.add_argument("-jas",
                    "--jasper",
                    action="store_true",
                    dest="jasper",
                    default=False,
                    help="Do you want to make a lammps file on JASPER? (default: %(default)s)")
    parser.add_argument("-gr",
                    "--grex",
                    action="store_true",
                    dest="grex",
                    default=False,
                    help="Do you want to make a lammps file on GREX? (default: %(default)s)")


    args = parser.parse_args()

    print args
    print colored('parameters have been red', 'green')

    return args


def read_createlammps():
    """
    read information about the lammps file you want to create
    --cooling 0.82 -> 0.72 in 60M steps
    --quench 0.82 -> 0.72 instant, and then 60M steps
    --underquench 0.82 -> 0.78 instant, and then in 60M steps
    --stretch stretching simulation 1M steps
    """
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.
                                     ArgumentDefaultsHelpFormatter)

    group = parser.add_mutually_exclusive_group()
    # group = parser.add_mutually_exclusive_group(help='what kind of\
    # sim do you need')
    group.add_argument('--cooling',
                       action='store_true',
                       help='cooling 0.82 -> 0.72 in 60M steps')
    group.add_argument('--quench',
                       action='store_true',
                       help='quench 0.82 -> 0.72 instant, and then 60M steps')
    group.add_argument('--underquench',
                       action='store_true',
                       help='underquench 0.82 -> 0.78 instant,\
                       and then in 60M steps')
    group.add_argument('--stretch',
                       action='store_true',
                       help='stretching in 1M steps')

    parser.add_argument("-n",
                        "--NumOfRuns",
                        dest="NumOfRuns",
                        default=20,
                        type=int,
                        help="Number of runs of 3M steps to do. Regular \
                        simulation require 20*3M steps = 60 M steps and it\
                        runs for 20 hours on jasper")
    args = parser.parse_args()

    print args
    print colored('parameters from read_createlammps have been red', 'green')

    return args


def read_from_file():
    """
    input: file.in, file.param
    read the ploting information
    output: the information for ploting
    """
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.
                                     ArgumentDefaultsHelpFormatter)
    parser.add_argument("-in",
                        "--input",
                        dest="input",
                        default='file.in',
                        type=str,
                        help="input file, which has information about\
                        all npz files (default: %(default)s)")
    parser.add_argument("-p",
                        "--param",
                        dest="param",
                        default='~/pythonfiles/file.param',
                        type=str,
                        help="parameters file that has information about\
                        how to plot (default: %(default)s)")
    parser.add_argument("-t",
                        "--title",
                        dest="titlename",
                        type=str,
                        help="Name of the result pdf and npz files type\
                        (default: %(default)s)")

    args = parser.parse_args()

    print args
    print colored('parameters have been red', 'green')

    return args


def read_plot():
    """read the ploting information
    output: the information for ploting
    """
    parser = argparse.ArgumentParser(description=__doc__,
                        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    group = parser.add_mutually_exclusive_group()
    # group.add_argument('--lines', action='store_true', help='Plot the data with the lines style')
    group.add_argument('--points', action='store_true', help='Plot the data with the points style')
    group.add_argument('--bars', action='store_true', help='Plot the data with the bars style')


    parser.add_argument('-a', action='append', dest='collection',
                        type=lambda x: is_valid_file(parser, x),
                        default=[],
                        help='Add repeated values to a list type (default: %(default)s)',
                        )
    parser.add_argument("-x",
                        "--name_x",
                        dest="name_x",
                        default='arr_0',
                        type=str,
                        help="x direction type (default: %(default)s)")
    parser.add_argument("-y",
                        "--name_y",
                        dest="name_y",
                        default='arr_1',
                        type=str,
                        help="y direction type (default: %(default)s)")
    parser.add_argument("-z",
                        "--name_z",
                        dest="name_z",
                        default='arr_2',
                        type=str,
                        help="z direction type (default: %(default)s)")
    parser.add_argument("-t",
                        "--title",
                        dest="titlename",
                        type=str,
                        help="Name of the result pdf and npz files type (default: %(default)s)")
    parser.add_argument("-xl",
                        "--xlabel",
                        dest="xlabel",
                        default='time,10^6\ lj\ units',
                        type=str,
                        help="xlabel of the plot (default: %(default)s)")
    parser.add_argument("-yl",
                        "--ylabel",
                        dest="ylabel",
                        default='crystalinity',
                        type=str,
                        help="ylabel of the plot (default: %(default)s)")
    parser.add_argument("-b",
                        "--boxes",
                        action="store_true",
                        dest="boxes",
                        default=False,
                        help="Do you need boxes? type (default: %(default)s)")
    parser.add_argument("-l",
                        "--legend",
                        action="store_true",
                        dest="legend",
                        default=True,
                        help="Do you need legend? type (default: %(default)s)")
    parser.add_argument("-lp",
                        "--logplot",
                        action="store_true",
                        dest="logplot",
                        default=False,
                        help="Do you want time to be plotted with log(t) scale? (default: %(default)s)")
    parser.add_argument("-ds",
                        "--dumptrajskip",
                        dest="dumpskip",
                        default=200000,
                        type=float,
                        help="This is a factor for transforming trajectory frames into lj units, numframes*dumpskip = simtime [lj] (default: %(default)s)")

    args = parser.parse_args()


    print args
    print colored('parameters have been red', 'green')

    return args

def read_create_mono():
    """
    create a monodisperse chains
    """
    parser = argparse.ArgumentParser(description=__doc__,
                        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("-nch",
                        "--nchains",
                        dest="Nchains",
                        default=500,
                        type=int,
                        help="Total number of chains type (default: %(default)s)")
    parser.add_argument("-l",
                        "--length",
                        dest="ChainLength",
                        default=100,
                        type=int,
                        help="Chain length type (default: %(default)s)")

    args = parser.parse_args()

    print args
    print colored('parameters have been red', 'green')

    return args

def read_create_pd():
    """
    create a polydisperse chains
    input:
    mu - average chain length
    sigma - chain length distribution width
    Nchains - number of chains
    Nbins - number of bins to discretize the distribution
    dist - distribution we want to use, gaussian or poisson
    """
    parser = argparse.ArgumentParser(description=__doc__,
                        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("-m",
                        "--mu",
                        dest="mu",
                        default=80,
                        type=float,
                        help="average chain length type (default: %(default)s)")
    parser.add_argument("-s",
                        "--sigma",
                        dest="sigma",
                        default=10,
                        type=float,
                        help="width of the chain length distribution type (default: %(default)s)")
    parser.add_argument("-nch",
                        "--nchains",
                        dest="Nchains",
                        default=500,
                        type=int,
                        help="Total number of chains type (default: %(default)s)")
    parser.add_argument("-Nb",
                        "--Nbins",
                        dest="Nbins",
                        default=10,
                        type=int,
                        help="Number of bins type (default: %(default)s)")

    # group = parser.add_mutually_exclusive_group(help='What do you want to do create melt or run sim')
    # group = parser.add_mutually_exclusive_group()
    # group.add_argument('--gaussian', action='store_true',help='Use gaussian distribution to create atoms')
    # group.add_argument('--poisson', action='store_true', help='Use poisson distribution to create atoms')



    args = parser.parse_args()

    print args
    print colored('parameters have been red', 'green')

    return args

def read_create_ndisp():
    """
    create a ndisperse chains
    input should be done in data tuples
    (#of chains, chain length) (#number of chains, chain length) ...

    example:
    for the case of a bidisperse system with 100 C100 and 50 C20 :
        (100,100), (50,20)
    """
    parser = argparse.ArgumentParser(description=__doc__,
                        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument ('-bi', '--bidisperse', dest='bichains', nargs=2, type=int, action='append',help='Import tuples of data \
    numberOfChains ChainLength')
    #create pd melt #create md melt
    args = parser.parse_args()

    print args
    print colored('parameters have been red', 'green')

    return args

def read_calc_pdi():
    """
    this will read .npz file
    calculate pdi index
    plot the pdf with the histogram
    and the PDI index information
    """
    parser = argparse.ArgumentParser(description=__doc__,
                        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-n", "--npz", dest="npzfile",
                default="./figures/polymer_0.8.data",
                # type=lambda x: is_valid_file(parser, x),
                help="read datafile and if exists then convert it to psf file by invoking a vmd script",
                metavar="FILE")
    parser.add_argument("-m",
                        "--mu",
                        dest="mu",
                        type=int,
                        help="average chain length type (default: %(default)s)")
    parser.add_argument("-s",
                        "--sigma",
                        dest="sigma",
                        type=int,
                        help="width of the chain length distribution type (default: %(default)s)")
    args = parser.parse_args()

    print args
    print colored('parameters have been red', 'green')

    return args


def is_valid_file(parser, arg):
    """
    Check if arg is a valid file
    """
    arg = os.path.abspath(arg)
    if not os.path.exists(arg):
        parser.error("The file %s doesn't exist " % arg)
    else:
        return arg

def create_wrapdcd(args):
    """
    given args - if the trajectory file doesn't exist it will be used
    trajfile = trajectory
    looks for wrap+trajfile
    if it can't find it , it will be created from trajfile, with coordinates wrapped.
    and trajectory skipped for enough steps (number of steps will be taken from trajskip)
    later the traj skip will be set to 1
    """
    psffile = os.path.splitext(args.psffile)[0]
    PSF = psffile+'.psf'
    DCD = os.path.splitext(args.traj)[0] + '.dcd'
    DATA = os.path.splitext(args.datafile)[0] + '.data'
    vmdscript = "create_wrapdcd.vmd"

    # write VMD loader script
    parameters = {'vmdfile': vmdscript,
                  'topology': PSF,
                  'datafile': DATA,
                  'trajectory': os.path.basename(DCD),
                  'trajskip':args.trajskip}

    script = """\
    mol new "{0[topology]}"
    animate read dcd "{0[trajectory]}" skip {0[trajskip]} waitfor all
    package require pbctools
    pbc wrap -all
    animate write dcd wrap{0[trajectory]}
    echo "writing dcd to wrap{0[trajectory]} "
    exit
    """.format(parameters)


    with open(vmdscript, 'w') as tcl:
        tcl.write(script+'\n')

    os.system("vmd -dispdev text -e {0[vmdfile]}".format(parameters))
    print colored("If there is an error with {0}: 'source {0}' to load everything manually, then repeat ".format(vmdscript),"blue")
    print colored("running the python script with explicict parameters that were generated".format(vmdscript),"blue")

    args.traj = "wrap{0[trajectory]}".format(parameters)
    args.trajskip = 1
    print colored("Wrote VMD script {0}  ".format(vmdscript),"cyan")
    print colored('run: AnalyzeCrystYamamoto.py -f {0[topology]} -t wraptrajectory_quench.dcd -s 1 -ns 10 -ds 2000000'.format(parameters), 'cyan')


    return args


def create_psf(args):
    """
    given data file produce psf file if it doesn't exist yet
    if it does then use it
    """

    if not os.path.exists(os.path.abspath(args.psffile)):
        psffile = os.path.splitext(args.psffile)[0]
        PSF = psffile+'.psf'
        DCD = os.path.splitext(args.traj)[0] + '.dcd'
        DATA = os.path.splitext(args.datafile)[0] + '.data'

        vmdscript = "create_psf.vmd"

        # write VMD loader script
        parameters = {'vmdfile': vmdscript,
                      'topology': PSF,
                      'datafile': DATA,
                      'trajectory': DCD,
                      'trajskip':args.trajskip}

        script = """\
            package require topotools
            topo readlammpsdata "{0[datafile]}" angle
            animate write psf "{0[topology]}"
        exit
        """.format(parameters)


        with open(vmdscript, 'w') as tcl:
            tcl.write(script+'\n')

        os.system("vmd -dispdev text -e {0[vmdfile]}".format(parameters))

        print colored("Wrote VMD script {0}  ".format(vmdscript),"blue")
        print colored("If there is an error with {0}: 'source {0}' to load everything manually, then repeat ".format(vmdscript),"blue")
        print colored("running the python script with explicict parameters that were generated".format(vmdscript),"blue")
    else:
        print "the file %s already exists" % args.psffile


    return args





def main():
    """main
    """
    # print 1
    # read_parameters()
    read_traj_vmd()
    # print a
    # b = read_create_ndisp()
    # print b
    # c = read_log()
    # print c
    # d = read_pbs()
    # print d
    # e = read_plot()
    # print e
    # f = read_traj()
    # print f

if __name__ == '__main__':
    main()
