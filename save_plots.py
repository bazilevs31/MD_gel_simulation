#!/usr/bin/env python
import os
import datetime
import numpy as np
import termcolor
import matplotlib
matplotlib.use('Agg')
import pylab
import matplotlib.pyplot as plt
import matplotlib.backends.backend_pdf

# def print_everything_dict(**kwargs):
#     for name, value in kwargs.items():
#         print '{0}={1}'.format(name,value)
# print_everything_dict(apple='fruit',cherry='berry')

# introduce special labelling system for AnalyzePlot.py
# so it goes like this : somestringINTnINTadditional_sting,
# this additional string will be also printed in the labeling
# def save_bar_plot():
#     """
#     saves a bar plot of data
#     useful for different CreatePDChain.py
#     """
#     plt.xlabel(r'$\mathrm{Chain\ length}$')
#     plt.ylabel(r'$\mathrm{Number\ of\ chains}$')
#     plt.title(r'$\mathrm{Polydisperse\ chains\ with }\ \mu=%d,\ \sigma=%d,\ PDI=%5.3f$' % (mu,sigma,PDI) )
#     plt.axis([bins.min()-10, bins.max()+10, 0, n.max() + int(0.05*n.max())])
#     plt.grid(True)
#     plt.savefig(figuresdir + '/pdmu%ds%dPDI%5.3f.pdf' % (mu,sigma,PDI))
#     np.savez(figuresdir + '/pdmu%ds%dPDI%5.3f' % (mu,sigma,PDI), bins, n)
#     plt.show()

def save_write_to_file(bins, n, PDI):
    """
    saves a .chain file for given n, bins
    input: bins, n, PDI
    output: None
    write created distribution to a file
    input mu, sigma, Nchains
    output ( creates a file pdMuSigma)
    """

    n = n.astype(int)
    bins = bins.astype(int)
    randseed = np.random.randint(int(102220*bins.max()))
    filename = 'bidisp'
    parameters = {'randseed': randseed,
                  'Nchains': 0,
                  'ChainLength': 0,
                  'ngroups': (len(bins)),
                  'filename': filename,
                  'PDI': PDI}

    script = """\
Polymer chain definition

0.8442          rhostar
{0[randseed]}   random # seed (8 digits or less)
{0[ngroups]}    # of sets of chains (blank line + 6 values for each set)
0               molecule tag rule: 0 = by mol, 1 = from 1 end, 2 = from 2 ends

""".format(parameters)

    for i in xrange(len(bins)):
        parameters['ChainLength'] = bins[i]
        parameters['Nchains'] = n[i]
        script += """\
{0[Nchains]}             number of chains
{0[ChainLength]}             monomers/chain
1               type of monomers (for output into LAMMPS file)
1               type of bonds (for output into LAMMPS file)
0.85            distance between monomers (in reduced units)
1.05            no distance less than this from site i-1 to i+1 (reduced unit)

""".format(parameters)

    with open('{0[filename]}{0[ChainLength]}n{0[Nchains]}.chain'.
              format(parameters), 'w') as f:
        f.write(script+'\n')
    return None


def save_create_hist(bins, n, **kwargs):
    """
    plots a histogram, prints info PDI, n, bins
    input: x, y, type('bi','pd'), mu, sigma, Nchains
    output: PDI
    this is primaraly for plotting .chain files
    for createting bidisp and Polydisperse melts
    prints parameters of the system
    plots the system built
    PDI = M_w/M_n
        bins: chain lengths, binned into Nbins
        nchains = n
        chainlength = bins

        n_i: population of each bin
        Mw = n_i*bins_{i+1}^2 / n_i*bins_{i+1}
        Mn = n_i*bins_{i+1} / n_i
    """
    # curdir, figuresdir=get_currdir()
    # parameters = kwargs.items()
    parameters = kwargs
    print bins, n
    Natoms = (n*bins).sum()
    N1 = 0
    N2 = 0
    N3 = 0

    for i in range(len(n)):
        N1 += n[i]*bins[i]**2
        N2 += n[i]*bins[i]
        N3 += n[i]
    M_w = float(N1)/float(N2)
    M_n = float(N2)/float(N3)
    PDI = M_w/M_n

    print ("number of chains are %d " % n.sum())
    print ("number of atoms is %d" % Natoms)
    print ("PDI is %f" % PDI)
    plt.xlabel(r'$\mathrm{Chain\ length}$')
    plt.ylabel(r'$\mathrm{Number\ of\ chains}$')
    plt.title('${systemtype}\ chains\$'.format(systemtype=parameters['type']))
    plt.axis([0.8*bins.min(), 1.2*bins.max(), 0, 1.1*n.max()])
    plt.grid(True)

    if ((len(parameters) == 0) or (parameters['type'] == 'bi')):
        filename = 'BIchain'
        for i in range(len(n)):
            filename += "C"+str(bins[i])+"l"+str(n[i])
    elif parameters['type'] == 'pd':
        filename = 'PD'
        filename += '\mu' + str(parameters['mu'])
        filename += '\sigma' + str(parameters['sigma'])
        filename += 'N_{chains}' + str(parameters['Nchains'])
    plt.bar(bins, n, width=(bins[1]-bins[0])/3.0,
            facecolor='green',
            alpha=0.75,
            label='$'+filename+'$')
    # plt.legend(loc='best')
    # name1 = 'l'.join(map(str, bins))
    # name2 = 'n'.join(map(str, n))
    # filename = name1+'c'+name2
    plt.savefig(filename + '.pdf')
    np.savez(filename, bins, n)
    # print "number   %d" % x.size
    print "number of chains are %d " % n.sum()
    print "number of atoms is %d" % Natoms
    print "PDI is %f" % PDI
    print "bins"
    print bins
    print "population"
    print n
    return PDI


def save_plot(x, y, **kwargs):
    """
    plots stuff
    input: x, y, **kwargs
    x,y : arrays of the same length to plot
    kwargs : additional parameters to make the plot better
        plotname : name of algorithm used for calculating x,y
        (for example ICC,YCryst)
            ICC : Individual Chain Crystallinity Parameter
            YC : Yamamoto Crystallinity Parameter
            Cos : Average cos2*8eta
            Rg : array of Rg crystallinity parameters
            Re : array of Re : end to end distance evolution
            Bead mean:square displacement
            Sq : static structure factor
            Lml : lamellae analysis
            bondacf: bond correlation function analysis for L_p, when this
                parameter is provided, then program analyzes first points
                and fits them, to understand persistance length.
                C(n) = <u_0 u_n >

        name : name of the psffile
        xlog, ylog : whether you want log scale or not
        title : name of the plot in the title
        xlabel, ylabel : labels
        duplicate: True, False
            whether you want to duplicate results of the simulations in
            ~/results_all/ folder

    TODO: bondacf to work properly

    """
    x = np.asarray(x, dtype='float32')
    y = np.asarray(y, dtype='float32')
    curdir = os.getcwd()
    figuresdir = curdir+'/figures'
    if not os.path.exists(figuresdir):
        os.mkdir(figuresdir)
    parameters = dict()
    input_dict = {key: value for (key, value) in kwargs.items()}
    parameters['plotname'] = input_dict.get('plotname', 'ICC')
    parameters['name'] = input_dict.get('name', 'poly')
    parameters['xlabel'] = input_dict.get('xlabel', 'time, 10^6 lj units')
    parameters['ylabel'] = input_dict.get('ylabel', 'crystallinity')
    parameters['title'] = input_dict.get('title', 'crystallinity')
    parameters['xlog'] = input_dict.get('xlog', False)
    parameters['ylog'] = input_dict.get('ylog', False)
    parameters['ymin'] = input_dict.get('ymin', 0.0)
    parameters['ymax'] = input_dict.get('ymax', 1.0)
    parameters['timestep'] = input_dict.get('timestep', 0.005)
    parameters['dumpskip'] = input_dict.get('dumpskip', 200000)
    parameters['duplicate'] = input_dict.get('duplicate', False)
    print parameters['timestep']
    parameters['figuresdir'] = figuresdir
    print termcolor.colored('parameters have been red', 'green')
    print parameters

    plt.ylabel('${0[ylabel]}$'.format(parameters))
    plt.xlabel('${0[xlabel]}$'.format(parameters))
    plt.title('{0[title]}'.format(parameters))
    plt.grid(True)
    # plt.ylim('{0[ymin]},{0[ymax]}'.format(parameters))

# here is the stuff that requires time as their second argument
# so if we need scaling of time steps it will do that correctly
    time_stuff = ['ICC', 'YC', 'Cos', 'Rg', 'Re', 'Lml']
    # factor = 1e6 # the time units
    # time = time*dumpskip*timestep/factor

    if parameters['plotname'] in time_stuff:
        factor = 1e6
        print parameters['dumpskip']
        print parameters['timestep']
        print parameters['dumpskip']*parameters['timestep']
        x *= float(parameters['dumpskip'])*float(parameters['timestep'])
        x /= factor

    if parameters['plotname'] is 'bondacf':
        print "ylogscale"

        # plt.xscale('log')
        plt.ylabel(r'$\mathrm{log(C(n))}$')
        plt.xlabel(r'$\mathrm{n}$')
        Nfirst = 40

        # plt.ylim(-0.1, 1.1)
        xfirst = x[:Nfirst]
        m, b = pylab.polyfit(xfirst, y[:Nfirst], 1)
        print m, b
        print x, y
        plt.plot(x, y, 'go', x, np.exp(m*x + b), 'r-')
        plt.yscale('log')
        plt.xscale('linear')
        # plt.plot(x, np.exp(m*x + b), 'r-')
        np.savez('{0[figuresdir]}/{0[plotname]}{0[name]}'
                 .format(parameters), x, y)
        plt.savefig('{0[figuresdir]}/{0[plotname]}{0[name]}.pdf'
                    .format(parameters))
        return parameters

    if parameters['xlog']:
        plt.xscale('log')
        parameters['xlabel'] = ''.join(('log', 'parameters[\'xlabel\']'))
        print "xlogscale"
    else:
        print "regular xscale"
    if parameters['ylog']:
        plt.yscale('log')
        parameters['ylabel'] = ''.join(('log', 'parameters[\'ylabel\']'))
        print "ylogscale"
    else:
        print "regular yscale"

    plt.plot(x, y, 'bo-', lw=1.5)
    np.savez('{0[figuresdir]}/{0[plotname]}{0[name]}'.format(parameters), x, y)
    plt.savefig('{0[figuresdir]}/{0[plotname]}{0[name]}.pdf'
                .format(parameters))

    # if parameters['duplicate']:
    #     # duplicate the results in the results_all folder
    #     np.savez('~/results_all/{0[plotname]}{0[name]}'
    #              .format(parameters), x, y)
    #     plt.savefig('~/results_all/{0[plotname]}{0[name]}.pdf'
    #                 .format(parameters))

    return parameters

# this file will have bunch of different save_plot
# functions which will plot information correcponing to the file

# def save_plot_log(x,thermodata,outfile,name,Nevery,logplot):


def save_plot_log(x, thermodata, **kwargs):
    """
    plots data from thermo log file
    """

    curdir = os.getcwd()
    figuresdir = curdir+'/figures'
    if not os.path.exists(figuresdir):
        os.mkdir(figuresdir)
    parameters = dict()
    input_dict = {key: value for (key, value) in kwargs.items()}
    Nevery = int(input_dict['Nevery'])
    parameters['plotname'] = input_dict.get('plotname', 'log')
    parameters['name'] = input_dict.get('name', 'poly')
    parameters['xlabel'] = input_dict.get('xlabel', 'time, 10^6 lj units')
    parameters['figuresdir'] = figuresdir
    print termcolor.colored('parameters have been red', 'green')
    print parameters
    with matplotlib.backends.backend_pdf.PdfPages('{0[plotname]}{0[name]}.pdf'.
                                                  format(parameters),
                                                  'w') as pdf:

        for i, element in enumerate(thermodata):
            print element
            plt.figure(figsize=(10, 6))
            plt.plot(x[::Nevery], thermodata[str(element)][::Nevery], 'b-')
            plt.title(str(element))
            plt.ylabel(str(element))
            pdf.savefig()  # saves the current figure into a pdf page
            plt.close()
        d = pdf.infodict()
        d['Title'] = 'equilibration of polymer thermodynamic log'
        d['Author'] = u'Vasiliy M. Triandafilidi'
        d['Subject'] = 'How to plot lammps log file thermodynamic entries'
        d['Keywords'] = 'equilibration lammps log pizza'
        d['CreationDate'] = datetime.datetime(2014, 01, 07)
        d['ModDate'] = datetime.datetime.today()

    plt.figure()
    time_array = x[::Nevery]
    Temp_array = thermodata['Temp'][::Nevery]
    # Density_array = thermodata['Density'][::Nevery]
    PotEng_array = thermodata['PotEng'][::Nevery]
    plt.ylim(0.3, 1.1)
    plt.title('Thermodynamic Parameters Evoution')
    plt.ylabel('Normalized Thermodynamic parameters')
    plt.xlabel('${0[xlabel]}$'.format(parameters))
    plt.plot(time_array, Temp_array, 'r-',
     label=r'$Temp$')
    # plt.plot(time_array, Density_array/Density_array.max(), 'g-',
     # label=r'$\frac{\rho}{\rho_{max}}$')
    plt.plot(time_array, PotEng_array/PotEng_array.max(), 'b-',
     label=r'$\frac{E_{pot}}{E_{potmax}}$')
    # plt.legend(loc='best')
    # plt.legend()
    plt.legend(loc='best')
    plt.savefig(r'publication{0[name]}.pdf'.format(parameters))
    # np.savez('{0[figuresdir]}/{0[name]}'.format(parameters), time_array, Temp_array, Density_array)
    return None


def main():

    save_plot(np.arange(100), np.arange(100),
              xlabel='lala', title='lalalala', plotname='bondacf')
    return None

if __name__ == '__main__':
    main()
