import numpy as np
import interpret_plots
import interpret_cdf
import os
from math import sqrt, ceil
from spacepy import pycdf
from datetime import datetime
import matplotlib.pyplot as plt
import json

dir_local_plot = os.path.join("plots_baspro")
dir_local_configs = os.path.join("basproplot_configurations")
default_config_fname = "default_template.json"

def identify_solution(dir_solution, dynamic, overwrite_cdf = True): #dir_local_solution #dir_local_output
    # #scan for solution directories:
    # dir_local_solution_dict = {}
    # dir_local_solution_found = None
    # dir_content = os.listdir(dir_local_output)
    # dir_content_ordered = sorted(dir_content,reverse=True)
    # for content in dir_content_ordered:
    #     if os.path.isdir(os.path.join(dir_local_output, content)):
    #         if content == dir_local_solution:
    #             dir_local_solution_found = dir_local_solution
    #             break

    # if type(dir_local_solution_found) == type(None):
    #     #print("Error: cannot find solution {} in {}".format(dir_local_solution, dir_local_output))
    #     return None, None

    # dir_solution = os.path.join(dir_local_output, dir_local_solution_found)

    #create a CDF of the solution, if it doesn't aleady exist:
    if not overwrite_cdf:
        print("","warning: solution CDF will NOT be updated, change to overwrite_cdf = True if needed")
    # fname_cdf = interpret_cdf.convert_to_cdf(dir_local_output, dir_local_solution_found, dynamic, overwrite = overwrite_cdf)
    fname_cdf = interpret_cdf.convert_to_cdf(dir_solution, dynamic=dynamic, overwrite = overwrite_cdf)
    print()
    #

    #print a summary of the CDF:
    print("{} overview:".format(fname_cdf))
    cdf = pycdf.CDF(fname_cdf)
    keys = list(cdf.keys())
    for key in keys:
        attrs = list(cdf[key].attrs)
        print("", key)
        print("", "", "shape:", np.shape(cdf[key]))
        for attr in attrs:
            print("", "", "{}:".format(attr), cdf[key].attrs[attr])
        print()
    print()

    return fname_cdf #dir_solution

def mkdir_new_plot_subdir(dir_local_plot, dir_local_solution):
    # create a new subdir for the plots:
    print("creating a new directory for plot output")
    epochstr = datetime.now().strftime("_%Y%m%d")
    wd_local_plots = os.path.join(dir_local_plot, dir_local_solution + epochstr)
    plotver = 0
    while(os.path.isdir("{}_{}".format(wd_local_plots,plotver))):
        plotver+=1
    wd_local_plots = "{}_{}".format(wd_local_plots,plotver)
    
    wd_local_plots = os.path.join(wd_local_plots, "")
    print("","creating",wd_local_plots)
    os.mkdir(wd_local_plots)
    print("","done")
    print()
    return wd_local_plots



#
#make colourbar plots of time-dependent solution:
#
def plot_dynamic_flux(fname_cdf, wd_local_plots, energies, Llims, epochrange, aeqs, save=True):
    if not len(energies):
        print("","no energies specified, skipping...")
        return None, None

    print("plotting flux as colour through time...")
    figures = {}
    dataset = {}

    #daterange = [datetime.fromtimestamp(x) for x in epochrange]
    #daterangestr = interpret_tools.daterange_to_str(daterange)

    for aeq in aeqs:

        fname_ = "plot_dynamic_flux_E_"+"_".join(["{:.2f}".format(e_) for e_ in energies])+"__{:.1f}d".format(aeq)
        fname_ = fname_.replace(".","-")
        figure, L_data_all, epoch_data_all, flux_data_all = interpret_plots.make_overviewplot(fname_cdf, Llims, epochrange, energies,
                                                                                              aeq = aeq, jlimits = [], anisotropy=False)

        if type(figure) == type(None): return figures, dataset

        fname_plot=os.path.join(wd_local_plots, fname_)
        figures[fname_plot] = figure
        dataset[fname_plot] = [L_data_all, epoch_data_all, flux_data_all]

        if save:
            print("saving",fname_plot)
            figure.savefig(fname_plot+".png",dpi=300, bbox_inches='tight')
            print('',"saved",fname_plot)
    print()

    return figures, dataset

def plot_dynamic_anisotropy(fname_cdf, wd_local_plots, energies, Llims, epochrange, save=True):
    if not len(energies):
        print("","no energies specified, skipping...")
        return None, None

    print("plotting anisotropy as colour through time...")
    figures = {}
    dataset = {}

    #daterange = [datetime.fromtimestamp(x) for x in epochrange]
    #daterangestr = interpret_tools.daterange_to_str(daterange)

    fname_ = "plot_dynamic_anisotropy_E_"+"_".join(["{:.2f}".format(e_) for e_ in energies])+"__anis"
    fname_ = fname_.replace(".","-")

    figure, L_data_all, epoch_data_all, flux_data_all = interpret_plots.make_overviewplot(fname_cdf, Llims, epochrange, energies,
                                                                                          aeq = 90, jlimits = [], anisotropy=True)

    if type(figure) == type(None): return figures, dataset

    fname_plot=os.path.join(wd_local_plots, fname_)
    figures[fname_plot] = figure
    dataset[fname_plot] = [L_data_all, epoch_data_all, flux_data_all]

    if save:
        figures[fname_plot].savefig(fname_plot+".png",dpi=300, bbox_inches='tight')
        print('',"saved",fname_plot)
    print()

    return figures, dataset


#
#make plots of spectrum:
#
def plot_spectrums(fname_cdf, wd_local_plots, Lshells, nt, epochrange, save=True):
    if not len(Lshells):
        print("","no Lshells specified, skipping...")
        return None, None
    if len(epochrange) == 2:
        epochs = np.linspace(epochrange[0], epochrange[1], nt)
    else:
        epochs = []

    print("plotting spectrums...")
    figures = {}
    dataset = {}

    for L_ in Lshells:
        figures[L_] = {}

        iK = 0 #plot the spectrum at every K from this index and up
        inrange = True
        while inrange:
            figure, data, inrange = interpret_plots.plot_spectrum(fname_cdf, iK, L_, nt = nt, epochs = epochs, universal_axes_limits = True)
            
            if type(figures) == type(None): return figures, dataset

            if inrange:
                figures[L_][iK] = figure
                dataset[iK] = data
            iK += 1
    
    if save:
        for L_ in Lshells:
            #make a new directory for each L:
            dir_L = os.path.join(wd_local_plots, "spectrum_L{:.2f}".format(L_).replace(".","-"))
            os.mkdir(dir_L)
            for iK in figures[L_].keys():
                fname_plot = os.path.join(dir_L, 'plot_spectrums_iK{:03d}'.format(iK).replace(".","-"))
                figures[L_][iK].savefig(fname_plot + ".png")
                print("","saved figure to", fname_plot)
    print()

    return figures, dataset


#
#make plots of psd/flux profiles at fixed mu/energy:
#
def arrange_in_table(lst):
    #reshape list into a regular table:
    lst = list(lst)
    if len(lst) < 3:
        arranged = [lst]
    else:
        ncols = ceil(sqrt(1.0*len(lst)))
        arranged = []
        nrows = len(lst) // ncols + ceil((len(lst) % ncols)/ncols)

        nrectangle = ncols * nrows
        lst += [None]*(nrectangle-len(lst)) #add None elements to make rectangular

        for irow in range(nrows):
            arranged.append(lst[irow*ncols:irow*ncols+ncols])
    return arranged

def plot_psd(fname_cdf, wd_local_plots, mus, Llims, nt, epochrange, save=True):
    if not len(mus):
        print("","no mus specified, skipping...")
        return None, None
    if len(epochrange) == 2:
        epochs = np.linspace(epochrange[0], epochrange[1], nt)
    else:
        epochs = []

    iKs = [0, 3]
    print("plotting psd profiles at K indices {}...".format(iKs))
    
    figures = []
    dataset = {}

    fname_pre = "plot_psd"
    fname1_plot = os.path.join(wd_local_plots, fname_pre)

    mus_arranged = arrange_in_table(mus)

    figure1, L_data_all, aeq_data_all, epoch_data_all, sol_data_all = interpret_plots.plot_f_vs_L_mu_panels(
        fname_cdf, mus_arranged, iKs, nt, epochs, Llims)
    dataset[fname1_plot] = [L_data_all, aeq_data_all, epoch_data_all, sol_data_all]

    if type(figure1) == type(None): return figures, dataset

    figures.append(figure1)

    if save:
        figure1.savefig(fname1_plot + ".png")
        print("","saved figure")
    print()

    return figures, dataset

def plot_j_aeq(fname_cdf, wd_local_plots, energies, Llims, aeqs, nt, epochrange, save=True):
    if not len(energies):
        print("","no energies specified, skipping...")
        return None, None
    if len(epochrange) == 2:
        epochs = np.linspace(epochrange[0], epochrange[1], nt)
    else:
        epochs = []

    print("plotting flux at fixed aeq...")

    figures = []
    dataset = {}
    for aeq in aeqs:
        fname_pre = "plot_j_aeq" + "__{:.1f}d".format(aeq)
        fname_pre = fname_pre.replace('.','-')

        energies_arranged = arrange_in_table(energies)

        fname3_plot = os.path.join(wd_local_plots, fname_pre)
        figure3, L_data_all, epoch_data_all, sol_data_all = interpret_plots.plot_f_vs_L_E_fixedaeq_panels(
            fname_cdf, energies_arranged, Llims, aeq, nt, epochs, False, ratio=False)
        dataset[fname3_plot] = [L_data_all, aeq, epoch_data_all, sol_data_all]
        
        if type(figure3) == type(None): return figures, dataset

        figures.append(figure3)

        if save:
            #for name, figure in zip(list(dataset.keys()), figures):
            figure3.savefig(fname3_plot + ".png")
            print("","saved figure")
        print()

    return figures, dataset

def plot_j_ratio(fname_cdf, wd_local_plots, energies, Llims, nt, epochrange, save=True):
    if not len(energies):
        print("","no energies specified, skipping...")
        return None, None
    if len(epochrange) == 2:
        epochs = np.linspace(epochrange[0], epochrange[1], nt)
    else:
        epochs = []

    print("plotting ratio of flux at 90 to 45 degrees...")

    figures = []
    dataset = {}

    fname_pre = "plot_j_ratio"
    fname3_plot = os.path.join(wd_local_plots, fname_pre)

    energies_arranged = arrange_in_table(energies)

    figure3, L_data_all, epoch_data_all, sol_data_all = interpret_plots.plot_f_vs_L_E_fixedaeq_panels(
        fname_cdf, energies_arranged, Llims, -1, nt, epochs, False, ratio=True)
    dataset[fname3_plot] = [L_data_all, -1, epoch_data_all, sol_data_all]
    
    if type(figure3) == type(None): return figures, dataset

    figures.append(figure3)

    if save:
        #for name, figure in zip(list(dataset.keys()), figures):
        figure3.savefig(fname3_plot + ".png")
        print("","saved figure")
    print()

    return figures, dataset

def plot_j_omni(fname_cdf, wd_local_plots, energies, Llims, nt, epochrange, save=True):
    if not len(energies):
        print("","no energies specified, skipping...")
        return None, None
    if len(epochrange) == 2:
        epochs = np.linspace(epochrange[0], epochrange[1], nt)
    else:
        epochs = []

    print("plotting flux at fixed aeq...")

    figures = []
    dataset = {}
    
    fname_pre = "plot_j_omni"
    fname4_plot = os.path.join(wd_local_plots, fname_pre)

    energies_arranged = arrange_in_table(energies)

    figure4, L_data_all, epoch_data_all, sol_data_all = interpret_plots.plot_f_vs_L_E_fixedaeq_panels(
        fname_cdf, energies_arranged, Llims, -1, nt, epochs, True, False)
    dataset[fname4_plot] = [L_data_all, -1, epoch_data_all, sol_data_all]
    
    if type(figure4) == type(None): return figures, dataset

    figures.append(figure4)

    if save:
        #for name, figure in zip(list(dataset.keys()), figures):
        figure4.savefig(fname4_plot + ".png")
        print("","saved figure")
    print()

    return figures, dataset


#
#make plots of pitch angle distributions at constant energies:
#
def plot_pad(fname_cdf, wd_local_plots, Lshells, energies, nt, epochrange, save=True):
    if not len(energies) or not len(Lshells):
        print("","no energies and/or Lshells specified, skipping...")
        return None, None
    if len(epochrange) == 2:
        epochs = np.linspace(epochrange[0], epochrange[1], nt)
    else:
        epochs = []
        
    print("plotting PADs...")    

    figures = []
    dataset = {}

    figname = os.path.join(wd_local_plots, 'plot_pad')
    figure, aeq_data_all, sol_data_all, epoch_data_all = interpret_plots.plot_padist_panels(fname_cdf, energies, Lshells, nt, epochs, plotflux = True)
    dataset[figname] = [aeq_data_all, sol_data_all, epoch_data_all]

    if type(figure) == type(None): return figures, dataset

    figures.append(figure)
    if save:
        figure.savefig(figname + ".png")
        print("","saved figure to", figname)
    print()

    return figures, dataset


def plot_timeseries(fname_cdf, wd_local_plots, energies, Lshells, epochrange, aeqs, save=True):
    if not len(energies) or not len(Lshells):
        print("","no energies and/or Lshells specified, skipping...")
        return None, None
    if type(epochrange) == type(None) or len(epochrange) != 2:
        print("", "no epoch range specified, skipping...")
        return None, None

    print("plotting flux through time at fixed E, L, aeq...")

    figures = []
    dataset = {}
    for aeq in aeqs:
        fname_pre = "plot_timeseries"+"__{:.1f}d".format(aeq)
        fname_pre = fname_pre.replace('.','-')


        fname_plot = os.path.join(wd_local_plots, fname_pre)
        figure, epoch_data_all, sol_data_all = interpret_plots.plot_j_timeseries_fixedaeq_panels(
            fname_cdf, Lshells, energies, aeq, epochrange)
        dataset[fname_plot] = [energies, Lshells, aeq, epoch_data_all, sol_data_all]
        
        if len(sol_data_all):

            figures.append(figure)

            if save:
                figure.savefig(fname_plot + ".png")
                print("","saved figure")
            print()

    return figures, dataset

def plot_map_flux(fname_cdf, wd_local_plots, altitudes_km, energies, aeqs, dt_max_days, save=True):
    if not len(energies):
        print("","no energies and/or Lshells specified, skipping...")
        return None, None

    figures = []
    dataset = {}
    for alt_km in altitudes_km:
        for energy in energies:

            figures_alt_energy, meshgrid_lonlat, data = interpret_plots.plot_map(
                fname_cdf, energy, alt_km, dt_max_days)

            for idxf, figure in enumerate(figures_alt_energy):
                fname_pre = "plot_map_{}_{:.1f}km_{:.2f}MeV".format(idxf, alt_km, energy)
                fname_pre = fname_pre.replace('.','-')
                fname_plot = os.path.join(wd_local_plots, fname_pre)

                dataset[fname_plot] = [energy, alt_km, meshgrid_lonlat, data[idxf]]

                figures.append(figure)

                if save:
                    figure.savefig(fname_plot + ".png")
                    print("","saved figure")
                print()

    return figures, dataset

# #
# #make plots of pitch angle n at constant energies:
# #
# def plot_n_profile(fname_cdf, wd_local_plots, energies, Llims, nt, epochs, save=True):
#     if not len(energies):
#         print("","no energies specified, skipping...")
#         return None, None

#     print("plotting anistropy... (may take a long time)")
    
#     figure = interpret_plots.plot_n_vs_L_panels(fname_cdf, energies, nt, epochs, Llims) #L_axis_plot

#     if save:
#         fname_plot = os.path.join(wd_local_plots, "n_profile")
#         figure.savefig(fname_plot + ".png")
#         print("","saved figure to", fname_plot)
#     print()

#     dataset = {}

#     return figure, dataset



functions_plotting = [
    plot_dynamic_flux,
    plot_dynamic_anisotropy,
    plot_spectrums,
    plot_psd,
    plot_j_aeq,
    plot_j_ratio,
    plot_j_omni,
    plot_pad,
    plot_timeseries,
    plot_map_flux]#,
    #plot_n_profile]

def regenerate_template_config():
    if not os.path.isdir(dir_local_configs):
        print("", "creating {}".format(dir_local_configs))
        os.mkdir(dir_local_configs)

    config = {"use_dynamic_solution" : True}

    for function in functions_plotting:
        function_name = function.__name__

        print("plot function key {} has 'enabled':True, expecting arguments:".format(function_name))
        argument_keys = function.__code__.co_varnames[:function.__code__.co_argcount]
        argument_dict = {}
        for argument in argument_keys:
            if argument == 'fname_cdf':
                continue
            elif argument == 'wd_local_plots':
                continue
            elif argument == 'save':
                continue
            elif argument == "nt":
                argument_dict[argument] = 2
            else:
                argument_dict[argument] = []
        argument_dict["enabled"] = True

        config[function_name] = argument_dict

    fname = os.path.join(dir_local_configs, default_config_fname)
    print("regenerating default config file at {}".format(fname))
    with open(fname, 'w') as fo:
        json.dump(config, fo, indent=4)


def check_key(dictionary, keytocheck, return_if_absent = None):
    try:
        value = dictionary[keytocheck]
    except:
        #print("Error: key {} not found in configuration, regenerate the template to check syntax".format(keytocheck))
        value = return_if_absent
    return value

def generate_from_config(dir_solution, fname_config, epochrange_focus = [], overwrite_cdf = True):
    #check filesystem:
    if not os.path.isdir(dir_solution):
        print("error: expecting solution directory at {}".format(dir_solution))
        return 0

    if not os.path.isdir(dir_local_plot):
        print("creating plot directory at {}".format(dir_local_plot))
        os.mkdir(dir_local_plot)

    # #regenerate default config file if needed:
    # if not os.path.isfile(os.path.join(dir_local_plot, default_config_fname)):
    #     regenerate_template_config()

    #read config:
    print("reading config file at {}".format(fname_config))
    try:
        with open(fname_config) as fi:
            config = json.loads(fi.read())
    except:
        print("error: could not load configuration file {}".format(fname_config))
        return 0
    print("","got config file")

    key_dynamic = 'use_dynamic_solution'
    dynamic = check_key(config, key_dynamic, None)
    if type(dynamic) ==type(None):
        print("error: expected top level key '{}' was not found in {}".format(key_dynamic, fname_config))
        return 0
    elif dynamic:
        print("","using dynamic solution")
    else:
        print("","using final static solution")
    print()

    fname_cdf = identify_solution(dir_solution, dynamic, overwrite_cdf = overwrite_cdf)
    # if type(dir_solution) == type(None):
    #     print("error: could not find solution {}".format(solution_name))
    #     return 0

    #derive a solution name from the solution directory:
    solution_name = os.path.splitext(os.path.basename(dir_solution))[0]

    wd_local_plots = mkdir_new_plot_subdir(dir_local_plot, solution_name)

    #we need to pass the following to each plotting function:
    # fname_cdf, wd_local_plots


    for function in functions_plotting:
        function_name = function.__name__

        #check the config file for the plot function name as a top level key:
        function_config = check_key(config, function_name, return_if_absent = "")

        if function_config == "":
            print("plot function key '{}' not found, skipping...".format(function_name))
            continue
        else:
            enabled = check_key(function_config, 'enabled', return_if_absent = False)
            if enabled:
                print("plot function key '{}' has 'enabled':True, expecting arguments:".format(function_name))
                argument_keys = function.__code__.co_varnames[:function.__code__.co_argcount]
                

                #try to read each function argument from the configuration file, set the value = None if not found:
                attempt_plot = True
                argument_set = []
                for argument in argument_keys:
                    #some arguments are not provided in the configuration file:
                    if argument == 'fname_cdf':
                        argument_value = fname_cdf
                    elif argument == 'wd_local_plots':
                        argument_value = wd_local_plots
                    elif argument == 'save':
                        argument_value = True
                    else:
                        if argument == 'epochrange' and len(epochrange_focus) == 2:
                            #if the user specified some epochs to look at, ignore configuration file options
                            print("using user-specified epochs")
                            argument_value = epochrange_focus
                        else:
                            argument_value = check_key(function_config, argument, return_if_absent = None)
                            if type(argument_value) == type(None):
                                print("error: argument '{}' not provided, skipping this plot".format(argument))
                                attempt_plot = False
                            print("","{} = {}".format(argument, argument_value))
                    argument_set.append(argument_value)

                #call the plotting function:
                if attempt_plot:
                    function(*argument_set)
                    plt.close('all') #close plots so as not to use up too much memory
                    print()
            else:
                print("plot function key '{}' has 'enabled':False (or absent), skipping...".format(function_name))
                print()
                continue
        
    try:
        plt.close('all')
    except:
        pass

    return 1





# #parse arguments:
# parser = argparse.ArgumentParser(description='')
# parser.add_argument('--config', type=str, default = default_config_fname, required = False)
# parser.add_argument('--solution', type=str, default = "", required = False)

# args = parser.parse_args()
# config_file = args.config
# solution_name = args.solution

# plot(solution_name, config_file)
