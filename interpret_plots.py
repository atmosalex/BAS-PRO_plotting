import numpy as np
import interpret_tools
import interpret_cdf
from spacepy import pycdf
pycdf.lib.set_backward(False)
from datetime import datetime
from scipy import optimize as opti
import matplotlib
import matplotlib.pyplot as plt
from math import sqrt, ceil
from matplotlib.ticker import (MultipleLocator)
import matplotlib.patheffects as path_effects
from scipy.interpolate import RegularGridInterpolator
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.dates as mdates
from math import log10, inf, sin
from mpl_toolkits.axes_grid1.inset_locator import mark_inset, inset_axes
import matplotlib.patheffects as pe
from matplotlib.dates import AutoDateFormatter, AutoDateLocator

import geopandas as gpd
import spacepy.time as spt ###
import spacepy.coordinates as spc ###
import spacepy.irbempy as ib ###
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.dates as mdates
import matplotlib.patheffects as pe
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter, AutoMinorLocator)

# font sizes and styles:--------------------------------+
titleSize = 16
yaxisSize = 16
xaxisSize = 16
ytickSize = 14
xtickSize = 14
legSize = 14
labSize = 15

# set the default font family like so:
#matplotlib.rc('font', family='Arial')

# we can also define a new font
import matplotlib.font_manager as font_manager

#cmusfont = {'fontname': 'Arial'}
#cmusfontFM = font_manager.FontProperties(family='Arial',
#                                         style='normal',
#                                         size=11)
years = mdates.YearLocator()  # every year
months = mdates.MonthLocator()  # every month
bimonthly = mdates.MonthLocator(interval=2)
years_fmt = mdates.DateFormatter('%Y/%m')
# months_fmt = mdates.DateFormatter('%m')
# ----------------------------------------------------END


    
smallpve = 1e-20

def make_overviewplot(fname_cdf, Llims, epochrange, energies, aeq, jlimits, anisotropy):
    overviewplot_cmap = 'plasma'
    overviewplot_cbad = 'black'

    xtick_locator = AutoDateLocator()
    xtick_formatter = AutoDateFormatter(xtick_locator)

    nplots = len(energies)

    cdf = pycdf.CDF(fname_cdf)

    anisotropy_denominator = 60 #ratio is 90/this
    if anisotropy:
        overviewplot_cmap = 'viridis'
        aeq = 90
        jlimits = [[0.9, 1000] for x in range(nplots)]
        print("","plotting anisotropy as the 90:{}d flux ratio".format(anisotropy_denominator))
    


    dynamic = cdf.attrs[interpret_cdf.lab_dynamic][0]
    if dynamic > 0:
        pass
    else:
        print("", "warning: solution is not dynamic, skipping...")
        return 0

    #which times to plot (for a dynamic run):
    ax_t = cdf[interpret_cdf.lab_axt]
    #ax_mu = cdf[interpret_cdf.lab_axmu]
    #ax_K = cdf[interpret_cdf.lab_axK]
    ax_L = cdf[interpret_cdf.lab_axL]
    map_alpha = cdf[interpret_cdf.lab_map]
    #t_plot = np.array(ax_t)#np.linspace(ax_t[-1], ax_t[0], nt)
    print("", "{}r x 1c plot".format(nplots))

    if not len(Llims) == 2:
        Llims = [ax_L[0],ax_L[-1]]
    if type(jlimits) == type(None):
        jlimits = []


    if len(epochrange) == 2:
        date_start_ts = epochrange[0]
        date_finish_ts = epochrange[1]
    else:
        date_start_ts = ax_t[0]
        date_finish_ts = ax_t[-1]
    datelims = [datetime.fromtimestamp(date_start_ts), datetime.fromtimestamp(date_finish_ts)]

    fig = plt.figure()
    overviewplot_figsize_x = 5.5
    fig_toppad = 0.1
    fig_bottompad = 0.1
    fig_sidepad = 0.1
    fig_axhpad = 0.22

    axh = 1.2*(Llims[1]-Llims[0])/(4-2) #inches
    hspace_const = 0.1 #inches
    hspace = hspace_const/axh
    #percentage of figure taken by vertical padding:
    pcvpad = fig_toppad + fig_bottompad
    
    overviewplot_figsize_y = (nplots*axh + (nplots-1)*(hspace_const))* (1 + pcvpad)


    cmap = plt.cm.get_cmap(overviewplot_cmap)
    cmap.set_bad(color=overviewplot_cbad)

    data_L_all = []
    data_epoch_all = []
    data_flux_all = []

    axs = []
    for pidx in range(nplots):
        if pidx == 0:
            ax = fig.add_subplot(nplots, 1, pidx+1)
        else:
            ax = fig.add_subplot(nplots, 1, pidx+1, sharex=ax)
        axs.append(ax)



        divider = make_axes_locatable(ax)
        cax = divider.append_axes('right', size='5%', pad=0.075)


        if pidx == 0:
            if anisotropy:
                cax.set_title("Ratio of flux $j_{90}$ : $j_{60}$",
                    fontdict={'horizontalalignment': 'right'})
            else:
                cax.set_title("Flux j [cm-2 s-1 MeV-1 sr-1]",
                    fontdict={'horizontalalignment': 'right'})

        Eval = energies[pidx]

        Evalstr = "{:.2f}".format(Eval)
        elabel = Evalstr + "MeV"
        
        #print("","plotting {:.2f}MeV".format(Eval))
        flux_max = 0
        flux_min = float("inf")


        data_epoch = [] #x
        data_L = [] #y
        data_flux = [] #z
        ax_plot_L = []
        ax_plot_t = []

        print ("", "processing E = {:.2f}MeV...".format(Eval))
        for idx_L, Lval in enumerate(ax_L):
            if Lval < Llims[0] or Lval > Llims[1]: continue
            ax_plot_L.append(Lval)
            ax_plot_t = []
            #print ("", "processing L = {:.2f}".format(Lval))

            #data_L_L = []
            #data_epoch_L = []
            #data_flux_L = []

            for idx_t, time_plot in enumerate(ax_t):
                if time_plot < date_start_ts or time_plot > date_finish_ts: continue
                ax_plot_t.append(time_plot)

                #take the mu-K slice of the solution cdf at current epoch:
                sol_f = cdf[interpret_cdf.lab_f][idx_t, :, :, idx_L]
                sol_en = cdf[interpret_cdf.lab_en][0, :, :, idx_L]

                #find the surrounding K indicies for the current L, pitch angle:
                iK0 = np.argmax(aeq >= map_alpha[0,idx_L, :])
                aeq_iK0 = map_alpha[0,idx_L, iK0]
                if aeq == aeq_iK0: #for 90 deg
                    iK1 = iK0
                else:
                    iK1 = iK0 - 1
                aeq_iK1 = map_alpha[0,idx_L, iK1]

                if aeq_iK0 < 0:
                    #print("Warning: {:.2f}d is in numerical loss cone at L = {:.2f}".format(aeq, Lval))
                    plot_j = smallpve
                else:

                    if Eval < sol_en[0, iK0] or Eval > sol_en[-1, iK0] or Eval < sol_en[0, iK1] or Eval > sol_en[-1, iK1]:
                        print("error: energy {:.2f}MeV out of solution range between iKs {} and {}".format(Eval, iK0, iK1))
                        return None, None, None, None
                    f_iK0 = np.interp(Eval, sol_en[:, iK0], sol_f[:, iK0])
                    f_iK1 = np.interp(Eval, sol_en[:, iK1], sol_f[:, iK1])
                    plot_f = np.interp(aeq, [aeq_iK0, aeq_iK1], [f_iK0, f_iK1])
                    plot_j = max(smallpve, interpret_tools.f2j(Eval, plot_f))


                    if anisotropy: #then the above will correspond to aeq=90

                        #repeat to store the ratio 90/anisotropy_denominator:
                        iK0 = np.argmax(anisotropy_denominator >= map_alpha[0,idx_L, :])
                        aeq_iK0 = map_alpha[0,idx_L, iK0]
                        if anisotropy_denominator == aeq_iK0: #for 90 deg
                            iK1 = iK0
                        else:
                            iK1 = iK0 - 1
                        aeq_iK1 = map_alpha[0,idx_L, iK1]

                        if aeq_iK0 < 0:
                            #print("warning: {:.2f}d is in numerical loss cone at L = {:.2f}".format(anisotropy_denominator, Lval))
                            plot_j = np.nan
                        else:  

                            if Eval < sol_en[0, iK0] or Eval > sol_en[-1, iK0] or Eval < sol_en[0, iK1] or Eval > sol_en[-1, iK1]:
                                print("error: energy {:.2f}MeV out of solution range between iKs {} and {}".format(Eval, iK0, iK1))
                                return None, None, None, None
                            f_iK0 = np.interp(Eval, sol_en[:, iK0], sol_f[:, iK0])
                            f_iK1 = np.interp(Eval, sol_en[:, iK1], sol_f[:, iK1])
                            plot_f_45 = np.interp(anisotropy_denominator, [aeq_iK0, aeq_iK1], [f_iK0, f_iK1])
                            plot_j_45 = interpret_tools.f2j(Eval, plot_f_45)
                            if plot_j_45 >0: 
                                plot_j = plot_j/ plot_j_45
                            else:
                                plot_j = np.nan




                data_epoch.append(time_plot) #x
                data_L.append(Lval) #y
                data_flux.append(plot_j)

            # data_L.append(data_L_L)
            # data_epoch.append(data_epoch_L)
            # data_flux.append(data_flux_L)
        
        data_epoch = np.array(data_epoch)
        data_L = np.array(data_L)
        data_flux = np.array(data_flux)

        #set the limits:
        flux_min = min([np.min(x) for x in data_flux])
        flux_max = max([np.max(x) for x in data_flux])

        if flux_max < smallpve:
            flux_min = flux_max
        else:
            flux_min = max(flux_min, 10**(log10(flux_max)-7.5)) #make sure min flux is not less than flux_max - some limit

        
        if len(jlimits):
            if jlimits[pidx][0] > 0 and jlimits[pidx][1] > 0:
                flux_max = jlimits[pidx][1]
                flux_min = jlimits[pidx][0]
                #flux_min = 10**(log10(flux_max)-jlimits[pidx][0]) #jlimits[pidx][0]#
        print("","model solution flux limits: ({}, {})".format(flux_min, flux_max))
        #print("","newly applied limits: ({}, {})".format(flux_min, flux_max))
        if flux_min == flux_max:
            flux_min = flux_max+1 #so that flux < flux_min ==> cmap(0)


        # #draw:
        # Lbefore = 0
        # idx_L_start = 0
        # for idx_L, Lval in enumerate(ax_L):
        #     if Lval < Llims[0]:
        #         idx_L_start += 1
        #         continue
        #     elif Lval > Llims[1]:
        #         continue
        #     #CHANGE THE NEXT LINE TO USE A CUSTOM TIME ARRAY, SO THAT NOT EVERY TIMESTEP IS PLOTTED:
        #     for idx_t, ep in enumerate(data_epoch[idx_L-idx_L_start]):
                

        #         flux = data_flux[idx_L-idx_L_start][idx_t]
        #         if flux < flux_min:
        #             c = cmap(0)
        #         else:
        #             c = cmap((log10(flux)-log10(flux_min))/(log10(flux_max) - log10(flux_min)))

        #         ymin_vl = (Lbefore-Llims[0])/(Llims[1]-Llims[0])
        #         ymax_vl = (Lval-Llims[0])/(Llims[1]-Llims[0])

        #         #print(c, time_plot, ymin_vl, ymax_vl)

        #         ax.axvline(ep, solid_capstyle='butt',ymin=ymin_vl, ymax=ymax_vl, color = c,alpha=1)#,transform=ax.transAxes) 
        #         #ax.plot((time_plot,time_plot), (ymin_vl, ymax_vl), color = c,alpha=1)

        #     Lbefore = Lval



         

        #interpolate onto a regular grid and plot using imshow
        ext_t = [min(data_epoch), max(data_epoch)]
        ext_L = [ax_plot_L[0], ax_plot_L[-1]]

        grid_t = np.linspace(ext_t[0],ext_t[1],201)
        grid_L = np.linspace(ext_L[0],ext_L[1],len(ax_plot_L))
        

        data = np.log10(data_flux).reshape((len(ax_plot_L),len(ax_plot_t)))
        interp = RegularGridInterpolator((ax_plot_L, ax_plot_t), data)
        ee, LL = np.meshgrid(grid_t, grid_L, indexing='ij')

        flux_reg = np.zeros(np.shape(ee))
        for idx_L in range(len(grid_L)):
            for idx_t in range(len(grid_t)):
                flux_reg[idx_t][idx_L] = 10**interp((grid_L[idx_L], grid_t[idx_t]))


        norm = matplotlib.colors.LogNorm(vmin=flux_min,vmax=flux_max)
        #norm = matplotlib.colors.Normalize(vmin=log10(flux_min),vmax=log10(flux_max))

        x_lims_dt = [datetime.fromtimestamp(x) for x in ext_t]
        x_lims = mdates.date2num(x_lims_dt)
        imshow = ax.imshow(flux_reg.T, aspect='auto', extent=(x_lims[0],x_lims[1],ext_L[0],ext_L[1]), origin='lower', cmap=cmap, norm = norm)
        ax.xaxis_date()

        #format y axis labels:
        ax.yaxis.set_major_locator(MultipleLocator(0.5))

        # #format x axis labels:
        # if (ax_plot_t[-1] - ax_plot_t[0]) >= 370*24*3600:
        #     ax.xaxis.set_major_locator(mdates.YearLocator())
        #     ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y'))
        #     ax.xaxis.set_minor_locator(mdates.MonthLocator())
        # elif (ax_plot_t[-1] - ax_plot_t[0]) >= 31*24*3600:
        #     ax.xaxis.set_major_locator(mdates.MonthLocator())
        #     ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m'))
        #     ax.xaxis.set_minor_locator(mdates.DayLocator(interval = 10))
        # else:
        #     ax.xaxis.set_major_locator(mdates.DayLocator())
        #     ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d'))
        #     ax.xaxis.set_minor_locator(mdates.HourLocator(interval = 6))

        ax.xaxis.set_major_locator(xtick_locator)
        ax.xaxis.set_major_formatter(xtick_formatter)
        
        fig.autofmt_xdate()

        ax.set_xlim(datelims)
        ax.tick_params(axis='x', rotation=45)
        
        #norm = matplotlib.colors.Normalize(vmin=log10(flux_min),vmax=log10(flux_max))
        #plt.cm.ScalarMappable(cmap=cmap, norm=norm)
        cbar = fig.colorbar(imshow,cax=cax)


        ax.tick_params(which='both', right=True, top=True, labelright=False, labeltop=False, labelbottom=False)
        


        ax.xaxis.set_tick_params(which='both', width=1.15)
        ax.yaxis.set_tick_params(which='both', width=1.15)

        ax.text(0.97,0.01, elabel, va="bottom", ha="right", transform=ax.transAxes, fontsize = 14, color="white", path_effects=[pe.withStroke(linewidth=4, foreground="black")])
        #ax.text(0.97,0.01, elabel, va="bottom", ha="right", transform=ax.transAxes, fontsize = 14, color="black")
        #ax.text(0.97,0.01, elabel, va="bottom", ha="right", transform=ax.transAxes)
        ax.set_facecolor((0, 0, 0))


        data_L_all.append(data_L)
        data_epoch_all.append(data_epoch)
        data_flux_all.append(data_flux)


    ax.set_ylabel("L")
    ax.tick_params(labelbottom=True)

    fig.set_size_inches(overviewplot_figsize_x, overviewplot_figsize_y)
    plt.subplots_adjust(left=fig_sidepad, right=1-fig_sidepad, top=1-fig_toppad, bottom=fig_bottompad, hspace = hspace, wspace = 0)
    

    data_L_all = np.array(data_L_all)
    data_epoch_all = np.array(data_epoch_all)
    data_flux_all = np.array(data_flux_all)

    return fig, data_L_all, data_epoch_all, data_flux_all
    #plt.close()




def plot_spectrum(fname_cdf, iK, L_, nt, epochs, universal_axes_limits):

    cdf = pycdf.CDF(fname_cdf)
    dynamic = cdf.attrs[interpret_cdf.lab_dynamic][0]
    if dynamic > 0: 
        dynamic = True
    else:
        dynamic = False

    #which times to plot:
    ax_t = cdf[interpret_cdf.lab_axt][:]
    if dynamic:
        if len(epochs):
            t_plot = np.array(epochs)
            nt = len(epochs)
            if t_plot[-1] > ax_t[-1] or t_plot[0] < ax_t[0]:
                print("error: plot time out of range for solution axis")
                return None, None, 0
        else:
            t_plot = np.linspace(ax_t[-1], ax_t[0], nt)
    else:
        t_plot = [ax_t[0]]
        nt = 1


    #print("", "1r x 1c plot")

    #set up plot:
    fig, ax = plt.subplots(1, 1)
    usegrid = True
    cmap = matplotlib.cm.get_cmap('viridis')

    ax_mu = cdf[interpret_cdf.lab_axmu]
    ax_K = cdf[interpret_cdf.lab_axK]
    ax_L = cdf[interpret_cdf.lab_axL]
    map_alpha = cdf[interpret_cdf.lab_map]

    if L_< ax_L[0] or L_ > ax_L[-1]:
        print("Error: L = {:.2f} is out of range".format(L_))
        return None, None, 0
    #------------------------------------------------------------------------------+
    # make plot                                                                    |
    #------------------------------------------------------------------------------+
    ax.set_xscale('log')
    ax.set_yscale('log')

    idx_L = 0
    while ax_L[idx_L] < L_:
        idx_L += 1
    idx_L_1 = idx_L
    idx_L_0 = idx_L_1 - 1
    frac_L = (L_ - ax_L[idx_L_0])/(ax_L[idx_L_1] - ax_L[idx_L_0])
    if frac_L == 1.: idx_L_0 = idx_L_1

    if idx_L_0 < 0:
        print("Warning: L = {:.2f} out of range".format(L_))
        return None, None, 0

    for idx_K in range(np.size(ax_K)):
        map_alpha = cdf[interpret_cdf.lab_map]
        if map_alpha[0, idx_L_0, idx_K] <= 0: break
    idx_K_max = idx_K - 1
    if iK > idx_K_max:
        print("Warning: iK = {} out of range for L = {:.2f}".format(iK, L_))
        return None, None, 0

    aeq = (1 - frac_L) * map_alpha[0,idx_L_0, iK] + frac_L * map_alpha[0,idx_L_1, iK]
    #enonzero = cdf[interpret_cdf.lab_en][0, :, :, :] > 0
    #enonzero = cdf[interpret_cdf.lab_en][0, :, :, :][enonzero]
    emin = 1#np.min(enonzero)
    emax = 140#np.max(enonzero)

    #find axes limits in terms of flux:
    jmin = inf#cdf[interpret_cdf.lab_f][0,0,0,-1]
    jmax = 0
    for idx_t in range(len(ax_t)):
        fnonzero_idx = cdf[interpret_cdf.lab_f][idx_t, :, :, :] > 0
        fnonzero = cdf[interpret_cdf.lab_f][idx_t, :, :, :][fnonzero_idx]
        e_correspnging = cdf[interpret_cdf.lab_en][0, :, :, :][fnonzero_idx]

        fbelowemax = fnonzero[e_correspnging <= emax]
        ebelowemax = e_correspnging[e_correspnging <= emax]

        faboveemin_belowemax = fbelowemax[ebelowemax >= emin]
        eaboveemin_belowemax = ebelowemax[ebelowemax >= emin]

        if not len(faboveemin_belowemax):
            continue

        fmin_idx = np.argmin(faboveemin_belowemax)
        fmin = faboveemin_belowemax[fmin_idx]
        e_fmin = eaboveemin_belowemax[fmin_idx]

        fmax_idx = np.argmax(faboveemin_belowemax)
        fmax = faboveemin_belowemax[fmax_idx]
        e_fmax = eaboveemin_belowemax[fmax_idx]
        
        jmin_t = interpret_tools.f2j(e_fmin, fmin)
        jmax_t = interpret_tools.f2j(e_fmax, fmax)
        if jmin_t==jmin_t and jmin_t < jmin:
            jmin = jmin_t
        if jmax_t==jmax_t and jmax_t > jmax and jmax_t < np.Inf:
            jmax = jmax_t


    #find the x axis (energy):
    spec_en_L0 = cdf[interpret_cdf.lab_en][0, :, iK, idx_L_0]
    spec_en_L1 = cdf[interpret_cdf.lab_en][0, :, iK, idx_L_1]
    spec_en = (1 - frac_L) * spec_en_L0 + frac_L * spec_en_L1

    abovezero = False #for checking we don't output the f = 0 boundary condition in the loss cone

    data = []

    #interpolate to the required time
    for time_plot in t_plot:
        #get idx_t_0 and idx_t_1 surrounding the time we want to plot:
        idx_t_1 = 0
        while ax_t[idx_t_1] < time_plot:
            idx_t_1 += 1
        
        if not dynamic:
            idx_t_0 = idx_t_1
            frac_t = 1.
        else:
            idx_t_0 = idx_t_1 - 1
            frac_t = (time_plot - ax_t[idx_t_0])/(ax_t[idx_t_1] - ax_t[idx_t_0])

        spec_f_t0_L0 = cdf[interpret_cdf.lab_f][idx_t_0, :, iK, idx_L_0]
        spec_f_t1_L0 = cdf[interpret_cdf.lab_f][idx_t_1, :, iK, idx_L_0]
        spec_f_t0_L1 = cdf[interpret_cdf.lab_f][idx_t_0, :, iK, idx_L_1]
        spec_f_t1_L1 = cdf[interpret_cdf.lab_f][idx_t_1, :, iK, idx_L_1]

        spec_f_t0 = (1 - frac_L) * spec_f_t0_L0 + frac_L * spec_f_t0_L1
        spec_f_t1 = (1 - frac_L) * spec_f_t1_L0 + frac_L * spec_f_t1_L1

        spec_f = (1 - frac_t) * spec_f_t0 + frac_t * spec_f_t1
        spec_j = interpret_tools.f2j(spec_en, spec_f)

        colour = interpret_tools.get_colour_from_time(time_plot, ax_t, cmap)
        if not dynamic: color = "black"

        if np.sum(spec_j>0): abovezero = True
        ax.plot(spec_en, spec_j, color=colour, linewidth=0.8, alpha=1)
        data.append([spec_en, spec_j])

    #label K
    ax.text(0.05,0.05,"$K=$ " + "{:.3f}".format(ax_K[iK]) + "$G^{0.5} R_{E}$\n" + "$\\alpha_{\\mathrm{eq}}=$" + "{:.2f}".format(aeq) + "$^{\\circ}$\n" + "$L=$" + "{:.2f}".format(L_),rotation=0,
        color='black', size=9, ha="left", va="bottom", transform=ax.transAxes)


    if dynamic:
        #label each dynamic simulation:
        shadow = False
        textax = ax
        n = 1
        for idx_t, time_plot in enumerate(t_plot):
            colour = interpret_tools.get_colour_from_time(time_plot, ax_t, cmap)
            time_dt = datetime.fromtimestamp(time_plot)
            if len(t_plot) > n*12:
                labelheight = n*(idx_t * 1.0/(len(t_plot)-1))
            else:
                labelheight = n*(idx_t * 1.0/(n*12))

            text = textax.text(1.1, labelheight, time_dt.strftime("%j/%Y"), transform=textax.transAxes, color=colour,
                va='center',fontdict={'fontsize': yaxisSize-4})

            if shadow:
                text.set_path_effects([path_effects.Stroke(linewidth=0.8, foreground='black'),
                                       path_effects.Normal()])




    if usegrid:
        ax.grid(which='minor', alpha=0.2)
        ax.grid(which='major', alpha=0.5)

    
    ax.yaxis.set_ticks_position('both')
    ax.xaxis.set_tick_params(direction='out', length=5, labelsize=xtickSize)#, which='bottom'
    ax.yaxis.set_tick_params(direction='out', length=2, labelsize=ytickSize)#, which='bottom'
    ax.tick_params(labelsize=xtickSize)
    ax.set_ylabel('$j$ [cm$^{-2}$s$^{-1}$str$^{-1}$MeV$^{-1}$]', fontdict={'fontsize': yaxisSize-2})
    ax.set_xlabel('$E$ [MeV]', fontdict={'fontsize': yaxisSize})
    jmin = max(jmin, 10**(log10(jmax)-9.5))
    if universal_axes_limits:
        ax.set_xlim((emin, emax))
        ax.set_ylim((jmin, jmax))
    #ax.set_ylim((1e-5, 2e7)) ###########################################
    plt.tight_layout()

    if not abovezero:
        print("","data out of range")
        return fig, None, 0

    return fig, data, 1


def plot_padist_panels(fname_cdf, energies, Lshells, nt, epochs, plotflux):
    plotshape = np.zeros((len(energies),len(Lshells)))
    plotenergies = np.array([energies]).T
    plotenergies = np.repeat(plotenergies, len(Lshells), axis = 1)
    plotLshells = np.array([Lshells])
    plotLshells = np.repeat(plotLshells, len(energies), axis = 0)

    cdf_loaded = interpret_cdf.Preloaded(fname_cdf)
    dynamic = cdf_loaded.dynamic
    if dynamic > 0: #boolean is not preserved by CDF format (?)
        dynamic = True
    else:
        dynamic = False


    #which times to plot:
    ax_t = cdf_loaded.ax_t
    if dynamic:
        if len(epochs):
            t_plot = np.array(epochs)
            nt = len(epochs)
            if t_plot[-1] > ax_t[-1] or t_plot[0] < ax_t[0]:
                print("error: plot time out of range for solution axis")
                return None, None, None, None
        else:
            t_plot = np.linspace(ax_t[-1], ax_t[0], nt)
    else:
        t_plot = [ax_t[0]]
        nt = 1



    n = np.shape(plotshape)[0]
    m = np.shape(plotshape)[1]
    print("", "{}r x {}c plot".format(n, m))

    #fig, ax_array_all = plt.subplots(n, m, sharex=True)#, gridspec_kw={'width_ratios': [0.6]*m})
    fig, ax_array_all = plt.subplots(n, m+(m-1), sharex=True, sharey=False,
                                 gridspec_kw={'width_ratios': [1]+[0.075, 1]*(m-1)})
    cmap = matplotlib.cm.get_cmap('viridis')
    usegrid = True
    addlabel = True

    #arrange the axes:
    ax_array = []
    if n == 1:
        ax_array_all = [ax_array_all]
    for row in ax_array_all:
        ax_row = []
        if m == 1:
            ax_row.append(row)
        else:
            for count, col in enumerate(row):
                if (count % 2 == 0):
                    ax_row.append(col)
                else:
                    col.axis('off')
        ax_array.append(ax_row)


    aeq_data_all = [[[] for ic in range(len(ax_row))] for ax_row in ax_array]
    sodata_L_all = [[[] for ic in range(len(ax_row))] for ax_row in ax_array]
    data_epoch_all = [[t_plot for ic in range(len(ax_row))] for ax_row in ax_array]

    for rowidx, ax_row in enumerate(ax_array):
        print("","starting row #{}".format(rowidx+1))
        for colidx, ax_col in enumerate(ax_row):
            print("","starting col #{}".format(colidx+1))
            en = plotenergies[rowidx][colidx]
            L_ = plotLshells[rowidx][colidx]

            ax_target = ax_col

            #interpolate to the required time
            for time_plot in t_plot:

                #get PAD:
                plot_alpha, plot_f = interpret_tools.getpad_(cdf_loaded, L_, en, dynamic, time_plot)

                plot_j = interpret_tools.f2j(en, plot_f)

                colour = interpret_tools.get_colour_from_time(time_plot, ax_t, cmap)
                # if not dynamic: #don't want scatter marks in the dynamic plots
                #     ax_target.scatter(plot_alpha, plot_j, color=colour,# linewidth=1.4,
                #                label=str(round(en,1))+" MeV", alpha=0.75,marker='x',s=15)
                
                plot_alpha = plot_alpha[plot_f>=0]
                plot_j = plot_j[plot_f>=0]
                plot_f = plot_f[plot_f>=0]

                if not len(plot_f):
                    continue

                if plotflux:
                    ax_target.plot(plot_alpha, plot_j, color=colour, linewidth=1.4, label="{:.1f} MeV".format(en), alpha=1)
                    sodata_L_all[rowidx][colidx].append(plot_j)
                    miny = min(plot_j)
                else:
                    ax_target.plot(plot_alpha, plot_f, color=colour, linewidth=1.4, label="{:.1f} MeV".format(en), alpha=1)
                    sodata_L_all[rowidx][colidx].append(plot_f)
                    miny = min(plot_f)

                aeq_data_all[rowidx][colidx].append(plot_alpha)
                

            x1 = ax_target.get_xlim()[0]
            x2 = ax_target.get_xlim()[1]
            y1 = ax_target.get_ylim()[0]
            y2 = ax_target.get_ylim()[1]

            ax_target.set_ylim([0, y2])


            start, end = ax_target.get_xlim()
            da = 20
            xtickarray = np.arange(end, start - da, -da)
            ax_target.xaxis.set_ticks(xtickarray)
            #x axis tidy up:
            #ax_target.set_xlim([interpret_tools.get_lc(Lshells[-1]), 90])
            ax_target.set_xlim([0, 90])
            
            ax_target.yaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter(useMathText=True, useOffset=False))
            ax_target.yaxis.get_offset_text()
            exponent = ax_target.yaxis.get_offset_text()
            exponent.set_fontsize(14)
            ax_target.ticklabel_format(axis='y', style='sci', scilimits=(0,0))

            usegrid = True
            if usegrid:
                ax_target.grid(linestyle='--', color='grey')


            # axis ticks:
            ax_target.yaxis.set_ticks_position('default')
            ax_target.yaxis.get_offset_text().set_fontsize(xtickSize-2)
            #ax_target.xaxis.set_tick_params(direction='out', length=5, labelsize=xtickSize-2)  # , which='bottom'
            ax_target.yaxis.set_tick_params(direction='out', length=2, labelsize=ytickSize-2, right = False)  # , which='bottom'
            ax_target.tick_params(labelsize=xtickSize-2)
            ax_target.set_ylabel('')


            if addlabel:
                ax_row[colidx].text(0.05,0.85,"L = {:.2f}, {:.2f}MeV".format(L_, en), transform=ax_row[colidx].transAxes) 

            #ax_row[colidx].yaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter(useMathText=True, useOffset=False))
        #ax_row[0].set_ylabel('f [km$^{-6}$s$^{3}$]', fontdict={'fontsize': yaxisSize})


    # Additional plot information:
    #aeq on x axes:
    for colidx, ax_col in enumerate(ax_array[-1]):
        ax_col.set_xlabel('$\\alpha_{eq} [^{\\circ}]$', fontdict={'fontsize': xaxisSize})


    #E labels:
    for rowidx, ax_row in enumerate(ax_array):
        en = plotenergies[rowidx][0]
        if plotflux:
            qlabel = 'j [cm$^{-2}$s$^{-1}$str$^{-1}$MeV$^{-1}$]'
        else:
            qlabel = 'f [km-6 s-1]'
        ax_row[0].text(-0.45, 0.5, "{} MeV".format(en) + '\n'+qlabel, transform=ax_row[0].transAxes,
                va='center',ha='center',fontdict={'fontsize': yaxisSize}, rotation=90)
    #L labels:
    for colidx, ax_col in enumerate(ax_array[0]):
        L_ = plotLshells[0][colidx]
        ax_col.text(0.5, 1.22, "L = {}".format(L_), transform=ax_col.transAxes,
                va='center',ha='center',fontdict={'fontsize': yaxisSize}, rotation=0)


    if dynamic:
        #label each dynamic simulation:
        shadow = False
        textax = ax_array[-1][-1]
        for idx_t, time_plot in enumerate(t_plot):
            colour = interpret_tools.get_colour_from_time(time_plot, ax_t, cmap)
            time_dt = datetime.fromtimestamp(time_plot)
            if len(t_plot) > n*12:
                labelheight = n*(idx_t * 1.0/(len(t_plot)-1))
            else:
                labelheight = n*(idx_t * 1.0/(n*12))

            text = textax.text(1.1, labelheight, time_dt.strftime("%j/%Y"), transform=textax.transAxes, color=colour,
                va='center',fontdict={'fontsize': yaxisSize-4})

            if shadow:
                text.set_path_effects([path_effects.Stroke(linewidth=0.8, foreground='black'),
                                       path_effects.Normal()])


    fig.set_size_inches(m*2.1, n*3.2)
    fig.set_size_inches(m*4.1, n*3.2)
    plt.tight_layout()
    #plt.savefig(fname_plot, dpi = 400)

    return fig, aeq_data_all, sodata_L_all, data_epoch_all


def plot_f_vs_L_panels(fname_cdf, z_, iKs, epochs, fixedenergy, nt, Llims):
    zoombox_choice = False
    plotshape = np.zeros((len(z_), len(z_[0])))
    plot_max_order_yrange = 3 #energy
    #if not fixedenergy: plot_max_order_yrange = 9
    if not fixedenergy: plot_max_order_yrange = 6 #psd

    cdf = pycdf.CDF(fname_cdf)

    dynamic = cdf.attrs[interpret_cdf.lab_dynamic][0]
    if dynamic > 0: #boolean is not preserved by CDF format (?)
        dynamic = True
    else:
        dynamic = False

    #which times to plot:
    ax_t = cdf[interpret_cdf.lab_axt][:]
    if dynamic:
        if len(epochs):
            t_plot = np.array(epochs)
            nt = len(epochs)
            if t_plot[-1] > ax_t[-1] or t_plot[0] < ax_t[0]:
                print("error: plot time out of range for solution axis")
                return None, None, None, None, None
        else:
            #t_ignore = 1391212800
            #print("Warning: ignoring times before t=",t_ignore)
            t_plot = np.linspace(ax_t[-1], ax_t[0], nt)
    else:
        t_plot = [ax_t[0]]
        nt = 1


    n = np.shape(plotshape)[0]
    m = np.shape(plotshape)[1]
    print("", "{}r x {}c plot".format(n, m))

    #set up plot:
    fig, ax_array_all = plt.subplots(n, m+(m-1), sharex=True, sharey=False,
                                 gridspec_kw={'width_ratios': [1]+[0.0, 1]*(m-1)})
    scale_axes_factor = 3.35
    scale_axes_factor = 5.3
    width = 0
    for axwidth in [1]+[0.0, 1]*(m-1): width += axwidth
    width = width * scale_axes_factor #scale factor
    fig.set_size_inches(width, scale_axes_factor*n/1.6)
    usegrid = True
    cmap = matplotlib.cm.get_cmap('viridis')
    if nt <= 2: cmap = matplotlib.cm.get_cmap('bwr')
    plot_energy_bars = True
    #legendon = False


    #arrange the axes:
    ax_array = []
    if n == 1:
        ax_array_all = [ax_array_all]

    for row in ax_array_all:
        ax_row = []
        if m == 1:
            ax_row.append(row)
        else:
            for count, col in enumerate(row):
                if (count % 2 == 0):
                    ax_row.append(col)
                else:
                    col.axis('off')
        ax_array.append(ax_row)


    ax_t = cdf[interpret_cdf.lab_axt]
    ax_mu = cdf[interpret_cdf.lab_axmu]
    ax_K = cdf[interpret_cdf.lab_axK]
    ax_L = cdf[interpret_cdf.lab_axL]
    map_alpha = cdf[interpret_cdf.lab_map]

    data_L_all = [[{} for ic in range(len(ax_row))] for ax_row in ax_array]
    aeq_data_all = [[{} for ic in range(len(ax_row))] for ax_row in ax_array]
    data_epoch_all = [[t_plot for ic in range(len(ax_row))] for ax_row in ax_array]
    sodata_L_all = [[{} for ic in range(len(ax_row))] for ax_row in ax_array]

    zoombox = False
    if not fixedenergy:
        x1_zoom = 1.15
        x2_zoom = 1.35
        zoombox = zoombox_choice


    ylim_max = 0
    for rowidx, ax_row in enumerate(ax_array):
        print("","starting row #{}".format(rowidx+1))
        for colidx, ax_col in enumerate(ax_row):
            ax_target = ax_col
            if z_[rowidx][colidx] == None:
                ax_target.axis('off')
                continue
            print("","starting col #{}".format(colidx+1))
            en = z_[rowidx][colidx]
            mu = en #could be mu or energy supplied
            


            #------------------------------------------------------------------------------+
            # make plot                                                                    |
            #------------------------------------------------------------------------------+
            minys = []
            maxys = []
            minys_zoom = []
            maxys_zoom = []
            text_labels = []

            ax_target.set_yscale('log')

            if zoombox:
                #zoom box:
                #ax_target_zoom = zoomed_inset_axes(ax_target, 2, loc=4, # zoom = 6
                #    bbox_to_anchor = [0.6, 0.075, 0.35, 0.8], bbox_transform =ax_target.transAxes)
                ax_target_zoom = inset_axes(ax_target, 1.4, 1.6, loc=4, # zoom = 6
                    bbox_to_anchor = [0.98, 0.075], bbox_transform =ax_target.transAxes)
                #mark_inset(ax_target, ax_target_zoom, loc1=2, loc2=4, fc="none", ec="red", alpha=0.5, lw=1, zorder = 1)
                mark_inset(ax_target, ax_target_zoom, loc1=2, loc2=4, fc="none", ec="black", alpha=0.5, lw=1, ls="solid")


            #     x1_zoom = 1.19
            #     x2_zoom = 1.35
            idx_K_max = -1
            for idx_K in iKs:
                if idx_K > len(ax_K) - 1:
                    continue
                K_now = ax_K[idx_K]
                sol_en = cdf[interpret_cdf.lab_en][0, :, idx_K, :]
                sol_f = cdf[interpret_cdf.lab_f][:, :, idx_K, :]

                idxLvf = 0
                while min(sol_en[:, idxLvf]) < 0:
                    idxLvf+=1
                    if idxLvf==len(ax_L):
                        break
                if idxLvf>=len(ax_L)-2:
                    print("","warning: cannot plot at K={:.2f}".format(K_now))
                    break
                idx_K_max = idx_K

                if (fixedenergy):
                    for idxL in range(idxLvf, len(ax_L)):
                        if (en > sol_en[-1, idxL] or en < sol_en[0, idxL]):
                            print("error: energy {}MeV is out of bounds".format(en))
                            return None, None, None, None, None


                data_L_all[rowidx][colidx][K_now] = []
                aeq_data_all[rowidx][colidx][K_now] = []
                sodata_L_all[rowidx][colidx][K_now] = []

                for time_plot in t_plot: #do it backward so we can stick our K label to t0 solution plotted

                    if dynamic: #interpolate to the current t we need:
                        #get idx_t_0 and idx_t_1 surrounding the time we want to plot:
                        idx_t_0 = -1
                        idx_t_1 = idx_t_0
                        if time_plot < ax_t[0] or time_plot > ax_t[-1]:
                            print("error: time_plot is out of range on K idx", idx_K)
                            return None, None, None, None, None

                        for idx_t, time_sol in enumerate(ax_t):
                            if time_plot >= time_sol:
                                idx_t_0 = idx_t
                                if time_plot == time_sol:
                                    idx_t_1 = idx_t_0 #we 'interpolate' across 1 grid
                                else:
                                    idx_t_1 = idx_t + 1
                            else:
                                break

                        sol_f1d_t_0 = []
                        sol_f1d_t_1 = []
                        if (not fixedenergy):
                            sol_en1d = [] #used for displaying energy bars when inspectmu == True
                            #get f at every L at the energy under investigation:
                            # t0
                            for idxL in range(idxLvf, len(ax_L)):
                                #sol_f1d_t_0.append(float(np.interp(np.log10(mu),ax_mu[:], sol_f[idx_t_0*len(ax_mu):(1+idx_t_0)*len(ax_mu), idxL])))
                                sol_f1d_t_0.append(np.interp(np.log10(mu), ax_mu[:], sol_f[idx_t_0, :, idxL]))
                                sol_en1d.append(np.interp(np.log10(mu), ax_mu[:], sol_en[:, idxL]))
                            

                            # t1
                            for idxL in range(idxLvf, len(ax_L)):
                                #sol_f1d_t_1.append(float(np.interp(np.log10(mu),ax_mu[:], sol_f[idx_t_1*len(ax_mu):(1+idx_t_1)*len(ax_mu), idxL])))  
                                #sol_en1d.append(float(np.interp(np.log10(mu), ax_mu[:], sol_en[:, idxL])))    
                                sol_f1d_t_1.append(np.interp(np.log10(mu), ax_mu[:], sol_f[idx_t_1, :, idxL]))

                        else:
                            #get f at every L at the energy under investigation:
                            # t0
                            for idxL in range(idxLvf, len(ax_L)):
                                #sol_f1d_t_0.append(float(np.interp(en, sol_en[:, idxL], sol_f[idx_t_0*len(ax_mu):(1+idx_t_0)*len(ax_mu), idxL])))
                                sol_f1d_t_0.append(np.interp(en, sol_en[:, idxL], sol_f[idx_t_0, :, idxL]))

                            # t1
                            for idxL in range(idxLvf, len(ax_L)):
                                #sol_f1d_t_1.append(float(np.interp(en, sol_en[:, idxL], sol_f[idx_t_1*len(ax_mu):(1+idx_t_1)*len(ax_mu), idxL])))
                                sol_f1d_t_1.append(np.interp(en, sol_en[:, idxL], sol_f[idx_t_1, :, idxL]))


                        sol_f1d = []
                        ax_t_surround = [ax_t[idx_t_0], ax_t[idx_t_1]]
                        for idxL in range(len(sol_f1d_t_0)): #interpolate f at each L from surrounding times
                            f_surround = [sol_f1d_t_0[idxL], sol_f1d_t_1[idxL]]
                            sol_f1d.append(np.interp(time_plot, ax_t_surround, f_surround))

                        
                    else:
                        #ask user whether to plot f vs. mu or energy:
                        sol_f1d = []
                        if (not fixedenergy):
                            sol_en1d = [] #used for displaying energy bars when inspectmu == True
                            for idxL in range(idxLvf,len(ax_L)):
                                #sol_f1d.append(float(np.interp(np.log10(mu),ax_mu[:],sol_f[:,idxL])))
                                sol_f1d.append(np.interp(np.log10(mu), ax_mu[:], sol_f[0, :, idxL]))
                                sol_en1d.append(np.interp(np.log10(mu), ax_mu[:], sol_en[:, idxL]))
                        else:
                            for idxL in range(idxLvf,len(ax_L)):
                                #sol_f1d.append(float(np.interp(en,sol_en[:,idxL],sol_f[:,idxL])))
                                sol_f1d.append(np.interp(en,sol_en[:,idxL],sol_f[0, :, idxL]))
                        

                    sol_f1d = np.array(sol_f1d)
                    sol_j1d = interpret_tools.f2j(en, sol_f1d)

                    #if not (labelarr):
                    if (not fixedenergy):
                        sol_plot = sol_f1d
                    else:
                        sol_plot = sol_j1d

                    ax_L_plot = ax_L[idxLvf:]
                    aeq_L_plot = map_alpha[0,:,idx_K][idxLvf:]

                    #filter negative values
                    ax_L_plot = ax_L_plot[sol_plot>0]
                    aeq_L_plot = aeq_L_plot[sol_plot>0]
                    sol_plot = sol_plot[sol_plot>0]

                    #filter nans
                    ax_L_plot = ax_L_plot[~np.isnan(sol_plot)]
                    aeq_L_plot = aeq_L_plot[~np.isnan(sol_plot)]
                    sol_plot = sol_plot[~np.isnan(sol_plot)]

                    if not len(sol_plot): continue
                    
                    data_L_all[rowidx][colidx][K_now].append(ax_L_plot)
                    aeq_data_all[rowidx][colidx][K_now].append(aeq_L_plot)
                    sodata_L_all[rowidx][colidx][K_now].append(sol_plot)



            for idx_K in iKs[:idx_K_max+1]:
                if idx_K > len(ax_K) - 1:
                    continue
                K_now = ax_K[idx_K]

                minys_K = []
                maxys_K = []
                minys_K_x = []
                maxys_K_x = []

                for idx_t in range(len(t_plot)): #do it backward so we can stick our K label to t0 solution plotted
                    time_plot = t_plot[idx_t]

                    if idx_t > len(data_L_all[rowidx][colidx][K_now])-1:
                        break

                    #colour = interpret_tools.get_colour_from_time(time_plot, ax_t, cmap)
                    colour = interpret_tools.get_colour_from_time(time_plot, t_plot, cmap)
       
                    if (not fixedenergy):
                        label = "{:.0f}MeV/G".format(mu)
                    else:
                        label = "{:.2f}MeV".format(en)

                    ax_L_plot = data_L_all[rowidx][colidx][K_now][idx_t]
                    aeq_L_plot = aeq_data_all[rowidx][colidx][K_now][idx_t]
                    sol_plot = sodata_L_all[rowidx][colidx][K_now][idx_t]

                    #save the y axis limits:
                    if len(Llims) ==2: #zoom in to a specific region:
                        #if the plottable part of the solution is outside this L range,
                        # do not let it affect the plot's limits
                        if ax_L_plot[0] <= Llims[1] and ax_L_plot[-1] >= Llims[0]: 
                            sol_plot_inyrange = sol_plot[ax_L_plot <= Llims[1]]
                            ax_L_inyrange = ax_L_plot[ax_L_plot <= Llims[1]]
                            sol_plot_inyrange = sol_plot_inyrange[ax_L_inyrange >= Llims[0]]
                            ax_L_inyrange = ax_L_inyrange[ax_L_inyrange >= Llims[0]]

                            minys_K.append(np.min(sol_plot_inyrange))
                            maxys_K.append(np.max(sol_plot_inyrange))
                            minys_K_x.append(ax_L_inyrange[np.argmin(minys_K[-1] < sol_plot_inyrange[sol_plot_inyrange>0])])
                            maxys_K_x.append(ax_L_inyrange[np.argmin(maxys_K[-1] > sol_plot_inyrange)])

                    else:
                        minys_K.append(np.min(sol_plot))
                        maxys_K.append(np.max(sol_plot))
                        minys_K_x.append(ax_L_plot[np.argmin(minys_K[-1] < sol_plot[sol_plot>0])])
                        maxys_K_x.append(ax_L_plot[np.argmin(maxys_K[-1] > sol_plot)])


                    zorder = len(ax_K) - idx_K
                    ax_target.plot(ax_L_plot, sol_plot, color=colour, linewidth=0.8, label=label, alpha=1, zorder = zorder)
                    #ax_target.scatter(ax_L_plot, sol_plot, color=colour, label=label, alpha=1, zorder = zorder)


                    if zoombox:
                        ax_target_zoom.plot(ax_L_plot,sol_plot, color=colour, linewidth=1, label=label, alpha=1,zorder = zorder)
                        #save the y axis limits:
                        minys_zoom.append(np.min(sol_plot[ax_L_plot>=x1_zoom]))
                        maxys_zoom.append(np.max(sol_plot[ax_L_plot<=x2_zoom]))

                if len(minys_K):
                    minys.append(np.min(minys_K))
                    maxys.append(np.max(maxys_K))
                    #label K
                    if len(sol_plot):
                        K_str = "{:.3f}".format(K_now) + " G$^{0.5}R_{E}$"
                        if K_now == 0: K_str = "0"
                        text_labels.append(ax_target.text(maxys_K_x[0]+0.01,maxys_K[0],"$K=$" + K_str,
                            color='black', size=9,ha="left",va="center"))#,transform=ax_col.transAxes)

                        



            #------------------------------------------------------------------------------+
            # adjust axes, etc.                                                            |
            #------------------------------------------------------------------------------+


            #automatically detect limits:
            xmin_select = ax_target.get_xlim()[0]
            xmax_select = ax_target.get_xlim()[1]
            if len(minys): ymin_select = min(minys)
            if len(maxys): ymax_select = max(maxys)
            ymin_select = max(ymin_select, 10**(log10(ymax_select)-plot_max_order_yrange))
            
            #delete any labels out of range:
            for lab in text_labels:
                #print(lab.get_position(), ymin_select)
                if lab.get_position()[1] < ymin_select:
                    lab.set_text("")

            major_ticks = np.power(10.0,np.linspace(-10,10,21))

            ax_target.set_yticks(major_ticks)
            ax_target.minorticks_off()

            if len(Llims) == 2:
                ax_target.set_xlim(Llims)
            else:
                ax_target.set_xlim([ax_L[0],ax_L[-1]])

            ax_target.set_ylim([ymin_select,ymax_select])


            if usegrid:
                ax_target.grid(which='minor', alpha=0.2)
                ax_target.grid(which='major', alpha=0.5)
                # if zoombox:

            if zoombox:
                #adjust zoomed axes:
                ax_target_zoom.set_yscale('log')
                ax_target_zoom.set_xlim(x1_zoom, x2_zoom)
                ax_target_zoom.set_ylim(max(min(minys_zoom), ymin_select), min(max(maxys_zoom), ymax_select))
                if usegrid:
                    ax_target_zoom.grid(which='minor', alpha=0.2)
                    ax_target_zoom.grid(which='major', alpha=0.5)






            #------------------------------------------------------------------------------+
            # print a set of lines indicating energy for given mu                          |
            #------------------------------------------------------------------------------+

            #enplot = np.array([0.5, 0.75, 1, 1.5, 2, 5, 7.5, 10, 15, 20, 30, 45, 65, 80, 100])
            #enplot_L = np.array([1.2,1.4,1.6,1.8,2.0])
            #enplot_L = np.arange(int(ax_L[0]), ax_L[-1]+0.2, 0.2)
            enplot_L = ax_target.get_xticks()#np.linspace(1.0, 4.0, 13)
            if plot_energy_bars and not fixedenergy:
                sol_en = cdf[interpret_cdf.lab_en][0, :, 0, :]
                #sol_f = cdf[interpret_cdf.lab_f][0, :, 0, :] #t0, K = 0
                sol_f1d = sodata_L_all[rowidx][colidx][0][0]
                idxLvf = 0
                while min(sol_en[:, idxLvf]) < 0:
                    idxLvf+=1
                    if idxLvf==len(ax_L):
                        break
                if idxLvf>=len(ax_L)-2:
                    print("","warning: cannot plot energy bars at K=0")
                else:
                    #sol_f1d = []
                    sol_en1d = [] #used for displaying energy bars when inspectmu == True
                    for idxL in range(idxLvf,len(ax_L)):
                        #sol_f1d.append(np.interp(np.log10(mu),ax_mu[:],sol_f[:, idxL]))
                        sol_en1d.append(np.interp(np.log10(mu), ax_mu[:], sol_en[:, idxL]))
                    #sol_f1d = np.array(sol_f1d)
                    ax_en = np.array(sol_en1d)

                    ax_L_nonfill = ax_L[idxLvf:]#[ax_en>0]

                    enplot_E = np.interp(enplot_L,ax_L_nonfill,ax_en)

                    #ax_top = ax_target.twiny()
                    ax_top = ax_target.secondary_xaxis("top")
                    ax_top_ticks = []
                    ax_top_ticklabels = []
                    #ax_top.set_xticklabels(enplot_E, rotation = 45, va="center")#, position=(0,-0.28))
                    #print(enplot_L);sys.exit(1)
                    for idx in range(len(enplot_L)):
                        if enplot_L[idx] < ax_L[0] or enplot_L[idx] > ax_L[-1]:
                            continue
                        ax_top_ticks.append(enplot_L[idx])
                        ax_top_ticklabels.append(round(enplot_E[idx],2))
                    ax_top.set_xticks(ax_top_ticks)
                    ax_top.set_xticklabels(ax_top_ticklabels)#, rotation=45)



            # text:
            #ax_col.set_title("$\\mu$ = " + str(mu) + "MeV/G", fontdict={'fontsize': yaxisSize - 1}, rotation=0)
            if not fixedenergy:
                ax_target.text(0.035, 0.975,"$\\mu=$" + "{:.1f}MeV/G".format(mu), fontdict={'fontsize': labSize}, rotation=0,
                            horizontalalignment='left', verticalalignment='top',
                            transform=ax_col.transAxes)
            else:
                ax_target.text(0.035, 0.975,"$E=$" + "{:.1f}MeV".format(en), fontdict={'fontsize': labSize}, rotation=0,
                            horizontalalignment='left', verticalalignment='top',
                            transform=ax_col.transAxes)
            # axis ticks:
            ax_target.set_ylabel('')
            ax_target.set_xlabel('')
            ax_target.yaxis.set_ticks_position('both')
            ax_target.xaxis.set_tick_params(direction='out', length=5, labelsize=xtickSize)#, which='bottom'
            ax_target.yaxis.set_tick_params(direction='out', length=2, labelsize=ytickSize)#, which='bottom'
            ax_target.tick_params(labelsize=xtickSize)

            ylims = ax_col.get_ylim()

        if not fixedenergy:
            ax_row[0].set_ylabel('$f$ [km$^{-6}$s$^{3}$]', fontdict={'fontsize': yaxisSize-2}) #m_{0}^{3}
        else:
            ax_row[0].set_ylabel('$j$ [cm$^{-2}$s$^{-1}$str$^{-1}$MeV$^{-1}$]', fontdict={'fontsize': yaxisSize-2})


    for ax in ax_array[n-1]:
        ax.set_xlabel('$L$', fontdict={'fontsize': yaxisSize})

    if not fixedenergy:
        ax_array[0][0].text(-0.1,1.1,'$E$ (MeV)\nat $K=0:$',rotation=0,
            color='black', transform=ax_array[0][0].transAxes,size=9,ha="center",va="center")
        #ax.arrow(-0.122,1.035, 0.09, 0, transform=ax_array[0][0].transAxes, head_width=0.022, head_length=0.01,
        #    fc='dimgrey', ec='dimgrey', clip_on=False)


    #write K units:
    xlim=ax_array[-1][-1].get_xlim()
    #ax_array[-1][-1].text(1.02,0.,"K units\nG$^{0.5}R_{E}$",rotation=0,
    #    transform=ax_array[-1][-1].transAxes,
    #    color='black', size=9,ha="left",va="bottom")#,transform=ax_col.transAxes)


    if dynamic:
        #label each dynamic simulation:
        shadow = False
        textax = ax_array[-1][-1]
        for idx_t, time_plot in enumerate(t_plot):
            colour = interpret_tools.get_colour_from_time(time_plot, t_plot, cmap)
            time_dt = datetime.fromtimestamp(time_plot)
            if len(t_plot) > n*12:
                labelheight = n*(idx_t * 1.0/(len(t_plot)-1))
            else:
                labelheight = n*(idx_t * 1.0/(n*12))
            text = textax.text(1.03, labelheight, time_dt.strftime("%j/%Y"), transform=textax.transAxes, color=colour,
                va='center',fontdict={'fontsize': yaxisSize-4})
            if shadow:
                text.set_path_effects([path_effects.Stroke(linewidth=0.8, foreground='black'),
                                       path_effects.Normal()])


    #data_L_all = np.array(data_L_all)
    #aeq_data_all = np.array(aeq_data_all)
    #data_epoch_all = np.array(data_epoch_all)
    #sodata_L_all = np.array(sodata_L_all)
    #plt.tight_layout()
    return fig, data_L_all, aeq_data_all, data_epoch_all, sodata_L_all

def plot_f_vs_L_mu_panels(cdf, mus, iKs, nt, epochs, Llims):
    fixedenergy = False
    return plot_f_vs_L_panels(cdf, mus, iKs, epochs, fixedenergy, nt, Llims)
def plot_f_vs_L_E_panels(cdf, energies, iKs, nt, epochs, Llims):
    fixedenergy = True
    return plot_f_vs_L_panels(cdf, energies, iKs, epochs, fixedenergy, nt, Llims)
#def plot_f_vs_L_mu_fixedaeq_panels(cdf, mus, iKs, nt): return plot_f_vs_L_panels(cdf, mus, iKs, fixedenergy = False, nt = nt, plot_aeq = True)


def plot_f_vs_L_E_fixedaeq_panels(fname_cdf, energies, Llims, aeq, nt, epochs, omni, ratio, overlapcols=False):
    #reshape the energies array:
    plotshape = np.zeros((len(energies), len(energies[0])))
    cdf = pycdf.CDF(fname_cdf)
    plot_max_order_yrange = 7.5#2.5

    #neq_omni = 12


    dynamic = cdf.attrs[interpret_cdf.lab_dynamic][0]
    if dynamic > 0: #boolean is not preserved by CDF format (?)
        dynamic = True
    else:
        dynamic = False

    ax_mu = cdf[interpret_cdf.lab_axmu]
    ax_K = cdf[interpret_cdf.lab_axK]
    ax_L = cdf[interpret_cdf.lab_axL]
    map_alpha = cdf[interpret_cdf.lab_map]
    ax_t = cdf[interpret_cdf.lab_axt]

    if dynamic:
        if len(epochs):
            t_plot = np.array(epochs)
            nt = len(epochs)
            if t_plot[-1] > ax_t[-1] or t_plot[0] < ax_t[0]:
                print("error: plot time out of range for solution axis")
                return None, None, None, None
        else:
            t_plot = np.linspace(ax_t[-1], ax_t[0], nt)
    else:
        t_plot = np.array([ax_t[0]])
        nt = 1

    n = np.shape(plotshape)[0]
    m = np.shape(plotshape)[1]
    if overlapcols: m = 1
    print("", "{}r x {}c plot".format(n, m))
    sharey = False
    if n == 1: sharey = True

    if len(Llims) != 2:
        Llims = [ax_L[0], ax_L[-1]]


    usegrid = True

    fig, ax_array_all = plt.subplots(n, m, sharex=True, sharey=sharey)#, gridspec_kw={'width_ratios': [0.6]*m})
    cmap = matplotlib.cm.get_cmap('viridis')
    if nt <= 2: cmap = matplotlib.cm.get_cmap('bwr')


    #arrange the axes:
    ax_array = []
    if n == 1:
        ax_array_all = [ax_array_all]
    for row in ax_array_all:
        ax_row = []
        if m == 1:
            ax_row.append(row)
        else:
            for count, col in enumerate(row):
                ax_row.append(col)
        ax_array.append(ax_row)

    sodata_L_all = [[[] for ic in range(np.shape(plotshape)[1])] for ax_row in ax_array]
    data_L_all = [[[] for ic in range(np.shape(plotshape)[1])] for ax_row in ax_array]
    data_epoch_all = [[t_plot for ic in range(np.shape(plotshape)[1])] for ax_row in ax_array]


    ylim_max = 0
    for rowidx, ax_row in enumerate(ax_array):
        print("","starting row #{}".format(rowidx+1))
        for colidx in range(np.shape(plotshape)[1]):
            Eval = energies[rowidx][colidx]
            if overlapcols:
                ax_col = ax_row[0]
                if Eval == None:
                    continue
            else:
                ax_col = ax_row[colidx]
                if Eval == None:
                    #ax_col.axis('off')
                    continue

            

            print("","starting col #{}".format(colidx+1))


            #interpolate to the required time
            for time_plot in t_plot:
                #get surrounding solution times:
                idx_t_1 = np.argmin(time_plot > ax_t)
                t_1 = ax_t[idx_t_1]
                if time_plot == t_1:
                    idx_t_0 = idx_t_1
                    t_0 = t_1
                    frac_t = 0
                else:
                    idx_t_0 = idx_t_1 - 1
                    t_0 = ax_t[idx_t_0]
                    frac_t = 1 - (t_1 - time_plot)/(t_1 - t_0)


                data_L = []
                sodata_L = []

                for idx_L, Lval in enumerate(ax_L):
                    if Lval < Llims[0] or Lval > Llims[1]: continue

                    sol_en = cdf[interpret_cdf.lab_en][0, :, :, idx_L]

                    #take the mu-K slice of the solution cdf at each epoch:
                    sol_f_t_0 = cdf[interpret_cdf.lab_f][idx_t_0, :, :, idx_L]
                    sol_f_t_1 = cdf[interpret_cdf.lab_f][idx_t_1, :, :, idx_L]

                    
                    if omni:
                        aeq_sample = map_alpha[0,idx_L, :]
                        aeq_sample = aeq_sample[aeq_sample>0]
                    elif ratio:
                        aeq_sample = [45, 90]
                    else:
                        aeq_sample = [aeq]
                    j_sample = []


                    for aeq in aeq_sample:
                        #find the surrounding K indicies for the current L, pitch angle:
                        iK0 = np.argmax(aeq >= map_alpha[0,idx_L, :])
                        aeq_iK0 = map_alpha[0,idx_L, iK0]
                        if aeq_iK0 < 0:
                            #print("warning: {:.2f}d is in numerical loss cone at L = {:.2f}".format(aeq, Lval))
                            plot_j = 0
                        else:
                            #interpolate f at iK0 for t0, t1:
                            f_t_0_iK0 = np.interp(Eval, sol_en[:, iK0], sol_f_t_0[:, iK0])
                            f_t_1_iK0 = np.interp(Eval, sol_en[:, iK0], sol_f_t_1[:, iK0])

                            #find iK1, then f at iK1 for t0, t1:
                            if aeq == aeq_iK0: #this is the case when omni=True
                                iK1 = iK0
                                f_t_0_iK1 = f_t_0_iK0
                                f_t_1_iK1 = f_t_1_iK0
                                aeq_iK1 = aeq_iK0

                                plot_f_t_0 = f_t_0_iK1
                                plot_f_t_1 = f_t_1_iK1
                            else:
                                iK1 = iK0 - 1
                                f_t_0_iK1 = np.interp(Eval, sol_en[:, iK1], sol_f_t_0[:, iK1])
                                f_t_1_iK1 = np.interp(Eval, sol_en[:, iK1], sol_f_t_1[:, iK1])
                                aeq_iK1 = map_alpha[0,idx_L, iK1]

                                sinaeq_srnd = [sin(aeq_iK0*np.pi/180), sin(aeq_iK1*np.pi/180)]
                                plot_f_t_0 = np.interp(sin(aeq*np.pi/180), sinaeq_srnd, [f_t_0_iK0, f_t_0_iK1])
                                plot_f_t_1 = np.interp(sin(aeq*np.pi/180), sinaeq_srnd, [f_t_1_iK0, f_t_1_iK1])
                            
                            if Eval < sol_en[0, iK0] or Eval > sol_en[-1, iK0] or Eval < sol_en[0, iK1] or Eval > sol_en[-1, iK1]:
                                print("error: energy {:.2f}MeV out of solution range between iKs {} and {}".format(Eval, iK0, iK1))
                                return None, None, None, None

                            plot_f = plot_f_t_0 + frac_t*(plot_f_t_1 - plot_f_t_0)
                            plot_j = interpret_tools.f2j(Eval, plot_f)
                        j_sample.append(plot_j)
                    
                    if omni:
                        #integrate j across pitch angle using the trapezoid rule (underestimate):
                        #print(aeq_sample)
                        #print(j_sample)
                        aeq_sample_rad = np.radians(aeq_sample)

                        integral = 0
                        for aidx in range(len(aeq_sample)-1):
                            aeq0 = aeq_sample_rad[aidx+1]
                            j_aeq0 = j_sample[aidx+1]

                            aeq1 = aeq_sample_rad[aidx]
                            j_aeq1 = j_sample[aidx]

                            aeqm = (aeq1 + aeq0)/2
                            j_aeqm = (j_aeq1 + j_aeq0)/2

                            daeq = aeq1 - aeq0
                            integral += 2*np.pi*sin(aeqm)*j_aeqm*daeq
                        integral = integral *2 #from 0 to 180
                        plot_j = integral
                    elif ratio:
                        if j_sample[0] == 0:
                            plot_j = np.nan
                        else:
                            plot_j = j_sample[1]/j_sample[0]


                    data_L.append(Lval)
                    sodata_L.append(plot_j)


                data_L_all[rowidx][colidx].append(data_L)
                sodata_L_all[rowidx][colidx].append(sodata_L)



    ylim_max = 0

    for rowidx, ax_row in enumerate(ax_array):
        print("","plotting row #{}".format(rowidx+1))

        for colidx in range(np.shape(plotshape)[1]):
            Eval = energies[rowidx][colidx]
            if overlapcols:
                ax_col = ax_row[0]
            else:
                ax_col = ax_row[colidx]
            #format axis ticks:
            ax_target = ax_col
            ax_target.set_yscale('log')
            ax_target.yaxis.set_ticks_position('both')
            ax_target.xaxis.set_tick_params(direction='out', length=5, labelsize=xtickSize)  # , which='bottom'
            ax_target.yaxis.set_tick_params(direction='out', length=2, labelsize=ytickSize)  # , which='bottom'
            ax_target.tick_params(labelsize=xtickSize)
            ax_target.set_ylabel('')
            ax_target.set_xlim([Llims[0],Llims[1]])
            if Eval == None:
                if not overlapcols:
                    for spine in ['top', 'right', 'left']:
                        ax_target.spines[spine].set_visible(False)
                    ax_target.set_yticks([], [])
                    ax_target.yaxis.set_tick_params(which='both', left=False, right=False)  # , which='bottom'
                continue

            if usegrid: ax_target.grid(linestyle='--', color='grey')


            print("","plotting col #{}".format(colidx+1))

            data_L = data_L_all[rowidx][colidx]
            sodata_L = sodata_L_all[rowidx][colidx]
            
            #interpolate to the required time
            for idx, time_plot in enumerate(t_plot):
                #print(data_L[idx], sodata_L[idx])
                colour = interpret_tools.get_colour_from_time(time_plot, t_plot, cmap)
                ax_target.plot(data_L[idx], sodata_L[idx], color=colour, linewidth=1.,alpha=1)
                ylim_max = max(ylim_max, max(sodata_L[idx]))
            

            
            if ratio:
                #ax_target.set_xlim([2,4])
                ymin_select = 1
                ymax_select = 1000
            else:
                ymin_select = smallpve
                ymax_select = max(smallpve,ylim_max)
                ymin_select = max(ymin_select, 10**(log10(ymax_select)-plot_max_order_yrange))
            ax_target.set_ylim([ymin_select,ymax_select])
            

            axtext_v = 0.89
            if overlapcols:
                axtext_v = axtext_v - 0.1*colidx
            ax_col.text(0.99, axtext_v, "{} MeV".format(Eval), transform=ax_col.transAxes,
                    va='center',ha='right',fontdict={'fontsize': yaxisSize})
            #ax_col.yaxis.set_major_locator(MultipleLocator(20))
            #ax_col.yaxis.set_minor_locator(MultipleLocator(10))
            #if ylim_max <20:
            #    ax_col.yaxis.set_major_locator(MultipleLocator(10))
            ax_col.grid(linestyle='-', color='grey',which='both',alpha=0.6,linewidth=1.)
            #ax_col.set_ylim([0,ylim_max])

        if omni:
            ax_row[0].set_ylabel('$J$ [cm$^{-2}$s$^{-1}$MeV$^{-1}$]', fontdict={'fontsize': yaxisSize-2})
        elif ratio:
            ax_row[0].set_ylabel('$j_{90}$ : $j_{45}$', fontdict={'fontsize': yaxisSize-2})
        else:
            ax_row[0].set_ylabel('$j$ [cm$^{-2}$s$^{-1}$str$^{-1}$MeV$^{-1}$]', fontdict={'fontsize': yaxisSize-2})

    #L on x axes:
    for colidx, ax_col in enumerate(ax_array[-1]):
        ax_col.set_xlabel('$L$', fontdict={'fontsize': xaxisSize})

    if omni:
        ax_array[0][0].text(0, 1.02, 'Omnidirectional flux', transform=ax_array[0][0].transAxes, ha="left", va="bottom", fontdict={'fontsize': yaxisSize})
    elif ratio:
        ax_array[0][0].text(0, 1.02, 'Ratio of flux', transform=ax_array[0][0].transAxes, ha="left", va="bottom", fontdict={'fontsize': yaxisSize})
    else:
        ax_array[0][0].text(0, 1.02, 'Flux at $\\alpha_{\\mathrm{eq}}=$' + "{:.1f}".format(aeq) + '$^{\\circ}$', transform=ax_array[0][0].transAxes, ha="left", va="bottom", fontdict={'fontsize': yaxisSize})
    

    if dynamic:
        #label each dynamic simulation:
        shadow = False
        textax = ax_array[-1][-1]
        for idx_t, time_plot in enumerate(t_plot):
            colour = interpret_tools.get_colour_from_time(time_plot, t_plot, cmap)
            time_dt = datetime.fromtimestamp(time_plot)
            if len(t_plot) > n*12:
                labelheight = n*(idx_t * 1.0/(len(t_plot)-1))
            else:
                labelheight = n*(idx_t * 1.0/(n*12))
            text = textax.text(1.005, labelheight, time_dt.strftime("%j/%Y"), transform=textax.transAxes, color=colour,
                va='center',fontdict={'fontsize': yaxisSize-4})
            if shadow:
                text.set_path_effects([path_effects.Stroke(linewidth=0.8, foreground='black'),
                                       path_effects.Normal()])
    figsy = n*3
    if overlapcols:
        figsx = figsy
    else:
        figsx = m*3
    
    fig.set_size_inches(figsx, figsy)
    fig.tight_layout() 
    

    return fig, data_L_all, data_epoch_all, sodata_L_all


def plot_n_vs_L_panels(fname_cdf, energies, nt, epochs, Llims):
    #reshape the energies array:
    plotenergies = np.array(energies)
    plotshape = np.zeros(np.shape(plotenergies))

    cdf_loaded = interpret_cdf.Preloaded(fname_cdf)
    dynamic = cdf_loaded.dynamic
    if dynamic > 0: #boolean is not preserved by CDF format (?)
        dynamic = True
    else:
        dynamic = False

    Lmax_plot = 3
    L_axis_plot = cdf[interpret_cdf.lab_axL][:][::3]
    if len(Llims):
        Lmin_plot = Llims[0]
        Lmax_plot = Llims[1]
        L_axis_plot = L_axis_plot[L_axis_plot<=Lmax_plot]
        if L_axis_plot[-1] != Lmax_plot:
            L_axis_plot = list(L_axis_plot)
            L_axis_plot.append(Lmax_plot)
            L_axis_plot = np.array(L_axis_plot)
        L_axis_plot = L_axis_plot[L_axis_plot>=Lmin_plot]
        if L_axis_plot[0] != Lmin_plot:
            L_axis_plot = list(L_axis_plot)
            L_axis_plot.insert(0,Lmin_plot)
            L_axis_plot = np.array(L_axis_plot)
    #print(L_axis_plot);sys.exit()


    #which times to plot (for a dynamic run):
    ax_t = cdf_loaded.ax_t
    if dynamic:
        if len(epochs):
            t_plot = np.array(epochs)
            nt = len(epochs)
            if t_plot[-1] > ax_t[-1] or t_plot[0] < ax_t[0]:
                print("error: plot time out of range for solution axis")
                return None, None, None, None
        else:
            t_plot = np.linspace(ax_t[-1], ax_t[0], nt)
    else:
        t_plot = np.array([ax_t[0]])
        nt = 1

    n = np.shape(plotshape)[0]
    m = np.shape(plotshape)[1]
    print("", "{}r x {}c plot".format(n, m))
    sharey = False
    if n == 1: sharey = True
    fig, ax_array_all = plt.subplots(n, m, sharex=True, sharey=sharey)#, gridspec_kw={'width_ratios': [0.6]*m})
    cmap = matplotlib.cm.get_cmap('viridis')
    usegrid = True

    #arrange the axes:
    ax_array = []
    if n == 1:
        ax_array_all = [ax_array_all]
    for row in ax_array_all:
        ax_row = []
        if m == 1:
            ax_row.append(row)
        else:
            for count, col in enumerate(row):
                ax_row.append(col)
        ax_array.append(ax_row)

    ylim_max = 0
    #print(t_plot)
    for rowidx, ax_row in enumerate(ax_array):
        print("","starting row #{}".format(rowidx+1))
        for colidx, ax_col in enumerate(ax_row):
            en = float(plotenergies[rowidx][colidx])
            ax_target = ax_col
            print("","starting col #{}".format(colidx+1))

            #interpolate to the required time
            for time_plot in t_plot:
                #print(time_plot)
                n_ = []
                ratio = []
                for L_ in L_axis_plot:
                    alpha_lc = interpret_tools.get_lc(L_)

                    alpha_eqpad_deg, f_ = interpret_tools.getpad_(cdf_loaded, L_, en, dynamic, time_plot)

                    if not len(f_):
                        n_.append(np.nan)
                        continue

                    j_ = interpret_tools.f2j(en, f_)

                    if j_[-1] <= 0: 
                        n_.append(np.nan)
                        continue

                    j_tofit_sinnc = [alpha_eqpad_deg, j_]


                    #get the x axis in radians:
                    #alpha_eqpad_rad = np.radians(alpha_eqpad_deg)


                    p0 = [j_[-1], 4] #A, n
                    popt, pcov = opti.curve_fit(interpret_tools.f_sinn_simple, j_tofit_sinnc[0], j_tofit_sinnc[1], p0= p0,
                                                bounds=((0, 0), (2*j_[-1], 200)), maxfev=2000)


                    perr = np.sqrt(np.diag(pcov))
                    n_.append(popt[-1])
                    #a = f_sinn_simple(np.array([70.0,90.0]),*popt)
                    #ratio.append(a[1]/a[0])

                colour = interpret_tools.get_colour_from_time(time_plot, ax_t, cmap)
                ax_target.plot(L_axis_plot, n_, color=colour, linewidth=1.4,alpha=1)
                #ax_target.scatter(L_axis_plot, n_, color=colour,marker='o', alpha=1)
                #ax_target.plot(L_axis_plot, ratio, color=colour, linewidth=1.4,alpha=1,linestyle=":")
                #ax_target.scatter(L_axis_plot, ratio, color=colour,marker='o', alpha=1)

            ax_target.set_xlim([L_axis_plot[0],L_axis_plot[-1]])
            if usegrid: ax_target.grid(linestyle='--', color='grey')



            
            #format axis ticks:
            ax_target.yaxis.set_ticks_position('both')
            ax_target.xaxis.set_tick_params(direction='out', length=5, labelsize=xtickSize)  # , which='bottom'
            ax_target.yaxis.set_tick_params(direction='out', length=2, labelsize=ytickSize)  # , which='bottom'
            ax_target.tick_params(labelsize=xtickSize)
            ax_target.set_ylabel('')

            #keep a record of the max value:
            ylim = ax_target.get_ylim()
            if (ylim[1] > ylim_max):
                ylim_max = ylim[1]


            #if addlabel:
            #    ax_row[colidx].text(0.1,0.1,"L = {}, E = {}".format(Lplot, en), transform=ax_row[colidx].transAxes) #DELETE THIS 

            #ax_row[colidx].yaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter(useMathText=True, useOffset=False))

        #ax_row[0].set_ylabel('f [km$^{-6}$s$^{3}$]', fontdict={'fontsize': yaxisSize})


    #add additional plot information:
    #aeq on x axes:
    for colidx, ax_col in enumerate(ax_array[-1]):
        ax_col.set_xlabel('$L$', fontdict={'fontsize': xaxisSize})

    plt.tight_layout(pad=5, w_pad=0, h_pad=0)
    

    for rowidx, ax_row in enumerate(ax_array):
        for colidx, ax_col in enumerate(ax_row):
            en = plotenergies[rowidx][colidx]

            ax_col.text(0.99, 0.89, "{} MeV".format(en), transform=ax_col.transAxes,
                    va='center',ha='right',fontdict={'fontsize': yaxisSize})
            ax_col.yaxis.set_major_locator(MultipleLocator(20))
            ax_col.yaxis.set_minor_locator(MultipleLocator(10))
            if ylim_max <20:
                ax_col.yaxis.set_major_locator(MultipleLocator(10))
            ax_col.grid(linestyle='--', color='grey',which='both',alpha=0.6,linewidth=1.3)
            #ax_col.set_ylim([0,ylim_max])
            ax_col.set_ylim([0,80])

            #set_size(ax_col, 2*(figsx/m),2*(figsy/n)*(ylim_max/60))


    if dynamic:
        #label each dynamic simulation:
        shadow = False
        textax = ax_array[-1][-1]
        for idx_t, time_plot in enumerate(t_plot):
            colour = interpret_tools.get_colour_from_time(time_plot, ax_t, cmap)
            time_dt = datetime.fromtimestamp(time_plot)
            if len(t_plot) > n*12:
                labelheight = n*(idx_t * 1.0/(len(t_plot)-1))
            else:
                labelheight = n*(idx_t * 1.0/(n*12))
            text = textax.text(1.005, labelheight, time_dt.strftime("%j/%Y"), transform=textax.transAxes, color=colour,
                va='center',fontdict={'fontsize': yaxisSize-4})
            if shadow:
                text.set_path_effects([path_effects.Stroke(linewidth=0.8, foreground='black'),
                                       path_effects.Normal()])
    

    
    figsx = m*4
    figsy = n*2
    fig.set_size_inches(figsx, figsy)
    fig.tight_layout() 
    
    
    #plt.savefig(fname_plot, dpi = 400)

    return fig







def plot_j_timeseries_fixedaeq_panels(fname_cdf, Lshells, energies, aeq, epochrange):
    cdf = pycdf.CDF(fname_cdf)
    nplots = len(Lshells)

    #check solution is dynamic:
    dynamic = cdf.attrs[interpret_cdf.lab_dynamic][0]
    if dynamic > 0:
        pass
    else:
        print("", "warning: solution is not dynamic, skipping...")
        return 0

    #process epoch range:
    date_start_ts = epochrange[0]
    date_finish_ts = epochrange[1]

    datelims = [datetime.fromtimestamp(date_start_ts), datetime.fromtimestamp(date_finish_ts)]


    xtick_locator = AutoDateLocator()
    xtick_formatter = AutoDateFormatter(xtick_locator)
    plot_max_order_yrange = 5
    plot_min_order_yrange = 1

    #axes:
    ax_mu = cdf[interpret_cdf.lab_axmu]
    ax_K = cdf[interpret_cdf.lab_axK]
    ax_L = np.array(cdf[interpret_cdf.lab_axL])
    map_alpha = cdf[interpret_cdf.lab_map]
    ax_t = cdf[interpret_cdf.lab_axt]
    print("", "{}r x 1c plot".format(nplots))



    fig, axs = plt.subplots(len(Lshells), sharex=True)
    figsx = 5
    figsy = figsx/2*len(Lshells)

    data_t = []
    data_j_energies_time_L = []
    for idx_sp, L_plot in enumerate(Lshells):
        ax = axs[idx_sp]

        if L_plot < ax_L[0] or L_plot > ax_L[-1]:
            print("","warning, L={:.2f} is outside the range of the model, skipping...".format(L_plot))
            continue

        #interpolate to the required L:
        idx_L_1 = np.argmin(L_plot > ax_L)
        L_1 = ax_L[idx_L_1]
        if L_plot == L_1:
            idx_L_0 = idx_L_1
            L_0 = L_1
            frac_L = 0
        else:
            idx_L_0 = idx_L_1 - 1
            L_0 = ax_L[idx_L_0]
            frac_L = 1 - (L_1 - L_plot)/(L_1 - L_0)


        #for every available time:
        data_j_energies_time = []
        data_t = []
        for idx_t, time_plot in enumerate(ax_t):
            if time_plot < date_start_ts or time_plot > date_finish_ts: continue
            data_t.append(time_plot)

            #for each energy:
            data_j_energies = []
            for idx_E, Eval in enumerate(energies):

                #
                # idx_L_0:
                #
                #take the mu-K slice of the solution cdf at current epoch:
                sol_f = cdf[interpret_cdf.lab_f][idx_t, :, :, idx_L_0]
                sol_en = cdf[interpret_cdf.lab_en][0, :, :, idx_L_0]
                #
                #find the surrounding K indicies for the current L, pitch angle:
                iK0 = np.argmax(aeq >= map_alpha[0,idx_L_0, :])
                aeq_iK0 = map_alpha[0, idx_L_0, iK0]
                if aeq == aeq_iK0: #for 90 deg
                    iK1 = iK0
                else:
                    iK1 = iK0 - 1
                aeq_iK1 = map_alpha[0,idx_L_0, iK1]
                #
                if aeq_iK0 < 0:
                    plot_j_L_0 = smallpve
                else:
                    if Eval < sol_en[0, iK0] or Eval > sol_en[-1, iK0] or Eval < sol_en[0, iK1] or Eval > sol_en[-1, iK1]:
                        print("error: energy {:.2f}MeV out of solution range between iKs {} and {}".format(Eval, iK0, iK1))
                        return None, None, None
                    f_iK0 = np.interp(Eval, sol_en[:, iK0], sol_f[:, iK0])
                    f_iK1 = np.interp(Eval, sol_en[:, iK1], sol_f[:, iK1])
                    plot_f_L_0 = np.interp(aeq, [aeq_iK0, aeq_iK1], [f_iK0, f_iK1])
                    plot_j_L_0 = max(smallpve, interpret_tools.f2j(Eval, plot_f_L_0))

                #
                # idx_L_1:
                #
                #take the mu-K slice of the solution cdf at current epoch:
                sol_f = cdf[interpret_cdf.lab_f][idx_t, :, :, idx_L_1]
                sol_en = cdf[interpret_cdf.lab_en][0, :, :, idx_L_1]
                #
                #find the surrounding K indicies for the current L, pitch angle:
                iK0 = np.argmax(aeq >= map_alpha[0,idx_L_1, :])
                aeq_iK0 = map_alpha[0, idx_L_1, iK0]
                if aeq == aeq_iK0: #for 90 deg
                    iK1 = iK0
                else:
                    iK1 = iK0 - 1
                aeq_iK1 = map_alpha[0,idx_L_1, iK1]
                #
                if aeq_iK0 < 0:
                    plot_j_L_1 = smallpve
                else:
                    if Eval < sol_en[0, iK0] or Eval > sol_en[-1, iK0] or Eval < sol_en[0, iK1] or Eval > sol_en[-1, iK1]:
                        print("error: energy {:.2f}MeV out of solution range between iKs {} and {}".format(Eval, iK0, iK1))
                        return None, None, None
                    f_iK0 = np.interp(Eval, sol_en[:, iK0], sol_f[:, iK0])
                    f_iK1 = np.interp(Eval, sol_en[:, iK1], sol_f[:, iK1])
                    plot_f_L_1 = np.interp(aeq, [aeq_iK0, aeq_iK1], [f_iK0, f_iK1])
                    plot_j_L_1 = max(smallpve, interpret_tools.f2j(Eval, plot_f_L_1))


                #average j at L_0 and L_1:
                plot_j = 10**((1 - frac_L) * log10(plot_j_L_0) + frac_L * log10(plot_j_L_1))

                data_j_energies.append(plot_j)
            data_j_energies_time.append(data_j_energies)

        #make lists into arrays:
        data_j_energies_time = np.array(data_j_energies_time)
        data_t = np.array(data_t)
        data_j_energies_time_L.append(data_j_energies_time) #return variable

        #make timestamps into datetimes for x axis:
        x_axis_dates = [datetime.fromtimestamp(x) for x in data_t]

        #derive limits:
        flux_min = min([np.min(x) for x in data_j_energies_time.T])
        flux_max = max([np.max(x) for x in data_j_energies_time.T])
        if flux_max <= smallpve:
            flux_min_plot = flux_max
            ax.axis('off')
            continue
        else:
            flux_min_plot = max(flux_min, 10**(log10(flux_max)-plot_max_order_yrange)) #make sure min flux is not less than flux_max - some limit
            flux_min_plot = min(flux_min, 10**(log10(flux_max)-plot_min_order_yrange))

        #plot:
        for idx_E, Eval in enumerate(energies):
            j_timeseries = data_j_energies_time[:, idx_E]
            ax.scatter(x_axis_dates, j_timeseries, marker = '.', color=interpret_tools.getcol(idx_E, len(energies)))
            ax.plot(x_axis_dates, j_timeseries, color=interpret_tools.getcol(idx_E, len(energies)))

                
        #set limits:
        flux_min_plot = 0.8*flux_min_plot
        flux_max_plot = 1.2*flux_max
        ax.set_ylim([flux_min_plot,flux_max_plot])
        ax.set_ylabel('$j$ *')
        ax.grid(linestyle='--', color='grey')

        #y axis:        
        ax.set_yscale('log')

        #x axis:
        ax.xaxis_date()
        ax.xaxis.set_major_locator(xtick_locator)
        ax.xaxis.set_major_formatter(xtick_formatter)
        #ax.xaxis.set_major_formatter(mdates.DateFormatter("%Y/%m/%d"))
        fig.autofmt_xdate()
        #ax.set_xlim(datelims)
        

        #energy labels:
        label_xspan = 0.13
        for idx_E in range(len(energies)): 
            color = interpret_tools.getcol(idx_E, len(energies))
            elabel = "{:.2f}".format(energies[idx_E])
            ax_elabel = ax.text(0.97-(len(energies)-idx_E)*label_xspan,0.01, elabel, va="bottom", ha="right", transform=ax.transAxes, fontsize = 13, color = color)
            ax_elabel.set_bbox(dict(facecolor='white', alpha=0.5, edgecolor='none'))
        ax_elabel = ax.text(0.97,0.01, "MeV", va="bottom", ha="right", transform=ax.transAxes, fontsize = 13, color = "black")
        ax_elabel.set_bbox(dict(facecolor='white', alpha=0.5, edgecolor='none'))

        #L labels:
        ax_Llabel = ax.text(0.5,0.98, "L={:.2f}".format(L_plot), va="top", ha="center", transform=ax.transAxes, fontsize = 13, color = "black")
        ax_Llabel.set_bbox(dict(facecolor='white', alpha=0.5, edgecolor='none'))
    axs[0].text(1, 1.02, '* [cm$^{-2}$s$^{-1}$str$^{-1}$MeV$^{-1}$], $\\alpha_{\\mathrm{eq}}=$' + "{:.1f}".format(aeq) + '$^{\\circ}$', transform=axs[0].transAxes, ha="right", va="bottom", fontdict={'fontsize': 11})
      
    fig.set_size_inches(figsx, figsy)
    fig.tight_layout() 

    return fig, data_t, data_j_energies_time_L



def plot_map(fname_cdf, energy, alt_km, dt_max_days = 30):
    dt_max_seconds = dt_max_days * 24 * 60 * 60

    def generate_new_map_fig(x_, y_, L_):
        #return a figure and axes with a world map and L contours drawn on it

        #make figure:
        figsize = (12,10)
        fig_map, ax_map = plt.subplots(figsize=figsize)
        ax_map.set_facecolor("white")

        #make cbar axis:
        divider = make_axes_locatable(ax_map)
        cax = divider.append_axes("right", size="5%", pad=0.04)

        #generate world map:
        path_map = gpd.datasets.get_path("naturalearth_lowres")
        #print("","loading map from",path_map)
        world = gpd.read_file(path_map)
        world.plot(ax=ax_map, color='none', edgecolor='black', linewidths=0.7, zorder = 3)
        world.plot(ax=ax_map, color='none', edgecolor='white', linewidths=1.5, zorder = 2)

        #plot L contours:
        CS = ax_map.contour(x_, y_, L_,zorder=5, colors='black',linewidths=1.5,
            levels=np.array([1.1,1.2,1.3,1.4,1.6,2.0,2.4,3.0, 3.9]))
        clbls = ax_map.clabel(CS, inline=True, fontsize=10, zorder=5)
        [txt.set_bbox(dict(facecolor='white', edgecolor='none', pad=0)) for txt in clbls]
        CS = ax_map.contour(x_, y_, L_,zorder=4, colors='white',linewidths=2.5,
            levels=np.array([1.1,1.2,1.3,1.4,1.6,2.0,2.4,3.0, 3.9]))

        #ticks and labels:
        ax_map.set(xlabel='Longitude [$^{\circ}$]', ylabel='Latitude [$^{\circ}$]')
        ax_map.set_ylim((-85,85))
        tickcolor = 'black'
        labelcolor = 'black'
        spinecolor = 'black'
        ax_map.tick_params(axis="x", which='both', top=True, colors=tickcolor, direction='in', labelcolor=labelcolor)
        ax_map.tick_params(axis="y", which='both', right=True, colors=tickcolor, direction='in', labelcolor=labelcolor)
        ax_map.tick_params(axis="x", which='major',length=6, width=1.2)
        ax_map.tick_params(axis="x", which='minor',length=4, width=1)
        ax_map.tick_params(axis="y", which='major',length=7, width=1.2)
        ax_map.tick_params(axis="y", which='minor',length=5, width=1)
        ax_map.xaxis.set_major_locator(MultipleLocator(50))
        ax_map.xaxis.set_major_formatter('{x:.0f}')
        ax_map.xaxis.set_minor_locator(MultipleLocator(10))
        ax_map.yaxis.set_major_locator(MultipleLocator(20))
        ax_map.yaxis.set_major_formatter('{x:.0f}')
        ax_map.yaxis.set_minor_locator(MultipleLocator(10))
        ax_map.spines['bottom'].set_color(spinecolor)
        ax_map.spines['top'].set_color(spinecolor) 
        ax_map.spines['right'].set_color(spinecolor)
        ax_map.spines['left'].set_color(spinecolor)

        return fig_map, ax_map, cax


    cdf_loaded = interpret_cdf.Preloaded(fname_cdf)
    dynamic = cdf_loaded.dynamic

    if dynamic > 0:
        #plot:
        # map at epochs[1]
        # map at epochs[0]
        # the ratio between both
        dynamic = True
    else:
        #plot:
        # map at most recent available epoch
        dynamic = False

    #which times to plot:
    ax_t = cdf_loaded.ax_t
    if dynamic:
        t_plot = np.array([max(ax_t[0], ax_t[-1] - dt_max_seconds), ax_t[-1]])
    else:
        t_plot = np.array([ax_t[0]])

    ax_L = cdf_loaded.ax_L



    #make lat, lon, alt grid:
    lon_deg = np.linspace(-180,180,num=61)
    lat_deg = np.linspace(-75, 75, num=36)
    alt_m = alt_km * 1e3

    x_, y_ = np.meshgrid(lon_deg, lat_deg)


    #get magnetic coordinates at each grid point: ------------------------------+
    Lmin_calc = 1 #L limits outside which computation is skipped
    Lmax_calc = 5

    sc_time0 = datetime.fromtimestamp(t_plot[0])
    dtdates0 = spt.Ticktock([sc_time0]*x_.size, 'UTC')

    sc_R_RE_av = np.repeat([1 + alt_m/interpret_tools.RE], x_.size)

    x_1d = np.reshape(x_,(1,np.size(x_)))
    y_1d = np.reshape(y_,(1,np.size(y_)))
    sc_xsc = np.vstack((sc_R_RE_av, y_1d, x_1d)).T #R, lat, lon

    #create spacepy coords:
    spxsc = spc.Coords(sc_xsc, 'GEO', 'sph') #the units must be re
    spxsc.ticks = dtdates0

    IRBEM_global_option_extMag = '0'
    IRBEM_global_option_intMag = 'IGRF'
    L_ = []
    Be_ = []
    Bm_ = []
    jomnis = [[], []] #list of values at t0, t1
    energy_ = np.array([energy]) #get_jomni only accepts np array of energy
    #iterate over coordinate grid:
    print("progress extracting model omnidirectional flux at geographical coordinates:")
    for idx in range(len(dtdates0)):

        #get magnetic coordinates from IRBEM:
        ibdata = ib.get_Lm(dtdates0[idx], spxsc[idx], [90], extMag=IRBEM_global_option_extMag, intMag=IRBEM_global_option_intMag)
        Lm = ibdata['Lm'][0][0] #third argument allows to take account of L shell splitting
        
        Be = ibdata['Bmin'][0]
        Bm = ibdata['Bmirr'][0]

        #ibdata = ib.find_magequator(dtdates0[idx], spxsc[idx], extMag=IRBEM_global_option_extMag)
        #Be = ibdata['Bmin'][0] 
        #ibdata = ib.get_Bfield(dtdates0[idx], spxsc[idx], extMag=IRBEM_global_option_extMag)
        #Bm = ibdata['Blocal'][0] 

        L_.append(Lm)
        Be_.append(Be)
        Bm_.append(Bm)
        

        #check L is in range
        if Lm == Lm and Lm >= Lmin_calc and Lm <= Lmax_calc and Lm <= ax_L[-1]:
            #extract observable omnidirectional flux from the model results:
            for idx_t, time in enumerate(t_plot):
                jomnis[idx_t].append(interpret_tools.get_jomni(cdf_loaded, time, Lm, Be/Bm, energy_)[0])
        else:
            for idx_t, time in enumerate(t_plot):
                jomnis[idx_t].append(0)

        if idx%200 == 0:
            print(" {:.2f}%".format(100*(idx)/len(dtdates0)))


    L_ = np.array(L_)
    Be_ = np.array(Be_)
    Bm_ = np.array(Bm_)

    L_ = np.reshape(L_, (np.shape(x_)[0],np.shape(x_)[1]))
    Be_ = np.reshape(Be_, (np.shape(x_)[0],np.shape(x_)[1]))
    Bm_ = np.reshape(Bm_, (np.shape(x_)[0],np.shape(x_)[1]))

    L_[L_<Lmin_calc] = np.nan
    L_[L_>Lmax_calc] = np.nan


    #find the min, max flux limits and convert flux to a numpy array
    jomni_max = 0
    jomni_min = np.inf
    for idx_t, time in enumerate(t_plot):
        jomnis_t = np.array(jomnis[idx_t])
        jomnis_t = np.reshape(jomnis_t, (np.shape(x_)[0],np.shape(x_)[1]))
        jomnis[idx_t] = jomnis_t

        #find min and max flux:
        jomnis_t_valid = jomnis_t[jomnis_t == jomnis_t]
        if not len(jomnis_t_valid):
            print("error: no valid model results at time index {idx_t} for this coordinate range".format(idx_t))
            return None, None, 0

        jomni_max = max(np.max(jomnis_t_valid), jomni_max)
        jomni_min = min(np.min(jomnis_t_valid), jomni_min)

    #check the min flux is in range:
    max_dlog10flux = 4 #######
    jomni_min = max([jomni_min, smallpve, 10**(log10(jomni_max) - max_dlog10flux)])

    # --------------------------------------------------------------------------+


    figures = []

    #plot flux:
    for idx_t in range(len(t_plot)):
        fig_map, ax_map, cax = generate_new_map_fig(x_, y_, L_)

        #plot flux colour:
        cmap_flux = 'rainbow'
        norm_flux = matplotlib.colors.LogNorm(vmin=jomni_min, vmax=jomni_max)
        pc = ax_map.pcolormesh(x_, y_, jomnis[idx_t], norm=norm_flux, cmap=cmap_flux)


        #flux colourbar: -----------------------------------------------------------+
        m_flux = cm.ScalarMappable(norm=norm_flux, cmap=cmap_flux)

        cbar_flux = plt.colorbar(m_flux, cax = cax)

        cax.set_ylabel("$j $ [MeV$^{-1}$cm$^{-2}$s$^{-1}$]")
        # --------------------------------------------------------------------------+

        #title:
        timestr = datetime.fromtimestamp(t_plot[idx_t]).strftime("%Y/%m/%d %H:%M:%S")
        ax_map.set_title('Modelled {:.1f}MeV Omnidirectional Flux at {:.0f}km, {}'.format(energy, alt_km, timestr),
           color='black', size=17, loc='left', ha="left", va="center")

        plt.tight_layout()

        figures.append(fig_map)

    #if we have flux at 2 different epochs, make an extra plot of the ratio:
    if len(t_plot) == 2:
        fig_map, ax_map, cax = generate_new_map_fig(x_, y_, L_)

        #plot flux ratio:
        cmap_flux = 'bwr'
        norm_flux = matplotlib.colors.LogNorm(vmin=1./20, vmax=20.)
        jomni_ratio = jomnis[1]/jomnis[0]
        pc = ax_map.pcolormesh(x_, y_, jomni_ratio, norm=norm_flux, cmap=cmap_flux)

        #flux colourbar: -----------------------------------------------------------+
        m_flux = cm.ScalarMappable(norm=norm_flux, cmap=cmap_flux)

        cbar_flux = plt.colorbar(m_flux, cax = cax)

        cax.set_ylabel("flux ratio (after:before)")
        # --------------------------------------------------------------------------+

        #title:
        timestr0 = datetime.fromtimestamp(t_plot[0]).strftime("%Y/%m/%d %H:%M:%S")
        timestr1 = datetime.fromtimestamp(t_plot[1]).strftime("%Y/%m/%d %H:%M:%S")
        ax_map.set_title('Change in {:.1f}MeV Omnidirectional Flux, {} to {}'.format(energy, timestr0, timestr1),
           color='black', size=14.5, loc='left', ha="left", va="center")

        plt.tight_layout()

        figures.append(fig_map)
    else:
        jomni_ratio = []

    return figures, [x_, y_], [jomnis[0], jomnis[1], jomni_ratio]
