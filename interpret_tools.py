import interpret_cdf
import colorsys
import sys
import numpy as np
from math import sqrt, asin, pi, sin
from scipy.interpolate import interp1d

mass0_proton = 1.6726219e-27
MeV2J = 1.60218e-13
c_ = 299792458
G2T = 1.e-4
RE = 6.3712e6

# def ticks_to_datetimeUTCaware(ticks):
#     datetimes = ticks.data
#     datetimes_aware = []
#     for idx in range(len(datetimes)):
#         datetimes_aware.append(datetimes[idx].replace(tzinfo = timezone.utc))
#     return np.array(datetimes_aware)
    
def daterange_to_str(datelims):
    if len(datelims):
        datestr = "_"+datelims[0].strftime("%Y%m%d-")+datelims[1].strftime("%Y%m%d")
    else:
        datestr = ""

def userselectkey(filedict, allowmulti=False):
    #user selects keys from a dictionary of items supplied as argument
    #input is sanitised and returned as a list of keys
    
    for key in filedict.keys():
        print(key, '...', filedict[key])

    filedict_selectkeys = [] #save the selection as keys
    
    #ask the user to selet results files to plot:
    decided = False
    more = False
    while not decided:
        choice = input("> ")
        #sanitise the input and re-ask if necessary:
        try:
            if choice[-1] == ",":
                more = True * allowmulti
                choice = choice[:-1]
            else:
                more = False
            choice = int(choice)
            if choice in filedict.keys():
                if choice in filedict_selectkeys:
                    print ("  already selected")
                else:
                    filedict_selectkeys.append(choice)
                if not more:
                    decided = True
            else:
                print("  out of range")
        except:
            print("  invalid")
            pass
            
    return filedict_selectkeys

def userinputfloat(allowblank=False, allowmulti=False):
    floatselectkeys = [] #save the selection as keys
    
    decided = False
    more = False
    while not decided:
        choice = input("> ")
                
        if choice == "":
            if allowblank:
                decided = True
                floatselectkeys.append(False)
            else:
                print("  a value must be specified")
        else:
            #sanitise the input and re-ask if necessary:
            try:
                if choice[-1] == ",":
                    more = True * allowmulti
                    choice = choice[:-1]
                else:
                    more = False

                choice = float(choice)
                
                if choice in floatselectkeys:
                    print ("  already selected")
                else:
                    floatselectkeys.append(choice)
                    
                if not more:
                    decided = True
            except:
                print("  invalid")
                pass
    return floatselectkeys


def get_N_HexCol(N):
    HSV_tuples = [(x*1.0/N, 0.7, 0.7) for x in range(N)]
    RGB_tuples = map(lambda x: colorsys.hsv_to_rgb(*x), HSV_tuples)
    return RGB_tuples

def getlastline(fname):
    with open(fname, 'r') as f:
        lines = f.read().splitlines()
    last_line = lines[-1]
    return last_line.strip("\n")

# def bytes2float(bytestring):
#     return float(((bytestring).decode("utf-8")).strip(' ').strip('\n'))
# def bytes2str(bytestring):
#     return str(((bytestring).decode("utf-8")).strip(' ').strip('\n'))

# def eqL2B(Larray):
#     B0 = 3.12e-5
#     T2G = 1.0e4
#     return T2G*B0*(np.power(np.reciprocal(Larray),3))

def getcol(i, n):
    ncols = list(get_N_HexCol(n))
    return ncols[i]

def get_axis_file(fname):
    with open(fname) as fi:
        ax_ = fi.readlines()
        for idxL in range(0, len(ax_)):
            ax_[idxL] = float(ax_[idxL].strip('/n'))
    ax_ = np.array(ax_)
    return ax_

def get_sol_file(fname):
    lines = []
    with open(fname) as fo:
        lines = fo.readlines()

    sol_f = []
    for idx in range(0, len(lines)):
        sol_f.append([float(x) for x in lines[idx].strip('\n').split()])

    sol_f = np.array(sol_f)
    return sol_f

def get_gamma(ke_MeV):
    #calculate gamma of a proton given its kinetic energy in MeV

    ke_J = ke_MeV * MeV2J

    erest_J = mass0_proton*c_*c_

    gamma = 1 + ke_J/(erest_J)
    return gamma

def get_pr(ke_MeV):
    # !calculate relativistic momentum in SI units (kg m s-1) of a proton given its kinetic energy in MeV


    erest_J = mass0_proton*c_*c_
    
    gamma = get_gamma(ke_MeV)

    mr = mass0_proton * gamma
    etot_J = mr*c_*c_
    
    p_ = np.sqrt((etot_J**2) - (erest_J**2))/c_ #relativistic momentum kg m/s
    return p_

def f2j(energy, psd):
    # !
    # ! input units of energy: MeV
    # ! input units of phase space density: m-6 s3 kg-3
    # !
    # !

    #change units to m-6 s3: 
    temp = psd / 1e18
    # units: m-6 s3
    
    
    temp = temp / (mass0_proton**3.)
    # units are m-6 s3 kg-3

    #get momentum squared:
    p2 = (get_pr(energy) **2.)
    # units: kg2 m2 s-2, or: ** kg J **

    flux = temp * p2
    # units: m-2 s-1 str-1 J-1

    flux = flux / 1e4
    # units: cm-2 s-1 str-1 J-1

    flux = flux * MeV2J
    # units: cm-2 s-1 str-1 MeV-1
    return flux

def j2f(energy, flux):
    # !
    # ! input units of energy: MeV
    # ! input units of flux: cm-2 s-1 str-1 MeV-1
    # !
    # !


    #get momentum squared:
    p2 = (get_pr(energy) **2.)
    # units: kg2 m2 s-2, or: ** kg J **


    flux = flux / MeV2J
    # units: cm-2 s-1 str-1 J-1
    
    flux = flux * 1e4
    # units: m-2 s-1 str-1 J-1
    

    temp = flux / p2
    
    temp = temp * (mass0_proton**3.)

    psd = temp * 1e18
    # units: m-6 s3 kg-3
    
    return psd


def get_colour_from_time(time, ax_t, cmap):
    if (len(ax_t) > 1):
        frac = (time-ax_t[0])/(ax_t[-1]-ax_t[0])
    else:
        frac=0
    return cmap(frac)

def f_sinn(x, A, b, c, n):
    # baspropy does not like a negative number to a decimal power
    # p0 should be something like [10, 0, 25, 4] in practise
    d2r = np.pi / 180.
    sinn = np.abs(np.sin((x+b)*d2r))
    return A * np.power(sinn,n) + c

def f_sinn_simple(x, A, n):
    # baspropy does not like a negative number to a decimal power
    # p0 should be something like [10, 0, 25, 4] in practise
    d2r = np.pi / 180.
    sinn = np.abs(np.sin((x)*d2r))
    return A * np.power(sinn,n)

def get_lc(Lb): #centred dipole loss cone approximation for 2015
    RE = 6.3712e6
    atm_height_std = 100000
    B0 = 2.986731323946967e-05

    ra = (RE + atm_height_std)/RE #~Earth's surface + atm_height_dipolelc m
    
    if ra >= Lb:
        return 90
    else:
        Ba = (B0/(ra**3)) * (4 - 3*ra/Lb)**(0.5)
        dipole_lc = asin(sqrt((B0 / Lb**3)/Ba)) * 180 / pi
        return dipole_lc

def getpad_(cdf_loaded, L_extract, en_all, dynamic, time_plot):
    #
    #
    #   WARNING: when interpolating between L output by the model, we can't determine where the loss cone is!
    #    so the PAD will be returned with a point at (0,0) but then shoot up to the first point outside the lc
    #
    #
    # returns pitch angle distribution alpha and corresponding f, straight from the solution grid
    # check energy was supplied as a list:
    if isinstance(en_all, list) or isinstance(en_all, (np.ndarray)):
        multiple = True
    else:
        multiple = False
        en_all = [en_all]
    # check energies are in ascending order:
    for idx in range(len(en_all) - 1):
        if en_all[idx + 1] < en_all[idx]:
            print("Error: energies must be supplied in ascending order")
            sys.exit(1)

    ax_t = cdf_loaded.ax_t
    ax_K = cdf_loaded.ax_K
    ax_L = cdf_loaded.ax_L
    map_alpha = cdf_loaded.map_alpha

    sol_f1d_allK = []
    sol_alpha_allK = []

    for idx_K in range(len(ax_K)):
        if not sum(map_alpha[0, :, idx_K] > 0):
            continue
        sol_en = cdf_loaded.sol_en[0, :, idx_K, :]
        sol_f = cdf_loaded.sol_f[:, :, idx_K, :]

        # find the minimum L that is outside the loss cone at the current K:
        idxL_outsidelc = np.argwhere(map_alpha[0, :, idx_K] > 0)[0][0]

        # check if L is out of range for interpolation at this K:
        if L_extract < ax_L[idxL_outsidelc] or L_extract > ax_L[-1]:
            continue

        # get alpha from the map file, but ignore fill values:
        sol_alpha = np.interp(L_extract, ax_L[idxL_outsidelc:], map_alpha[0, :, idx_K][idxL_outsidelc:])

        # check the energy array
        if np.sum(sol_en[:, idxL_outsidelc:len(ax_L)] < 0.) > 0:
            # if there are any elements below 0 in sol_en:
            print("Error: energy arrays contain fill values")
            sys.exit(1)

        # error if energy is out of range for interpolation:
        for idxL in range(idxL_outsidelc, len(ax_L)):
            if en_all[-1] > sol_en[-1, idxL]:
                print("Error: energy {:.2f}MeV is out of bounds at alpha={:.2f} (iK={})".format(en_all[-1], sol_alpha,
                                                                                                idx_K + 1))
                sys.exit(1)
            elif en_all[0] < sol_en[0, idxL]:
                print("Error: energy {:.2f}MeV is out of bounds at alpha={:.2f} (iK={})".format(en_all[0], sol_alpha,
                                                                                                idx_K + 1))
                sys.exit(1)

        # L interpolation:
        L_idx1 = np.argmin(L_extract > ax_L[idxL_outsidelc:])
        if ax_L[idxL_outsidelc:][L_idx1] == L_extract:
            L_idx0 = L_idx1
            L_frac = 0
        else:
            L_idx0 = L_idx1 - 1
            L_frac = (L_extract - ax_L[idxL_outsidelc:][L_idx0]) / (
                        ax_L[idxL_outsidelc:][L_idx1] - ax_L[idxL_outsidelc:][L_idx0])

        if dynamic:  # interpolate to the current t we need:

            if time_plot < ax_t[0] or time_plot > ax_t[-1]:
                print("", "Error: time_plot is out of range on K idx", idx_K)
                sys.exit(1)

            idx_t_1 = np.argmin(time_plot > ax_t)
            if ax_t[idx_t_1] == time_plot:
                idx_t_0 = idx_t_1
                time_frac = 0
            else:
                idx_t_0 = idx_t_1 - 1
                time_frac = (time_plot - ax_t[idx_t_0]) / (ax_t[idx_t_1] - ax_t[idx_t_0])

            # get f at every L at the energy under investigation:
            # t0
            sol_f1d_t_0 = np.zeros((len(ax_L) - idxL_outsidelc, len(en_all)))  # preallocate empty array
            sol_f1d_t_1 = np.zeros((len(ax_L) - idxL_outsidelc, len(en_all)))  # preallocate empty array

            # for idxL in range(idxL_outsidelc, len(ax_L)):
            for idxL in [L_idx0 + idxL_outsidelc, L_idx1 + idxL_outsidelc]:
                sol_en_idx1_last = 0

                for idxe, en in enumerate(en_all):
                    # automatically interpolate in energy:
                    # sol_f1d_t_0.append(np.interp(en, sol_en[:, idxL], sol_f[idx_t_0, :, idxL]))

                    # manually interpolate in energy:
                    sol_en_idx1 = sol_en_idx1_last + np.argmin(en > sol_en[sol_en_idx1_last:, idxL])

                    if sol_en[sol_en_idx1, idxL] == en:
                        sol_en_idx0 = sol_en_idx1
                        frac = 0
                    else:
                        sol_en_idx0 = sol_en_idx1 - 1
                        frac = (en - sol_en[sol_en_idx0, idxL]) / (
                                    sol_en[sol_en_idx1, idxL] - sol_en[sol_en_idx0, idxL])

                    sol_f1d_t_0[idxL - idxL_outsidelc, idxe] = (sol_f[idx_t_0, sol_en_idx0, idxL] + frac * (
                                sol_f[idx_t_0, sol_en_idx1, idxL] - sol_f[idx_t_0, sol_en_idx0, idxL]))
                    sol_f1d_t_1[idxL - idxL_outsidelc, idxe] = (sol_f[idx_t_1, sol_en_idx0, idxL] + frac * (
                                sol_f[idx_t_1, sol_en_idx1, idxL] - sol_f[idx_t_1, sol_en_idx0, idxL]))

                    sol_en_idx1_last = sol_en_idx1  # save time on the next energy

            # get f at the L under investigation, interpolating for each energy at the same time:

            sol_f1d_t_0 = sol_f1d_t_0[L_idx0] + L_frac * (sol_f1d_t_0[L_idx1] - sol_f1d_t_0[L_idx0])

            sol_f1d_t_1 = sol_f1d_t_1[L_idx0] + L_frac * (sol_f1d_t_1[L_idx1] - sol_f1d_t_1[L_idx0])

            sol_f1d = sol_f1d_t_0 + time_frac * (sol_f1d_t_1 - sol_f1d_t_0)

        else:
            idx_t_0 = 0

            # t0
            sol_f1d_t_0 = np.zeros((len(ax_L) - idxL_outsidelc, len(en_all)))  # preallocate empty array

            # for idxL in range(idxL_outsidelc, len(ax_L)):
            for idxL in [L_idx0 + idxL_outsidelc, L_idx1 + idxL_outsidelc]:
                sol_en_idx1_last = 0

                for idxe, en in enumerate(en_all):
                    # automatically interpolate in energy:
                    # sol_f1d_t_0.append(np.interp(en, sol_en[:, idxL], sol_f[idx_t_0, :, idxL]))

                    # manually interpolate in energy:
                    sol_en_idx1 = sol_en_idx1_last + np.argmin(en > sol_en[sol_en_idx1_last:, idxL])

                    if sol_en[sol_en_idx1, idxL] == en:
                        sol_en_idx0 = sol_en_idx1
                        frac = 0
                    else:
                        sol_en_idx0 = sol_en_idx1 - 1
                        frac = (en - sol_en[sol_en_idx0, idxL]) / (
                                    sol_en[sol_en_idx1, idxL] - sol_en[sol_en_idx0, idxL])

                    sol_f1d_t_0[idxL - idxL_outsidelc, idxe] = (sol_f[idx_t_0, sol_en_idx0, idxL] + frac * (
                                sol_f[idx_t_0, sol_en_idx1, idxL] - sol_f[idx_t_0, sol_en_idx0, idxL]))

                    sol_en_idx1_last = sol_en_idx1  # save time on the next energy

            sol_f1d = sol_f1d_t_0[L_idx0] + L_frac * (sol_f1d_t_0[L_idx1] - sol_f1d_t_0[L_idx0])

        sol_f1d_allK.append(sol_f1d)
        sol_alpha_allK.append(sol_alpha)


    sol_f1d_allK = np.array(sol_f1d_allK)
    sol_f = np.array([row[::-1] for row in sol_f1d_allK.T])
    if not len(sol_f):
        sol_f = [[] for _ in range(len(en_all))]
    sol_alpha_allK = np.array(sol_alpha_allK[::-1])

    if not multiple:
        sol_f = sol_f[0]

    return sol_alpha_allK, sol_f


import random
def get_jomni_dummy(cdf, time, L_, BBe, E_all):
    #returns random number for testing purposes (much faster to execute during test)
    if L_ != L_: return [np.nan]
    return [10**random.randint(1,5)]

def get_jomni(cdf_loaded, time, L_, BBe, E_all):
    #returns several estimate of local omnidirectional flux from an equatorial pitch angle distribution
    if BBe < 1: BBe = 1.0


    #get equatorial pitch angle distribution:
    #
    #   we are interpolating at L between model outputs
    #   so we may not be able to determine the loss cone from getpad(...)
    #
    # #getpad_(cdf, L_extract, en, dynamic, time_plot)
    aeq_PAD_d, f_all = getpad_(cdf_loaded, L_, E_all, True, time)
    for f_ in f_all:
        if not len(f_):
            #print("","warning: satellite outside numerical boundary, L={:.2f}, BBe={:.2f}".format(L_, BBe))
            return [0]*len(E_all)

    aeq_PAD_r = np.radians(aeq_PAD_d)
    j_all = [f2j(E_,f_) for E_, f_ in zip(E_all, f_all)] #directional flux in standard units
    # E_, f_ is a scalar, list pair
    


    #solve for the equatorial pitch angle of local 90 at the spacecraft:
    aeq_local90 = asin(1/sqrt(BBe))*180/np.pi

    #solve for the local pitch angle at the spacecraft of the smallest equatorial pitch angle outside the loss cone:
    y0eq_rBBe = sin(aeq_PAD_r[0]) * sqrt(BBe)
    if y0eq_rBBe > 1: #satellite is in loss cone
        print("","warning: satellite in numerical loss cone, L={:.2f}, BBe={:.2f}".format(L_, BBe))
        return [0]*len(E_all)
    al_minaeq = asin(y0eq_rBBe)*180/np.pi


    #define the optimal integration range in terms of local pitch angle:
    al_n = 30
    al_range_d = np.linspace(al_minaeq, 90, al_n)
    al_range_r = np.radians(al_range_d)



    # find aeq of this range:
    aeq_range_r = np.arcsin(np.sin(al_range_r)/sqrt(BBe))
    aeq_range_d = np.degrees(aeq_range_r)
    # alpha_lc_local = np.arcsin(np.sin(alpha_lc*pi/180)*sqrt(BBe))*180/pi

    #just to make sure numerical error doesn't cause us to extrapolate:
    aeq_range_r[0] = aeq_PAD_r[0]
    aeq_range_d[0] = aeq_PAD_d[0]

    #interpolate j at each aeq:
    pa_dist_r = interp1d(aeq_PAD_r, j_all, axis = 1)
    j_local_raw_all = pa_dist_r(aeq_range_r)



    
    #integrate directional flux to get omnidirectional flux:
    # only include flux above the loss cone
    # trapezoid rule:
    j_integrate = j_local_raw_all
    y_range = np.expand_dims(np.sin(al_range_r), axis=0)

    integrand = j_integrate*2*pi*y_range.repeat(len(E_all),axis=0)
    
    jomni = 2 * np.trapz(integrand, al_range_r, axis = 1)


    return jomni# [jomni_fit1, jomni_fit2, jomni_fit3, jomni, jomni_taper]
