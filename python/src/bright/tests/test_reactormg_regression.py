import os
import re

import numpy as np
import tables as tb

from bright import *
import isoname
from mass_stream import MassStream

import serpent
import matplotlib.pyplot as plt
from matplotlib import rc

rc('text', usetex=True)
rc('font', family='roman')

from char.m2py import convert_res, convert_dep
from metasci.graph import StairStepEnergy
import metasci.nuke as msn

# Hack to use origen as burnup calculator
from origen_reactormg import OrigenReactorMG
ReactorMG = OrigenReactorMG

def run_serpent():
    #bu_file = 'serp_bu_10'
    #xs_file = 'serp_xs_10'

    bu_file = 'serp_bu_19'
    xs_file = 'serp_xs_19'

    if (not os.path.exists(bu_file + '_res.py')) or (not os.path.exists(bu_file + '_dep.py')):
        serpent.main(bu_file)
        convert_res(bu_file + '_res.m')
        convert_dep(bu_file + '_dep.m')

    if (not os.path.exists(xs_file + '_res.py')):
        serpent.main(xs_file)
        convert_res(xs_file + '_res.m')

    res_bu = {}
    execfile(bu_file + '_res.py', {}, res_bu)

    dep_bu = {}
    execfile(bu_file + '_dep.py', {}, dep_bu)

    res_xs = {}
    execfile(xs_file + '_res.py', {}, res_xs)

    #
    # Renomralize Mass
    #
    mw_conversion = 270.0 / (237.0 * dep_bu['TOT_VOLUME'] * 10.7)
    mw = dep_bu['TOT_MASS'] * mw_conversion

    iso_LL = {}
    iso_index = {}
    for iso_zz in dep_bu['ZAI']:
        # Find valid isotope indeces
        try:
            iso_LL[iso_zz] = isoname.mixed_2_LLAAAM(int(iso_zz))
        except:
            continue
        iso_index[iso_zz] = dep_bu['i{0}'.format(iso_zz)] - 1

    # Caclulate actual mass of isotopes present
    mass = mw[iso_index.values()].sum(axis=0)
    mass = mass / mass[0]

    dep_bu['mw'] = mw
    dep_bu['mass'] = mass
    dep_bu['iso_index'] = iso_index

    return res_bu, dep_bu, res_xs


def run_reactormg():
    # Init bright
    libfile = os.getenv("BRIGHT_DATA") + 'lwr_mg.h5'
    load_track_isos_hdf5(libfile)
    bright_config.write_text = False
    bright_config.write_hdf5 = False
    bright_config.verbosity = 100

    # Init reactor paramters
    rp = lwr_defaults()
    rp.batches = 3
    rp.flux = 4*(10**14)

    rp.fuel_form = {"IHM": 1.0, "O16": 2.0}
    rp.cladding_form = {"ZR93": 0.5, "ZR95": 0.5}
    rp.coolant_form = {"H1": 2.0, "O16": 1.0}

    rp.fuel_density = 10.7
    rp.cladding_density = 5.87
    rp.coolant_density = 0.73

    rp.pnl = 0.92
    #rp.pnl = 0.98
    #rp.pnl = 1.0
    rp.BUt = 50.0
    rp.use_disadvantage_factor = True
    #rp.lattice_type = 'Spherical'
    rp.lattice_type = 'Cylindrical'
    #rp.lattice_type = 'Planar'
    rp.rescale_hydrogen = True
    rp.burnup_via_constant = 'power'

    rp.fuel_radius = 0.412
    rp.void_radius = 0.4205
    rp.clad_radius = 0.475
    rp.unit_cell_pitch = 1.33

    rp.open_slots = 25
    rp.total_slots = 289

    bt = list(np.linspace(0.0, 365.0, 10))
    #bt.insert(1, 1.0)
    bt.insert(1, 4.0)
    bt = np.array(bt)
    rp.burn_times = bt

    # Init mass stream
    leu = MassStream({922340: 0.01, 922350: 0.05, 922380: 0.94})

    # Init ReactorMG
    rmg = ReactorMG(reactor_parameters=rp, name="rmg")
    rmg.loadlib(libfile)

    # Run the reactor
    #rmg.calc(leu)
    rmg.ms_feed = leu
    rmg.burnup_core()
    return rmg


def test_regression():
    res_bu, dep_bu, res_xs = run_serpent()
    rmg = run_reactormg()
    return rmg, res_bu, dep_bu, res_xs


def make_graphs(r, s, diff, serr=None, name=""):
    global burn_times
    plt.clf()

    # make the compare plot
    if serr is None:
        plt.plot(burn_times, s, 'k-', label="Serpent")
    else:
        plt.errorbar(burn_times, s, serr, color='k', label="Serpent")
    
    plt.errorbar(burn_times, r, diff*r, color='r', label="RMG")

    plt.xlabel("burn time [days]")
    plt.ylabel(name)
    plt.ticklabel_format(scilimits=(-5, 5))
    plt.legend(loc=0)
    plt.savefig(name.split('[')[0].replace(' ', '_') + '.png')
    plt.savefig(name.split('[')[0].replace(' ', '_') + '.eps')
    plt.clf()

    plt.hist(burn_times, weights=diff, align='mid', color='g', rwidth=0.95)
    plt.xlabel("burn time [days]")
    plt.ylabel(name + " Relative Error")
    plt.savefig(name.split('[')[0].replace(' ', '_') + '_rel_err.png')
    plt.savefig(name.split('[')[0].replace(' ', '_') + '_rel_err.eps')
    plt.clf()


def make_flux_graphs(r, s, diff, serr=None, name=""):
    global E_g
    plt.clf()

    u_g = np.log(E_g[0] / E_g)
    u_diff = u_g[1:] - u_g[:-1]
    u_norm = u_diff[0] / u_diff

    # make the compare plot
    gkw = {'datalabel': "Serpent",
           'colorline': 'k-',
           'scale': 'log',
           'ylabel': name,
           'write': False, 
           'show': False, 
           'ymin': 0.0,
           'ymax': 0.25,
           }
    StairStepEnergy((s * u_norm)[::-1], E_g[::-1], **gkw)

    gkw['datalabel'] = "RMG"
    gkw['colorline'] = 'r-'
    gkw['write'] = True 
    gkw['name'] = name.split('[')[0].replace(' ', '_')
    StairStepEnergy((r * u_norm)[::-1], E_g[::-1], **gkw)    

    """plt.hist(-1 * (E_g[1:] - E_g[:-1])[::-1], weights=diff[::-1], align='mid', color='g', rwidth=0.95)
    plt.xscale('log')
    plt.xlabel("Energy [MeV]")
    plt.ylabel(name + " Relative Error")
    plt.savefig(name.split('[')[0].replace(' ', '_') + '_rel_err.png')
    plt.savefig(name.split('[')[0].replace(' ', '_') + '_rel_err.eps')
    plt.clf()

    plt.hist(-1 * (E_g[1:] - E_g[:-1])[::-1], weights=diff[::-1] * r[::-1], align='mid', color='g', rwidth=0.95)
    plt.xscale('log')
    plt.xlabel("Energy [MeV]")
    plt.ylabel(name + " Flux Weighted Relative Error")
    plt.savefig(name.split('[')[0].replace(' ', '_') + '_fw_rel_err.png')
    plt.savefig(name.split('[')[0].replace(' ', '_') + '_fw_rel_err.eps')
    plt.clf()"""


def calc_diff(r, s, serr=None, name=""):
    print "Summary for {0}:".format(name)
    print "Reactor: "
    print repr(r)
    print "Serpent: "
    print repr(s)
    print "Fractional Diff: "
    diff = r / s - 1.0
    diff[np.isnan(diff)] = 0.0
    print repr(diff)
    print

    try:
        make_graphs(r, s, diff, serr, name)
    except Exception as e:
        print e

    return r, s, diff

def compare_1g_xs(rmg):
    global time_range
    f = tb.openFile('benchmark_xs.h5', 'r')

    isos_LL = f.root.transmute_isos_LL.read()
    isos_zz = f.root.transmute_isos_zz.read()
    isos = zip(isos_LL, isos_zz)

    r_norm_phi = rmg.phi_tg / rmg.phi_t[:, np.newaxis]
    r_norm_phi = r_norm_phi[1:]

    reactions = ['sigma_t', 'sigma_f', 'sigma_gamma', 'sigma_2n', 'sigma_a', 'sigma_s', 
                 'sigma_alpha', 'sigma_gamma_x', 'sigma_2n_x']

    sig = {}
    for rx in reactions:
        if rx == 'sigma_s':
            r_sig = getattr(rmg, rx + '_itgh')
            #r_sig = {key: value.sum(axis=-1) for key, value in r_sig.items()}
            r_sig = {key: value.sum(axis=-2) for key, value in r_sig.items()}
        else:
            r_sig = getattr(rmg, rx + '_itg')
        s_sig = getattr(f.root, rx)

        for iso_LL, iso_zz in isos:
            r_xs = (r_sig[iso_zz][1:] * r_norm_phi).sum(axis=1)
            s_xs = (getattr(s_sig, iso_LL) * r_norm_phi).sum(axis=1)
            diff = r_xs / s_xs - 1.0
            diff[np.isnan(diff)] = 0.0
            sig[rx, iso_LL] = (r_xs, s_xs, diff)

    f.close()

    return sig

def sort_sig(sig):
    key_func = lambda x: abs(x[1][2]).max()
    s = sorted(sig.items(), key=key_func, reverse=True)
    return s


def make_1g_xs_graphs(nuc, sig):
    global burn_times
    nuc_zz = isoname.LLAAAM_2_zzaaam(nuc)
    plt.clf()

    reactions = ['sigma_t', 'sigma_s', 'sigma_a', 'sigma_f']
    markers = ['-', 'o', 's', 'x']
    ncol = 2

    if nuc_zz < 860000:
        reactions = reactions[:-1]
        markers = markers[:-1]
        ncol = 3

    for reaction, marker in zip(reactions, markers):
        r, s, diff = sig[reaction, nuc]
        plt.plot(burn_times, s, color='k', marker=marker, linestyle='-', label="Serpent $\\{0}$".format(reaction))
        plt.errorbar(burn_times, r, diff*r, color='r', marker=marker, linestyle='-', label="RMG $\\{0}$".format(reaction))

    plt.xlabel("burn time [days]")
    plt.ylabel(nuc + " One-Group Cross Sections [barns]")
    plt.ticklabel_format(scilimits=(-5, 5))
    plt.legend(loc=0, ncol=ncol)
    plt.savefig(nuc + '_1g_xs.png')
    plt.savefig(nuc + '_1g_xs.eps')
    plt.clf()


def make_rank_table(reaction, nuc_class='Actinide', nrows=20, hl_cutoff=86400.0):
    global hl, sig, u

    imp_acts = set(['U234',  'U235',  'U236',   'U238',  'NP237', 'PU238', 'PU239', 'PU240', 'PU241', 'PU242', 
                    'AM241', 'AM242', 'AM242M', 'AM243', 'AM244', 'CM242', 'CM243', 'CM244', 'CM245', 'CM246', ])

    if nuc_class == 'Actinide':
        #filt_u = [(key, value) for key, value in u if key[0] == reaction and 890000 < isoname.LLAAAM_2_zzaaam(key[1]) and hl_cutoff < hl[key[1]]]
        filt_u = [(key, value) for key, value in u if key[0] == reaction and key[1] in imp_acts and hl_cutoff < hl[key[1]]]
    elif nuc_class == 'Fission Product':
        filt_u = [(key, value) for key, value in u if key[0] == reaction and 200000 < isoname.LLAAAM_2_zzaaam(key[1]) < 730000 and hl_cutoff < hl[key[1]]]
    else:
        raise ValueError

    latex_table = ("\\begin{table}[htbp]\n"
                   "\\begin{center}\n")
    latex_table += "\\caption{{Maximum {0} $\\{1}$ Relative Error}}\n".format(nuc_class, reaction)
    latex_table += "\\label{{rank_{0}_{1}_table}}\n".format(nuc_class.replace(" ", "_"), reaction)
    latex_table +=("\\begin{tabular}{|l|c|}\n"
                   "\\hline\n"
                   "\\textbf{Nuclide} & \\textbf{$\\varepsilon$} \\\\\n"
                   "\\hline\n")

    nuc_latex = "\\nuc{{{0}}}{{{1}}}"
    nuc_pattern = re.compile("([A-Z]{1,2})(\d{1,3})(M?)")
    for key, value in filt_u[:nrows]:
        rx, nuc_LL = key
        r, s, diff = value

        z, a, m = nuc_pattern.match(nuc_LL).groups()
        nl = nuc_latex.format(z.capitalize(), a)
        if 0 < len(m):
            nl += "\\superscript{*}"
        latex_table += nl + " & {0:+.4F} \\\\\n".format(value[2][abs(value[2]).argmax()])

    latex_table += ("\\hline\n"
                    "\\end{tabular}\n"
                    "\\end{center}\n"
                    "\\end{table}\n")

    fname = "rank_table_{0}_{1}.tex".format(nuc_class.replace(" ", "_"), reaction)
    with open(fname, 'w') as f:
        f.write(latex_table)


def make_nn_table(nnt, burn_times):
    len_p = len(nnt[0])

    latex_table = ("\\begin{table}[htbp]\n"
                   "\\begin{center}\n")
    latex_table += "\\caption{Nearest Neighbors over Burn}\n"
    latex_table += "\\label{nn_table}\n"
    latex_table += "\\tiny\n"
    latex_table += "\\begin{tabular}{|l||" + ("c"*len_p) + "|}\n"
    latex_table += "\\hline\n"
    latex_table += "\\textbf{days} & \\multicolumn{" + str(len_p) + "}{|c|}{\\textbf{$p^*$}} \\\\\n"
    latex_table += "\\hline\n"

    for s, bt in enumerate(burn_times):
        latex_table += "{0:.3G}".format(bt) + " & "
        latex_table += " & ".join([str(p) for p in nnt[s]])
        latex_table += " \\\\\n"

    latex_table += ("\\hline\n"
                    "\\end{tabular}\n"
                    "\\end{center}\n"
                    "\\end{table}\n")

    with open("nearest_neighbor_table.tex", 'w') as f:
        f.write(latex_table)
    
if __name__ == "__main__":
    rmg, res_bu, dep_bu, res_xs = test_regression()

    burn_times = rmg.burn_times
    time_range = [0] + range(2, len(burn_times))
    burn_times = burn_times[time_range]

    E_g = rmg.E_g

    """\
    """
    r_k, s_k, diff_k = calc_diff(rmg.k_t[1:], res_bu['SIX_FF_KEFF'][1:, 0], res_bu['SIX_FF_KEFF'][1:, 1], "k")

    r_norm_phi = rmg.phi_tg / rmg.phi_t[:, np.newaxis]
    s_norm_phi = res_bu['FLUX'][:, 2::2] / res_bu['FLUX'][:, np.newaxis, 0]
    serr_phi = res_bu['FLUX'][:, 3::2] / res_bu['FLUX'][:, np.newaxis, 0]
    r_phi, s_phi, diff_phi = calc_diff(r_norm_phi[1:], s_norm_phi[1:], "Normalized Flux")

    for t in range(len(burn_times)):
        print "Making flux figure for time t = {0}".format(t)
        make_flux_graphs(r_phi[t], s_phi[t], diff_phi[t], serr=serr_phi[t], name="Normalized Flux at {0} days".format(int(burn_times[t])))


    r_total, s_total, diff_total = calc_diff(rmg.Sigma_t_fuel_tg[0], res_xs['TOTXS'][0, 2::2], res_xs['TOTXS'][0, 3::2], name="Total XS")
    make_flux_graphs(r_total, s_total, diff_total, name="Total XS")

    r_fiss, s_fiss, diff_fiss = calc_diff(rmg.Sigma_f_fuel_tg[0], res_xs['FISSXS'][0, 2::2], res_xs['FISSXS'][0, 3::2], "Fission XS")
    make_flux_graphs(r_fiss, s_fiss, diff_fiss, name="Fission XS")

    r_abs, s_abs, diff_abs = calc_diff(rmg.Sigma_a_fuel_tg[0], res_xs['ABSXS'][0, 2::2], res_xs['ABSXS'][0, 3::2],"Absorption XS")
    make_flux_graphs(r_abs, s_abs, diff_abs, name="Absorption XS")

    r_gamma, s_gamma, diff_gamma = calc_diff(rmg.Sigma_gamma_fuel_tg[0], res_xs['CAPTXS'][0, 2::2], res_xs['CAPTXS'][0, 3::2], "Capture XS")
    make_flux_graphs(r_gamma, s_gamma, diff_gamma, name="Capture XS")

    T_it = rmg.T_it

    bu_mass_norm = np.array([sum([dep_bu['mw'][dep_bu['iso_index'][key]][n] for key in dep_bu['iso_index'].keys() if 860000 < key]) for n in range(dep_bu['mw'].shape[1])])
    dep_bu['mw'] = dep_bu['mw'] / bu_mass_norm[0]
    bu_mass_norm = bu_mass_norm / bu_mass_norm[0]

    # Patch known isos
    dep_bu['mw'][dep_bu['iso_index'][922340]] *= (0.01 / dep_bu['mw'][dep_bu['iso_index'][922340]][0])
    dep_bu['mw'][dep_bu['iso_index'][922380]] *= (0.94 / dep_bu['mw'][dep_bu['iso_index'][922380]][0])

    nuclides = ['U234',  'U235',  'U236',  'U238',  'NP237', 'PU238', 'PU239', 'PU240', 'PU241', 'PU242', 
                'AM241', 'AM242', 'AM243', 'CM242', 'CM243', 'CM244', 'CM245', 'CM246', 'SE79',  'KR85',  
                'SR90',  'ZR93',  'TC99',  'I129',  'PD107', 'CS134', 'CS135', 'CS137', 'SM151', 'EU155']

    # Make mass fraction figures
    for nuc_LL in nuclides:
        nuc_zz = isoname.LLAAAM_2_zzaaam(nuc_LL)
        nuc_ind = dep_bu['iso_index'][nuc_zz]
        r_i, s_i, diff_i = calc_diff(T_it[nuc_zz][time_range], dep_bu['mw'][nuc_ind][time_range], name=nuc_LL + " Mass Fraction [kg/kgIHM]")

    mss = [MassStream({i: T_it[i][t] for i in T_it.keys()}) for t in time_range]
    r_mass, s_mass, diff_mass = calc_diff(np.array([ms.mass for ms in mss]), bu_mass_norm[time_range], name="Mass")

    sig = compare_1g_xs(rmg)
    u = sort_sig(sig)
    
    print "Max fractional deviations"
    for key, value in u[:20]:
        rx, iso = key
        r, s, diff = value
        print iso, rx, max(abs(diff))

    print
    print

    # Make 1g xs figs
    for nuc in nuclides:
        print "Making cross section figures for {0}".format(nuc)
        make_1g_xs_graphs(nuc, sig)
    

    print 
    print

    # Load in half-lives
    with tb.openFile(msn.nuc_data, 'r') as nd:
        hl = {row['from_iso_LL']: row['half_life'] for row in nd.root.decay}

    # make rank tables    
    make_rank_table('sigma_s', 'Actinide')
    make_rank_table('sigma_s', 'Fission Product')
    make_rank_table('sigma_a', 'Actinide')
    make_rank_table('sigma_a', 'Fission Product')
    make_rank_table('sigma_f', 'Actinide')


    # Nearest neighbor calc
    nnt = np.array(rmg._nearest_neighbors[1:]) + 1

    nnt[5 < nnt] += 10
    nnt[nnt < 6] += 5
    nnt = np.append(nnt, nnt-5, axis=1)

    make_nn_table(nnt, burn_times)

