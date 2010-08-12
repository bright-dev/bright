######################
### CHAR Libraries ###
######################
from char import *

def Make_Figs_MCNP(q=False):
    "Make MCNP Graphs, for fun and learning"
    t1 = time.time()
    if not ( 'figs' in os.listdir('.') ):
        os.mkdir('figs/')
    os.chdir('figs/')
    for fig in os.listdir('.'):
        metasci.SafeRemove(fig)

    libfile = tb.openFile('../libs/' + reactor + '.h5', 'r')
    lfr = libfile.root

    #Flux Graphs
    QPrint("\033[1;32mNow making Flux graphs...\033[0m")
    msg.SimpleGraph(lfr.Coarse.time, lfr.Coarse.flux, xlabel="Burn Time [days]", 
        ylabel=r"Flux [n/s/cm$^2$]", scale='linear', name='Flux' )

    msg.SimpleGraph(lfr.CINDER.time, lfr.CINDER.flux, xlabel="CINDER Burn Time [days]", 
        ylabel=r"CINDER Flux [n/s/cm$^2$]", scale='linear', name='Flux_Cinder' )

    #Group Flux Graphs
    QPrint("\033[1;32mNow making Group Flux graphs...\033[0m")
    group_flux_kwg = {
        'scale':  'log',
        'name':   'Group_Flux',
        'ylabel': r"Group Flux [n/s/cm$^2$]",
        'G':      lfr.E_up.nrows-1,
        'write':  False,
        'ymin':  10.0**8,
        }

    for t in range(len(lfr.Coarse.flux_g)):
        if t == len(lfr.Coarse.flux_g)-1:
            group_flux_kwg['write'] = True
        msg.StairStepEnergy(lfr.Coarse.flux_g[t], lfr.E_up, **group_flux_kwg)

    group_flux_kwg['write']  = False
    group_flux_kwg['G']      = lfr.CINDER.E_up.nrows-1
    group_flux_kwg['name']   = 'Group_Flux_Cinder'
    group_flux_kwg['ylabel'] = r"CINDER Group Flux [n/s/cm$^2$]"
    for t in range(len(lfr.CINDER.flux_g)):
        if t == len(lfr.CINDER.flux_g)-1:
            group_flux_kwg['write'] = True
        msg.StairStepEnergy(lfr.CINDER.flux_g[t], lfr.CINDER.E_up, **group_flux_kwg)

    #Cross-Section Graphs
    xs_kwg = {
        'name':   'XS',
        'ylabel': "XS [barns]",
        'scale':  'log',
        'G':      lfr.E_up.nrows-1,
        'write':  False,
        'energy_bins': lfr.E_up,
        }
    RootGroups = lfr.Coarse._v_groups
    for cgkey in RootGroups.keys():
        if not (cgkey[:5] == "sigma"):
            continue
        QPrint("\033[1;32mNow making %s graphs...\033[0m"%cgkey)
        sigma = RootGroups[cgkey]._v_leaves
        for iso in sigma.keys():
            xs_kwg['name']   = '{0}_{1}'.format(iso, cgkey)
            xs_kwg['ylabel'] = r'{0} $\{1}}}}}$ [barns]'.format(iso, cgkey).replace("_", "_{\mbox{")
            for t in range(len(sigma[iso])):
                if t == len(sigma[iso])-1:
                    xs_kwg['write'] = True
                msg.StairStepEnergy(sigma[iso][t], **xs_kwg)

    #Chi Graphs
    QPrint("\033[1;32mNow making Chi graphs...\033[0m")
    chi_kwg = {
        'name':   'chi',
        'ylabel': r"$\chi$",
        'scale':  'logx',
        'G':      lfr.E_up.nrows-1,
        'write':  False,
        'energy_bins': lfr.E_up,
        }
    chi = lfr.Coarse.chi._v_leaves
    for iso in chi.keys():
        chi_kwg['name']   = '{0}_{1}'.format(iso, "chi")
        chi_kwg['ylabel'] = r'{0} $\chi$'.format(iso)
        for t in range(len(chi[iso])):
            if t == len(chi[iso])-1:
                chi_kwg['write'] = True
            msg.StairStepEnergy(chi[iso][t], **chi_kwg)

    #Nubar Graphs
    QPrint("\033[1;32mNow making NuBar graphs...\033[0m")
    nubar_kwg = {
        'name':   'nubar',
        'ylabel': r"$\bar{\nu}$",
        'scale':  'log',
        'G':      lfr.E_up.nrows-1,
        'write':  False,
        'energy_bins': lfr.E_up,
        }
    nubar = lfr.Coarse.nubar._v_leaves
    for iso in nubar.keys():
        nubar_kwg['name']   = '{0}_{1}'.format(iso, "nubar")
        nubar_kwg['ylabel'] = r'{0} $\bar{{\nu}}$'.format(iso)
        for t in range(len(nubar[iso])):
            if t == len(nubar[iso])-1:
                nubar_kwg['write'] = True
            msg.StairStepEnergy(nubar[iso][t], **nubar_kwg)

    libfile.close()
    os.chdir('../') #Back to 'reactor' root directory
    t2 = time.time()
    QPrint( "\n\033[1;32mMCNP Graphs generated in %G minutes.\033[0m"%( (t2-t1)/60.0), q )

    return

def Make_Figs_ORIGEN(q=False):
    """Make ORIGEN Graphs, for more fun and now with more learning"""
    t1 = time.time()
    t1 = time.time()
    if not ( 'figs' in os.listdir('.') ):
        os.mkdir('figs/')
    os.chdir('figs/')

    #Try opening the HDF5 File
    try:
        libfile = tb.openFile('../libs/' + reactor + '.h5', 'r')
    except:
        QPrint("\033[1;31mCould not find a valid HDF5 file!  Looking for: %s.h5\033[0m"%reactor, q)
        return
    lfr  = libfile.root
    lfrf = lfr.Fine

    G = lfr.E_up.nrows-1 

    #Make Burnup Graphs
    try:
        BU = lfrf.Burnup._v_leaves
    except:
        QPrint("\033[1;31mHDF5 file %s.h5 does not contain ORIGEN data!\033[0m"%reactor, q)
        return
    QPrint("\033[1;32mNow making Burnup graphs...\033[0m")
    for iso in BU.keys():
        graphlib.GraphForTime(lfrf.time, BU[iso], G=G, units = "[MWd/kgIHM]", figname='%s Burnup'%iso ) 

    #Make k graphs
    k = lfrf.k._v_leaves
    QPrint("\033[1;32mNow making Multiplication Factor graphs...\033[0m")
    for iso in k.keys():
        graphlib.GraphForTime(lfrf.time[1:], k[iso][1:], G=G, figname='%s Multiplication Factor'%iso ) 

    #Make Production Rate grpahs
    Pro = lfrf.Production._v_leaves
    QPrint("\033[1;32mNow making Neutron Production Rate graphs...\033[0m")
    for iso in Pro.keys():
        graphlib.GraphForTime(lfrf.time[1:], Pro[iso][1:], G=G, figname='%s Neutron Production Rate'%iso ) 

    #Make Destruction Rate grpahs
    Des = lfrf.Destruction._v_leaves
    QPrint("\033[1;32mNow making Neutron Destruction Rate graphs...\033[0m")
    for iso in Des.keys():
        graphlib.GraphForTime(lfrf.time[1:], Des[iso][1:], G=G, figname='%s Neutron Destruction Rate'%iso ) 

    #Make Transmutation rate Graphs
    QPrint("\033[1;32mNow making Transmutation graphs...\033[0m")
    CutOff = 10.0**-3	#Concentration Cutoff
    T = lfrf.Transmutation
    for g in range(G, 0, -1):
        for i in lfr.CoreLoad_LLAAAM:
            Ti = {}
            for j in lfr.CoreTran_LLAAAM:
                try:
                    Ti[j] = T._f_getChild('%s/%s'%(i,j))
                except:
                    continue
            graphlib.TransformForTime(lfrf.time, Ti, min_lim=CutOff, g=g, figname=i)

    libfile.close()
    os.chdir('../') #Back to 'reactor' root directory
    t2 = time.time()
    QPrint( "\n\033[1;32mORIGEN Graphs generated in %G minutes.\033[0m"%( (t2-t1)/60.0), q )
    return
