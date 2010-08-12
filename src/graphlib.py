from pylab import *
from math import log

GoodBW = ['k-', 'r-', 'k--', 'r--', 'k-.', 'r-.', 'k:', 'r:', 'k.', 'r.', 'k+', 'r+', 'kx', 'rx']
GoodBWlen = len(GoodBW)


#####################
### ORIGEN Graphs ###
#####################

def GraphForTime(t, Data, G=0, units="[unitless]", figname=""):
    "Plots the multigroup data as a function of time t."

    Data = libcom.MatrixSwitchDim(Data)
    if G == 0:
        G = len(Data)
    gmod = G / GoodBWlen 

    for n_g in range(G):
        if G <= GoodBWlen:
            plot(t, Data[-(n_g+1)], GoodBW[n_g], label='Group %s'%(n_g+1))
        elif (GoodBWlen < G) and (0 == n_g%gmod) and (n_g < (GoodBWlen*gmod)):
            plot(t, Data[-(n_g+1)], GoodBW[n_g], label='Group %s'%(n_g+1))
        else:
            continue

    xlabel('Time [days]')
    ylabel( "%s %s"%(figname, units) )

    legend(loc = 0)

    filename = figname.replace(' ', '_').replace('/', '')
    savefig(filename + ".eps")
    savefig(filename + ".png")

    clf()
    return filename

def TransformForTime(t, Ti, min_lim=1.0, g=0, figname=""):
    "Plots a single group of a transmutation matrix."

    #Pick out the gth Group
    E_Ti = {}
    for j in Ti.keys():
        E_Ti[j] = libcom.MatrixSwitchDim(Ti[j])[-g]

    max_dic = {}
    for j in E_Ti.keys():
        max_dic[j] = max(E_Ti[j])

    graph_order = []
    for v in reversed(sorted(max_dic.values())):
        if v < min_lim:
            break
        if len(graph_order) == GoodBWlen:
            break

        for j in max_dic.keys():
            if (max_dic[j] == v):
                graph_order.append(j)
                break

    to_graph = []
    for j in graph_order:
        to_graph.append( (j, []) )
        for n_t in range(len(t)):
            if E_Ti[j][n_t] < min_lim:
                to_graph[-1][1].append( min_lim * 0.01 )
            else:
                to_graph[-1][1].append( E_Ti[j][n_t] )

    for n in range(len(to_graph)):
        semilogy(t, to_graph[n][1], GoodBW[n], label=to_graph[n][0])
                
    #Set Axes Bounds...if log'd
    ax = axis()
    axis([ax[0], ax[1], min_lim * 0.1, ax[3]]) 

    xlabel('Time [days]')
    ylabel( "%s Mass for Group %s [kg/kgIHM]"%(figname, g) )

    legend(loc = 0)

    filename = "%s Mass E%s"%(figname, g)
    filename = filename.replace(' ', '_').replace('/', '')
    savefig(filename + ".eps")
    savefig(filename + ".png")

    clf()
    return filename
