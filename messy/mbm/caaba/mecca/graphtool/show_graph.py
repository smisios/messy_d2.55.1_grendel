#!/usr/bin/env python3
# -*- coding: utf-8 -*- Time-stamp: <2019-06-03 13:57:52 sander>

##############################################################################

# import os
# import subprocess
# import datetime
import sys
# import numpy as np
import re
import graph_tool.all as gt
from utils_graph import n2v
import matplotlib.pyplot as plt

HLINE =  "-" * 78
DEBUG = 1

##############################################################################

def create_graphviz_draw(g):
    gt.graphviz_draw(g, # http://www.graphviz.org/doc/info/attrs.html
                     #pos=pos,
                     gprops={#'concentrate':'true',
                         #'overlap':'prism0',
                         #'overlap':'scalexy',
                         'overlap':'voronoi', 
                         #'rankdir':'LR',
                         'splines':'curved',
                         #'splines':'ortho',
                         #'splines':'polyline',
                         #'splines':'spline',
                         'size':'40,20'},
                     vprops={'shape':'oval',
                             'fontsize':20,#gt.prop_to_size(x, mi=1, ma=30),
                             #'color':'Black', # circle around vertex
                             #'colorscheme':'accent8', # doesn't work?
                             #'style':'filled',
                             #'label':g.vp.name,
                             'label':g.vp.label},
                     eprops={'label':g.ep.eqntag,
                             'dir':'forward', # forward arrows
                             'arrowhead':'normal',
                             'arrowsize':1,
                             'fontsize':15,
                             #'style':'dashed',
                             #'style':'dotted',
                             #'color':'black;0.25:red;0.5:yellow;0.25',
                             'color':'blue',
                             'penwidth':1},
                     overlap='compress',
                     vcolor=g.vp.fillcolor, # vertex fill color
                     output='mecca_graph.pdf')

##############################################################################

def create_max_flow(g):
    
    g.ep.cap = g.new_edge_property("double")
    rxnrates = read_rxn_rates('caaba_mecca_rr.nc')
    for e in g.edges():
        eqntag = g.ep.eqntag[e]
        rxnrate = rxnrates[g.ep.eqntag[e]]
        prodstoic = g.ep.prodstoic[e] # stoic factor of product
        print('%9s: %7g %g' % (eqntag, prodstoic, rxnrate))
        g.ep.cap[e] = prodstoic*rxnrate

    # choose one:
    res = gt.edmonds_karp_max_flow(g, src, tgt, g.ep.cap)
    #res = gt.push_relabel_max_flow(g, src, tgt, g.ep.cap)
    #res = gt.boykov_kolmogorov_max_flow(g, src, tgt, g.ep.cap)
    res.a = g.ep.cap.a - res.a  # the actual flow
    #print res.a
    max_flow = sum(res[e] for e in tgt.in_edges())
    print('max flow: %g' % (max_flow))
    edgewidth = gt.prop_to_size(res, mi=0.1, ma=5, power=0.3)
    #-------------------------------------------------------------------------
    print() ; print() ; print('Sorted listing of important reactions from %s to %s:' % (
        g.vp.name[src], g.vp.name[tgt]))
    # put all info into mylist:
    mylist = []
    for e in g.edges():
        if (res[e]>0):
            mylist.append(
                [g.ep.cap[e], res[e], edgewidth[e], g.ep.eqntag[e],
                 g.vp.name[e.source()], g.vp.name[e.target()]])
    # sort and print mylist:
    print('            rxn rate                 flow            edgewidth     eqntag')
    #for myitem in sorted(mylist, reverse=True, key=lambda tup: tup[0]):
    for myitem in sorted(mylist, reverse=True):
        print('%20s %20s %20s %10s %s -> %s' % (
            myitem[0], myitem[1], myitem[2], myitem[3], myitem[4], myitem[5]))
    #-------------------------------------------------------------------------
    # # filter out all edges that contribute < 1 %
    # for e in g.edges():
    #     if (res[e] < max_flow/100.):
    #         g.ep.myfilter[e] = False
    #     else:
    #         g.ep.myfilter[e] = True
    #     print g.ep.myfilter[e], res[e], g.ep.eqntag[e], g.ep.reaction[e]
    # g.set_edge_filter(g.ep.myfilter) # keep if True
    #-------------------------------------------------------------------------
    # filter out all vertices that contribute < 1 %
    keep = []
    delete = []
    for v in g.vertices():
        rxnrate = (sum(res[e] for e in v.in_edges()) +
                    sum(res[e] for e in v.out_edges())) / 2.
        if (rxnrate < max_flow/50.):
            g.vp.myfilter[v] = False
            delete.append(g.vp.name[v])
        else:
            g.vp.myfilter[v] = True
            keep.append(g.vp.name[v])
    print('Keep: ', sorted(keep))
    print('Delete: ', sorted(delete))
    g.set_vertex_filter(g.vp.myfilter) # keep if True
    #-------------------------------------------------------------------------
    # vertex colors:
    myvcolor = g.new_vertex_property('string')
    for v in g.vertices():
        myvcolor[v] = 'yellow'
    myvcolor[src] = 'red'
    myvcolor[tgt] = 'green'
    #-------------------------------------------------------------------------

    gt.graphviz_draw(g, # http://www.graphviz.org/doc/info/attrs.html
                     #pos=pos,
                     gprops={#'concentrate':'true',
                         #'overlap':'prism0',
                         'overlap':'voronoi', 
                         #'rankdir':'LR',
                         'splines':'curved',
                         #'splines':'ortho',
                         #'splines':'polyline',
                         #'splines':'spline',
                         'size':'40,20'},
                     vprops={'shape':'oval',
                             'fontsize':20,#gt.prop_to_size(x, mi=1, ma=30),
                             #'color':'Black', # circle around vertex
                             'style':'filled',
                             #'label':g.vp.name,
                             'label':'%s (%g)' % (g.vp.name,g.vp.oic)},
                     eprops={'label':g.ep.eqntag,
                             'dir':'forward', # forward arrows
                             'arrowhead':'normal',
                             'arrowsize':1,
                             'fontsize':20,
                             #'style':'dashed',
                             #'style':'dotted',
                             #'color':'black;0.25:red;0.5:yellow;0.25',
                             'color':'blue',
                             'penwidth':edgewidth},
                     overlap='compress',
                     vcolor=myvcolor,
                     #vcolor='yellow', # vertex fill color
                     output='mecca_edmonds-karp.pdf')

def read_rxn_rates(rxnfilename):
    from netCDF4 import Dataset, num2date
    ncid = Dataset(rxnfilename)
    time = ncid.variables['time']
    mytime = len(time)-36 # last day at noon if delta_t = 20 min
    print('Selecting time', num2date(time[mytime],time.units), 'from', rxnfilename)
    if DEBUG: print(ncid.variables['RRJ41000'][:,0,0,0]) # show rxn rates
    rxnrates = {}
    for rxn in ncid.variables:
        if (rxn[0:2]=='RR'):
            mydata = ncid.variables[rxn][mytime,0,0,0]
            rxnrates[rxn[2:]] = mydata
    ncid.close()
    return rxnrates # dict with eqntags -> rxnrates

##############################################################################

def create_interactive_window(g):
    # middle mouse button: move
    # mouse wheel: zoom
    # shift + mouse wheel: zoom (including vertex + edge sizes)
    # control + mouse wheel: rotate
    # 'a': autozoom
    # 'r': resize and center
    # 's': spring-block layout for all non-selected vertices
    # 'z': zoom to selected vertices
    # left mouse button: select vertex
    # right mouse button: unselect all
    # shift + left mouse button drag: select several vertices
    # scroll middle mouse button: zoom in/out
    
    # https://graph-tool.skewed.de/static/doc/draw.html#graph_tool.draw.interactive_window
    # for options, see "List of vertex properties" at:
    # https://graph-tool.skewed.de/static/doc/draw.html
    gt.interactive_window(g,
                          geometry=(1000, 800), # initial window size
                          edge_text=g.ep.eqntag,
                          vertex_text=g.vp.name,
                          #layout_callback=layoutchanged,
                          key_press_callback=keypressed,
                          vertex_fill_color='yellow',
                          #vertex_aspect=2, # very slow!!!
                          #display_props=g.vp.name,
                          display_props=[g.vp.name,g.vp.atoms],
                          display_props_size=20
    )

def keypressed(gtk, g, keyval, picked, pos, vprops, eprops):
    print('keypressed')
    #print 'gtk='    ,gtk    # gtk_draw.GraphWidget
    #print 'g='      ,g      # graph being drawn
    print('keyval=' ,keyval) # key id (a=97, z=122)
    if picked:
        print('picked=' ,g.vp.name[picked]) # vertex or boolean vertex property map for selected vertices
    #print 'pos='    ,pos    # vertex positions
    print('vprops=')         # vertex property dictionary
    print(vprops['text'][picked])
    for vprop in vprops:
        print('%s = %s' % (vprop, vprops[vprop]))
    print('eprops=')         # edge property dictionary
    for eprop in eprops:
        print('%s = %s' % (eprop, eprops[eprop]))
    # for eqntag in eprops['text']:
    #     print eqntag

def layoutchanged(gtk, g, picked, pos, vprops, eprops):
    print('layoutchanged')
    print('gtk='    ,gtk)    # gtk_draw.GraphWidget
    print('g='      ,g)      # graph being drawn
    print('picked=' ,picked) # vertex or boolean vertex property map for selected vertices
    print('pos='    ,pos)    # 
    print('vprops=' ,vprops) # vertex property dictionary 
    print('eprops=' ,eprops) # edge property dictionary

##############################################################################

def list_vertices(g):
    for v in g.vertices():
        print(g.vp.name[v], end=' ')
    print('--> %d species' % (g.num_vertices()))

##############################################################################

def list_edges(g):
    for e in g.edges():
        print(g.ep.eqntag[e], g.ep.reaction[e])
        print('%s (%d) -> %s (%d)' % (
            g.vp.name[e.source()], elem_count[g.vp.name[e.source()]]['C'],
            g.vp.name[e.target()], elem_count[g.vp.name[e.target()]]['C']))
    print('--> %d reactions' % (g.num_edges()))

##############################################################################

def elemental_composition(g):
    elem_count = {} # define dictionary for elemental composition
    for v in g.vertices(): # loop over all species
        # create temporary dictionary with element count set to zero:
        tmp_dict = {'Ignore':0,'Pls':0,'Min':0,'O':0,'H':0,'N':0,'C':0,'F':0,
                    'Cl':0,'Br':0,'I':0,'S':0,'Hg':0,'Fe':0}
        # loop over all elements of current species:
        for element in g.vp.atoms[v].split('+'):
            # search for count and element symbol:
            search_result = re.search('([0-9]*)([A-Za-z]+)', element)
            if (search_result.group(1)==''):
                count = 1
            else:
                count = int(search_result.group(1))
            # update current element in temporary dictionary:
            tmp_dict[search_result.group(2)] = count
        # add elemental composition of current species to dictionary:
        elem_count[g.vp.name[v]] = tmp_dict
        if (DEBUG>1):
            print('%-15s %-15s' % (g.vp.name[v],g.vp.atoms[v]), end=' ')
            print(tmp_dict)
    # example usage: HCHO has elem_count['HCHO']['H'] H atoms
    return elem_count

##############################################################################

def set_vcolor_oic():

    g.vp.oic = g.new_vertex_property('double') # oic
    oicfilename = 'OIC.dat'
    OICFILE = open(oicfilename)
    for line in iter(OICFILE):
        search_result = re.search( r'^ *([-+.Ee0-9]+) *([A-Za-z0-9_]+) *$', line)
        oic     = float(search_result.group(1))
        species = search_result.group(2)
        v = n2v(g,species) # vertex number of the species in the graph
        if (v<0):
            #print 'not in current mechanism'
            continue
        else:
            g.vp.oic[v] = oic
            # cd ~/messy/mecca/mechanism_reduction/skeleton_simple_organic
            # nvar = 462/663, nreact = 1444/2091, eps = 0.0007
            # nvar = 429/663, nreact = 1320/2091, eps = 0.00098
            # nvar = 411/663, nreact = 1262/2091, eps = 0.001372
            g.vp.fillcolor[v] = '#FF5555' # red = default
            if (oic>0.0007):   g.vp.fillcolor[v] = '#F1A65B'
            if (oic>0.00098):  g.vp.fillcolor[v] = 'yellow'
            if (oic>0.001372): g.vp.fillcolor[v] = 'green'
    OICFILE.close()

# ##############################################################################

def set_vcolor_Fe():
    for v in g.vertices():
        # yelow for Fe(II):
        if (g.vp.name[v]=='Fepp_a01'):  g.vp.fillcolor[v] = 'yellow'
        if (g.vp.name[v]=='FeClp_a01'): g.vp.fillcolor[v] = 'yellow'
        if (g.vp.name[v]=='FeOHp_a01'): g.vp.fillcolor[v] = 'yellow'
        # red for Fe(IV):
        if (g.vp.name[v]=='FeOpp_a01'): g.vp.fillcolor[v] = 'red'

##############################################################################

def set_labels_oic():
    # define vertex label:
    g.vp.label = g.new_vertex_property('string')
    for v in g.vertices():
        # g.vp.label[v] = '%s\n%dC, %8.2E' % (
        #     g.vp.name[v], elem_count[g.vp.name[v]]['C'], g.vp.oic[v])
        g.vp.label[v] = '%s\n%8.2E' % (g.vp.name[v], g.vp.oic[v])

##############################################################################

def set_filter_organic(nC):
    # only organic species (nC = number of carbon atoms):
    g.vp.myfilter.a = [elem_count[g.vp.name[v]]['C']==nC for v in g.vertices()]

def set_filter_Fe():
    g.vp.myfilter.a = [elem_count[g.vp.name[v]]['Fe']>0 for v in g.vertices()]

def set_filter_src_tgt():

    # remove edges where number of C atoms increases:
    g.ep.myfilter.a = [elem_count[g.vp.name[e.source()]]['C'] >=
                       elem_count[g.vp.name[e.target()]]['C'] for e in g.edges()]
    g.set_edge_filter(g.ep.myfilter) # keep if True

    #-------------------------------------------------------------------------
    # define sources and targets:
    #src = n2v(g,'MACR')
    #src = n2v(g,'APINENE')
    #src = n2v(g,'CH3COCH3')
    #src = n2v(g,'HCHO')
    #src = n2v(g,'C5H8')
    src = n2v(g,'CH4')  
    #tgt = n2v(g,'HCOOH')
    tgt = n2v(g,'CO2')
    #tgt = n2v(g,'HCHO')
    #-------------------------------------------------------------------------
    # remove unreachable species:
    s_dist = gt.shortest_distance(g, source=src)
    g.set_reversed(True)
    t_dist = gt.shortest_distance(g, source=tgt)
    g.set_reversed(False)
    g.vp.myfilter.a = [s_dist[v]<N_v and t_dist[v]<N_v for v in g.vertices()]

    for v in g.vertices():
        g.vp.myfilter.a[v] = g.vp.myfilter.a[v] and elem_count[g.vp.name[v]]['C']==5
        #g.vp.myfilter.a[v] = g.vp.myfilter.a[v] and elem_count[g.vp.name[v]]['C']>0
        if (DEBUG>1):
            print(g.vp.name[v], s_dist[v], t_dist[v], g.vp.myfilter.a[v])
    
    # show shortest path from source to target:
    print('from %s to %s' % (g.vp.name[src], g.vp.name[tgt]))
    vlist, elist = gt.shortest_path(g, src, tgt)
    for v in vlist:
        print(g.vp.name[v], ' ', end=' ')
    print()
    for e in elist:
        print(g.ep.eqntag[e], g.ep.reaction[e])
    return src, tgt

##############################################################################

def misc(src, tgt):

    # create graph view gC with only organic species (C>=1):
    gC = gt.GraphView(g, vfilt=lambda v: elem_count[g.vp.name[v]]['C']>0)
    
    # for v in gC.vertices():
    #     print '%-15s %d %d %d %d' % (
    #         gC.vp.name[v], int(gC.vp.myfilter[v]),
    #         elem_count[gC.vp.name[v]]['C'], s_dist[v], t_dist[v])

    # remove edges where number of C atoms increases:
    #print HLINE ; list_edges(gC)
    gC = gt.GraphView(gC, efilt=lambda e:
                      elem_count[gC.vp.name[e.source()]]['C'] >=
                      elem_count[gC.vp.name[e.target()]]['C'])
    #print HLINE ; list_edges(gC)

    s_dist = gt.shortest_distance(gC, source=src)
    gC.set_reversed(True)
    t_dist = gt.shortest_distance(gC, source=tgt)
    gC.set_reversed(False)
    fraction = g.new_vertex_property('float')
    list_vertices(gC)
    # remove unreachable species from graph:
    gC = gt.GraphView(gC, vfilt=lambda v: s_dist[v]<N_v and t_dist[v]<N_v)
    list_vertices(gC)
    print('name, #C, %s->, ->%s, fraction' % (gC.vp.name[src], gC.vp.name[tgt]))
    for v in gC.vertices():
        fraction[v] = float(s_dist[v])/(s_dist[v]+t_dist[v])
        print('%-15s %d %d %d %s' % (
            gC.vp.name[v], elem_count[gC.vp.name[v]]['C'],
            s_dist[v], t_dist[v], fraction[v]))

    #-------------------------------------------------------------------------

    # centrality:
    # x = gt.katz(g)
    # print x.a

    #-------------------------------------------------------------------------

    # correlations:
    # https://graph-tool.skewed.de/static/doc/correlations.html
    # h = gt.corr_hist(g, "out", "out")
    # plt.clf()
    # plt.xlabel("Source out-degree")
    # plt.ylabel("Target out-degree")
    # plt.imshow(h[0].T, interpolation="nearest", origin="lower")
    # plt.colorbar()
    # plt.savefig("corr.pdf")

##############################################################################

if __name__ == '__main__':

    # load graph produced by define_graph.py:
    g = gt.load_graph("mecca_graph.xml.gz")
    N_v = g.num_vertices() # number of species

    # define dictionary with elemental composition:
    elem_count = elemental_composition(g)

    g.vp.myfilter = g.new_vertex_property('bool') # vertex filter (internal property map)
    g.ep.myfilter = g.new_edge_property('bool')   # edge filter   (internal property map)

    ##########################################################################
    # choose a filter:
    #set_filter_organic(10)
    #set_filter_Fe()
    src, tgt = set_filter_src_tgt()

    g.set_vertex_filter(g.vp.myfilter) # keep if True
    # g.clear_filters()

    ##########################################################################
    # choose fill colors for the species:
    #set_vcolor_Fe()
    set_vcolor_oic()
    ##########################################################################
    # choose labels for vertices:
    set_labels_oic()

    ##########################################################################
    
    misc(src, tgt)
    
    #-------------------------------------------------------------------------

    #create_interactive_window(g)
    #create_graphviz_draw(g)
    create_max_flow(gC)
    #list_edges(g)

    #-------------------------------------------------------------------------

    #sys.exit('END') #qqq
