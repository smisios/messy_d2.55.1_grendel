#!/usr/bin/env python3
# -*- coding: utf-8 -*- Time-stamp: <2019-06-03 13:44:40 sander>

##############################################################################

import re
import graph_tool.all as gt
from utils_graph import n2v

HLINE = '-' * 78
DEBUG = 0
BIPART = False # create also bipartite graph for Gupta analysis?

paramfilename = 'messy_mecca_kpp_parameters.f90'
spcfilename   = 'mecca.spc'
eqnfilename   = 'mecca.eqn'
graphfilename = 'mecca_graph.xml.gz'

##############################################################################

def init_graph():
    g = gt.Graph()
    # make internal property maps:
    g.vp.name      = g.new_vertex_property('string')
    g.vp.atoms     = g.new_vertex_property('string') 
    g.vp.fillcolor = g.new_vertex_property('string') 
    g.ep.reaction  = g.new_edge_property('string')
    g.ep.eqntag    = g.new_edge_property('string')
    g.ep.prodstoic = g.new_edge_property('double') # stoic factor of product
    if DEBUG: g.list_properties()
    return g

##############################################################################

def generate_spclist(paramfilename):
    # create a list of species that actually occur in the mechanism,
    # i.e. for which ind_XYZ in the paramfile is not zero:
    spclist = []
    PARAMFILE = open(paramfilename)
    for line in iter(PARAMFILE):
        # is current line something like 'INTEGER ... ind_XYZ = nnn' with nnn != 0 ?
        search_result = re.search( r'INTEGER.*ind_([A-Za-z0-9_]+) = [^0]', line)
        if (search_result): # skip lines that do not define a species
            species = search_result.group(1)
            if (species[0:2]!='RR'): # ignore RR* rxn rate pseudo-species:
                spclist.append(species)
        #if DEBUG: print '|%s|, |%s|' % (search_result.group(1), search_result.group())
    PARAMFILE.close()
    if DEBUG: print(spclist)
    return spclist

##############################################################################

def define_vertices(spcfilename):
    # get elemental composition from *.spc file:
    SPCFILE = open(spcfilename)
    for line in iter(SPCFILE):
        search_result = re.search( r'^ *([A-z0-9_]+) *=([A-z0-9+ ]+);', line)
        if (not search_result): # skip lines that do not define a species
            if (DEBUG>2): print('NO:  |%s|' % (line))
            continue
        spc_name = search_result.group(1)
        spc_composition = search_result.group(2).replace(' ','')
        # split composition into atoms and remove whitespace:
        atoms = [x.strip() for x in spc_composition.split('+')]
        if (DEBUG>2):
            print('%s\nALL: |%s|' % (HLINE,line))
            print('ALL: |%s|' % (search_result.group()))
            print('1:   |%s|' % (spc_name))
            print('2:   |%s|' % (spc_composition))
            for atom in atoms:
                print('|%s|' % (atom))
        #---------------------------------------------------------------------
        # add vertex to graph describing this species
        # but ignore RR* rxn rate pseudo-species:
        if ((spc_name in spclist) and (spc_name[0:2]!='RR')):
            #if DEBUG: print 'adding:      %s' % (spc_name)
            newvertex = g.add_vertex()
            g.vp.name[newvertex]  = spc_name
            g.vp.atoms[newvertex] = spc_composition
            g.vp.fillcolor[newvertex]  = 'white'
            if (BIPART): bipart_insertSpeciesName(spc_name)
        #---------------------------------------------------------------------
    SPCFILE.close()

##############################################################################

def define_edges(eqnfilename):
    # get reactions from *.eqn file:
    rxnlist = []
    EQNFILE = open(eqnfilename)
    for line in iter(EQNFILE):
        line = re.sub('{[^}]*}', '', line) # delete all comments {...} from line
        search_result = re.search( r'<([A-Za-z_0-9]+)> *(.*)=(.*):.*;', line)
        if (not search_result): # skip lines that do not define a reaction
            if (DEBUG>2): print('NO:  |%s|' % (line))
            continue
        if (DEBUG>2): print(HLINE)
        eqntag        = search_result.group(1)
        rxnlist.append(eqntag)
        reactants_str = search_result.group(2).replace(' ','')
        products_str  = search_result.group(3).replace(' ','')
        # remove 1st species from product list if it is a reaction rate dummy RR*+:
        if (DEBUG>2): print('BEFORE: %s' % (products_str))
        products_str = re.sub('^RR[^+]*\+', '', products_str)
        if (DEBUG>2): print('AFTER:  %s' % (products_str))
        fullrxn = reactants_str + '->' + products_str
        if (DEBUG>2): print('FULL: %9s %s' % ('<'+eqntag+'>', fullrxn))
        reactants = reactants_str.split('+')
        products  = products_str.split('+')
        if DEBUG: print('%10s %s -> %s' % ('<'+eqntag+'>', reactants_str, products_str))
        for reactant in reactants:
            if (reactant=='hv'): continue # ignore pseudo-reactant hv
            for product in products:
                # separate stoichiometric factor:
                search_result = re.search( '^([0-9.]*)(.*)', product)
                stoic = search_result.group(1)
                stoic = 1 if (stoic=='') else float(stoic)
                prod = search_result.group(2)
                if DEBUG: print(eqntag, prod, stoic)
                # add edge to graph describing this reaction:
                newedge = g.add_edge(g.vertex(n2v(g,reactant)), g.vertex(n2v(g,prod)))
                g.ep.reaction[newedge]  = fullrxn
                g.ep.prodstoic[newedge] = stoic
                g.ep.eqntag[newedge]    = eqntag
                if DEBUG: print('edge: %9s %s -> %s' % ('<'+eqntag+'>', reactant, product))
        if (BIPART): bipart_insertEdges(reactants, products, eqntag)
    EQNFILE.close()

##############################################################################

def load_and_plot(graphfilename):
    g = gt.load_graph(graphfilename)
    gt.graphviz_draw(g,
                     gprops={'concentrate':'true', 'rankdir':'LR',
                             'splines':'spline',
                             'size':'40,20'},
                     vprops={'shape':'oval',
                             'fontsize':10,
                             'style':'filled',
                             'label':g.vp.name},
                     eprops={'label':g.ep.eqntag,
                             'dir':'forward', # forward arrows
                             'arrowhead':'normal',
                             'arrowsize':1,
                             'fontsize':10,
                             'color':'blue',
                             'penwidth':1},
                     overlap='compress',
                     vcolor='yellow', # vertex fill color
                     output='define_graph.pdf')

##############################################################################

# some functions for bipartite graphs for Gupta analysis:

def bipart_read_rxn_rates(rxnfilename):
    from netCDF4 import Dataset, num2date
    ncid = Dataset(rxnfilename)
    time = ncid.variables['time']
    mytime = len(time)-36 # last day at noon if delta_t = 20 min
    print('Selecting time', num2date(time[mytime],time.units), 'from', rxnfilename)
    rxnrates = {}
    for rxn in ncid.variables:
        if (rxn[0:2]=='RR'):
            mydata = ncid.variables[rxn][mytime,0,0,0]
            rxnrates[rxn[2:]] = mydata
    ncid.close()
    return rxnrates # dict with eqntags -> rxnrates

def bipart_insertEdges(reactants, products, eqntag):
    rxnnum = -1-rxnlist.index(eqntag)
    for reactant in reactants:
        if (not reactant=='hv'):
            print('G.insertEdges(%4d, %4d, 1); // %s: %s ->' % (
                spclist.index(reactant)+1, rxnnum, eqntag, reactant), file=HPPFILE)
    for product in products:
        # remove stoichiometric factor:
        prod  = re.sub('^[0-9.]+', '', product)
        print('G.insertEdges(%4d, %4d, 2); // %s: -> %s' % (
            spclist.index(prod)+1, rxnnum, eqntag, prod), file=HPPFILE)
    # categorize reactions as slow and fast:
    if (rxnrates[eqntag]>1e-14):
        speed = 'FAST'
    else:
        speed = 'SLOW'
    print('G.insertIrreversibleRxns(%4d, %s); // %s' % (
        rxnnum, speed, eqntag), file=HPPFILE)        

def bipart_insertSpeciesName(spc_name):
    # if not hasattr(bipart_insertSpeciesName, 'counter'):
    #     bipart_insertSpeciesName.counter = 1
    # else:
    #     bipart_insertSpeciesName.counter += 1
    # spc_num = bipart_insertSpeciesName.counter
    # print >> HPPFILE, 'G.insertSpeciesName(%4d, "%s")' % (
    #     spc_num, spc_name)
    print('G.insertSpeciesName(%4d, "%s");' % (
        spclist.index(spc_name)+1, spc_name), file=HPPFILE)
                
##############################################################################

if __name__ == '__main__':

    # create a graph that will contain the reaction mechanism:
    g = init_graph()
    
    if (BIPART): HPPFILE = open('mecca.hpp','w')
    if (BIPART): rxnrates = bipart_read_rxn_rates('caaba_mecca_rr.nc')

    # create a list of species:
    spclist = generate_spclist(paramfilename)

    # define chemical species and their composition as vertices of the graph:
    define_vertices(spcfilename)

    # define reactions and their stoichiometry as edges of the graph:
    define_edges(eqnfilename)

    if (BIPART): HPPFILE.close()

    # save graph to file:
    g.save(graphfilename)

    # check that saved graph is okay:
    load_and_plot(graphfilename)

##############################################################################
