#!/usr/bin/env python
import sys, os
import argparse
import networkx as nx
import numpy as np
import scipy.io as scio
from useful import system, which, round_to_n

def grangerNetwork(matfile, fn_out, pvalname='GCpvalB', weightname='GCdevB', clustername='clustersB', bciunits = 'brainunits'):
    #Load mat file
    a = scio.loadmat(matfile)
    weights = a[weightname]
    pvals = a[pvalname].astype(float)
    names = [str(name[0]) for name in a['unitnames'][0]]
    #Add info about whether node was used in the BCI or not
    bci = [str(name[0]) for name in a[bciunits][0]]
    usedinbci = [int(u in bci) for u in names]
    clusters = a[clustername][0]
    nodeclusters = {}
    for i in range(len(clusters)):
        for j in range(np.size(clusters[i][0])):
            node = clusters[i][0][j]
            nodeclusters[node] = i

    nU = np.size(weights,1)
    #Set pvals of zero to something very small
    pvals[pvals==0] = 1e-20
    #Set diagonal elements to 0
    for i in range(nU):
        weights[i,i]=0
    #weights[weights==0] = 0
    weights[pvals>0.05] = 0
    maxw = np.max(weights)
    #Create network
    g = nx.MultiDiGraph()
    #Color nodes according to sum of input and output change in deviance
    directionality = np.zeros((nU,1))
    sumcols = weights.sum(0)
    sumrows = weights.sum(1)
    for i in range(nU):
        directionality[i] = sumcols[i]-sumrows[i]
    [g.add_node(i, cluster = float(nodeclusters[i+1]), unitname = float(names[i]), dirs=float(directionality[i]), bci=usedinbci[i]) for i in range(nU)]
    for i in range(nU):
        for j in range(nU):
            if abs(weights[i,j])>5:
                g.add_edge(i,j, weight=round_to_n(float(weights[i,j]),2), pval=round_to_n(float(pvals[i,j]),2))
    #Write to graphml file
    nx.write_graphml(g, fn_out, prettyprint = False);

def runGranger():
    matfile = './GLMGrangerB.mat'; fn_out = './GLMGrangerB.xml';
    grangerNetwork(matfile, fn_out, pvalname='GCpvalB', weightname='GCdevB', clustername='clustersB', bciunits = 'brainunits')
    matfile = './GLMGrangerBP.mat'; fn_out = './GLMGrangerBP.xml';
    grangerNetwork(matfile, fn_out, pvalname='GCpvalBP', weightname='GCdevBP', clustername='clustersBP', bciunits = 'brainunits')
    matfile = './GLMGrangerM.mat'; fn_out = './GLMGrangerM.xml';
    grangerNetwork(matfile, fn_out, pvalname='GCpvalM', weightname='GCdevM', clustername='clustersM', bciunits = 'manualunits')
    matfile = './GLMGrangerMP.mat'; fn_out = './GLMGrangerMP.xml';
    grangerNetwork(matfile, fn_out, pvalname='GCpvalMP', weightname='GCdevMP', clustername='clustersMP', bciunits = 'manualunits')

def runCytoscape(fn_graphml, fn_image, fn_style = '/home/lansdell/projects/bci/matlab/eval/GrangerStyle.xml'):
    #Check if cytoscape.sh is present, if not then exit
    if not which('cytoscape.sh'):
        return None 

    fn_out = './tmp.cy'
    fn_attr = fn_image + 'attr'
    fn_circ = fn_image + 'circ'

    #Generate script 
    script = """#Import network
network import file indexColumnTargetInteraction=1 indexColumnSourceInteraction=2 file="%s"
#Import and set style
vizmap load file file="%s"
vizmap apply styles=DirectedGranger
#Set layout to attributes
layout attributes-layout NodeAttribute=cluster maxwidth=400
#Set view to fit display
view fit content
#Save 
view export OutputFile="%s" options=PDF
#Set layout to cirlce
layout attribute-circle NodeAttribute="unit name"
#Set view to fit display
view fit content
#Save 
view export OutputFile="%s" options=PDF
#Quit
command quit""" % (fn_graphml, fn_style, fn_attr, fn_circ)

    print script

    with open(fn_out, 'w') as f_out:
        f_out.write(script)

    #Delete old images before making new ones
    cmd = 'rm ' + fn_attr + '.pdf'
    system(cmd)
    cmd = 'rm ' + fn_circ + '.pdf'
    system(cmd)
    #Run cytoscape
    cmd = 'cytoscape.sh -S ' + fn_out
    system(cmd)

def runDiffCytoscape(fn_graphml, fn_image, fn_style = '/home/lansdell/projects/bci/matlab/eval/GrangerStyle.xml'):
    #Check if cytoscape.sh is present, if not then exit
    if not which('cytoscape.sh'):
        return None 

    fn_out = './tmp.cy'
    fn_attr = fn_image + 'attr'
    fn_circ = fn_image + 'circ'

    #Generate script 
    script = """#Import network
network import file indexColumnTargetInteraction=1 indexColumnSourceInteraction=2 file="%s"
#Import and set style
vizmap load file file="%s"
vizmap apply styles=DirectedGrangerDiff
#Set layout to attributes
layout attributes-layout NodeAttribute=cluster maxwidth=400
#Set view to fit display
view fit content
#Save 
view export OutputFile="%s" options=PDF
#Set layout to cirlce
layout attribute-circle NodeAttribute="unitname"
#Set view to fit display
view fit content
#Save 
view export OutputFile="%s" options=PDF
#Quit
#command quit""" % (fn_graphml, fn_style, fn_attr, fn_circ)

    print script

    with open(fn_out, 'w') as f_out:
        f_out.write(script)

    #Delete old images before making new ones
    cmd = 'rm ' + fn_attr + '.pdf'
    system(cmd)
    cmd = 'rm ' + fn_circ + '.pdf'
    system(cmd)
    #Run cytoscape
    cmd = 'cytoscape.sh -S ' + fn_out
    system(cmd)

def main(argv):
    usage = """networkx_granger.py <matfile> <pdffile> --weightname [name] 
                --pvalname [name] --clustername [name] --bciunits [name] --style [xmlfile]

    Create networkx object from clustered Granger causality code

    Requires a .mat file containing:

    - A square NxN matrix which contains the weights for each edge. The entry in row i
     and column j indicates unit i's "G-causal effect" on unit j.
    - A NxN matrix with the p-values testing for Granger causality
    - A cell array containing lists of integers specifying units within the same cluster
    - A cell array containing strings specifying names of units were used in BCI control

    Variables names are specified by the command line arguments.

    Writes to graphml file for importing into cytoscape.
   
    If cytoscape.sh is in path will generate and run a cytoscape script to plot the network.

    Ben Lansdell: ben.lansdell@gmail.com
    2014"""
    parser = argparse.ArgumentParser(description=usage)
    parser.add_argument('matfile', type=str,
        help='matfile containing Granger causality data')
    parser.add_argument('fn_out', type=str,
        help='filename for output pdf file')
    parser.add_argument('--weightname', type=str, default='GCdevB',
        help='name of weight matrix within mat file to plot')
    parser.add_argument('--pvalname', type=str, default='GCpvalB',
        help='name of pval matrix within mat file to plot')
    parser.add_argument('--clustername', type=str, default='clustersB',
        help='name of cluster cell array within mat file to plot')
    parser.add_argument('--bciunits', type=str, default='brainunits',
        help='name of cell array listing BCI units within mat file to plot')
    parser.add_argument('--style', type=str, default='/home/lansdell/projects/bci/matlab/eval/GrangerStyle.xml',
        help='path to style file for cytoscape image')
    parser.add_argument('--diff', action='store_true',
        help='Whether matrices indicate difference in Granger causality, or absolute change')
    args = parser.parse_args()
    xml_out_rel = './tmp.xml'
    xml_out_abs = os.path.abspath(xml_out_rel)
    fn_style = os.path.abspath(args.style)
    grangerNetwork(args.matfile, xml_out_abs, args.pvalname, args.weightname, args.clustername, args.bciunits)
    if (args.diff == False):
        runCytoscape(xml_out_abs, args.fn_out, fn_style)
    else:
        runDiffCytoscape(xml_out_abs, args.fn_out, fn_style)

if __name__ == "__main__":
    sys.exit(main(sys.argv))