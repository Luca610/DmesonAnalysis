from ROOT import *
import sys
import argparse
import yaml
import os

'''
Simple script to project THnSparse filled with resonance task and save it as output in a root file.
Author: Luca Aglietta
'''

# main
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Arguments to pass')
    parser.add_argument('cfgFileName', metavar='text', default='cfgFileName.yml',
                    help='config file name with root input files')
    parser.add_argument('--outputDirectory','-d', default='.',
                    help='Directory for output plots. If ".", save in the same directory of the input file.')
    parser.add_argument('--plot', '-p', action='store_true',
                    help='Boolean to produce very simple plot')
    args = parser.parse_args()

    # acces config
    with open(args.cfgFileName, 'r') as ymlCfgFile:
        inputCfg = yaml.load(ymlCfgFile, yaml.FullLoader)
    conf = inputCfg['proj']
    inputFile = conf['filename']
    dirName = conf['dirname']
    histName = conf['histname']
    massAxis = conf['massAxis']
    massLims = conf['masslims']
    ptAxis = conf['ptAxis']
    ptLims = conf['ptLims']

    # open root file and get Sparse
    file = TFile.Open(inputFile)
    taskDir = file.Get(dirName)
    sparse = taskDir.Get(histName)

    # Limit mass range and projection
    minBin = sparse.GetAxis(massAxis).FindBin(massLims[0] - 1e-6*massLims[0])
    maxBin = sparse.GetAxis(massAxis).FindBin(massLims[1] + 1e-6*massLims[1])
    sparse.GetAxis(massAxis).SetRange(minBin,maxBin)
    ptBinMin = sparse.GetAxis(ptAxis).FindBin(ptLims[0] + 1e-6*ptLims[0])
    ptBinMax = sparse.GetAxis(ptAxis).FindBin(ptLims[-1] - 1e-6*ptLims[-1])
    sparse.GetAxis(ptAxis).SetRange(ptBinMin, ptBinMax)
    massProj = sparse.Projection(massAxis)
    # pt differential
    ptDiffProj = []
    for iPt, pt in enumerate(ptLims[:-1]):
        binMin = sparse.GetAxis(ptAxis).FindBin(pt + 1e-6*pt)
        binMax = sparse.GetAxis(ptAxis).FindBin(ptLims[iPt + 1] - 1e-6*ptLims[iPt + 1])
        sparse.GetAxis(ptAxis).SetRange(binMin, binMax)
        ptDiffProj.append(sparse.Projection(massAxis))
        ptDiffProj[iPt].SetName(f'hSparse_proj_{ptLims[iPt]}-{ptLims[iPt+1]}')

    # Optional produce a QA plot
    if args.plot:
        c1 = TCanvas("c1", "Projection with Range", 800, 600)
        massProj.Draw()
        c1.SaveAs("projectionProva.png")  # Save the plot as a PNG file

    # Save projection in an output file
    inFileName = os.path.basename(inputFile)
    outFileName = 'Proj_' + inFileName 
    outPath = os.path.join(args.outputDirectory, outFileName)
    outputFile = TFile(outPath, "RECREATE")
    massProj.Write()
    for h in ptDiffProj:
        h.Write()
    outputFile.Close()

    file.Close()

