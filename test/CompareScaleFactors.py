#!/usr/bin/env python
import os
import sys
import ROOT
import math
import optparse
from array import *

# From https://github.com/cms-nanoAOD/nanoAOD-tools/blob/master/python/postprocessing/modules/btv/btagSFProducer.py
supported_btagSF = {
    'csvv2' : {
        '2016' : {
            'inputFileName' : "btagSF_CSVv2_ichep2016.csv",
            'measurement_types' : {
                0 : "comb",  # b
                1 : "comb",  # c
                2 : "incl"   # light
            },
            'supported_wp' : [ "L", "M", "T", "shape_corr" ]
        },
        '2017' : {
            'inputFileName' : "CSVv2_94XSF_V2_B_F.csv",
            'measurement_types' : {
                0 : "comb",  # b
                1 : "comb",  # c
                2 : "incl"   # light
            },
            'supported_wp' : [ "L", "M", "T", "shape_corr"]
        }
    },
    'deepcsv' : {
        'Legacy2016' : {
            'inputFileName' : "DeepCSV_2016LegacySF_V1.csv",
            'measurement_types' : {
                0 : "comb",  # b
                1 : "comb",  # c
                2 : "incl"   # light
            },
            'supported_wp' : [ "L", "M", "T", "shape_corr"]
        },
        '2017' : {
            'inputFileName' : "DeepCSV_94XSF_V4_B_F.csv",
            'measurement_types' : {
                0 : "comb",  # b
                1 : "comb",  # c
                2 : "incl"   # light
            },
            'supported_wp' : [ "L", "M", "T", "shape_corr"]
        },
        '2018' : {
            'inputFileName' : "DeepCSV_102XSF_V1.csv",
            'measurement_types' : {
                0 : "comb",  # b
                1 : "comb",  # c
                2 : "incl"   # light
            },
            'supported_wp' : [ "L", "M", "T", "shape_corr"] 
        }    
    },
    'deepjet' : {
        'Legacy2016' : {
            'inputFileName' : "DeepJet_2016LegacySF_V1.csv",
            'measurement_types' : {
                0 : "comb",  # b
                1 : "comb",  # c
                2 : "incl"   # light
            },
            'supported_wp' : [ "L", "M", "T", "shape_corr"]
        },
        '2017' : {
            'inputFileName' : "DeepFlavour_94XSF_V3_B_F.csv",
            'measurement_types' : {
                0 : "comb",  # b
                1 : "comb",  # c
                2 : "incl"   # light
            },
            'supported_wp' : [ "L", "M", "T", "shape_corr"]
        },
        '2018' : {
            'inputFileName' : "DeepJet_102XSF_V1.csv",
            'measurement_types' : {
                0 : "comb",  # b
                1 : "comb",  # c
                2 : "incl"   # light
            },
            'supported_wp' : [ "L", "M", "T", "shape_corr"]
        }    
    },
    'cmva' : {
        '2016' : {
            'inputFileName' : "btagSF_cMVAv2_ichep2016.csv",
            'measurement_types' : {
                0 : "ttbar", # b
                1 : "ttbar", # c
                2 : "incl"   # light
            },
            'supported_wp' : [ "L", "M", "T", "shape_corr" ]
        }
    }
}

def plotScaleFactors(plottitle, scaleFactorHistograms, plotformat):

    yoff = 0.85 if ('_T_' in plottitle) else 0.35
    leg = ROOT.TLegend(0.25, yoff, 0.50, yoff-0.2)
    leg.SetFillColor(ROOT.kWhite) 
    leg.SetBorderSize(0)
    leg.SetTextColor(1) 
    leg.SetTextSize(0.04)
    leg.SetTextFont(62) 
    
    histoiter = 1 
    minY, maxY = 999., -1.

    for histo in scaleFactorHistograms:

        histo.SetXTitle("p_{T} [GeV]")
        histo.SetYTitle("SF_{b}")
        histo.SetLineColor(histoiter)
        histo.SetLineWidth(3)
        histo.SetFillColor(histoiter)
        histo.SetFillStyle(0)
        
        histo.SetLabelSize(0.05, "XYZ")
        histo.SetTitleSize(0.06, "XYZ") 
        histo.SetLabelFont(42, "XYZ") 
        histo.SetTitleFont(42, "XYZ")
        histo.SetTitleOffset(0.95,"X")
        histo.SetTitleOffset(0.8,"Y")
        histo.SetTickLength(0.06,"X")
        histo.SetNdivisions(509, "XYZ")
        histo.GetXaxis().SetMoreLogLabels()
        histo.GetXaxis().SetNoExponent()
        
        leg.AddEntry(histo, histo.GetName(), "l")

        if histo.GetMaximum()>maxY: 
            maxY = histo.GetMaximum()

        if histo.GetMinimum()<minY: 
            minY = histo.GetMinimum()

        histoiter += 1
      
    ROOT.gStyle.SetOptFit(ROOT.kFALSE)
    ROOT.gStyle.SetOptStat(ROOT.kFALSE)
    ROOT.gStyle.SetOptTitle(ROOT.kFALSE)
    
    CC = ROOT.TCanvas("CC", "", 1200, 800)

    CC.SetFillColor(10)
    CC.SetFillStyle(4000)
    CC.SetBorderSize(2)
    CC.Divide(1, 1)
    CC.Range(0,0,1,1)
    CC.SetFillColor(10)
    CC.SetBorderMode(0)
    CC.SetBorderSize(2)
    CC.SetTickx(1)
    CC.SetTicky(1)
    CC.SetLeftMargin(0.16)
    CC.SetRightMargin(0.02)
    CC.SetTopMargin(0.05)
    CC.SetBottomMargin(0.13)
    CC.SetFrameFillColor(0)
    CC.SetFrameFillStyle(0)
    CC.SetFrameBorderMode(0)
    
    pad1 = CC.GetPad(1)

    pad1.SetFillColor(0)
    pad1.SetBorderMode(0)
    pad1.SetBorderSize(2)
    pad1.SetLogx()
    pad1.SetTickx(1)
    pad1.SetTicky(1)
    pad1.SetLeftMargin(0.16)
    pad1.SetRightMargin(0.02)
    pad1.SetTopMargin(0.065)
    pad1.SetBottomMargin(0.13)
    pad1.SetFrameFillStyle(0)
    pad1.SetFrameBorderMode(0)
    pad1.SetFrameFillStyle(0)
    pad1.SetFrameBorderMode(0)
    
    pad1.cd()
    
    pad1.SetLogx()
    
    Option = ""
    for histo in scaleFactorHistograms:

        histo.SetMaximum(maxY+0.1)
        histo.SetMinimum(minY-0.1)
        histo.DrawCopy("e3" + Option)
        Option = "same"
        histo.DrawCopy("lhistsame")

    leg.Draw()

    for ext in plotformat.split('-'):
        CC.Print("./Plots/" + plottitle + '.' + ext)

def getScaleFactorHistogram(title, reader, flavour):
    
    #ptEdges = [20., 20.1, 30., 50., 70., 100., 140., 200., 300., 600., 999., 1000.]
    ptEdges = [20.]
    for xp in range(300, 691):
        ptEdges.append(math.exp(xp/100.))
    ptEdges.append(1000.)

    histo = ROOT.TH1F(title, '', len(ptEdges)-1, array('d',ptEdges))

    for xb in range(1, len(ptEdges)):
        
        pt = histo.GetXaxis().GetBinCenter(xb)
        sf  = reader.eval_auto_bounds('central', flavour, 0., pt)
        sfu = reader.eval_auto_bounds('up', flavour, 0., pt) - sf

        histo.SetBinContent(xb, sf)
        histo.SetBinError(xb, sfu)
        
    return histo 

if __name__ == '__main__':

    # Input parameters
    usage = 'usage: %prog [options]'
    parser = optparse.OptionParser(usage)
    
    parser.add_option('--taggers'     , dest='taggers'     , help='Tagger(s) to be compared'            , default='deepcsv')
    parser.add_option('--years'       , dest='years'       , help='Year(s) to be compared'              , default='2016-2017-2018')
    parser.add_option('--inputpath'   , dest='inputpath'   , help='Path where csv files are stored '    , default='')
    parser.add_option('--inputfiles'  , dest='inputfiles'  , help='Additional csv files to be compared' , default='')
    parser.add_option('--customfiles' , dest='customfiles' , help='Only use custom csv files'           , default=False, action='store_true')
    parser.add_option('--wps'         , dest='wps'         , help='Working point(s) to be compared'     , default='M')
    parser.add_option('--meastypes'   , dest='meastypes'   , help='Measurement type(s) to be compared'  , default='default')
    parser.add_option('--flavour'     , dest='flavour'     , help='Flavour to be studied'               , default='b')
    parser.add_option('--plotformat'  , dest='plotformat'  , help='Formats of the plot, e.g.: png-pdf'  , default='png')
    (opt, args) = parser.parse_args()
    
    flavour = { 'b' : 0, 'c' : 1, 'l' : 2, '0' : 0, '1' : 1, '2' : 2 }.get(opt.flavour, None)
    if flavour==None:
        print 'Wrong choice of jet flavour:', opt.flavour, '-> exiting'
        exit()

    inputFilePath = opt.inputpath if (opt.inputpath!='') else os.environ['CMSSW_BASE'] + '/src/PhysicsTools/NanoAODTools/data/btagSF/'

    systs = ROOT.vector('string')()
    systs.push_back('up')
    systs.push_back('down')

    """
    # 
    if opt.customfiles:
        supported_btagSF.clear()
                     
    if opt.inputfiles!='':
        for inputfile in opt.inputfiles.split(','):

            csvfile = inputfile.split('/')[0]
            tagger = inputfile.split('/')[1]

            #supported_btagSFsupported_btagSF
    """

    # Fill histograms
    scaleFactorHistograms = [ ]

    for tagger in opt.taggers.split('-'):
        if tagger in supported_btagSF:
            
            for year in opt.years.split('-'):
                for campaign in supported_btagSF[tagger]:
                    if year in campaign:
                            
                        calibration = ROOT.BTagCalibration(tagger, os.path.join(inputFilePath, supported_btagSF[tagger][campaign]['inputFileName']))

                        for wp in supported_btagSF[tagger][campaign]['supported_wp']:
                            if wp in opt.wps:

                                wp_btv = { "l" : 0, "m" : 1, "t" : 2 }.get(wp.lower(), None)
                                if wp_btv==None:
                                    print 'Working point', wp, 'not supported -> skipping'
                                    continue
                                    
                                reader = ROOT.BTagCalibrationReader(wp_btv, 'central', systs)

                                meastypes = opt.meastypes if (opt.meastypes!='default') else supported_btagSF[tagger][campaign]['measurement_types'][flavour]

                                for meastype in meastypes.split('-'):
                                    if meastype in supported_btagSF[tagger][campaign]['measurement_types'][flavour]:

                                        reader.load(calibration, flavour, meastype)
                                        title = tagger + '_' + year + '_' + wp + '_' + meastype
                                        scaleFactorHistograms.append(getScaleFactorHistogram(title, reader, flavour))
                                        
    # Make plot
    plottitle = opt.taggers + '_' + opt.years + '_' + opt.wps + '_' + opt.flavour
    if opt.meastypes!='default':
        plottitle += '_' + opt.meastypes 
    plotScaleFactors(plottitle, scaleFactorHistograms, opt.plotformat)
