import re
import copy
import os
import time
from ROOT import gPad, TCanvas, TPad, TPaveText, TBox, gStyle
# from CMGTools.H2TauTau.proto.plotter.officialStyle import CMSPrelim as CMSPrelimOfficial
from CMGTools.RootTools.DataMC.Stack import Stack

can = None,
pad = None
padr = None
ratio = None
ocan = None
save = []

xtitles = {
    'l1_pt':'p_{T,#mu 1} [GeV]',
    'l1Jet_pt':'p_{T,#tau jet} [GeV]',
    'l2_pt':'p_{T,#mu 2} [GeV]',
    'l2Jet_pt':'p_{T,#mu jet} [GeV]',
    'l1_eta':'#eta_{#mu 1}',
    'l2_eta':'#eta_{#mu 2}',
    'l1_rawMvaIso':'MVA #tau iso',
    'l1_relIso05':'#mu 1 iso',
    'l2_relIso05':'#mu 2 iso',
    'mt':'m_{T} [GeV]',
    'visMass':'m_{vis} [GeV]',
    'svfitMass':'m_{#tau#tau} [GeV]',
    'nJets':'N_{jets}',
    'nVert':'N_{vertices}',
    'jet1_pt':'p_{T,jet1} [GeV]',
    'jet2_pt':'p_{T,jet2} [GeV]',
    'jet1_eta':'#eta_{jet1}',
    'jet2_eta':'#eta_{jet2}',
    'diL_mass':'m_{ll} [GeV]',
    'diL_pt':'p_{T}^{ll} [GeV]',
    'diL_eta':'#eta^{ll} [GeV]',
    'VBF_deta':'#Delta #eta jets',
    'jet1_puMvaFull':'PU full MVA (jet 1)',
    'jet1_puMvaFull53X':'PU full MVA 53X (jet 1)',
    'jet1_puMvaSimple':'PU simple MVA (jet 1)',
    'jet1_puJetId':'PU jet ID, old WP (jet 1)',
    'jet1_puJetId53X':'PU jet ID, 53X WP (jet 1)',
    'leg1_relIso05':'Relative isolation 0.5 (muon 1)',
    'leg2_relIso05':'Relative isolation 0.5 (muon 2)',
    'jet1_diL_ptBalance':'p_{T} balance (jet 1, ll) [GeV]',
    'jet1_diL_dphi':'#Delta#phi (jet 1, ll) [GeV]',
    'puJetID53XnewWP':'PU full MVA 53X new WP',
    'puJetIDOldMVA':'PU previous full MVA previous WP',
    }

xtitles_TauEle = {
    'l1_pt'       :'p_{T,#tau} [GeV]',
    'l1Jet_pt'    :'p_{T,#tau jet} [GeV]',
    'l2_pt'       :'p_{T,e} [GeV]',
    'l2Jet_pt'    :'p_{T,e jet} [GeV]',
    'l1_eta'      :'#eta_{#tau}',
    'l2_eta'      :'#eta_{e}',
    'l1_phi'      :'#phi_{#tau}',
    'l2_phi'      :'#phi_{e}',
    'l1_charge'   :'charge_{#tau}',
    'l2_charge'   :'charge_{e}',
    'l1_rawMvaIso':'MVA #tau iso',
    'l1_relIso05' :'#tau iso',
    'l2_relIso05' :'e iso',
    'mt'          :'m_{T} [GeV]',
    'met'         :'MET [GeV]',
    'pfmet'       :'PF MET [GeV]',
    'visMass'     :'m_{vis} [GeV]',
    'svfitMass'   :'m_{#tau#tau} [GeV]',
    'nJets'       :'N_{jets}',
    'nVert'       :'N_{vertices}',
    'jet1_pt'     :'p_{T,jet1} [GeV]',
    'jet2_pt'     :'p_{T,jet1} [GeV]',
    'jet1_eta'    :'#eta_{#tau}',
    'jet2_eta'    :'#eta_{#tau}',
    'l1_EOverp'   :'#tau E/p',
    }

translatedVars = {
    'jet1_diL_ptBalance' : "sqrt( (jet1_pt * cos(jet1_phi) + diL_pt * cos(diL_phi) )**2 + (jet1_pt * sin(jet1_phi) + diL_pt * sin(diL_phi) )**2 )",
    'jet1_diL_dphi':'abs((jet1_phi - diL_phi > TMath::Pi()) * (jet1_phi - diL_phi - 2 * TMath::Pi()) + (jet1_phi - diL_phi < - TMath::Pi()) * (jet1_phi - diL_phi + 2 * TMath::Pi()) + (abs(jet1_phi - diL_phi) <= TMath::Pi()) * (jet1_phi - diL_phi))',
    'puJetID53XnewWP':'( (abs(jet1_eta) < 2.5 && jet1_puMvaFull53X > -0.63) || (abs(jet1_eta) < 2.75 && jet1_puMvaFull53X > -0.60) || (abs(jet1_eta) < 3.0 && jet1_puMvaFull53X > -0.55) || (abs(jet1_eta) >= 3.0 && jet1_puMvaFull53X > -0.45) )',
    'puJetIDOldMVA':'( (abs(jet1_eta) < 2.5 && jet1_puMvaFull > -0.8) || (abs(jet1_eta) < 2.75 && jet1_puMvaFull > -0.74) || (abs(jet1_eta) < 3.0 && jet1_puMvaFull > -0.68) || (abs(jet1_eta) >= 3.0 && jet1_puMvaFull > -0.77) )'
}

def savePlot(name):
    if gPad is None:
        print 'no active canvas'
        return
    fileName = '%s/%s' % (anaDir, name)
    print 'pad', gPad.GetName(), 'saved to', fileName    
    gPad.SaveAs( fileName )   


def buildCanvas():
    global can, pad, padr
    can = TCanvas('can','',800,800)
    can.cd()
    can.Draw()
    sep = 0.35
    pad = TPad('pad','',0.01,sep,0.99, 0.99)
    pad.SetBottomMargin(0.04)
    padr = TPad('padr','',0.01, 0.01, 0.99, sep)
    padr.SetTopMargin(0.04)
    padr.SetBottomMargin(0.3)
    pad.Draw()
    padr.Draw()
    return can, pad, padr


def datasetInfo(plot):
    year = ''
    if plot.dataComponents[0].find('2012')!=-1:
        year = '2012'
        energy = 8
    elif plot.dataComponents[0].find('2011')!=-1:
        year = '2011'
        energy = 7       
    lumi = plot.intLumi/1e3
    return year, lumi, energy 

def insertText(plot, pad, insetText, lowY=0.77):
    pad.cd()
    theStr = insetText
    lowX = 0.15

    if not hasattr(plot, 'insertTexts'):
        plot.insertTexts = []

    plot.insertTexts.append(TPaveText(lowX, lowY, lowX+0.8, lowY+0.16, "NDC"))
    txt = plot.insertTexts[-1]
    txt.SetBorderSize(   0 )
    txt.SetFillStyle(    0 )
    txt.SetTextAlign(   12 )
    txt.SetTextSize ( 0.05 )
    txt.SetTextFont (   62 )
    txt.AddText(theStr)
    txt.Draw('same')

def CMSPrelim(plot, pad, channel ):
    pad.cd()
    year, lumi, energy = datasetInfo( plot )
    theStr = 'CMS Prelim. {year}, {lumi:3.3} fb^{{-1}}, #sqrt{{s}} = {energy:d} TeV'.format( year=year, lumi=lumi, energy=energy)
    lowX = 0.21
    lowY = 0.87
    plot.cmsprel = TPaveText(lowX, lowY, lowX+0.8, lowY+0.16, "NDC")
    plot.cmsprel.SetBorderSize(   0 )
    plot.cmsprel.SetFillStyle(    0 )
    plot.cmsprel.SetTextAlign(   12 )
    plot.cmsprel.SetTextSize ( 0.05 )
    plot.cmsprel.SetTextFont (   62 )
    plot.cmsprel.AddText(theStr)
    plot.cmsprel.Draw('same')
    
    plot.chan = TPaveText(0.85, lowY, 0.90, lowY+0.18, "NDC")
    plot.chan.SetBorderSize(   0 )
    plot.chan.SetFillStyle(    0 )
    plot.chan.SetTextAlign(   12 )
    plot.chan.SetTextSize ( 0.1 )
    plot.chan.SetTextFont (   62 )
    plot.chan.AddText(channel)
    plot.chan.Draw('same')


unitpat = re.compile('.*\((.*)\)\s*$')

keeper = []


def draw(plot, channel='MuMu', plotprefix=None, SetLogy=0, inset=None, insets=[], directory=''):
    # print plot
    Stack.STAT_ERRORS = True
    blindxmin = None
    blindxmax = None
    titles = xtitles
    if channel=='TauEle':
        titles = xtitles_TauEle
    xtitle = titles.get( plot.varName, None )
    if xtitle is None:
        xtitle = ''
    global can, pad, padr, ratio
    if pad is None:
        can, pad, padr = buildCanvas()
    pad.cd()
    pad.SetLogy (SetLogy)
    plot.DrawStack('HIST')
    h = plot.supportHist
    h.GetXaxis().SetLabelColor(0)
    h.GetXaxis().SetLabelSize(0)
    gevperbin = h.GetXaxis().GetBinWidth(1)
    h.GetYaxis().SetTitle('Events')
    h.GetYaxis().SetTitleOffset(1.4)
    padr.cd()
    ratio = copy.deepcopy(plot)
    ratio.legendOn = False
    ratio.DrawRatioStack('HIST', ymin=0.4, ymax=1.6)
    hr = ratio.stack.totalHist
    # hr.weighted.Fit('pol1')
    hr.GetYaxis().SetNdivisions(4)
    hr.GetYaxis().SetTitle('Exp./Obs.')
    hr.GetYaxis().SetTitleSize(0.1)
    hr.GetYaxis().SetTitleOffset(0.5)
    hr.GetXaxis().SetTitle('{xtitle}'.format(xtitle=xtitle))
    hr.GetXaxis().SetTitleSize(0.13)
    hr.GetXaxis().SetTitleOffset(0.9)
    rls = 0.075
    hr.GetYaxis().SetLabelSize(rls)
    hr.GetXaxis().SetLabelSize(rls)
    hr.GetYaxis().SetRangeUser(0.5, 1.5)
    padr.Update()    
    # blinding
    if plot.blindminx:
        pad.cd()
        max = plot.stack.totalHist.GetMaximum()
        box = TBox( plot.blindminx, 0,  plot.blindmaxx, max )
        box.SetFillColor(1)
        box.SetFillStyle(3004)
        box.Draw()
        # import pdb; pdb.set_trace()
        keeper.append(box)
    print channel

    if inset: insertText(plot, pad, inset)
    lowY = 0.67
    for inset2 in insets:
        insertText(plot, pad, inset2, lowY=lowY)
        lowY -= 0.10

    if channel == 'MuMu' : CMSPrelim( plot, pad, '#mu#mu')
    can.cd()
    if plotprefix == None : plotname = plot.varName
    else : plotname = plotprefix + '_' + plot.varName
    if directory:
        if not os.path.exists(os.path.join(os.getcwd(), directory)):
            os.mkdir(os.path.join(os.getcwd(), directory))
    can.SaveAs( os.path.join(directory, plotname + '.png'))
    pad.SetLogy(0)


def buildCanvasOfficial():
    global ocan
    ocan = TCanvas('ocan','',600,600)
    ocan.cd()
    ocan.Draw()
    return ocan 



## def drawOfficial(plot, doBlind=False, channel='TauMu', plotprefix = None, ymin = 0.1):
##     global ocan
##     print plot
##     Stack.STAT_ERRORS = False
##     blindxmin = None
##     blindxmax = None
##     doBlind = (plot.varName == 'svfitMass') and doBlind
##     if doBlind:
##         blindxmin = 100
##         blindxmax = 160
##         plot.Blind(blindxmin, blindxmax, False)
##     titles = xtitles
##     if channel=='TauEle':
##         titles = xtitles_TauEle
##     xtitle = titles.get( plot.varName, None )
##     if xtitle is None:
##         xtitle = ''
##     global ocan
##     if ocan is None:
##         ocan = buildCanvasOfficial()
##     ocan.cd()
##     plot.DrawStack('HIST', ymin=ymin)
##     h = plot.supportHist
##     h.GetXaxis().SetTitle('{xtitle}'.format(xtitle=xtitle))
##     # blinding
##     if plot.blindminx:
##         ocan.cd()
##         max = plot.stack.totalHist.GetMaximum()
##         box = TBox( plot.blindminx, 0,  plot.blindmaxx, max )
##         box.SetFillColor(1)
##         box.SetFillStyle(3004)
##         box.Draw()
##         # import pdb; pdb.set_trace()
##         keeper.append(box)
##     year, lumi, energy = datasetInfo( plot )
##     datasetStr = "CMS Preliminary, #sqrt{{s}} = {energy} TeV, L = {lumi:3.2f} fb^{{-1}}".format(energy=energy, lumi=lumi)
##     if channel == 'TauMu' : a,b = CMSPrelimOfficial( datasetStr, '#tau_{#mu}#tau_{h}',0.15,0.835)
##     elif channel == 'TauEle' : a,b = CMSPrelimOfficial( datasetStr, '#tau_{e}#tau_{h}', 0.15, 0.835)
##     a.Draw()
##     b.Draw()
##     save.extend([a,b])
##     ocan.Modified()    
##     ocan.Update()    
##     ocan.cd()
##     if plotprefix == None : plotname = plot.varName
##     else : plotname = plotprefix + '_' + plot.varName
##     ocan.SaveAs( plotname + '.png')
##     ocan.SetLogy()
##     ocan.SaveAs( plotname + '_log.png')
##     ocan.SetLogy(0)


cantemp = None

def drawTemplates(plot):
    global cantemp
    cantemp = TCanvas('cantemp', 'templates', 800,800)
    cantemp.cd()
    for hist in plot.histosDict.values():
        if hist.name not in ['TTJets','Ztt','WJets', 'QCD']:
            continue
        hist.Draw()
        gPad.Update()
        cantemp.SaveAs(hist.name + '.png')
        time.sleep(1)
