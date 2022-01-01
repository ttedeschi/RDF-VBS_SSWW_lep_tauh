import ROOT as rt
#from ROOT 

# CMS_lumi
#   Initiated by: Gautier Hamel de Monchenault (Saclay)
#   Translated in Python by: Joshua Hardenbrook (Princeton)
#   Updated by:   Dinko Ferencek (Rutgers)
#

cmsText     = "CMS";
#cmsText     = "an LHC";
cmsTextFont   = 61  

writeExtraText = True
#extraText   = "Experiment"
extraText   = ""
extraText   = "Work in progress"
extraTextFont = 52 

#lumiTextSize     = 0.6
lumiTextSize     = 0.55
lumiTextOffset   = 0.2

#cmsTextSize      = 0.75
cmsTextSize      = 0.6
cmsTextOffset    = 0.2

relPosX    = 0.045 + 0.07
relPosY    = 0.035
relExtraDY = 1.2

extraOverCmsTextSize  = 0.7

lumi_13TeV = "20.1 fb^{-1}"
lumi_8TeV  = "19.7 fb^{-1}" 
lumi_7TeV  = "5.1 fb^{-1}"
lumi_sqrtS = ""

drawLogo   = False

def CMS_lumi(pad,  lumi_sqrtS,  iPosX , lepText):
    outOfFrame    = False
    if(iPosX/10==0 ): outOfFrame = True

    alignY_=3
    alignX_=2
    if( iPosX/10==0 ): alignX_=1
    if( iPosX==0    ): alignY_=1
    if( iPosX/10==1 ): alignX_=1
    if( iPosX/10==2 ): alignX_=2
    if( iPosX/10==3 ): alignX_=3
    align_ = 10*alignX_ + alignY_

    H = pad.GetWh()
    W = pad.GetWw()
    l = pad.GetLeftMargin()
    t = pad.GetTopMargin()
    r = pad.GetRightMargin()
    b = pad.GetBottomMargin()
    e = 0.025

    pad.cd()

    lumiText = ""
    lumiText += lumi_sqrtS


    latex = rt.TLatex()
    latex.SetNDC()
    latex.SetTextAngle(0)
    latex.SetTextColor(rt.kBlack)    
    
    extraTextSize = extraOverCmsTextSize*cmsTextSize
    
    latex.SetTextFont(42)
    latex.SetTextAlign(31) 
    latex.SetTextSize(lumiTextSize*t)    

    latex.DrawLatex(1-r,1-t+lumiTextOffset*t,lumiText)

    if( outOfFrame ):
        latex.SetTextFont(cmsTextFont)
        latex.SetTextAlign(11) 
        latex.SetTextSize(cmsTextSize*t)    
        latex.DrawLatex(l,1-t+lumiTextOffset*t,cmsText)
  
    pad.cd()

    posX_ = 0
    if( iPosX%10<=1 ):
        posX_ =   l + relPosX*(1-l-r)
    elif( iPosX%10==2 ):
        posX_ =  l + 0.5*(1-l-r)
    elif( iPosX%10==3 ):
        posX_ =  1-r - relPosX*(1-l-r)

    posY_ = 1-t - relPosY*(1-t-b)

#    print lumiText," 2 "

    if( not outOfFrame ):
        if( drawLogo ):
            posX_ =   l + 0.045*(1-l-r)*W/H
            posY_ = 1-t - 0.045*(1-t-b)
            xl_0 = posX_
            yl_0 = posY_ - 0.15
            xl_1 = posX_ + 0.15*H/W
            yl_1 = posY_
            CMS_logo = rt.TASImage("CMS-BW-label.png")
            pad_logo =  rt.TPad("logo","logo", xl_0, yl_0, xl_1, yl_1 )
            pad_logo.Draw()
            pad_logo.cd()
            CMS_logo.Draw("X")
            pad_logo.Modified()
            pad.cd()          
        else:
            latex.SetTextFont(cmsTextFont)
            latex.SetTextSize(cmsTextSize*t)
            latex.SetTextAlign(align_)
            latex.DrawLatex(posX_, posY_, cmsText)
            if( writeExtraText ) :
                latex.SetTextFont(extraTextFont)
                latex.SetTextAlign(align_)
                latex.SetTextSize(extraTextSize*t)
                latex.DrawLatex(posX_, posY_- relExtraDY*cmsTextSize*t, extraText)
                #aggiunta per la stampa di lepText all'interno della Canvas
                latex.SetTextFont(42)
                #latex.SetTextAlign(31)
                latex.SetTextSize(lumiTextSize*t)
                latex.DrawLatex(posX_, posY_- relExtraDY*cmsTextSize*t - relExtraDY*extraTextSize*t , lepText)      
    elif( writeExtraText ):
        if( iPosX==0):
            posX_ =   l +  relPosX*(1-l-r)
            posY_ =   1-t+lumiTextOffset*t

        latex.SetTextFont(extraTextFont)
        latex.SetTextSize(extraTextSize*t)
        latex.SetTextAlign(align_)
        latex.DrawLatex(posX_, posY_, extraText)      

        #    print lumiText," 3 "
    #pad.Update()
    #    print lumiText," 4 "

import ROOT
import os
    
def stackplot(region, feature, final_state, folder, h, blinded = False):
    
    print(region)


    #if blinded == False:
    hdata = h[region]['Data'][feature._name][final_state]
    hsig = h[region]['VBS_ssWW'][feature._name][final_state]

    if final_state == "etau":
        lep_tag = "e+"
    else:
        lep_tag = "#mu+"
    cmsreg = "#tau_{h}" 
    cmsreg = lep_tag + cmsreg 

    stackname = "stack"
    canvasname = "canvas"
    blind = False

    h_err = ROOT.TH1F()

    ROOT.gROOT.SetStyle('Plain')
    ROOT.gStyle.SetPalette(1)
    ROOT.gStyle.SetOptStat(0)
    ROOT.TH1.SetDefaultSumw2()

    # Draw stack with MC contributions
    stack = ROOT.THStack(stackname, feature._name)
    leg_stack = ROOT.TLegend(0.32,0.58,0.93,0.87)

    colors = [(222, 90, 106), (155, 152, 204), (208, 240, 193), (122, 130, 106), (200, 131, 274), (218, 190, 193), (222, 10, 106), (122, 90, 106), (22,10,67), (102,100,67), (2,10,167), (22,10,67), (22,10,67), (32,121,100), (65,54,63), (132,121,100), (4,11,100), (100,1,10), (50,79,88)]

    i = 0
    for v in h[region].keys():
        if v == 'Data' or v == 'VBS_ssWW':
            continue
        h_ = h[region][v][feature._name][final_state]
        h_.SetLineWidth(1)
        h_.SetLineColor(1)
        color = colors[i]
        h_.SetFillColor(ROOT.TColor.GetColor(*color))
        stack.Add(h_.GetValue())
        leg_stack.AddEntry(h_.GetValue(), v, "f")
        i +=1

    leg_stack.SetNColumns(2)
    leg_stack.SetFillColor(0)
    leg_stack.SetFillStyle(0)
    leg_stack.SetTextFont(42)
    leg_stack.SetBorderSize(0)
    leg_stack.SetTextSize(0.035)

    # Create canvas
    c1 = ROOT.TCanvas(canvasname,"c1",50,50,700,600)
    c1.SetFillColor(0)
    c1.SetBorderMode(0)
    c1.SetFrameFillStyle(0)
    c1.SetFrameBorderMode(0)
    c1.SetLeftMargin( 0.12 )
    c1.SetRightMargin( 0.9 )
    c1.SetTopMargin( 1 )
    c1.SetBottomMargin(-1)
    c1.SetTickx(1)
    c1.SetTicky(1)
    c1.cd()

    pad1= ROOT.TPad("pad1", "pad1", 0, 0.31 , 1, 1)
    pad1.SetTopMargin(0.1)
    pad1.SetBottomMargin(0.02)
    pad1.SetLeftMargin(0.12)
    pad1.SetRightMargin(0.05)
    pad1.SetBorderMode(0)
    pad1.SetTickx(1)
    pad1.SetTicky(1)
    pad1.Draw()
    pad1.cd()

    if blinded == False:
        maximum = max(stack.GetMaximum(),hdata.GetMaximum())
    else:
        maximum = stack.GetMaximum()
    logscale = True # False #
    if(logscale) and stack.GetStack().Last().Integral()>0.:
        stack.SetMinimum(0.01)
        pad1.SetLogy()
        stack.SetMaximum(maximum*10000)
    else:
        stack.SetMaximum(maximum*1.6)

    stack.Draw("HIST")

    if not feature._iscustom:
        step = float(feature._xmax - feature._xmin)/float(feature._nbins)
        #print(str(step))
        if "GeV" in feature._title:
            if step.is_integer():
                ytitle = "Events/ %.0f GeV" %step
            else:
                ytitle = "Events / %.2f GeV" %step
        else:
           if step.is_integer():
               ytitle = "Events / %.0f units" %step
           else:
               ytitle = "Events / %.2f units" %step
    else:
        if "GeV" in feature._title:
            ytitle = "Events / GeV"
        else:
            ytitle = "Events / a.u"

    print(stack)
    stack.GetYaxis().SetTitle(ytitle)
    stack.GetYaxis().SetTitleFont(42)
    stack.GetXaxis().SetLabelOffset(1.8)
    stack.GetYaxis().SetTitleOffset(0.85)
    stack.GetXaxis().SetLabelSize(0.15)
    stack.GetYaxis().SetLabelSize(0.07)
    stack.GetYaxis().SetTitleSize(0.07)
    stack.SetTitle("")

    #hsig.Scale(1000)
    print("VBS_ssWW")
    hsig.Draw("hist same")
    leg_stack.AddEntry(hsig.GetValue(), "VBS_ssWW", "l")

    h_err = stack.GetStack().Last().Clone("h_err")
    h_err.SetLineWidth(100)
    h_err.SetFillStyle(3154)
    h_err.SetMarkerSize(0)
    h_err.SetFillColor(ROOT.kGray+2)
    h_err.Draw("e2same0")
    leg_stack.AddEntry(h_err, "Stat. Unc.", "f")

    if blinded == False:
        print("Data")
        leg_stack.AddEntry(hdata.GetValue(), "Data", "ep")
    leg_stack.Draw("same")

    if blinded == False:
        # Draw data
        hdata.SetMarkerStyle(20)
        hdata.SetMarkerSize(1.2)
        hdata.SetLineWidth(2)
        hdata.SetLineColor(ROOT.kBlack)
        #hdata.Draw("E SAME")
        hdata.Draw("eSAMEpx0")
    else:
        hdata = stack.GetStack().Last().Clone("h_data")
    
    
    lumi = 41.53
    lumi_scale = 21

    CMS_lumi.writeExtraText = 1
    CMS_lumi.extraText = ""

    print("lep_tag: ", lep_tag)
    lumi_sqrtS = "%s fb^{-1}  (13 TeV)"%(lumi)

    iPeriod = 0
    iPos = 11
    CMS_lumi(pad1, lumi_sqrtS, iPos, str(cmsreg))

    hratio = stack.GetStack().Last()

    c1.cd()
    pad2= ROOT.TPad("pad2", "pad2", 0, 0.01 , 1, 0.30)
    pad2.SetTopMargin(0.05)
    pad2.SetBottomMargin(0.45)
    pad2.SetLeftMargin(0.12)
    pad2.SetRightMargin(0.05)
    ROOT.gStyle.SetHatchesSpacing(2)
    ROOT.gStyle.SetHatchesLineWidth(2)
    c1.cd()
    pad2.Draw()
    pad2.cd()
    ratio = hdata.Clone("ratio")
    ratio.SetLineColor(ROOT.kBlack)
    ratio.SetMaximum(10)
    ratio.SetMinimum(0)
    ratio.Sumw2()
    ratio.SetStats(0)

    ratio.Divide(hratio)
    ratio.SetMarkerStyle(20)
    ratio.SetMarkerSize(0.9)
    ratio.Draw("epx0e0")
    ratio.SetTitle("")

    h_bkg_err = hratio.Clone("h_err")
    h_bkg_err.Reset()
    h_bkg_err.Sumw2()
    for i in range(1,hratio.GetNbinsX()+1):
        h_bkg_err.SetBinContent(i,1)
        if(hratio.GetBinContent(i)):
            h_bkg_err.SetBinError(i, (hratio.GetBinError(i)/hratio.GetBinContent(i)))
        else:
            h_bkg_err.SetBinError(i, 10^(-99))
    h_bkg_err.SetLineWidth(100)

    h_bkg_err.SetMarkerSize(0)
    h_bkg_err.SetFillColor(ROOT.kGray+1)
    h_bkg_err.Draw("e20same")

    if not feature._iscustom:
        xmin = feature._xmin
    else:
        xmin = feature._xmin[0]
    f1 = ROOT.TLine(xmin, 1., feature._xmax,1.)
    xmin = 0
    xmax = 2000
    f1 = ROOT.TLine(xmin, 1., xmax,1.)
    f1.SetLineColor(ROOT.kBlack)
    f1.SetLineStyle(ROOT.kDashed)
    f1.Draw("same")

    ratio.GetYaxis().SetTitle("Data / Bkg")
    ratio.GetYaxis().SetNdivisions(503)
    ratio.GetXaxis().SetLabelFont(42)
    ratio.GetYaxis().SetLabelFont(42)
    ratio.GetXaxis().SetTitleFont(42)
    ratio.GetYaxis().SetTitleFont(42)
    ratio.GetXaxis().SetTitleOffset(1.1)
    ratio.GetYaxis().SetTitleOffset(0.35)
    ratio.GetXaxis().SetLabelSize(0.15)
    ratio.GetYaxis().SetLabelSize(0.15)
    ratio.GetXaxis().SetTitleSize(0.16)
    ratio.GetYaxis().SetTitleSize(0.16)
    ratio.GetYaxis().SetRangeUser(0,1.5)
    ratio.GetXaxis().SetTitle(feature._title)
    ratio.GetXaxis().SetLabelOffset(0.04)
    ratio.GetYaxis().SetLabelOffset(0.02)
    ratio.Draw("epx0e0same")

    c1.cd()
    c1.RedrawAxis()
    pad2.RedrawAxis()
    c1.Update()
    #c1.Draw()
    if region not in os.listdir(folder):
        os.mkdir(folder + "/" + region)
    if final_state not in os.listdir(folder + "/" + region):
        os.mkdir(folder + "/" + region + "/" + final_state)
    c1.SaveAs("{}.png".format(folder + "/" + region + "/" + final_state + "/" + feature._name))