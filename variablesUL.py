from array import array

'''
class variable(object):
    def __init__(self, name, title, nbins, xmin, xmax=None):
        self._name=name
        self._title=title
        self._nbins=nbins
        self._xmin=xmin
        if xmax==None:
            self._xmax=xmin[nbins]
            self._iscustom = True
        else:
            self._xmax=xmax
            self._iscustom = False

    def __str__(self):
        return  '\"'+str(self._name)+'\",\"'+str(self._title)+'\",\"'+str(self._taglio)+'\",'+str(self._nbins)+','+str(self._xmin)+','+str(self._xmax)
'''

class variable(object):
    def __init__(self, name, title, nbins, bins):
        self._name=name
        self._title=title
        self._nbins=nbins
        self._bins=bins
        self._iscustom = True
        self._xmin=bins
        self._xmax=bins[nbins]

variables = {}

variables["SR"] = []
#variables["SR"].append(variable('Leadingjet_pt',  'Lead jet p_{T} [GeV]', 9 ,array("f", [0., 50., 100., 150., 200., 250., 300., 400., 500., 600.])))
#variables["SR"].append(variable('SubLeadingjet_pt',  'Sublead jet p_{T} [GeV]', 5 , array("f", [0., 50., 100., 150., 250., 500.])))
#variables["SR"].append(variable('Lepton_pt',  'Lepton p_{T} [GeV]', 9 , array("f", [0., 30., 45., 60., 80., 100., 125., 150, 200., 250.]) ))
#variables["SR"].append(variable('SelectedTau_pt',  'Tau p_{T} [GeV]', 8 , array("f", [30., 45., 60., 80., 100., 125., 150, 200., 250.])))
#variables["SR"].append(variable('SMbdt', 'SM BDT output',  5, array("f", [0., 0.2, 0.4, 0.6, 0.8, 1.])))
variables["SR"].append(variable('mjj', 'Invarant mass JJ',  5, array("f", [500., 700., 1000., 1500., 1800., 2000.])))

variables["CR_fakes"] = []
#variables["CR_fakes"].append(variable('Leadingjet_pt',  'Lead jet p_{T} [GeV]', 5 ,array("f", [0., 50., 100., 150., 250., 400.])))
#variables["CR_fakes"].append(variable('SubLeadingjet_pt',  'Sublead jet p_{T} [GeV]', 3 , array("f", [0., 50., 100., 200.])))
#variables["CR_fakes"].append(variable('Lepton_pt',  'Lepton p_{T} [GeV]', 7 , array("f", [0., 30., 45., 60., 80., 100., 150, 250.]) ))
#variables["CR_fakes"].append(variable('SelectedTau_pt',  'Tau p_{T} [GeV]', 4 , array("f", [30., 45., 60., 100., 200.])))
#variables["CR_fakes"].append(variable('SMbdt', 'SM BDT output',  5, array("f", [0., 0.2, 0.4, 0.6, 0.8, 1.])))
variables["CR_fakes"].append(variable('mjj', 'Invarant mass JJ',  6, array("f", [0., 300., 500., 700., 1000., 1500., 2000.])))

#variables["CR_QCD"] = []
#variables["CR_QCD"].append(variable('Leadingjet_pt',  'Lead jet p_{T} [GeV]', 5 ,array("f", [0., 50., 100., 150., 250., 400.])))
#variables["CR_QCD"].append(variable('SubLeadingjet_pt',  'Sublead jet p_{T} [GeV]', 3 ,array("f", [0., 50., 100., 200.])))
#variables["CR_QCD"].append(variable('Lepton_pt',  'Lepton p_{T} [GeV]', 7 , array("f", [0., 30., 45., 60., 80., 100., 150, 250.]) ))
#variables["CR_QCD"].append(variable('SelectedTau_pt',  'Tau p_{T} [GeV]', 4 , array("f", [30., 45., 60., 100., 200.])))
#variables["CR_QCD"].append(variable('SMbdt', 'SM BDT output',  5, array("f", [0., 0.2, 0.4, 0.6, 0.8, 1.])))
#variables["CR_QCD"].append(variable('mjj', 'Invarant mass JJ',  8, array("f", [0., 150., 300., 500., 700., 1000., 1400., 1800., 3000.])))

variables["CR_ttbar"] = []
#variables["CR_ttbar"].append(variable('Leadingjet_pt',  'Lead jet p_{T} [GeV]', 9 ,array("f", [0., 50., 100., 150., 200., 250., 300., 400., 500., 600.])))
#variables["CR_ttbar"].append(variable('SubLeadingjet_pt',  'Sublead jet p_{T} [GeV]', 5 , array("f", [0., 50., 100., 150., 250., 500.])))
#variables["CR_ttbar"].append(variable('Lepton_pt',  'Lepton p_{T} [GeV]', 9 , array("f", [0., 30., 45., 60., 80., 100., 125., 150, 200., 250.]) ))
#variables["CR_ttbar"].append(variable('SelectedTau_pt',  'Tau p_{T} [GeV]', 8 , array("f", [30., 45., 60., 80., 100., 125., 150, 200., 250.])))
#variables["CR_ttbar"].append(variable('SMbdt', 'SM BDT output',  5, array("f", [0., 0.2, 0.4, 0.6, 0.8, 1.])))
variables["CR_ttbar"].append(variable('mjj', 'Invarant mass JJ',  6, array("f", [0., 300., 500., 700., 1000., 1500., 2000.])))

variables["CR_opposite_sign"] = []
#variables["CR_opposite_sign"].append(variable('Leadingjet_pt',  'Lead jet p_{T} [GeV]', 9 ,array("f", [0., 50., 100., 150., 200., 250., 300., 400., 500., 600.])))
#variables["CR_opposite_sign"].append(variable('SubLeadingjet_pt',  'Sublead jet p_{T} [GeV]', 5 , array("f", [0., 50., 100., 150., 250., 500.])))
#variables["CR_opposite_sign"].append(variable('Lepton_pt',  'Lepton p_{T} [GeV]', 9 , array("f", [0., 30., 45., 60., 80., 100., 125., 150, 200., 250.]) ))
#variables["CR_opposite_sign"].append(variable('SelectedTau_pt',  'Tau p_{T} [GeV]', 8 , array("f", [30., 45., 60., 80., 100., 125., 150, 200., 250.])))
#variables["CR_opposite_sign"].append(variable('SMbdt', 'SM BDT output',  5, array("f", [0., 0.2, 0.4, 0.6, 0.8, 1.])))
variables["CR_opposite_sign"].append(variable('mjj', 'Invarant mass JJ',  6, array("f", [0., 300., 500., 700., 1000., 1500., 2000.])))

#variables["CR_Zjets"] = []
#variables["CR_Zjets"].append(variable('Leadingjet_pt',  'Lead jet p_{T} [GeV]', 5 ,array("f", [0., 50., 100., 150., 250., 400.])))
#variables["CR_Zjets"].append(variable('SubLeadingjet_pt',  'Sublead jet p_{T} [GeV]', 3 ,array("f", [0., 50., 100., 200.])))
#variables["CR_Zjets"].append(variable('Lepton_pt',  'Lepton p_{T} [GeV]', 7 , array("f", [0., 30., 45., 60., 80., 100., 150, 250.]) ))
#variables["CR_Zjets"].append(variable('SelectedTau_pt',  'Tau p_{T} [GeV]', 4 , array("f", [30., 45., 60., 100., 200.])))
#variables["CR_Zjets"].append(variable('SMbdt', 'SM BDT output',  5, array("f", [0., 0.2, 0.4, 0.6, 0.8, 1.])))
#variables["CR_Zjets"].append(variable('mjj', 'Invarant mass JJ',  8, array("f", [0., 150., 300., 500., 700., 1000., 1400., 1800., 3000.])))