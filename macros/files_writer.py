# prende i path da path_writer (crab_paths.txt) e usando gfal-ls scorre su tutti i file e li salva su un file .txt
from PhysicsTools.NanoAODTools.postprocessing.samples.samples import *
import os
import optparse

usage = 'python files_writer.py -d sample_name'
parser = optparse.OptionParser(usage)
parser.add_option('-d', '--dat', dest='dat', type=str, default = '', help='Please enter a dataset name')
(opt, args) = parser.parse_args()

if not(opt.dat in sample_dict.keys()):
    print sample_dict.keys()
dataset = sample_dict[opt.dat]
samples = []

if hasattr(dataset, 'components'): # How to check whether this exists or not
    samples = [sample for sample in dataset.components]# Method exists and was used.
else:
    print "You are launching a single sample and not an entire bunch of samples"
    samples.append(dataset)

path = ".."

for sample in samples:
    if not os.path.exists("./files/"):
        os.makedirs("./files/")
    f = open("./files/"+str(sample.label)+"_bu.txt", "w")
    #url = os.popen('crab getoutput --xrootd --quantity="all" -d ' + path + '/crab_' + str(sample.label) + '/').readlines()#[0]
    url = os.popen('crab getoutput --xrootd --jobids=1 -d ' + path + '/crab_' + str(sample.label) + '/').readlines()[0]
    '''
    print("Printing out crabbed files for "+str(sample.label))
    for u in url:
        #f.write
        if 'root' not in u:
            npaths = int(u.split("for ")[-1].split(" ")[0])

            finished = False
            rang = 500
            t = 0
            while not finished:
                intmin = int(t*500+1)
                intmax = int(min((t+1)*500, npaths))
                crabgo = str(intmin)+'-'+str(intmax)
                curl = os.popen('crab getoutput --xrootd --jobids=' + str(crabgo) + ' -d ' + path + '/crab_' + str(sample.label) + '/').readlines()
                for cu in curl:
                    f.write(cu)
                if intmax == npaths:
                    finished = True
                t+=1
                
            break
        else:
            f.write(u)
f.close()

    '''
    s1=url.split(str(os.environ.get('USER')))
    print(s1[1])
    s2=s1[1].split('0000')
    
    path_xrd = 'root://cms-xrd-global.cern.ch//store/user/' + str(os.environ.get('USER')) + s2[0]
    newurl = 'srm://stormfe1.pi.infn.it/cms/store/user/' + str(os.environ.get('USER')) + s2[0]
    #path_xrd = 'root://cms-xrd-global.cern.ch//store/user/adeiorio/OutDir/JetHT/DataHTB_2016/200511_122826/'
    #newurl = 'srm://stormfe1.pi.infn.it/cms/store/user/adeiorio/OutDir/JetHT/DataHTB_2016/200511_122826/'

    #newurl = path_xrd
    #print(newurl)

    i=0
    print('\nChecking files in the folder '+newurl.strip('\n')+'\n')
    while True:
        folder = os.popen('eval `scram unsetenv -sh`; gfal-ls '+ newurl.strip('\n')+'000'+str(i)).readlines()
        print(folder)
        newpath_xrd = path_xrd.strip('\n')+'000'+str(i)

        if(len(folder)==0):
            print("The folder does not exist: "+ str(folder))
            break
        print('sottocartella: '+'000'+str(i))
        for file in range(len(folder)):
            if(folder[file].strip('\n') == 'log'):
                continue
            f.write(newpath_xrd+'/'+folder[file]) 
        i+=1

    f.close()
    print(" ")
    print("The file files/" + str(sample.label) + ".txt has been created.")
    
    #print('Setting cmsenv again')
    #os.popen("eval `scram runtime -csh`")
    print('Goodbye')

