# VBS Same Sign WW analysis with a hadronic tau and a e/mu in the final state
This analysis is performed with the RDataFrame tool.


## How to access analysis facility
- Go to [https://jhub.90.147.75.109.myip.cloud.infn.it/](https://jhub.90.147.75.109.myip.cloud.infn.it/)
- Login with your Dodas credentials (to get them, register at https://dodas-iam.cloud.cnaf.infn.it/start-registration)
- Choose a JupyterLab image and set memory and CPU
- Once deployed, open a new terminal and type 
  ``` 
  oidc-token infncloud > ~/.token 
  ```
- To test that you can reach the HTCondor deployment: 
  ```
  condor_q
  ```

## Usage of RDataFrame distributed on Dask, on top of HTCondor
- Open a a new Python3 notebook
- Deploy a Dask cluster on HTCondor: 
  ```
  from dask_remote_jobqueue import RemoteHTCondor
  cluster = RemoteHTCondor()
  ```
- Once deployed, initialize the Dask client:
  ```
  from dask.distributed import Client
  client = Client(address="tcp://127.0.0.1:"+str(cluster.sched_port))
  ```
- Insert the declaration of your custom functions inside an initialization function:
  ```
  import ROOT
  
  text_file = open("postselection.h", "r")
  data = text_file.read()
  distributed = ROOT.RDF.Experimental.Distributed

  def my_initialization_function():
      ROOT.gInterpreter.Declare('{}'.format(data))
    
  distributed.initialize(my_initialization_function)
  ```
- Create a distributed RDataFrame reading a list of samples:
  ```
  chain = [<path to 1st .root file>, <path to 2nd .root file>, ...]

  df = ROOT.RDF.Experimental.Distributed.Dask.RDataFrame("<name of tree>",
                chain,
                npartitions=<number of partitions>,
                daskclient=client)
  ```
  
## The analysis
Post selection steps of the analysis can be found [here](postselectionDaskHTC.ipynb)
