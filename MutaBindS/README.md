# MutaBindS
<font size=7> MutaBindS.py: source code

inputfiles: parameter files that will be called in MutaBindS.py.

Example: an example for running MutaBindS.py, which jobid is sexample. When you need to submit a job, you need to make a directory with .input(the format should be consistent with the jobid.input) and .pdb in it.</font>

## Installation
<font size=7> MutaBind2 requires the following software and packages.</font>

<font size=5> FoldX: http://foldxsuite.crg.eu/
  
CHARMM: https://www.charmm.org/charmm/

VMD：https://www.ks.uiuc.edu/Research/vmd/

NAMD：https://www.ks.uiuc.edu/Research/namd/

DSSP: https://swift.cmbi.umcn.nl/gv/dssp/

PROVEAN: http://provean.jcvi.org/index.php/ 

Python packages: pandas,collections and rpy2.</font>

R packages: randomForest and forestFloor.</font>

&nbsp; &nbsp; The FoldX software needs to be installed in the working directory.

## RUNNING MutaBind2
<font size=7> python MutaBindS.py -i jobid</font>

## Platform

<font size=4>

PremPS is only intended to run on *linux* operating systems.

</font>

## Issues

<font size=4>

You will need to have Python 2 (or 3) and R 3.4.0 (or higher) installed.

</font>
