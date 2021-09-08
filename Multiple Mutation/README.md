# MutaBindM
<font size=7> MutaBindM.py: source code

inputfiles: parameter files that will be called in MutaBindS.py.

Example: an example for running MutaBindM.py, which jobid is mexample. When you need to submit a job, you need to make a directory with .input(the format should be consistent with the jobid.input) and .pdb in it.</font>

## Installation
<font size=7> 1.Download this directory and put it in linux.
  
2.Download softwares needed in source code, including FoldX,CHARMM,VMD,NAMD,DSSP,PROVEAN.</font>

<font size=5> FoldX: http://foldxsuite.crg.eu/
  
CHARMM: https://www.charmm.org/charmm/

VMD：https://www.ks.uiuc.edu/Research/vmd/

NAMD：https://www.ks.uiuc.edu/Research/namd/

DSSP: https://swift.cmbi.umcn.nl/gv/dssp/

PROVEAN: http://provean.jcvi.org/index.php/ </font>

<font size=7> 3.Download python 2.7 and some python packages, including pandas,collections,rpy2.</font>

<font size=7> 4.Download R 3.4.0 or newer and some R packages, including randomForest,forestFloor.</font>


## Command
<font size=7> python MutaBindM.py -i jobid</font>
