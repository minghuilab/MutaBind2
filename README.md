# MutaBind2
## About
<font size=4> 
  
MutaBind2 evaluates the effects of mutations on protein-protein interactions for soluble complexes. It calculates changes in binding affinity upon single or multiple mutations and provides a structural model of mutated complex. The structure of a protein-protein complex is required for this method. MutaBind2 is the second version of MutaBind method. The method can be used for SARS-CoV-2 modeling.
  
</font>

## Scoring mutations with MutaBind2
<font size=4> 

We recommend that most users who just want to obtain MutaBind2 predictions use [MutaBind2 website](https://lilab.jysw.suda.edu.cn/research/mutabind2/).

</font>

## Source code releases
<font size=4> 
  
You can download [releases](https://github.com/minghuilab/MutaBind2/releases) on github.

</font>

## Installation

#### I. PREREQUISITES

<font size=4>
 
MutaBind2 requires the following software and packages.

1. DSSP

   This is available at the DSSP website.

   https://swift.cmbi.umcn.nl/gv/dssp/

2. PROVEAN

   This is available at the PROVEAN website.

   http://provean.jcvi.org/index.php/

3. FoldX

   This is available at the FoldX website.

   http://foldxsuite.crg.eu/

4. VMD

   This is available at the VMD website.

   https://www.ks.uiuc.edu/Research/vmd/

5. CHARMM

   This is available at the CHARMM website.

    https://www.charmm.org/charmm/

6. NAMD

   This is available at the NAMD website.

   https://www.ks.uiuc.edu/Research/namd/

7. Python packages: pandas, collections and rpy2

8. R packages: randomForest and forestFloor.

</font>

#### II. INSTALLATION INSTRUCTIONS

<font size=4>

1. Download and/or install prerequisites described above.

2. Download and unpack the distribution:

</font>

<font size=4>
	
	Single Mutation Model
	$ wget https://github.com/minghuilab/PremPS/archive/MutaBindS.tar.gz
	$ tar -zxvf MutaBindS.tar.gz

	Multiple Mutation Model 
	$ wget https://github.com/minghuilab/PremPS/archive/MutaBindM.tar.gz
	$ tar -zxvf MutaBindM.tar.gz

</font> 

<font size=4>

3. Change to the source directory:

</font>

<font size=4>

	$ cd MutaBindS/MutaBindM

</font> 

<font size=4>

4. Change the path parameters in MutaBindS.py/MutaBindM.py (line 15-20)/: the software configuration path and the user's working directory (workdir) need to be changed. The working directory must be lowercase, otherwise CHARMM cannot recognize it.

<font size=4>
	
	workdir = Your working directory
	pathvmd = path for running VMD software  
	pathcharmm = path for running CHARMM software
  	pathnamd = path for running NAMD software
	pathdssp = path for running DSSP software 
 	pathprovean = path for running provean software
	pathrscript = path for running Rscript
	 
</font> 


</font>

&nbsp; &nbsp; The FoldX software needs to be installed in the working directory.

#### III. RUNNING MutaBind2

<font size=4>

	$ python MutaBindS.py -i jobid

</font> 

## Platform

<font size=4>

MutaBind2 is only intended to run on *linux* operating systems.

</font>

## Issues

<font size=4>

You will need to have Python 2 and R 3.4.0 (or higher) installed.

</font>

