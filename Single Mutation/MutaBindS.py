#!/usr/bin/python
# coding=utf-8
import os, sys, os.path, math, commands, re, getopt, datetime,time
import pandas as pd
import rpy2.robjects as robjects
from retrying import retry
from string import ascii_uppercase
from string import ascii_lowercase
from collections import defaultdict
from multiprocessing import Pool

ascii_cases = ascii_uppercase + ascii_lowercase + ''.join([str(i) for i in range(0,10)]) # if this case, the maximal inputchains should be 52.
r = robjects.r
r('''library(randomForest)''')

# set up path for each job, which need to be changed
workdir = "/data/webservice/mutabind2/"  
pathoutput = workdir + jobid + "out"  

# set up path for softwares
pathvmd = '/usr/local/bin/'
pathcharmm = '/usr/local/bin/'
pathnamd = '/home/lim9/softwares/NAMD_2.12_Source/Linux-x86_64-g++/'
pathdssp = '/usr/local/bin/'
pathprovean = ''
pathrscript = '/usr/local/bin/'

# set up path for jobpath
jobpath = workdir + jobid

# mkdir for joboutput
pathoutput = jobpath + "out"
os.system("mkdir %s" % pathoutput)

# set up path for inputfiles
pathpara = workdir + "sinputfiles"
pathinput = workdir + "sinputfiles"

# set up input and output filename
in_file = jobpath + '/' + jobid + '.input'  # inputfile name
out_file = jobpath + '/' + jobid + '.sunddg'  # outputfile name

# input jobid
jobid = ''

#####################################
# Step-1: Parse parameters
#####################################
myopts, args = getopt.getopt(sys.argv[1:], "i:z")

for o, a in myopts:
    if o == '-i':
        jobid = a
    else:
        print("Usage: %s -i jobid" % sys.argv[0])
	
normal_format_pro = ['CYS','GLN','ILE','SER','VAL','MET','ASN','PRO','LYS','THR','PHE','ALA','HIS','GLY','ASP','LEU','ARG','TRP','GLU','TYR']

# map residue name three letters to one
map = {"GLY": "G", "ALA": "A", "SER": "S", "THR": "T", "CYS": "C",
       "VAL": "V", "LEU": "L", "ILE": "I", "MET": "M", "PRO": "P",
       "PHE": "F", "TYR": "Y", "TRP": "W", "ASP": "D", "GLU": "E",
       "ASN": "N", "GLN": "Q", "HIS": "H", "LYS": "K", "ARG": "R",
       "ASX": "X", "GLX": "X", "CSO": "X", "HIP": "X", "MSE": "X",
       "UNK": "X", "SEC": "X", "PYL": "X", "SEP": "X", "TPO": "X",
       "PTR": "X", "XLE": "X", "XAA": "X", "HSD": "H", "HID": "H",
       "HSE": "H"}
# for pro and NA
map1 = {"GLY": "G", "ALA": "A", "SER": "S", "THR": "T", "CYS": "C",
        "VAL": "V", "LEU": "L", "ILE": "I", "MET": "M", "PRO": "P",
        "PHE": "F", "TYR": "Y", "TRP": "W", "ASP": "D", "GLU": "E",
        "ASN": "N", "GLN": "Q", "HSD": "H", "LYS": "K", "ARG": "R",
        "ASX": "X", "GLX": "X", "CSO": "X", "HIP": "X", "MSE": "X",
        "UNK": "X", "SEC": "X", "PYL": "X", "SEP": "X", "TPO": "X",
        "PTR": "X", "XLE": "X", "XAA": "X", }

# one letter to three !! H to HSD otherwise charmm cannot recognize HIS
mapr = {"G": "GLY", "A": "ALA", "S": "SER", "T": "THR", "C": "CYS",
        "V": "VAL", "L": "LEU", "I": "ILE", "M": "MET", "P": "PRO",
        "F": "PHE", "Y": "TYR", "W": "TRP", "D": "ASP", "E": "GLU",
        "N": "ASN", "Q": "GLN", "H": "HSD", "K": "LYS", "R": "ARG",
        "X": "ASX", "X": "GLX", "X": "CSO", "X": "HIP", "X": "MSE",
        "X": "UNK", "X": "SEC", "X": "PYL", "X": "SEP", "X": "TPO",
        "X": "PTR", "X": "XLE", "X": "XAA", }

# process biounit model and write another colummn for indicating chains information from inputfiles
def ProPDB1():
    pdball = []
    f = open(in_file, 'r')
    _unused = f.next()
    for line in f:
        ff = line.split("\t")
        pdb = ff[0]
        if pdb not in pdball:
            ffpdb = open(pathoutput + '/' + pdb.split(".")[0].lower() + '_p.pdb', 'w')
            fpdb = open(jobpath + '/' + pdb, 'r')
            if pdb.split('.')[1].strip() != 'pdb':
                ST_play = False
                for linepdb in fpdb:
                    if linepdb[0:5] == "MODEL":
                        CountModel = linepdb.split()[1]
                        ST_play = True
                        continue
                    if ST_play:
                        if linepdb[:4] == "ATOM" and linepdb[17:20].strip() in normal_format_pro:
                            ffpdb.write("%s                  %s\n" % (linepdb[0:54].strip('\r\n'), str(linepdb[21:22]) + '_' + str(CountModel)))
            else:
                ST_play = False
                for linepdb in fpdb:
                    if linepdb[:4] == "ATOM" and linepdb[17:20].strip() in normal_format_pro:
                        ffpdb.write("%s                  %s\n" % (linepdb[0:54].strip('\r\n'), str(linepdb[21:22]) + '_' + str(1)))
                        ST_play = True
                    if ('MODEL' in linepdb) and ST_play:
                        if (linepdb.strip().split()[0] == 'MODEL') and (linepdb.strip().split()[1] == '2'):
                            break
            pdball.append(pdb)
            fpdb.close()
            ffpdb.close()
        else:
            continue
    f.close()


def del_unknown_incomplete():
    residues = [k for k, v in map.iteritems() if v!= 'X']
    f = open(in_file, 'r')
    f.next()
    for line in f:
        ff = line.split("\t")
        pdbid = ff[0].split('.')[0].lower()
        pdbfile = open('{}/{}_p.pdb'.format(pathoutput,pdbid)).readlines()
        # delete unknow resdues
        pdbfile_del_unknown = [i for i in pdbfile if i[17:20] in residues]
        # delete incomplete residues
        final_row = pdbfile_del_unknown[-1]  
        last = ''          
        above = []        
        allresidues = [] 
        for row in pdbfile_del_unknown:     
            if row[17:26] == last and row == final_row:  # when read final row，append it if equal to last
                above.append(row)
                atoms = [i[13:16].strip() for i in above]
                if set(['C','N','O','CA']).issubset(set(atoms)):
                    allresidues.append(above)
            elif row[17:26] == last and row != final_row:  # when read same residue, but not last row
                above.append(row)
            else:  # when read different residue
                if len(above)>=4:
                    atoms = [i[13:16].strip() for i in above]  
                    if set(['C','N','O','CA']).issubset(set(atoms)):
                        allresidues.append(above)
                above = [row] 
            last = row[17:26]
        # write out
        with open('{}/{}_p_test.pdb'.format(pathoutput,pdbid),'w') as fw:
            fw.write(''.join([y for x in allresidues for y in x]))
        break
    os.system('mv {}/{}_p_test.pdb {}/{}_p.pdb'.format(pathoutput,pdbid, pathoutput,pdbid))
    f.close()


# split chains and produce pdb files for each chain of wild type pdb.
def splitchain():
    pdball = []
    f = open(in_file, 'r')
    _unused = f.next()
    for line in f:
        ff = line.split("\t")
        pdbfile = ff[0]
        pdb = pdbfile.split(".")[0].upper()
        partner1 = ff[1].split(".")
        partner2 = ff[2].split(".")
        if pdbfile not in pdball:
            for chains in list(partner1 + partner2):
                #os.system('grep "%s" %s/%s_p.pdb > %s/%s_%s.pdb' % (chains, pathoutput, pdb.lower(), pathoutput, pdb, chains))
                f1 = open(pathoutput + '/' + pdb.lower() + '_p.pdb', 'r')
                fw = open(pathoutput + '/' + pdb + '_' + chains + '.pdb', 'w')
                for line1 in f1:
                    if chains == line1[72:].strip():
                        fw.write(line1)
                f1.close()
                fw.close()
            pdball.append(pdbfile)
        else:
            continue

    f.close()


# processing pdb files and prepair input files for provean
def CleanPdb():
    first_line = open(in_file, 'r').readlines()[0][:-1]
    fw = open(in_file + ".cleaned", "w")
    fw.write("%s\t%s\t%s\t%s\t%s\n" % (first_line, "PDBid", "NewPartner1", "NewPartner2", "Mutation_cleaned"))

    second_line = file(in_file, 'r').readlines()[1]
    ff = second_line.split("\t")
    pdb = ff[0].split(".")[0].upper()
    partner1 = ff[1].split(".")
    partner2 = ff[2].split(".")

    mapchainarray = []
    counti = 0
    for chains in list(partner1 + partner2):
        cc = (chains, ascii_cases[counti])
        mapchainarray.append(cc)
        counti += 1
        mapchaindict = dict(iter(mapchainarray))

    newpartner1 = ''
    for chains in list(partner1):
        newpartner1 += mapchaindict[chains]
    newpartner2 = ''
    for chains in list(partner2):
        newpartner2 += mapchaindict[chains]

    countchain = 1
    for chains in list(partner1 + partner2):
        fwpdb = open(pathoutput + "/" + pdb + "_CH" + str(countchain) + ".pdb", "w")
        fvar = open(pathoutput + "/" + pdb + "_" + chains + ".var", "w")
        countchain += 1
        count = 1
        fpdb = open(pathoutput + "/" + pdb + "_" + chains + ".pdb", "r")
        resname = fpdb.readlines()[0][17:20].strip()
        fpdb = open(pathoutput + "/" + pdb + "_" + chains + ".pdb", "r")
        resnum = fpdb.readlines()[0][22:27].strip()
        line1 = open(pathoutput + "/" + pdb + "_" + chains + ".pdb", "r").readlines()[0]
        fpdb = open(pathoutput + "/" + pdb + "_" + chains + ".pdb", "r")
        atomname = fpdb.readlines()[0][13:16].strip()
        fwpdb.write("%s %s%s   %s %s" % (line1[0:16], line1[17:21], mapchaindict[chains], str(count), line1[27:]))

        mutchainall = []
        f = open(in_file, 'r')
        _unused = f.next()
        for line in f:
            ff = line.split("\t")
            mut = ff[4].upper()
            mutchain = ff[3]
            if str(chains) == str(mutchain) and mutchain not in mutchainall:
                fseq = open(pathoutput + "/" + pdb + "_" + mutchain + ".seq", "w")
                fseq.write("%s %s\n" % (">", pdb + mutchain))
                fseq.write("%s" % (map[resname]))
            mutchainall.append(mutchain)
            if str(resnum) == str(mut[1:-1]) and str(map[resname]) == str(mut[0:1]) and str(chains) == str(mutchain):
                fw.write("%s\t%s\t%s\t%s\t%s\n" % (line[:-1], pdb, newpartner1, newpartner2, str(map[resname]) + mapchaindict[chains] + str(count) + mut[-1:]))
                fvar.write("%s\n" % (str(map[resname]) + str(count) + mut[-1:]))
        f.close()

        fpdb = open(pathoutput + "/" + pdb + "_" + chains + ".pdb", "r")
        for linepdb in fpdb.readlines()[1:]:
            resnamepdb = linepdb[17:20].strip()
            resnumpdb = linepdb[22:27].strip()
            atomnamepdb = linepdb[13:16].strip()
            if resnamepdb == resname and resnumpdb == resnum:
                if atomnamepdb != atomname:
                    if count < 10:
                        fwpdb.write("%s %s%s   %s %s" % (linepdb[0:16], linepdb[17:21], mapchaindict[chains], str(count), linepdb[27:]))
                    if count >= 10 and count < 100:
                        fwpdb.write("%s %s%s  %s %s" % (linepdb[0:16], linepdb[17:21], mapchaindict[chains], str(count), linepdb[27:]))
                    if count >= 100 and count < 1000:
                        fwpdb.write("%s %s%s %s %s" % (linepdb[0:16], linepdb[17:21], mapchaindict[chains], str(count), linepdb[27:]))
                    if count >= 1000 and count < 10000:
                        fwpdb.write("%s %s%s%s %s" % (linepdb[0:16], linepdb[17:21], mapchaindict[chains], str(count), linepdb[27:]))
                else:
                    continue
            else:
                if atomnamepdb != atomname:
                    count += 1
                    if count < 10:
                        fwpdb.write("%s %s%s   %s %s" % (linepdb[0:16], linepdb[17:21], mapchaindict[chains], str(count), linepdb[27:]))
                    if count >= 10 and count < 100:
                        fwpdb.write("%s %s%s  %s %s" % (linepdb[0:16], linepdb[17:21], mapchaindict[chains], str(count), linepdb[27:]))
                    if count >= 100 and count < 1000:
                        fwpdb.write("%s %s%s %s %s" % (linepdb[0:16], linepdb[17:21], mapchaindict[chains], str(count), linepdb[27:]))
                    if count >= 1000 and count < 10000:
                        fwpdb.write("%s %s%s%s %s" % (linepdb[0:16], linepdb[17:21], mapchaindict[chains], str(count), linepdb[27:]))
                    #print("DEBUG: how often we go in this branch", resnamepdb, resname, resnumpdb, resnum, atomnamepdb, atomname)

                    mutchainall = []
                    with open(in_file, 'r') as f:
                        _unused = f.next()
                        for line in f:
                            ff = line.split("\t")
                            mut = ff[4]
                            mutchain = ff[3]
                            if str(chains) == str(mutchain) and mutchain not in mutchainall:
                                fseq.write("%s" % (map[resnamepdb]))
                            mutchainall.append(mutchain)
                            if str(resnumpdb) == str(mut[1:-1]) and str(map[resnamepdb]) == str(mut[0:1]) and str(chains) == str(mutchain):
                                fw.write("%s\t%s\t%s\t%s\t%s\n" % (
                                line[:-1], pdb, newpartner1, newpartner2, str(map[resnamepdb]) + mapchaindict[chains] + str(count) + mut[-1:]))
                                fvar.write("%s\n" % (str(map[resnamepdb]) + str(count) + mut[-1:]))
                else:
                    continue

            resname = linepdb[17:20].strip()
            resnum = linepdb[22:27].strip()
            atomname = linepdb[13:16].strip()

        fpdb.close()
        fwpdb.close()
        fvar.close()
    fseq.close()
    fw.close()


# produce wild type pdb structure for DSSP running.
def wtpdb():
    pdball = []
    f = open(in_file + ".cleaned", 'r')
    _unused = f.next()
    for line in f:
        ff = line.split("\t")
        pdb = ff[7]
        if pdb not in pdball:
            commands.getoutput('cat %s/%s_CH*.pdb > %s/%s.pdb' % (pathoutput, pdb, pathoutput, pdb))
            pdball.append(pdb)
        else:
            continue
    f.close()


# produce psf and pdb files of wild-type with vmd.
def vmd_wt():
    pdball = []
    f = open(in_file + ".cleaned", 'r')
    _unused = f.next()
    template = open(pathpara + '/vmd.pgn').read()
    for line in f:
        ff = line.split("\t")
        pdb = ff[7]
        partner1 = ff[8]
        partner2 = ff[9]
        NumChain = int(len(partner1 + partner2))
        if pdb not in pdball:
            vmd_pdb = template.replace('protname', pdb).replace('NumChain', str(NumChain)).replace('pathinput', pathinput).replace('pathoutput', pathoutput)
            with open(pathoutput + '/vmd_' + pdb + '.pgn', 'w') as fw:
                fw.write(vmd_pdb)
            commands.getoutput('%svmd -dispdev text -e %s/vmd_%s.pgn' % (pathvmd, pathoutput, pdb))
            pdball.append(pdb)
        else:
            continue
    f.close()
    # check error in *_vmd.pdb
    global vmderror
    fvmdpdb = ''.join(open(pathoutput+'/'+pdb+'_vmd.pdb').readlines())
    vmderror =  'ATOM  *****' in fvmdpdb


# produce charmm psf and crd files for individual chains
def charmmfile_wt():
    pdball = []
    f = open(in_file + ".cleaned", 'r')
    _unused = f.next()
    for line in f:
        ff = line.split("\t")
        pdb = ff[7].lower()
        partner1 = ff[8]
        partner2 = ff[9]
        NumP = int(len(partner1 + partner2))
        if pdb not in pdball:
            countp = 1
            for chains in (list(partner1) + list(partner2)):
                os.system('grep "^.\{21\}%s" %s/%s_vmd.pdb > %s/%s_%s.pdb' % (chains, pathoutput, pdb.upper(), pathoutput, pdb, 'ch' + str(countp)))
                ffpdb = open(pathoutput + '/' + pdb + '_ch' + str(countp) + '.pdb', 'a')
                ffpdb.write("%s" % ('END'))
                ffpdb.close()
                os.system('%scharmm <%s/setup.inp path=%s protname=%s I=%s pathpara=%s> %s/setup_%s.out' % (pathcharmm, pathinput, pathoutput, pdb, countp, pathpara, pathoutput, pdb + '_ch' + str(countp)))
                countp += 1
            os.system('%scharmm <%s/append.inp path=%s protname=%s NumP=%s pathpara=%s> %s/append_%s.out' % (pathcharmm, pathinput, pathoutput, pdb, NumP, pathpara, pathoutput, pdb))
            os.system('%scharmm <%s/charmm_to_namd.inp path=%s protname=%s pathpara=%s> %s/charmm_to_namd_%s.out' % (pathcharmm, pathinput, pathoutput, pdb, pathpara, pathoutput, pdb))
            pdball.append(pdb)
        else:
            continue
    f.close()


# get interface for wild type. interface control.
def interface():
    pdball = []
    f = open(in_file + ".cleaned", 'r')
    _unused = f.next()
    for line in f:
        ff = line.split("\t")
        pdb = ff[7].lower()
        partner1 = ff[8]
        partner2 = ff[9]
        NumP1 = int(len(partner1))
        NumP2 = int(len(partner2))
        if pdb not in pdball:
            p1 = 'ch1'
            p2 = 'ch' + str(NumP1 + 1)
            for countp1 in range(2, NumP1 + 1):
                p1 += ' .or. segid ch' + str(countp1)
            for countp2 in range(NumP1 + 2, NumP1 + NumP2 + 1):
                p2 += ' .or. segid ch' + str(countp2)
            os.system('sed -e "s/p1/%s/g" -e "s/p2/%s/g" %s/interface.inp > %s/interface_%s.inp' % (p1, p2, pathinput, pathoutput, pdb))
            os.system('%scharmm <%s/interface_%s.inp path=%s protname=%s pathpara=%s> %s/interface_%s.out' % (pathcharmm, pathoutput, pdb, pathoutput, pdb, pathpara, pathoutput, pdb))
            os.system('grep \'DELTAMS <\' %s/interface_%s.out | cut -c 24-100 > %s/dsasa_%s.out' % (pathoutput, pdb, pathoutput, pdb))
            fdsasa = open(pathoutput + '/dsasa_' + pdb + '.out', 'r')
            global dsasa
            dsasa = float(fdsasa.readline()[1:-2])
            pdball.append(pdb)
    f.close()


# produce input files for FoldX buildmodel and energy of the molecule (stability).
def inputfoldx():
    """
	command=BuildModel
	pdb=4P23_test7.pdb
	mutant-file=individual_list_4P23_B_KB5A.txt
 	"""
    f = open(in_file + ".cleaned", 'r')
    _unused = f.next()
    for line in f:
        ff = line.split("\t")
        pdb = ff[7]
        partner1 = ff[8]
        partner2 = ff[9]
        mut = ff[10][:-1]
        fi = open('individual_list_' + jobid + '_' + mut + '.txt', 'w')
        fi.write('%s' % (mut + ';'))
        fi.close()
        fs = open('foldx_buildmodel_' + jobid + '_' + mut + '.txt', 'w')
        fs.write('command=BuildModel\npdb=%s\nmutant-file=%s' % (jobid + '_' + mut + '.pdb', 'individual_list_' + jobid + '_' + mut + '.txt'))
        fs.close()
        os.system('cp %s/%s.pdb %s_%s.pdb' % (pathoutput, pdb, jobid, mut))
    f.close()


# submit foldx jobs for mutations, command: ./foldx -f foldx_buildmodel_4P23_B_KB5A.txt
def runfoldx_mut():
    f = open(in_file + ".cleaned", 'r')
    _unused = f.next()
    for line in f:
        ff = line.split("\t")
        pdb = ff[7]
        mut = ff[10][:-1]
        os.system('./foldx -f foldx_buildmodel_%s_%s.txt' % (jobid, mut))
        # delete files
        commands.getoutput("rm %s/%s_%s_temp.pdb" % (pathoutput, jobid, mut))
        commands.getoutput("rm %s_%s.pdb" % (jobid, mut))
        commands.getoutput("rm WT_%s_%s_1.pdb" % (jobid, mut))
        commands.getoutput("rm individual_list_%s_%s.txt" % (jobid, mut))
        commands.getoutput("rm foldx_buildmodel_%s_%s.txt" % (jobid, mut))
        commands.getoutput("rm Average_%s_%s.fxout" % (jobid, mut))
        commands.getoutput("rm Raw_%s_%s.fxout" % (jobid, mut))
        commands.getoutput("rm PdbList_%s_%s.fxout" % (jobid, mut))
        commands.getoutput("rm %s_%s.fxout" % (jobid, mut))
        commands.getoutput("mv Dif_%s_%s.fxout %s" % (jobid, mut, pathoutput))
        commands.getoutput("mv %s_%s_1.pdb %s/%s_%s.pdb" % (jobid, mut, pathoutput, pdb, mut))
    f.close()


def runfoldx_mut_threads(i):
    pdb = i.split('_')[0]
    mut = i.split('_')[1]
    try:
        os.system('./foldx -f foldx_buildmodel_%s_%s.txt' % (jobid, mut))
        # delete files
        commands.getoutput("rm %s/%s_%s_temp.pdb" % (pathoutput, jobid, mut))
        commands.getoutput("rm %s_%s.pdb" % (jobid, mut))
        commands.getoutput("rm WT_%s_%s_1.pdb" % (jobid, mut))
        commands.getoutput("rm individual_list_%s_%s.txt" % (jobid, mut))
        commands.getoutput("rm foldx_buildmodel_%s_%s.txt" % (jobid, mut))
        commands.getoutput("rm Average_%s_%s.fxout" % (jobid, mut))
        commands.getoutput("rm Raw_%s_%s.fxout" % (jobid, mut))
        commands.getoutput("rm PdbList_%s_%s.fxout" % (jobid, mut))
        commands.getoutput("rm %s_%s.fxout" % (jobid, mut))
        commands.getoutput("mv Dif_%s_%s.fxout %s" % (jobid, mut, pathoutput))
        commands.getoutput("mv %s_%s_1.pdb %s/%s_%s.pdb" % (jobid, mut, pathoutput, pdb, mut))

    except Exception as e:
        print e.args[0], e.args[1]


# split chains and produce pdb files for each chain of mutant pdb, which are used for VMD.
# produce 2FTL_KI15G_CH1.pdb and 2FTL_KI15G_CH2.pdb by 2FTL_KI15G.pdb
def splitchain_mut():
    f = open(in_file + ".cleaned", 'r')
    _unused = f.next()
    for line in f:
        ff = line.split("\t")
        pdb = ff[7]
        partner1 = ff[8]
        partner2 = ff[9]
        mut = ff[10][:-1]
        pdb = pdb + '_' + mut
        countchain = 1
        for chains in (list(partner1) + list(partner2)):
            os.system('grep "^.\{21\}%s" %s/%s.pdb > %s/%s_%s.pdb' % (chains, pathoutput, pdb, pathoutput, pdb, 'CH' + str(countchain)))
            countchain += 1
    f.close()


# produce psf and pdb files with vmd. must do this for charmm input, Changing the file of vmd.pgn's path as you need and force field.
def vmd_mut():
    f = open(in_file + ".cleaned", 'r')
    _unused = f.next()
    template = open(pathpara + '/vmd.pgn').read()
    for line in f:
        ff = line.split("\t")
        pdb = ff[7]
        partner1 = ff[8]
        partner2 = ff[9]
        mut = ff[10][:-1]
        protname = pdb + '_' + mut
        NumChain = int(len(partner1 + partner2))
        vmd_pdb = template.replace('protname', protname).replace('NumChain', str(NumChain)).replace('pathinput', pathinput).replace('pathoutput', pathoutput)
        with open(pathoutput + '/vmd_' + protname + '.pgn', 'w') as fw:
            fw.write(vmd_pdb)
        commands.getoutput('%svmd -dispdev text -e %s/vmd_%s.pgn' % (pathvmd, pathoutput, protname))
    f.close()


# produce charmm psf and crd files for individual chains
def charmmfile_mut():
    f = open(in_file + ".cleaned", 'r')
    _unused = f.next()
    for line in f:
        ff = line.split("\t")
        pdb = ff[7].lower()
        partner1 = ff[8]
        partner2 = ff[9]
        mut = ff[10][:-1].lower()
        protname = pdb + '_' + mut
        NumP = int(len(partner1 + partner2))
        countp = 1
        for chains in (list(partner1) + list(partner2)):
            os.system('grep "^.\{21\}%s" %s/%s_vmd.pdb > %s/%s_%s.pdb' % (chains, pathoutput, protname.upper(), pathoutput, protname, 'ch' + str(countp)))
            ffpdb = open(pathoutput + '/' + protname + '_ch' + str(countp) + '.pdb', 'a')
            ffpdb.write("%s" % ('END'))
            ffpdb.close()
            os.system('%scharmm <%s/setup.inp path=%s protname=%s I=%s pathpara=%s> %s/setup_%s.out' % (pathcharmm, pathinput, pathoutput, protname, countp, pathpara, pathoutput, protname + '_ch' + str(countp)))
            countp += 1
        os.system('%scharmm <%s/append.inp path=%s protname=%s NumP=%s pathpara=%s> %s/append_%s.out' % (pathcharmm, pathinput, pathoutput, protname, NumP, pathpara, pathoutput, protname))
        os.system('%scharmm <%s/charmm_to_namd.inp path=%s protname=%s pathpara=%s> %s/charmm_to_namd_%s.out' % (pathcharmm, pathinput, pathoutput, protname, pathpara, pathoutput, protname))
    f.close()


# run rest.png (change directory) for preparing the restrain pdb file for minimization of NAMD running.
def rest():
    pdball = []
    f = open(in_file + ".cleaned", 'r')
    _unused = f.next()
    template = open(pathpara + '/rest.pgn').read()
    for line in f:
        ff = line.split("\t")
        pdb = ff[7].lower()
        mut = ff[10][:-1].lower()
        protname = pdb + '_' + mut
        rest_mut = template.replace('protname', protname).replace('pathoutput', pathoutput)
        with open(pathoutput + '/rest_' + protname + '.pgn', 'w') as fw:
            fw.write(rest_mut)
        commands.getoutput('%svmd -dispdev text -e %s/rest_%s.pgn' % (pathvmd, pathoutput, protname))
        if pdb not in pdball:
            rest_wild = template.replace('protname', pdb).replace('pathoutput', pathoutput)
            with open(pathoutput + '/rest_' + pdb + '.pgn', 'w') as fw:
                fw.write(rest_wild)
            commands.getoutput('%svmd -dispdev text -e %s/rest_%s.pgn' % (pathvmd, pathoutput, pdb))
            pdball.append(pdb)
        else:
            continue
    f.close()


# run minimization.
def minimization():
    pdball = []
    runstep = "100"
    f = open(in_file + ".cleaned", 'r')
    _unused = f.next()
    template = open(pathpara + '/min.conf').read()
    for line in f:
        ff = line.split("\t")
        pdb = ff[7].lower()
        mut = ff[10][:-1].lower()
        protname = pdb + '_' + mut
        min_mut = template.replace('protname', protname).replace('pathoutput', pathoutput).replace('runstep', runstep).replace('pathinput', pathinput)
        with open(pathoutput + '/min_' + protname + '.conf', 'w') as fw:
            fw.write(min_mut)
        commands.getoutput('%snamd2 %s/min_%s.conf' % (pathnamd, pathoutput, protname))
        if pdb not in pdball:
            min_wild = template.replace('protname', pdb).replace('pathoutput', pathoutput).replace('runstep', runstep).replace('pathinput', pathinput)
            with open(pathoutput + '/min_' + pdb + '.conf', 'w') as fw:
                fw.write(min_wild)
            commands.getoutput('%snamd2 %s/min_%s.conf' % (pathnamd, pathoutput, pdb))
        else:
            continue
    f.close()


def frames():
    """
        get frame from minimization. 
        produce 2ftl_min.crd 
        """
    pdball = []
    startframe = 100
    f = open(in_file + ".cleaned", 'r')
    _unused = f.next()
    for line in f:
        ff = line.split("\t")
        pdb = ff[7].lower()
        mut = ff[10][:-1].lower()
        protname = pdb + '_' + mut
        os.system('%scharmm <%s/frame.inp path=%s protname=%s startframe=%s pathinput=%s> %s/frame_%s.out' % (pathcharmm, pathinput, pathoutput, protname, startframe, pathpara, pathoutput, protname))
        if pdb not in pdball:
            os.system('%scharmm <%s/frame.inp path=%s protname=%s startframe=%s pathinput=%s> %s/frame_%s.out' % (pathcharmm, pathinput, pathoutput, pdb, startframe, pathpara, pathoutput, pdb))
            pdball.append(pdb)
        else:
            continue
    f.close()


def energyfile():
    """
    produce energy inputfile
        """
    f = open(in_file + ".cleaned", 'r')
    _unused = f.next()
    template = open(pathpara + '/energy.inp').read()
    for line in f:
        ff = line.split("\t")
        pdb = ff[7].lower()
        mut = ff[10][:-1].lower()
        protname = pdb + '_' + mut
        partner1 = ff[8]
        partner2 = ff[9]
        NumP1 = int(len(partner1))
        NumP2 = int(len(partner2))
        p1 = 'ch1'
        p2 = 'ch' + str(NumP1 + 1)
        for countp1 in range(2, NumP1 + 1):
            p1 += ' -\n .or. segid ch' + str(countp1)
        for countp2 in range(NumP1 + 2, NumP1 + NumP2 + 1):
            p2 += ' -\n .or. segid ch' + str(countp2)
        if mut[1].upper() in partner1:
            mp = p1
            op = p2
            partner = 'p1'
        else:
            mp = p2
            op = p1
            partner = 'p2'
        partnerall = []
        # wild, energy_3se8_a.inp
        if partner not in partnerall:
            energy_wild = template.replace('pathinput', pathinput).replace('p1', mp).replace('p2', op)
            fwild = open(pathoutput + "/energy_" + pdb+'_'+partner +".inp", "w")
            fwild.write(energy_wild)
            partnerall.append(partner)
            fwild.close()
        # mut, energy_3se8_ra30a.inp
        energy_mut = template.replace('pathinput', pathinput).replace('p1', mp).replace('p2', op)
        fmut = open(pathoutput + "/energy_" + protname + ".inp", "w")
        fmut.write(energy_mut)
        fmut.close()
    f.close()


def energy():
    """
        energy calculation
        """
    partnerall = []
    f = open(in_file + ".cleaned", 'r')
    _unused = f.next()
    for line in f:
        ff = line.split("\t")
        pdb = ff[7].lower()
        mut = ff[10][:-1].lower()
        mutchain = mut[1]   # cleaned
        protname = pdb + '_' + mut
        partner1 = ff[8]
        if mutchain.upper() in partner1:
            partner = 'p1'
        else:
            partner = 'p2'
        # wild, energy_3se8_a.out
        if partner not in partnerall:
            os.system('%scharmm <%s/energy_%s.inp path=%s protname=%s > %s/energy_%s.out' % (pathcharmm, pathoutput, pdb+'_'+partner, pathoutput, pdb, pathoutput, pdb+'_'+partner))
            partnerall.append(partner)
        # mut, energy_3se8_ra30a.inp
        #os.system('%scharmm <%s/energy_%s.inp path=%s protname=%s > %s/energy_%s.out' % (pathcharmm, pathoutput, protname, pathoutput, protname, pathoutput, protname))
    f.close()


@retry(stop_max_attempt_number=3)
def energy_mut_threads(i):
    protname = i # pdb + '_' + mut
    #try:
    os.system('%scharmm <%s/energy_%s.inp path=%s protname=%s > %s/energy_%s.out' % (pathcharmm, pathoutput, i, pathoutput, i, pathoutput, i))        
    #except Exception as e:
    #    print e.args[0], e.args[1]


def RunProvean():
        """
        run provean.
        provean_2FTL_I_1.out 
        """
        pdbchain = []
        f = open (in_file+".cleaned",'r')
        _unused = f.next()
        for line in f:
            ff = line.split("\t")
            pdb = ff[7]
            chain = ff[3]
            if pdb+chain not in pdbchain:
                os.system('%sprovean.sh -q %s/%s_%s.seq -v %s/%s_%s.var > %s/provean_%s_%s.out --num_threads 20' % (pathprovean, pathoutput,pdb,chain,pathoutput,pdb,chain,pathoutput,pdb,chain))
            pdbchain.append(pdb+chain)
        f.close()

def DSSP():
    """
        produce dssp files showing secondary structure infromation for crystal structure of wild-type PDB.
        produce 2FTL.dssp, 2FTL_P1.dssp and 2FTL_P1.dssp by 2FTL.pdb from last step
        """
    pdball = []
    f = open(in_file + ".cleaned", 'r')
    _unused = f.next()
    for line in f:
        ff = line.split("\t")
        pdb = ff[7]
        partner1 = list(ff[8])  # [A, B]
        partner2 = list(ff[9])  # [C]

        # wild
        if pdb not in pdball:
            countp1 = 1
            countp2 = 1
            # COMPLEX
            os.system('%smkdssp %s/%s.pdb %s/%s.dssp' % (pathdssp, pathoutput, pdb, pathoutput, pdb))
            # P1
            for chains in (partner1):
                os.system('grep "^.\{21\}%s" %s/%s.pdb > %s/%s_%s_%s.pdb' % (chains, pathoutput, pdb, pathoutput, pdb, 'P1', 'C' + str(countp1)))
                countp1 += 1
            # P2
            for chains in (partner2):
                os.system('grep "^.\{21\}%s" %s/%s.pdb > %s/%s_%s_%s.pdb' % (chains, pathoutput, pdb, pathoutput, pdb, 'P2', 'C' + str(countp2)))
                countp2 += 1
            commands.getoutput('cat %s/%s_P1_C*.pdb > %s/%s_P1.pdb' % (pathoutput, pdb, pathoutput, pdb))
            commands.getoutput('cat %s/%s_P2_C*.pdb > %s/%s_P2.pdb' % (pathoutput, pdb, pathoutput, pdb))
            os.system('%smkdssp %s/%s_P1.pdb %s/%s_P1.dssp' % (pathdssp, pathoutput, pdb, pathoutput, pdb))
            os.system('%smkdssp %s/%s_P2.pdb %s/%s_P2.dssp' % (pathdssp, pathoutput, pdb, pathoutput, pdb))
            os.system('rm %s/%s_P[1-9]_C[1-9].pdb' % (pathoutput, pdb))   
            os.system('rm %s/%s_P[1-9].pdb' % (pathoutput, pdb))
            pdball.append(pdb)


    f.close()


# calculate sasa for each residue
def surfacefile():
    pdball = []
    f = open(in_file + ".cleaned", 'r')
    f.next()
    template = open(pathpara + '/surface.inp').read()
    for line in f:
        ff = line.split("\t")
        pdb = ff[7].lower()
        mut = ff[10][:-1].lower()
        protname = pdb + '_' + mut
        partner1 = ff[8]
        partner2 = ff[9]
        NumP1 = int(len(partner1))
        NumP2 = int(len(partner2))
        p1 = 'ch1'
        p2 = 'ch' + str(NumP1 + 1)
        for countp1 in range(2, NumP1 + 1):
            p1 += ' -\n .or. segid ch' + str(countp1)
        for countp2 in range(NumP1 + 2, NumP1 + NumP2 + 1):
            p2 += ' -\n .or. segid ch' + str(countp2)

        # wild
        if pdb not in pdball:
            fwt = open(pathoutput + "/surface_" + pdb + ".inp", "w")
            fwt_surface = template.replace('p1', p1).replace('p2', p2)
            fwt.write(fwt_surface)
            fwt.close()
        # mutant
        fmut = open(pathoutput + "/surface_" + protname + ".inp", "w")
        fmut_surface = template.replace('p1', p1).replace('p2', p2)
        fmut.write(fmut_surface)
        fmut.close()
		
    f.close()


def surface():
    pdball = []
    f = open(in_file + ".cleaned", 'r')
    _unused = f.next()
    for line in f:
        ff = line.split("\t")
        pdb = ff[7].lower()
        mut = ff[10][:-1].lower()
        protname = pdb + '_' + mut
        chain = ff[10][1:2]
        resnum = ff[10][2:-2]
        partner1 = ff[8]
        partner2 = ff[9]
        mapchainarray = []
        count = 1
        for chains in (list(partner1) + list(partner2)):
            cc = (chains, 'ch' + str(count))
            mapchainarray.append(cc)
            count += 1
        mapchaindict = dict(iter(mapchainarray))

        # mutant
        if chain in list(partner1):
            os.system('%scharmm <%s/surface_%s.inp pathinput=%s path=%s protname=%s chainm=%s resnumber=%s partner=%s> %s/surface_%s_%s_%s.out' % (pathcharmm, pathoutput, protname, pathinput, pathoutput, protname, mapchaindict[chain], resnum, "a", pathoutput, protname, chain.lower(), resnum))
        else:
            os.system('%scharmm <%s/surface_%s.inp pathinput=%s path=%s protname=%s chainm=%s resnumber=%s partner=%s> %s/surface_%s_%s_%s.out' % (pathcharmm, pathoutput, protname, pathinput, pathoutput, protname, mapchaindict[chain], resnum, "b", pathoutput, protname, chain.lower(), resnum))
        # wild
        if pdb not in pdball:
            if chain in list(partner1):
                os.system('%scharmm <%s/surface_%s.inp pathinput=%s path=%s protname=%s chainm=%s resnumber=%s partner=%s> %s/surface_%s.out' % (pathcharmm, pathoutput, pdb, pathinput, pathoutput, pdb, mapchaindict[chain], resnum, "a", pathoutput, pdb))
            else:
                os.system('%scharmm <%s/surface_%s.inp pathinput=%s path=%s protname=%s chainm=%s resnumber=%s partner=%s> %s/surface_%s.out' % (pathcharmm, pathoutput, pdb, pathinput, pathoutput, pdb, mapchaindict[chain], resnum, "b", pathoutput, pdb))

        pdball.append(pdb)
    f.close()


# calculate contact residues within cutoff with charmm
def distance_contact_resid():
    with open(in_file + ".cleaned", 'r') as f:
        f.next()
        for line in f:
            ff = line.split("\t")
            pdb = ff[7].lower()  # 1aay
            mut = ff[10][:-1].lower()  # ra16a
            resnum = ff[10][2:-2]  # 16
            Mutation_cleaned = ff[10]  # RA16A
            protname = pdb + '_' + mut  # 1aay_ra16a
            wild_resid = protname + '_wt'  # 1aay_ra16a_wt
            partner1 = ff[8]  # cleaned
            partner2 = ff[9]
            NumP1 = int(len(partner1))
            NumP2 = int(len(partner2))
            for chain in list(partner1 + partner2):
                if chain == Mutation_cleaned[1]:
                    chmut = 'ch' + str(list(partner1 + partner2).index(chain) + 1)
            p1 = ['ch{}'.format(i) for i in range(1, NumP1 + 1)]
            p1 = ' -\n .or. segid '.join(p1)
            p2 = ['ch{}'.format(i) for i in range(NumP1 + 1, NumP1 + NumP2 + 1)]
            p2 = ' -\n .or. segid '.join(p2)
            # mut partner:p1 or p2
            if mut[1].upper() in partner1:
                mp = p1
                op = p2
            else:
                mp = p2
                op = p1

            cuthb = '10' 
            # wild type
            with open('{}/pdb_distance_withH.inp'.format(pathinput), 'r') as f1:
                template = f1.read()
            result = template.replace('p1', mp).replace('p2', op).replace('chmut', chmut).replace('resnum',resnum).replace('cutoff', cuthb)
            with open('{}/distance_{}_{}.inp'.format(pathoutput, wild_resid, cuthb), 'w') as fwild:
                fwild.write(result)
            os.system('%scharmm <%s/distance_%s_%s.inp path=%s protname=%s pathpara=%s > %s/distance_%s_%s.out' % (pathcharmm, pathoutput, wild_resid, cuthb, pathoutput, pdb, pathpara, pathoutput, wild_resid, cuthb))


# ave(distance), sum(vdw), sum(elec)
def distance_contact_extract():
    cutoff = '10'
    with open(in_file + ".cleaned", 'r') as f1:
        f1.next()
        for line1 in f1:
            ff = line1.split("\t")
            pdb = ff[7].lower()  # 1aay
            wt_resid = ff[10][:-1].lower() + '_wt'  # ra22a_wt
            count = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))
            # wild_type
            with open(pathoutput + '/' + 'distance_' + pdb + '_' + wt_resid + '_' + cutoff + '.out', 'r') as f:
                for line in reversed(f.readlines()):
                    if re.match(r'\s+\d+ CH', line) or re.match(r'\d+ CH', line):
                        if (line[21:23] == 'C ' and line[46:48] == 'N ') or (line[21:23] == 'N ' and line[46:48] == 'C '):
                            continue
                        else:
                            data = re.split('\s+', line)
                            data_p1 = [line[6:9], line[11:14], line[14:21]]
                            data_p2 = [line[31:34], line[36:39], line[39:46]]
                            vdw = float(data[-3])
                            count[tuple(data_p1)][tuple(data_p2)]['vdw'].append(vdw)

                    elif line.startswith('            DISTANCES FOR SELECTED ATOMS'):
                        with open(pathoutput + '/distance_'+pdb+'_p1_p2_'+wt_resid+'_withoutH_' + cutoff + '.txt', 'w') as fw:
                            for key1, value1 in count.iteritems():
                                for key2, value2 in value1.iteritems():
                                    vdw_sum = sum(value2['vdw'])
                                    fw.write('{}\t{}\t{}\n'.format('\t'.join(key1), '\t'.join(key2),vdw_sum))
                        break


# select cutoff 6/10, calculate the length of mutpartner in contact with otherpartner .
# read distance_1aay_p1_p2_ra22a_wt/mut_6/10.txt, produce distance_1aay_len_mp_contact_op_6/10.txt
def distance_protein_contact_protein():
    cutoff = '10'
    out = pathoutput + '/distance_' + jobid + '_len_mp_contact_op_' + cutoff + '.txt'
    title = ['PDB_ID', 'mutation_cleaned', 'p1_p2_wt']
    with open(in_file + ".cleaned", 'r') as f, open(out, 'w') as fw:
        fw.write('%s\n' % '\t'.join(title))
        f.next()
        for line in f:
            ff = line.split('\t')
            pdb = ff[7].lower()  # 1aay
            mutation_cleaned = ff[10][:-1].lower()  # ra22a
            wt_resid = mutation_cleaned + '_wt'  # ra22a_wt
            in_wt = pathoutput + '/distance_' + pdb + '_p1_p2_' + wt_resid + '_withoutH_' + cutoff + '.txt'
            # wild type
            if os.path.isfile(in_wt):
                with open(in_wt, 'r') as fwt:
                    count_wt = set()
                    for ffwt in fwt:
                        ffwt_list = re.split(r'\s+', ffwt)
                        chain = ffwt_list[0]
                        resid = ffwt_list[1]
                        location = ffwt_list[2]
                        count_wt.add((chain, resid, location))
                len_wt = len(count_wt)
            else:
                len_wt = '0'
            fw.write('%s\t%s\t%s\n' % (pdb, mutation_cleaned, len_wt))


# get all features to 1aay.input.cleaned.outdata
def getenergy():

    fw = open(in_file + ".cleaned.outdata", 'w')
    first_line = file(in_file + ".cleaned", 'r').readlines()[0][:-1]
    fw.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' %
             (first_line, 'ddVdw', 'ddPb', 'provean','Foldx_fold', 'ACC_wt','ACCmp_wt','len_mp_contact_op_wt10','interface'))

    f = open(in_file + ".cleaned", 'r')
    f.next()
    for line in f:
        ff = line.strip().split("\t")
        pdb = ff[7].lower()  # 1aay
        partner1 = ff[8]
        partner2 = ff[9]
        mut = ff[10].lower()  # ra16a
        mutchain = mut[1]  # a
        protname = pdb + '_' + mut  # 1aay_ra16a
        count = 1
        mapchains = defaultdict(str)
        for i in list(partner1+partner2):
            mapchains[i] = 'CH'+str(count)
            count+=1

        if mutchain.upper() in partner1:
            partner = 'p1'
        else:
            partner = 'p2'

        fdebug = [i.strip() for i in open(jobpath+'/debug.txt').readlines()]
        if 'Operating system error: Cannot allocate memory' in fdebug:
            vdw_mut,pb_mut,vdw_wild,pb_wild,ddvdw,ddpb='None','None','None','None','None','None'
        else:
            # add vdw ,pb , mmpbsa, gasa, ms msac, msa, msc of mutant type
            vdw_mut = float(os.popen('grep \'VDWENER <\' %s/energy_%s.out | cut -c 24-100' % (pathoutput, protname)).read().strip()[1:-1])
            pb_mut = float(os.popen('grep \'DELTAGELEC <\' %s/energy_%s.out | cut -c 27-100' % (pathoutput, protname)).read().strip()[1:-1])
            # add vdw ,pb , mmpbsa, gasa, ms energy of wild type
            vdw_wild = float(os.popen('grep \'VDWENER <\' %s/energy_%s.out | cut -c 24-100' % (pathoutput, pdb+'_'+partner)).read().strip()[1:-1])
            pb_wild = float(os.popen('grep \'DELTAGELEC <\' %s/energy_%s.out | cut -c 27-100' % (pathoutput, pdb+'_'+partner)).read().strip()[1:-1])
            # diff
            ddvdw = vdw_mut-vdw_wild
            ddpb = pb_mut-pb_wild

        # add foldx folding free energies
        fmut = open(pathoutput + '/Dif_' + jobid + '_' + mut.upper() + '.fxout', 'r')
        FoldX_fold = fmut.readlines()[9].split("\t")[1]

        # add provean score
        ST_play = False
        fp = open (pathoutput+"/provean_"+pdb.upper()+"_"+ff[3]+".out","r")
        for provean in fp:
            ffp = provean.split("\t")
            if provean[2:11] == "VARIATION":
                ST_play = True
                continue
            if ST_play:
                if ffp[0] == (mut[0]+mut[2:]).upper():
                        provean_score = ffp[1][:-1]

        # add location
        resnum = mut[2:-1]
        ## mut
        drasap_mut = float(os.popen('grep \' DRASAP <\' %s/surface_%s_%s_%s.out | cut -c 23-100' % (pathoutput, protname, mutchain, resnum)).read().strip()[1:-1])
        ## if mutation on interface
        locationx = ''
        if drasap_mut > 0.0:
            locationx = "yes"
        if drasap_mut == 0.0:
            locationx = "no"

        # add ACC_wild
        ST_play = False
        fd = open(pathoutput + "/" + pdb.upper() + ".dssp", "r")
        for dssp in fd:
            if dssp[2:3] == "#":
                ST_play = True
                continue
            if ST_play:
                resnum = dssp[5:10].strip()
                chain = dssp[11:12]
                resname = dssp[13:14]
                ACC = dssp[35:38]
                if (resname + chain + resnum).lower() == mut[0:-1]:
                    ACCi_wild = ACC
                    break
        fd.close()

        # add ACCmp_wild
        chainm = mut[1:2].upper()
        if chainm in list(partner1):
            ST_play = False
            fd = open(pathoutput + "/" + pdb.upper() + "_P1.dssp", "r")
            for dssp in fd:
                if dssp[2:3] == "#":
                    ST_play = True
                    continue
                if ST_play:
                    resnum = dssp[5:10].strip()
                    chain = dssp[11:12]
                    resname = dssp[13:14]
                    ACCmp = dssp[35:38]
                    if (resname + chain + resnum).lower() == mut[0:-1]:
                        ACCmpi_wild = ACCmp
                        break
            fd.close()
        if chainm in list(partner2):
            ST_play = False
            fd = open(pathoutput + "/" + pdb.upper() + "_P2.dssp", "r")
            for dssp in fd:
                if dssp[2:3] == "#":
                    ST_play = True
                    continue
                if ST_play:
                    resnum = dssp[5:10].strip()
                    chain = dssp[11:12]
                    resname = dssp[13:14]
                    ACCmp = dssp[35:38]
                    if (resname + chain + resnum).lower() == mut[0:-1]:
                        ACCmpi_wild = ACCmp
                        break
            fd.close()

        # add distance_1bpx_len_p2_contact_p1_10.txt
        with open(pathoutput + '/distance_' + jobid + '_len_mp_contact_op_10.txt', 'r') as fp1p2_10:
            fp1p2_10.next()
            dp1p2_10 = defaultdict(lambda: defaultdict(str))
            for ffp1p2_10 in fp1p2_10:
                ffp1p2_10_list = ffp1p2_10.strip().split('\t')
                pdbid = ffp1p2_10_list[0]
                mutation = ffp1p2_10_list[1]
                len_p1_p2_10 = ffp1p2_10_list[2]
                dp1p2_10[pdbid][mutation] = len_p1_p2_10


        fw.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % 
                 (line.strip(), ddvdw,ddpb,provean_score,FoldX_fold,ACCi_wild,ACCmpi_wild,dp1p2_10[pdb][mut],locationx))

    f.close()
    fw.close()


# using our fitting models.
def Prediction():
    outdata = in_file + '.cleaned.outdata'
    robjects.globalenv["outdata"] = outdata
    robjects.globalenv["workdir"] = workdir
    r('''test = read.table(outdata,header=T,sep="\t")''')
    r('''filename_rf = paste(workdir, 'sinputfiles/mutabinds.RData',sep = '')''')
    r('''load(file = filename_rf)''')
    PredR = r('''predict(mutabinds.rf,test)''')
    robjects.globalenv["PredR"] = PredR
    SunDDG = r('''PredR''')

    first_line = file(in_file, 'r').readlines()[0][:-1]
    fw = open(out_file, "w")
    fw.write("%s\t%s\t%s\t%s\t%s\n" % (first_line, "SunDDG", "Interface", "Deleterious", "Confidence"))

    f = open(in_file + '.cleaned.outdata', 'r')
    _unused = f.next()
    count = 0
    for line in f:
        ff = line.strip().split("\t")
        if ((float(SunDDG[count])) >= 1.5) or ((float(SunDDG[count])) <= (-1.5)):
            Deleterious = "yes"
            Confidence = "High"
        else:
            Deleterious = "no"
            Confidence = "High"
        if ff[-1] == "no":  # for all non-interfacial mutations, the confidence is low.
            Confidence = "Low"
        fw.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%3.2f\t%s\t%s\t%s\n" % (
        ff[0], ff[1], ff[2], ff[3], ff[4], ff[5], ff[6], SunDDG[count], ff[-1], Deleterious, Confidence))
        count += 1
    f.close()
    fw.close()


# Add RF contribution
def rf_contribution():
    outdata = in_file + '.cleaned.outdata'
    template = open(pathpara+'/MutaBindS_forestFloor.R').read()
    result = template.replace('test_outfeature',outdata).replace('test_outcontribution',outdata+'.contribution').replace('test_sunddg',out_file)
    with open(pathoutput+'/'+jobid+'_forestFloor.R','w') as fw:
        fw.write(result)
    os.system('%sRscript %s' % (pathrscript, pathoutput+'/'+jobid+'_forestFloor.R'))
    finalout = open(outdata+'.contribution').read().replace('"','')
    with open(out_file,'w') as fw:
        fw.write(finalout)


# change mutant pdb file format and get output.
def mutpdb1():
    pdball = []
    f = open(in_file + ".cleaned", 'r')
    _unused = f.next()
    for line in f:
        ff = line.split("\t")
        pdb = ff[7].lower()
        mut = ff[10][:-1].lower()
        protname = pdb + '_' + mut
        Partner1 = ff[1]
        Partner2 = ff[2]
        NewPartner1 = ff[8]
        NewPartner2 = ff[9]
        NumP1 = int(len(NewPartner1))
        NumP2 = int(len(NewPartner2))
        if pdb not in pdball:
            mapchainarray = []
            for countp1 in range(1, NumP1 + 1):
                cc = ('CH' + str(countp1), 'P1')
                mapchainarray.append(cc)
                mapchaindict = dict(iter(mapchainarray))
            for countp2 in range(NumP1 + 1, NumP1 + NumP2 + 1):
                cc = ('CH' + str(countp2), 'P2')
                mapchainarray.append(cc)
                mapchaindict = dict(iter(mapchainarray))

        if pdb not in pdball:
            chainresnumall = []
            for chains in (Partner1.split(".") + Partner2.split(".")):
                fp = open(pathoutput + "/" + pdb.upper() + '_' + chains + ".pdb", "r")
                for linep in fp:
                    chainresnum = linep[21:27] + '.' + linep[72:-1]
                    if chainresnum not in chainresnumall:
                        chainresnumall.append(chainresnum)
                    else:
                        continue
            pdball.append(pdb)
        os.system('sed -e "s/HSD/HIS/g" %s/%s_min.pdb > %s/%s_min_temp.pdb' % (pathoutput, protname, pathoutput, protname))
        fw = open(pathoutput + "/" + pdb + "_" + ff[5] + "_min.pdb", "w")
        chainresnumml = []
        count = 0
        fm = open(pathoutput + "/" + protname + "_min_temp.pdb", "r")
        for linem in fm.readlines()[2:-2]:
            chainresnumm = linem[17:26]
            if chainresnumm != chainresnumml:
                count += 1
                fw.write("%s%s%s%s %s\n" % (
                linem[0:21], chainresnumall[count - 1].split(".")[0], linem[27:72], mapchaindict[linem[72:-2]], chainresnumall[count - 1].split(".")[1][2:]))
            else:
                fw.write("%s%s%s%s %s\n" % (
                linem[0:21], chainresnumall[count - 1].split(".")[0], linem[27:72], mapchaindict[linem[72:-2]], chainresnumall[count - 1].split(".")[1][2:]))
            chainresnumml = linem[17:26]
    f.close()


# Put "TER" at the end of each chain.
def mutpdb2():
    pdball = []
    f = open(in_file + ".cleaned", 'r')
    _unused = f.next()
    for line in f:
        ff = line.split("\t")
        pdb = ff[7]
        if pdb not in pdball:
            fw = open(jobpath + "/" + pdb + "_" + ff[5] + "_min.pdb", "w")
            fw.write("REMARK  THIS MUTANT(" + pdb + '-' + ff[3] + '-' + ff[4] + ") STRUCTURE IS PRODUCED BY 100-STEP ENERGY MINIMIZATION\n")
            fw.write("REMARK   DATE:" + "  " + time.strftime("%x") + "  " + time.strftime("%X") + "   CREATED BY SERVER: MUTABIND2\n")
            fw.write("REMARK   REFERENCE: PLEASE CITE *** \n")
            fm = open(pathoutput + "/" + pdb.lower() + "_" + ff[5] + "_min.pdb", "r")
            chainl = file(pathoutput + "/" + pdb.lower() + "_" + ff[5] + "_min.pdb", "r").readlines()[0][21:22] + \
                     file(pathoutput + "/" + pdb.lower() + "_" + ff[5] + "_min.pdb", "r").readlines()[0][75:-1]
            for linem in fm:
                chain = linem[21:22] + linem[75:-1]
                if chain != chainl:
                    fw.write("%s\n" % ("TER"))
                    fw.write("%s" % (linem))
                else:
                    fw.write("%s" % (linem))
                chainl = chain
            fw.write("END")
    f.close()
    fw.close()
    fm.close()


#def change_status(jobid, newstatus, progress=0):
#    con = None
#    try:
#        con = mdb.connect(host="10.20.212.172", user="web", passwd="Web2017mutabind-pniServer", db="mutabindDB")
#        cur = con.cursor()
#        cur.execute("UPDATE sunddg_job set status=%s, progress=%s where job LIKE %s ", (newstatus, progress, jobid))
#	con.commit()
#    except mdb.Error, e:
#        print("Error %d: %s" % (e.args[0], e.args[1]))
        # sys.exit(1)
        # Continue on error
#    finally:
#        if con:
#            con.close()


def main1():
    global dsasa
    global jobid
    global vmderror 

    ProPDB1()
    del_unknown_incomplete()
    splitchain()
    CleanPdb()
    wtpdb()
    vmd_wt()
    if vmderror:
        f = open(out_file + '.error2', 'w')
        f.write("%s" % "Your interaction parters contain more than 99999 atoms after adding hydrogen atoms by VMD program. We cannot calculate that :(.")
    else:
        charmmfile_wt()
        interface()
        if dsasa < -100:
            inputfoldx()
        else:
            f = open(out_file + '.error1', 'w')
            f.write("%s" % "Your interaction parters do not contact with each other.")

main1()


# 如果没有error1和error2文件,代表用户上传的相互作用蛋白之间有相互作用,且在进行VMD加氢时,总原子数小于99999.
# FoldX生成突变体结构
if os.path.isfile(out_file + '.error1') or os.path.isfile(out_file + '.error2'):
    pass
else:
    f = pd.read_table(in_file+'.cleaned',dtype=str,sep='\t',header=0,index_col=False)
    last_mut = list(f['Mutation_cleaned'])[-1]
    if os.path.isfile(workdir + jobid + '_'+last_mut+'.pdb'):
        f['protname'] = f['PDBid']+'_'+f['Mutation_cleaned']
        pool = Pool(processes=5)
        pool.map(runfoldx_mut_threads,list(f['protname']))


def main2():
    if os.path.isfile(out_file + '.error1') or os.path.isfile(out_file + '.error2'):
        pass
    else:
        splitchain_mut()
        vmd_mut()
        charmmfile_mut()
        rest()
        minimization()
        frames()
        energyfile()
        energy()
        RunProvean()
        DSSP()
        surfacefile()
        surface()
        distance_contact_resid()
        distance_contact_extract()
        distance_protein_contact_protein()

main2()


# 如果存在distance_*_len_mp_contact_op_10.txt文件,表示除了突变体energy外,其他features均已计算完成.
# CHARMM计算VMD和Pb
if os.path.isfile(pathoutput + '/distance_' + jobid + '_len_mp_contact_op_10.txt'):
    f = pd.read_table(in_file+'.cleaned',dtype=str,sep='\t',header=0,index_col=False)
    f['protname'] = f['PDBid'].str.lower()+'_'+f['Mutation_cleaned'].str.lower()
    pool = Pool(processes=5)
    pool.map(energy_mut_threads,list(f['protname']))


def main3():
    if os.path.isfile(out_file + '.error1') or os.path.isfile(out_file + '.error2') or os.path.isfile(out_file + '.error3'):
        pass
    else:
        getenergy()
        f = pd.read_csv(in_file + ".cleaned.outdata", header=0, index_col=None,sep='\t')
        if f[pd.isnull(f['ddPb'])].shape[0]>0:
            ferror = open(out_file + '.error3', 'w')
            ferror.write("%s" % "Your job encountered a problem when calculating PB. Please check your structure.")
        elif 'None' in f['ddPb'].tolist():
            ferror = open(out_file + '.error3', 'w')
            ferror.write("%s" % "Your job encountered a problem when calculating VDW and PB because of memory limit. We cannot calculate that :(.")
        else:
            Prediction()
            rf_contribution()
            mutpdb1()
            mutpdb2()

main3()
