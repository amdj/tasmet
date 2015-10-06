#!/usr/bin/python
# Author: J.A. de Jong

# This script can be used to install TaSMET on a user-basis (a local
# install). The script will create a directory $HOME/bin/TaSMET, copy/link the
# installation to this directory and finally update the PYTHONPATH variable in
# .bashrc.

# We assume the following is present:
# - A .bashrc file
# cmake and make
# A valid compiler, and all other prerequisites of TaSMET (see documentation).


import sys,os,subprocess as sp,multiprocessing as mp
# ####################  Configuration    

# Find home directory
home=os.getenv("HOME")

# Fill in the path to TaSMET
TaSMETpath=home+'/bin/TaSMET'

# For editing the code, it might me advisable to link the code, in stead of copying it
link=True

# First make the files
if sp.call(['cmake','.'])!=0:
    exit(1)

ncpus=mp.cpu_count()
if sp.call(['make','-j%g' %ncpus if (ncpus>0 and ncpus<9) else 1])!=0:
    exit(1)

# ####################  Install the files 

# Recursive way of creating a directory
def mkdir_recursive(path):
    sub_path = os.path.dirname(path)
    if not os.path.exists(sub_path):
        mkdir_recursive(sub_path)
    if not os.path.exists(path):
        os.mkdir(path)


# Create the directories
mkdir_recursive(TaSMETpath+'/common')
print('Paths created succesfully')

# Find current dir
curdir=os.getcwd()

swigfiles=[['common/_common.so','common/_common.so'],\
       ['common/common.py','common/common.py'],\
       ['src/_TaSMET.so','_TaSMET.so'],\
       ['src/TaSMET.py','TaSMET.py']]
if link is True:
    cparg='-sf'
else:
    cparg='p'
for swigfile in swigfiles:
    if sp.call(['cp',cparg,curdir+'/'+swigfile[0],TaSMETpath+'/'+swigfile[1]])!=0:
        print("Error while copying ", swigfile[0], " to its destination.")
# And copy the python dir contents
if sp.call(['cp',cparg,'-r',curdir+'/src/python/.',TaSMETpath])!=0:
        print("Error while copying python files to their destination.")    

# Add TaSMET to PYTHONPATH, if not already.
if not TaSMETpath in sys.path:
    print(TaSMETpath + ' not yet in PYTHONPATH. Appending it.')
    with open(home+'/.bashrc', 'a') as bashrc:
        bashrc.write('export PYTHONPATH=$PYTHONPATH:'+TaSMETpath+'\n')
    os.system('export PYTHONPATH=$PYTHONPATH:'+TaSMETpath)
    print('PYTHONPATH updated. Please close this session to reload the updated PYTHONPATH.')

print("Installation completed.")
