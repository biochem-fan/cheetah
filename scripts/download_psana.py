#!/usr/bin/env python

import sys
import os

if(len(sys.argv) < 3):
    print "Usage: download_psana.py <ana-version> <output-directory>"
    print ""
    print "Example: download_psana.py ana-0.9.1 /opt/psana/g/psdm/portable/sw/releases/"
    print ""
    print "Will download the components of ana necessary to build cheetah into the output directory"
    sys.exit(0)

ana_version = sys.argv[1] 
output_dir = sys.argv[2]


                    

os.system('ssh psexport.slac.stanford.edu "cd /reg/g/psdm/sw/releases/ && tar czf - '+ana_version+'" | tar xzvf - -C '+output_dir)
os.system('ssh psexport.slac.stanford.edu "cd /reg/g/psdm/sw/external/pdsdata/6.1.4/x86_64-rhel6-gcc44-opt && tar czf - pdsdata" | tar xzvf - -C '+output_dir+'/'+ana_version+'/include)
os.system('ssh psexport.slac.stanford.edu "cd /reg/g/psdm/sw/external/root/5.30.06-python2.7/x86_64-rhel6-gcc44-opt/include && tar czf - root" | tar xzvf - -C '+output_dir+'/'+ana_version+'/include')

# Fix broken symlinks
indir = output_dir+'/'+ana_version+'/arch'
for root, dirs, filenames in os.walk(indir):
    for f in filenames:
        fpath = os.path.join(root, f)
        if(os.path.islink(fpath) and not os.path.exists(os.readlink(fpath))):
            orig_path = os.readlink(fpath)
            if(orig_path.find('external')):
                end_path = orig_path[orig_path.find('external/')+9:]
                new_path = '../../../../../external'+end_path
                if(os.path.exists(new_path)):
                    os.unlink(fpath)
                    os.symlink(new_path,fpath)
                else:
                    new_path = output_dir+"../../../../../common/package/"+end_path
                    if(os.path.exists(new_path)):
                        os.unlink(fpath)
                        os.symlink(new_path,fpath)