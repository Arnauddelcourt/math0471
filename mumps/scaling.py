#!/usr/bin/env python
# -*- coding: utf-8 -*-


import os


def runtest(ne, nproc, omp):

    exename = 'bin/test_mumps'

    if not os.path.isfile(exename):
        raise Exception("%s not found" % exename)

    inpfile = \
    """
    {
        "grid.o" : [10.0, 10.0, 10.0],
        "grid.L" : [10.0, 10.0, 10.0],
        "grid.np2" : [%d, %d, %d],
        "matlab" : false,
        "mumps.verb" : true
    }
    """



    #print inpfile % (np,np,np)

    f = open('input.json','w')
    f.write(inpfile % (ne-1,ne-1,ne-1))
    f.close()



    os.environ['OMP_NUM_THREADS'] = '%d' % omp
    cmd = 'mpiexec -np %d %s input.json >out.txt' % (nproc, exename)
    #print "running '%s'" % cmd
    print "running ne=%d, nproc=%d, omp=%d" % (ne, nproc, omp)
    os.system(cmd)

if __name__=="__main__":
    ne = 10    
    nproc = 2
    omp = 3
    runtest(ne, nproc, omp)
