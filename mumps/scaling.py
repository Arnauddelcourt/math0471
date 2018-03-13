#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os, time

def runtest(ne, nproc, omp):


    if os.sep=='\\':
        exename = 'bin\\test_mumps.exe'
    else:
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
        "mumps.verb" : true,
        "save.vti" : false
    }
    """

    # write input file
    f = open('input.json','w')
    f.write(inpfile % (ne-1,ne-1,ne-1))
    f.close()

    # run test
    os.environ['OMP_NUM_THREADS'] = '%d' % omp
    cmd = 'mpiexec -np %d %s input.json >out.txt' % (nproc, exename)
    #print "running '%s'" % cmd
    #print "running ne=%d, nproc=%d, omp=%d" % (ne, nproc, omp)
    os.system(cmd)

if __name__=="__main__":
    ne = 10    
    nproc = 6
    omp = 2

    print "ne\tnproc\tomp\telapsed"
    for ne in range(10, 80, 10):
        start = time.time()
        runtest(ne, nproc, omp)
        elapsed = time.time()-start
        print ne, nproc, omp, elapsed
