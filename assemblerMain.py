import dataGen

import cluster
import graphForm
import branchClear

import bridgeResolve
import alignmentBridge
import eulerCycle
import readAns

import compare
import logging
import time


def runAssembler(snpRate, typeOfGen, detail,parameterRobot):
    #N, G, L, p,K,snpRate, typeOfGen, detail = 100, 100, 10000, 0.01,30,0.001, 'r', "1000"

    
    motherGen, reads, noisyReads = dataGen.generateData(typeOfGen, detail,parameterRobot)
    motherGen, reads, noisyReads = logging.rawDataLoad(parameterRobot.defaultFolder+"UnitTest",parameterRobot.G,parameterRobot.N,parameterRobot.L, "dn")

    f1 = cluster.groupIndelNoisyKmers(noisyReads,parameterRobot)
    
    G1,startList, f1 = graphForm.getSeqGraph(f1,noisyReads, parameterRobot)
    f2, G2 = branchClear.clearResidual(f1, G1,parameterRobot)
    
    G3 = bridgeResolve.resolveRepeats(f2, G2,parameterRobot)  

    G4 = alignmentBridge.MSAresolve(f2, G3, noisyReads, snpRate,parameterRobot )
    #G4 = G3 
    recovSeq = eulerCycle.findEC(G4)
    
    recovGen = readAns.reportRecovSeq(recovSeq, f2, noisyReads,parameterRobot)
    
    numMistakes, success = compare.subAlignCompare(recovGen, motherGen,parameterRobot)
    #numMistakes, success = 0 , 0 
    return numMistakes, success

# Target of the new code : Fast Assemble and assembly in optimal amount of information 
#t0 = time.time()
#N, G, L, p,snpRate, typeOfGen, detail = 1000, 10000,200, 0.015, 0.001 ,'m', "500-300-50" 
#runAssembler(N, G, L, p,snpRate, typeOfGen, detail)
#print "Time (sec) :", time.time() - t0