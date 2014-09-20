import dataGen
import cluster
import graphForm
import branchClear
import bridgeResolve
import alignmentBridge
import eulerCycle
import readAns
import compare

import assemblerMain

import logging
import numpy as np 
import random

import networkx as nx
import matplotlib.pyplot as plt 
import time
import cProfile
import threading
import multiprocessing as mp
from multiprocessing import Process
import multiprocessing
import ctypes

import os
#from alignment.sequence import Sequence
#from alignment.vocabulary import Vocabulary
#from alignment.sequencealigner import SimpleScoring, GlobalSequenceAligner
### Testing
def dataGenUnitTest():
    print "dataGenUnitTest"
    dummyParameters = logging.parameterObj()
    dummyParameters.N, dummyParameters.G, dummyParameters.L, dummyParameters.p, typeOfGen, detail = 100, 10000,100, 0.01, 't', "500"
    dummyParameters.indel = True
    
    motherGen, reads, noisyReads = dataGen.generateData(typeOfGen,detail,dummyParameters)
    
    
    
def dataGenUnitTest2():
    print "dataGenUnitTest2"
    dummyParameters = logging.parameterObj()
    dummyParameters.N, dummyParameters.G, dummyParameters.L, dummyParameters.p, typeOfGen, detail = 100, 10000,100, 0.01, 't', "500"    
    G,N,L = dummyParameters.G, dummyParameters.N, dummyParameters.L
    
    motherGen, reads, noisyReads = logging.rawDataLoad("UnitTest",G,N,L, 'd')
    
    print reads 
    print noisyReads
    
    baseList = [] 
    
    
    for indexN in range(N):
        for indexL in range(L):
            if reads[indexN][indexL] == noisyReads[indexN][indexL]:
                baseList.append(0)
            else:
                baseList.append(1)
                
    print "Error Probability", np.mean(baseList)
    print "numberOfBasesGet ", len(baseList)

def clusterUnitTest():
    

    
    typeLen = 30 
    readLen = 40
    copyLen = 10 
    p = 0.01 
 
 
    noisyReads  = np.zeros(typeLen*copyLen*2*readLen, dtype = np.int8).reshape(typeLen*copyLen,2*readLen)
    rawReads = np.zeros(typeLen*readLen, dtype = np.int8).reshape(typeLen, readLen)    
    
    for typeindex in range(typeLen):
        for eachbaseindex in range(readLen):
            rawReads[typeindex][eachbaseindex] = random.randint(1,4)
            
    
        
    
    for index in range(copyLen):
        startindex = index* typeLen
        endindex = (index +1 )*typeLen
        tempNoisyReads = dataGen.addIndelNoise(rawReads, p)
        for j in range(startindex, endindex):
            noisyReads[j][0:len(tempNoisyReads[j-startindex])] = tempNoisyReads[j-startindex][:]
        
    
    logging.rawDataSave("clusterReadsUnitTest", "", "", noisyReads, 'n')
    
def clusterUnitTest2():
   # Unit Test 1 : Generate 30 reads with 10 copies with noise , length being 20 
    
    motherGen, reads, noisyReads = logging.rawDataLoad("clusterReadsUnitTest",10,300,40,"dn")
    
    
    parameterRobot = logging.parameterObj()
    parameterRobot.N = 300 
    parameterRobot.L = 40
    parameterRobot.G = 10000
    parameterRobot.liid = 30
    parameterRobot.K = 30
    parameterRobot.threshold = 5
    parameterRobot.p = 0.01
    
    parameterRobot.indel = True

    
    cluster.groupIndelNoisyKmers(noisyReads, parameterRobot,  "fast")
    
             
            
def clusteringTestingOfK(G, N, L, folderName, K, threshold, liid):
    motherGen, reads, noisyReads = logging.rawDataLoad(folderName + "UnitTest",G,N,L,'a')
    dummyParameters = logging.parameterObj()
    dummyParameters.brachingDepth = 20
    dummyParameters.clusterRounds , dummyParameters.fingerPrint,  dummyParameters.clusterRatio = 2, 6, 2
    
    
    dummyParameters.N, dummyParameters.L, dummyParameters.K, dummyParameters.G, dummyParameters.threshold, dummyParameters.liid = N, L, K, G, threshold, liid
     
    cluster.groupNoisyKmers(noisyReads,dummyParameters, 'fast' )

def graphFormUnitTest(N, G, L, folderName =""):
    print "dataGenUnitTest"
    dummyParameters = logging.parameterObj()
    dummyParameters.defaultFolder = folderName
    dummyParameters.liid = 40
    dummyParameters.K = 40
    dummyParameters.threshold = 6
    dummyParameters.p = 0.015
    dummyParameters.indel = True


    dummyParameters.N, dummyParameters.G, dummyParameters.L, dummyParameters.p, typeOfGen, detail = N, G, L, 0.015, 'm', "500-200-50" 
    
   # motherGen, reads, noisyReads = dataGen.generateData( typeOfGen,detail,dummyParameters)    
    motherGen, reads, noisyReads = logging.rawDataLoad(dummyParameters.defaultFolder+"UnitTest",dummyParameters.G,dummyParameters.N,dummyParameters.L, "dn")
    
    returnfmapping= cluster.groupIndelNoisyKmers(noisyReads, dummyParameters,  "fast")

def graphFormUnitTest2(N, G, L, folderName = ""):   
    dummyParameters = logging.parameterObj()
    dummyParameters.defaultFolder = folderName
    returnfmapping= logging.fmappingLoad(dummyParameters.defaultFolder+'clusteredGroup.csv')

    dummyParameters.N, dummyParameters.G, dummyParameters.L, dummyParameters.p, typeOfGen, detail = N, G, L, 0.015, 'm', "500-200-50" 

    motherGen, reads, noisyReads = logging.rawDataLoad(dummyParameters.defaultFolder+"UnitTest",dummyParameters.G,dummyParameters.N,dummyParameters.L, "dn")
    G1,startList,fmapping = graphForm.getSeqGraph(returnfmapping,noisyReads, dummyParameters )
    #checkCondensingBasic2(G1,startList, "simple")
    #printGraphTrace(G1, G1[1016], 50)

def checkCondensingBasic2(G1,startList, typeOfGraph):
    print "Checking : ----- checkCondensingBasic(G1,startList)"
    print "len(G1), len(startList)", len(G1), len(startList)
    G = nx.MultiDiGraph()
    counter = 0 
    for eachnode in G1:
        if eachnode.nodeIndex == 8039:
            print "counter", counter
        else:
            counter += 1 
            
        if len(eachnode.listOfNextNodes) == 0:
            G.add_edges_from([(str(eachnode.nodeIndex)+" : "+str(len(eachnode.nodeIndexList)), str(eachnode.nodeIndex)+" : "+str(len(eachnode.nodeIndexList)))])

            
        for eachnextnode in eachnode.listOfNextNodes:
            if eachnode.nodeIndex == eachnextnode.nodeIndex:
                print "Self Loop", len(eachnode.nodeIndexList), eachnode.nodeIndex
            
            if typeOfGraph == "simple":
                G.add_edges_from([(str(eachnode.nodeIndex)+" : "+str(len(eachnode.nodeIndexList)), str(eachnextnode.nodeIndex)+" : "+ str(len(eachnextnode.nodeIndexList)))])
            elif typeOfGraph == "MB":
                G.add_edges_from([(str(eachnode.nodeIndex)+" : "+str(len(eachnode.nodeIndexList)), str(eachnextnode[0].nodeIndex)+" : "+ str(len(eachnextnode[0].nodeIndexList)))])
            
            
            for eachprevnode in eachnextnode.listOfPrevNodes:
                if eachprevnode.nodeIndex == eachnextnode.nodeIndex:
                    print "Self Loop", len(eachnode.nodeIndexList), eachnode.nodeIndex
                if typeOfGraph == "simple":
                    G.add_edges_from([(str(eachprevnode.nodeIndex)+" : "+str(len(eachprevnode.nodeIndexList)), str(eachnextnode.nodeIndex)+" : "+ str(len(eachnextnode.nodeIndexList)))])
                elif typeOfGraph == "MB":
                    G.add_edges_from([(str(eachprevnode.nodeIndex)+" : "+str(len(eachprevnode.nodeIndexList)), str(eachnextnode[0].nodeIndex)+" : "+ str(len(eachnextnode[0].nodeIndexList)))])
            
    plt.figure(figsize=(30,15))
    nx.draw(G,node_size=200,layout=nx.spring_layout(G))
    #nx.draw_circular(G)
    #nx.draw_spring(G) 
        
    plt.show()



def printGraphTrace(G1, startNode, Max):
    G = nx.MultiDiGraph()    
    
    for eachnode in G1:
        eachnode.visited = False
        
    stack = [startNode]
    counter = 0 
    
    while (len(stack)> 0 and counter < Max):
        currentNode = stack.pop(0)
        currentNode.visited = True
        
        for eachnextnode in currentNode.listOfNextNodes:
            G.add_edges_from([(str(currentNode.nodeIndex)+" : "+str(len(currentNode.nodeIndexList)), str(eachnextnode.nodeIndex)+" : "+ str(len(eachnextnode.nodeIndexList)))])
            if eachnextnode.visited == False:
                stack.append(eachnextnode)
                counter +=1 
     
        for eachprevnode in currentNode.listOfPrevNodes:
            G.add_edges_from([(str(eachprevnode.nodeIndex)+" : "+str(len(eachprevnode.nodeIndexList)), str(currentNode.nodeIndex)+" : "+ str(len(currentNode.nodeIndexList)))])
            if eachprevnode.visited == False:
                stack.append(eachprevnode)
                counter += 1 
            
    plt.figure(figsize=(30,15))
    #nx.draw(G,node_size=200,layout=nx.spring_layout(G))
    nx.draw(G,node_size=200,layout=nx.spring_layout(G))
    #nx.draw_circular(G)
    #nx.draw_spring(G) 
        
    plt.show()    
def checkCondensingBasic(G1,startList, typeOfGraph):
    print "Checking : ----- checkCondensingBasic(G1,startList)"
    print "len(G1), len(startList)", len(G1), len(startList)
    G = nx.MultiDiGraph()
    for eachnode in G1:
        #print G1
        if len(eachnode.listOfNextNodes) == 0:
            G.add_edges_from([(str(eachnode.nodeIndex)+" : "+str(len(eachnode.nodeIndexList)), str(eachnode.nodeIndex)+" : "+str(len(eachnode.nodeIndexList)))])
        for eachnextnode in eachnode.listOfNextNodes:
            if typeOfGraph == "simple":
                G.add_edges_from([(str(eachnode.nodeIndex)+" : "+str(len(eachnode.nodeIndexList)), str(eachnextnode.nodeIndex)+" : "+ str(len(eachnextnode.nodeIndexList)))])
            elif typeOfGraph == "MB":
                G.add_edges_from([(str(eachnode.nodeIndex)+" : "+str(len(eachnode.nodeIndexList)), str(eachnextnode[0].nodeIndex)+" : "+ str(len(eachnextnode[0].nodeIndexList)))])

    plt.figure(figsize=(30,15))
    nx.draw(G,node_size=200,layout=nx.spring_layout(G))
    #nx.draw_spring(G) 
    
    plt.show()
    
def branchClearingUnitTest(G, N, L,K, liid, threshold, foldername = "",branchingDepth= 20):
    dummyParameters = logging.parameterObj()
    
    dummyParameters.defaultFolder = foldername

    dummyParameters.brachingDepth = branchingDepth    
    
    dummyParameters.K, dummyParameters.liid, dummyParameters.threshold = K, liid, threshold
    dummyParameters.G, dummyParameters.N, dummyParameters.L = G, N, L 
    
    motherGen, reads, noisyReads = logging.rawDataLoad(foldername+"UnitTest",G,N,L, "dn")
    
    returnfmapping= logging.fmappingLoad(foldername+'clusteredGroup.csv')
    G1,startList, fmapping = graphForm.getSeqGraph(returnfmapping,noisyReads, dummyParameters)
    #checkCondensingBasic(G1, startList, "simple")


    returnfmapping, G1= branchClear.clearResidual(returnfmapping,G1,dummyParameters)
    
    G2 = G1

    graphForm.debugSeqGraph(G2)
    
    #G2 = G1 
    checkCondensingBasic2(G2, [G2[0]], "simple")
    
    

def resolveRepeatsUnitTest(foldername= ""):    
    var = raw_input("Enter something: ")
    print "you entered ", var
    f2= logging.fmapfusedLoad(foldername+'clusteredGroup2.csv') 

    G2 = logging.loadGraph(foldername+'basicMapping.csv', foldername+'seqMapping.txt', 'simple')
    
    checkCondensingBasic(G2, [G2[0]], "simple")
    dummyParameters = logging.parameterObj()
    dummyParameters.bridgingDepth = 10
    
    G3 = bridgeResolve.resolveRepeats(f2,G2,dummyParameters)
    
    
    checkCondensingBasic(G3, [G3[0]], "MB")
    
    
def MSAResolverEpeatUnitTest(Nin, Lin, foldername = ""):

    f2= logging.fmapfusedLoad('clusteredGroup2.csv') 
    G2 = logging.loadGraph(foldername+'basicMapping.csv', foldername+'seqMapping.txt', 'simple')
    #checkCondensingBasic(G2, [G2[0]], "simple")
    dummyParameters = logging.parameterObj()
    dummyParameters.bridgingDepth = 5
    dummyParameters.msaWidth = 20 
    
    G3 = bridgeResolve.resolveRepeats(f2,G2,dummyParameters)
    
    
    #checkCondensingBasic(G3, [G3[0]], "MB")

    N, G, L, p,snpRate, typeOfGen, detail = Nin,  10000,Lin, 0.015, 0.001 ,'m', "500-300-50" 
    motherGen, reads, noisyReads = logging.rawDataLoad(foldername+"UnitTest",G,N,L,"dn")

    alignmentBridge.MSAresolve(f2, G3, noisyReads, snpRate,dummyParameters)
    
def ECUnitTest(Nin,G, Lin, foldername = "", bridgingDepth = 20, msaWidth = 20 ):
    
    dummyParameters = logging.parameterObj()
    dummyParameters.bridgingDepth = bridgingDepth
    dummyParameters.msaWidth = msaWidth
    dummyParameters.N, dummyParameters.G, dummyParameters.L, dummyParameters.p  = Nin, G,Lin, 0.015
    dummyParameters.indel = True
    dummyParameters.defaultFolder = foldername
    dummyParameters.threshold = 5
    dummyParameters.liid = 48
    
    snpRate, typeOfGen, detail = 0.001 ,'m', "500-200-50" 
    G,N,L =  dummyParameters.G, dummyParameters.N, dummyParameters.L,
    motherGen, reads, noisyReads = logging.rawDataLoad(foldername+"UnitTest",G,N,L, "dn")  
  
    f2= logging.fmapfusedLoad(foldername+'clusteredGroup2.csv') 
    G2 = logging.loadGraph(foldername+'basicMapping.csv', foldername+'seqMapping.txt', 'simple')
    checkCondensingBasic(G2, [G2[0]], "simple")

    G3 = bridgeResolve.resolveRepeats(f2,G2,dummyParameters)
    
    checkCondensingBasic(G3, [G3[0]], "MB")

    G4 = alignmentBridge.MSAresolve(f2, G3, noisyReads, snpRate,dummyParameters)
    #G4 = G3

    #checkCondensingBasic(G4, [G4[0]], "MB")
    
    recovSeq = eulerCycle.findEC(G4)
    
    recovGen = readAns.reportRecovSeq(recovSeq, f2, noisyReads,dummyParameters)
    
    numMistakes, success = compare.subAlignCompare(recovGen, motherGen,dummyParameters)

    return numMistakes, success
    
def overallTestRun():
    #"oneLine.fasta-0-1440371" 
    G, N, L, p, epsilon, K, liid, threshold,NKcov, Nbridge, Ncov, ratio, numberOfClusterRounds,brachingDepth,bridgingDepth,msaWidth,folderName  =10000   ,83 ,   1189,  0.01 ,   0.1  ,  600 ,   32,    8  ,  508  ,  377 ,   417 ,   1.21822542 ,   0,    0,    0   , 0  , ""
    parameterRobot = logging.parameterObj(G, N, L, p, epsilon, K, liid, threshold,NKcov, Nbridge, Ncov, ratio, numberOfClusterRounds,brachingDepth,bridgingDepth,msaWidth,folderName )
    snpRate, typeOfGen, detail = 0.001 ,'d', "genome.fasta-50000-60000" 
    parameterRobot.clusterRounds , parameterRobot.fingerPrint, parameterRobot.clusterRatio = 2 , 6, 1
    #snpRate, typeOfGen , detail = 0.001,'m', "500-200-50" 
    numMistakes, success = assemblerMain.runAssembler(snpRate, typeOfGen, detail, parameterRobot)
    return numMistakes, success 

def tuneParamenters():
    numberOfPoints = 5
    numberOfRounds = 1
    branchDepth,bridgingDepth, msaWidth  = 1 , 10 , 10
    
    headerName = "synthetic_reads\\"
    
    listOfNLKDataPts = logging.loadingLNKFile(headerName)
    
    resultList  = []
    for index in range(numberOfPoints):
        N, L = listOfNLKDataPts[index][1], listOfNLKDataPts[index][2]
        print N, L 
        for roundNum in range(numberOfRounds):
            
            folderName = headerName + "sample_point_"+ str(index) +"\\round_"+str(roundNum) +"\\"
            
            branchClearingUnitTest(folderName, branchDepth)
            numMistakes, success = ECUnitTest(N, L, folderName,bridgingDepth, msaWidth  )
            
            resultList.append([N, L, roundNum , numMistakes, success])
            print "resultList", resultList
    
    
          
    logging.logBatch(resultList)
            



def multiplier(indexStart,indexEnd):
    def cmpFn(kmer1,kmer2):
        i =indexStart
        threshold = indexEnd
        
        while kmer1[0][i]==kmer2[0][i] and i< threshold:
                i= i+1
        
        return cmp(kmer1[0][i] ,kmer2[0][i])
 
    return cmpFn

def testMemoryUnit():
    print "hello"
    #kmerList = sorted(kmerList, cmp=multiplier(fingerprint*index,fingerprint*(index+1)) )
    
    
    #kmerList = sorted(kmerList, key = itemgetterkk(range(, )))

def compareAns(folderName = ""):

    f2 = open(folderName+"rec.txt", 'r')

    temp2 = f2.read()
    

    recov = np.zeros(len(temp2), dtype = np.int8)
    

    j =0 
    while j <len(temp2):
        if temp2[j] != '-':
            recov[j] = int(temp2[j])
            j = j+1 
        else:
            recov[j] = -1 
            j = j+2
            
            
            
    # Transform to FASTA
    fout= open(folderName+"recov.fasta", 'w')
    
    fout.write(">Seg1\n")
    for i in range(len(recov)):
        if recov[i]  ==1 :
            fout.write('A')
        elif recov[i]  ==2 :
            fout.write('C')
        elif recov[i]  ==3 :
            fout.write('G')
        elif recov[i]  ==4 :
            fout.write('T')
        else:
            fout.write('A')
        
        if np.mod(i,70) == 69 :
            fout.write("\n")

    fout.close()
    #compare.subAlignCompare(motherGen, recov,parameterRobot)
    
    
    

    f2.close()
    
    
def testBandedClustering():
    dummyParameters = logging.parameterObj()
    #dummyParameters.defaultFolder = ""
    dummyParameters.liid = 40
    dummyParameters.K = 40
    dummyParameters.threshold = 6
    dummyParameters.p = 0.015
    dummyParameters.indel = True
    
    dummyParameters.N, dummyParameters.G, dummyParameters.L, dummyParameters.p, typeOfGen, detail = 600, 10000,200, 0.015, 'm', "500-200-50" 
    
    
    startFingerPrint, endFingerPrint, read1, read2 = [] , [] , [] , []

    read1= np.zeros(200, dtype= np.int64)
    read2 = np.zeros(200, dtype= np.int64)
    
    for i in range(200):
        read1[i] = random.randint(1,4)
        read2[i] = random.randint(1,4)
    
    
    for i in range(50):
        read2[i] = 4
        read1[150+i] = 4
    
    
    print read1, read2
    startFingerPrint, endFingerPrint = [160,8] ,[ 190,40]  
    
    
    score , returnalignedSeq1, returnalignedSeq2 , starti, startj , endi, endj =cluster.SWAlignmentBanded(startFingerPrint, endFingerPrint,read1 ,read2 , dummyParameters)
    print "score", score 
    cluster.printSeq(returnalignedSeq1)
    cluster.printSeq(returnalignedSeq2)
    
    print "starti, startj, endi , endj : ", starti, startj, endi , endj


    score , returnalignedSeq1, returnalignedSeq2 , starti, startj , endi, endj =cluster.SWAlignment(read1, read2 , dummyParameters)
    print "score", score 
    cluster.printSeq(returnalignedSeq1)
    cluster.printSeq(returnalignedSeq2)
    
    print "starti, startj, endi , endj : ", starti, startj, endi , endj
    
    
def MUMMERBatch():
    G = 10000


def segmentChopAndTestTest(Nin,G, Lin, foldername = "", bridgingDepth = 20, msaWidth = 20):         
    dummyParameters = logging.parameterObj()
    dummyParameters.bridgingDepth = bridgingDepth
    dummyParameters.msaWidth = msaWidth
    dummyParameters.N, dummyParameters.G, dummyParameters.L, dummyParameters.p  = Nin, G,Lin, 0.015
    dummyParameters.indel = True
    dummyParameters.defaultFolder = foldername
    dummyParameters.threshold = 5
    dummyParameters.liid = 48
    
    snpRate, typeOfGen, detail = 0.001 ,'m', "500-200-50" 
    G,N,L =  dummyParameters.G, dummyParameters.N, dummyParameters.L,
    motherGen, reads, noisyReads = logging.rawDataLoad(foldername+"UnitTest",G,N,L, "dn")  
    
    
    #G3, f2 = interactiveNodeRemoval(folderName, "skip")
    G3, f2 = interactiveNodeRemoval(folderName, "desiredOpt")
    
    chopAndAlign(G3, f2 , noisyReads, motherGen, dummyParameters)



def chopAndAlign(G3, f2 , noisyReads, motherGenome, parameterRobot):
    print "Hi, I am working NOW"
    # Input : G3 ~ [nodeIndexList] f2 ~ [clusterMapping] noisyReads ~ [noisyReads] 
    # Output : .FASTA file containing the recovered Segments , dotplots of the segment against the real genome

    segmentList = []
    
    # Chop and fill in 
    #print len(G3)
    for eachitem in G3: 
        #print eachitem
        tempSeq = eachitem.nodeIndexList
        recovGen = readAns.reportRecovSeq(tempSeq, f2, noisyReads,parameterRobot)
        segmentList.append(recovGen)
    
    
    # Create dotpolts and fasta output
    #myFile = open(parameterRobot.defaultFolder, 'w')
    
    print "len(segmentList)" , len(segmentList)
    for eachSegment, index in zip(segmentList, range(len(segmentList))):
        #if len(eachSegment) == 1963 :
        frecov = open(parameterRobot.defaultFolder+"rec_"+str(index)+".txt", 'w')
        for eachbase in eachSegment:
            frecov.write(str(eachbase))
        frecov.close()
        
        print "index, eachSegment[0:10]",index, eachSegment[0:10]
        
        compare.outputToFastaFiles(eachSegment, motherGenome, parameterRobot, index)    
    
        
    #myFile.close()
    
   
def interactiveNodeRemoval(foldername, modeOfOpt = "skip"):
     
    '''
    delnode 1002
    deledge 1004, 1028
    addedge 1023 , 3434
    '''
    
    G2 = []
    if modeOfOpt == "skip":
        varList = ["start none", ""]
    elif modeOfOpt == "desiredOpt":
        varList = ["start none", "delnode 734288", "delnode 1716432", "fusenode 1438104", "fusenode 1166060" , "view", ""] 
    
    var = varList.pop(0)
    
    while len(var) > 0 :
        
        command = var.split()
        if command[0] == "start":   
            print "Start"
            G2 = logging.loadGraph(foldername+'basicMapping.csv', foldername+'seqMapping.txt', 'simple')
                
        elif command[0] == "delnode":
            print "To delete node"
            currentNodeIndex = int(command[1])
            print "currentNodeIndex", currentNodeIndex
            for eachnode in G2: 
                if eachnode.nodeIndex == currentNodeIndex:
                    currentNode = eachnode
                    
            for eachitem in currentNode.listOfPrevNodes:
                eachitem.listOfNextNodes.remove(currentNode)
            currentNode.listOfPrevNodes = []
            
            for eachitem in currentNode.listOfNextNodes:
                eachitem.listOfPrevNodes.remove(currentNode)
            currentNode.listOfNextNodes = []
            
            currentNode.nodeIndexList = [] 
            graphForm.condenseGraph(G2)    
                          
        elif command[0] == "deledge":
            node1 = int(command[1])
            node2 = int(command[2])
            startNode, endNode = [] , []
            for eachnode in G2: 
                if eachnode.nodeIndex == node1:
                    startNode = eachnode
                     
                if eachnode.nodeIndex == node2:
                    endNode = eachnode
                    
            startNode.listOfNextNodes.remove(endNode)
            endNode.listOfPrevNodes.remove(startNode)
                    
            print "To insert nodes"
        elif command[0] == "addedge" : 
            node1 = int(command[1])
            node2 = int(command[2])
            startNode, endNode = [] , []
            for eachnode in G2: 
                if eachnode.nodeIndex == node1:
                    startNode = eachnode
                if eachnode.nodeIndex == node2:
                    endNode = eachnode
            
            startNode.listOfNextNodes.append(endNode)
            endNode.listOfPrevNodes.append(startNode)
        elif command[0] == "fusenode":
            myNodeIndex = int(command[1])
            for eachnode in G2:
                if eachnode.nodeIndex == myNodeIndex:
                    currentNode = eachnode
            if len(currentNode.listOfPrevNodes) == 1:
                prevNode = currentNode.listOfPrevNodes[0]
                for eachnextnode in currentNode.listOfNextNodes:
                    eachnextnode.listOfPrevNodes.remove(currentNode)
                    eachnextnode.listOfPrevNodes.append(prevNode)
                    prevNode.listOfNextNodes.append(eachnextnode)
                
                prevNode.listOfNextNodes.remove(currentNode)
                currentNode.listOfNextNodes = [] 
                currentNode.listOfPrevNodes =[]
                currentNode.nodeIndexList = [] 
            
            elif len(currentNode.listOfNextNodes) ==1 :
                nextNode = currentNode.listOfNextNodes[0]
                for eachprevnode in currentNode.listOfPrevNodes:
                    eachprevnode.listOfNextNodes.remove(currentNode)
                    eachprevnode.listOfNextNodes.append(nextNode)
                    nextNode.listOfPrevNodes.append(eachprevnode)
                
                nextNode.listOfPrevNodes.remove(currentNode)
                currentNode.listOfNextNodes = []
                currentNode.listOfPrevNodes = []
                currentNode.nodeIndexList = [] 
            
        elif command[0] == "view":
            G2 = graphForm.newCondensingStep(G2)
            G2 = G2[0]
            checkCondensingBasic(G2, [G2[0]], "simple")
        #var = raw_input("Enter Operations: ")
        var = varList.pop(0)
    
        

    dummyParameters = logging.parameterObj()
    dummyParameters.bridgingDepth = 20
    
    print "Loading fmap and graph "
    f2= logging.fmapfusedLoad(foldername+'clusteredGroup2.csv') 
    
    G3 = bridgeResolve.resolveRepeats(f2,G2,dummyParameters)
    assert(1==2)
    
    G3 = G2
    print "Done Loading the fmap and graph"
    
    #checkCondensingBasic(G3, [G3[0]], "MB")
    
    return  G3,f2
    
    #dynamicAsking(G3, f2)

def dynamicAsking(G3, f2):
    searchDepth = 5
    
    var = "start none"
    
    while len(var) > 0 : 
        command = var.split()
        if command[0] == "printlist":
            myNodeIndex = int(command[1])
            myNode = []
            for eachnode in G3:
                if eachnode.nodeIndex == myNodeIndex:
                    myNode = eachnode
                    
            print "nodeIndex, nodeIndexList : ",  myNode.nodeIndex, myNode.nodeIndexList 
            
        elif command[0] == "printstart":
            myNodeIndex = int(command[1])
            frankingdepth = int(command[2])
            myNode = []
            
            for eachnode in G3:
                if eachnode.nodeIndex == myNodeIndex:
                    myNode  = eachnode 


            frankinginList = myNode.nodeIndexList[frankingdepth: frankingdepth+searchDepth ]         

            inList = bridgeResolve.findRangeList(frankinginList, f2 )
            print "inList : ", inList 
            
            
        elif command[0] == "printend":
            myNodeIndex = int(command[1])
            frankingdepth = int(command[2])
            myNode = []
            
            for eachnode in G3:
                if eachnode.nodeIndex == myNodeIndex:
                    myNode  = eachnode 
                    
            frankingoutList = myNode.nodeIndexList[-frankingdepth- searchDepth: -frankingdepth ]         

            outList = bridgeResolve.findRangeList(frankingoutList, f2 )
            print "outList : ", outList 

    
        
        var = raw_input("Enter Operations: ")
        
    
def testChecking(folderName, N, L, G ):
    #N, L, G =     1345, 1000, 50000
    f2 = open(folderName+"rec.txt", 'r')
    temp2 = f2.read()
    recov = np.zeros(len(temp2), dtype = np.int32)

    j =0 
    while j <len(temp2):
        if temp2[j] != '-':
            recov[j] = int(temp2[j])
            j = j+1 
        else:
            recov[j] = -1 
            j = j+2
    f2.close()
    dummyParameters = logging.parameterObj()
    dummyParameters.defaultFolder = folderName 
    dummyParameters.G = G
    motherGen, reads, noisyReads = logging.rawDataLoad(folderName+"UnitTest",G,N,L, 'd')
    
    print len(recov) , len(motherGen)
    numberMistakes, success = compare.subAlignCompare(recov, motherGen,dummyParameters)
    
def process(count, jobid, output):
    pom = []
    for i in range(count):
        pom.append(random.random())

    print "Job ", jobid, " finished!"

def process2(count, jobid, output):
    pom = []
    for i in range(count):
        pom.append(random.random())
    print "Job ", jobid, " finished!"
    output.append(1)

def debuggingParallel(myoption):
    print "Hello world"
    times = 100000
    if myoption == 1 :
        out1 = list()
        out2 = list()
        
        thread1 = threading.Thread(target=process(times, 'A', out1))
        thread2 = threading.Thread(target=process(times, 'B', out2))
        
        job = []
        job.append(thread1)
        job.append(thread2)
        
        for i in job:
            i.start()
        for i in job:
            i.join()
        
        print "Finished!"
    elif myoption == 2 :


        out1 = list()
        out2 = list()
        
        job = []
        job.append(Process(target=process2, args=(times, 'A', out1)))
        job.append(Process(target=process2, args=(times, 'B', out2)))
        
        print "out1", out1
        print"out2", out2
        
        for j in job:
            j.start()
        for j in job:
            j.join()
        
        print "Finished!"
        

  
shared_array_base = []
shared_array = []
shared_array = []
A = []

def my_func( i, def_param=[shared_array]):
    shared_array[:,i] = i

    return [i]

def mycallback(x):
    #print('mycallback is called with {}'.format(x))
    A.extend(x)
     
def sharedMemory(numProc):
    # No copy was made
    n =1000
    global shared_array_base 
    shared_array_base = multiprocessing.Array(ctypes.c_double, n*n) 
    global shared_array 
    shared_array = np.ctypeslib.as_array(shared_array_base.get_obj())
    shared_array = shared_array.reshape(n, n)   
    global A 
    assert shared_array.base.base is shared_array_base.get_obj()
    

    # Parallel processing
    pool = multiprocessing.Pool(processes=numProc)
    #pool.map(my_func,range(n))
    r = pool.map_async(my_func, [[1,1], [2,1], [3,1]], callback=mycallback)
    r.wait()

    print shared_array
    print A
    A = []

def transformReads(G,N,L,folderName,myIndex):

    motherGen, reads, noisyReads = logging.rawDataLoad(folderName+"UnitTest",G,N,L, "dn")
    
    myRead = reads[myIndex]
    
    fout= open(folderName + "bridgingReads_"+str(myIndex)+".fasta", 'w')
    
    fout.write(">Seg1\n")
    for i in range(len(myRead)):
        if myRead[i]  ==1 :
            fout.write('A')
        elif myRead[i]  ==2 :
            fout.write('C')
        elif myRead[i]  ==3 :
            fout.write('G')
        elif myRead[i]  ==4 :
            fout.write('T')
        #else:
        #    fout.write('A')
        
        if np.mod(i,70) == 69 :
            fout.write("\n")

    fout.close()
    
    print "len(motherGen) , len(reads), len(noisyReads)  : ",  len(motherGen) , len(reads), len(noisyReads)  
#dataGenUnitTest()

#dataGenUnitTest2()

#clusterUnitTest()
#clusterUnitTest2()

#MUMMERBatch()
#testBioLib()
#testBandedClustering()
#graphFormUnitTest(folderName)

#cProfile.run('graphFormUnitTest(folderName)')
#graphFormUnitTest2(folderName)

#testMemoryUnit()
#folderName = "real_genome_synthetic_reads_2\\sample_point_12\\round_2\\"
#folderName = "synthetic_reads\\sample_point_0\\round_0\\"

folderName = "download/methWhole/"
#folderName = "synthetic_reads/sample_point_0/round_1/"

t0 = time.time()

def speedTest():
            
    #N, L, G =     869, 180, 10000
    #N, L, G = 339, 2000, 50000
    #N,L, G =     927, 200,10000
    #N, L, G=     1345, 1000, 50000
    #G, N, L = 1440371,41097,1000
    #G, N, L = 1589953,11325,2900
    G, N, L = 1772693,9096,4279
    #G, N, L = 10000    ,927,    200
    #G, N, L = 50000, 6708 ,180
    K = 32
    threshold = 5
    liid = 48
    #K = 600
    
    
    
    result = []
    
    for i in [7]:
        #folderName = "synthetic_reads/sample_point_0/round_"+ str(i) +"/"
        #folderName = "download/round_"+ str(i)+"/"
        #graphFormUnitTest(N, G, L, folderName)
        #graphFormUnitTest2(N, G, L, folderName)
        #branchClearingUnitTest(G, N, L,K, liid, threshold, folderName,branchingDepth=20)

        #resolveRepeatsUnitTest(folderName)
        #numMistakes, success = ECUnitTest(N,G, L, folderName,30,30) 
        #testChecking(folderName, N, L , G )
        #print len(motherGen), len(reads), len(noisyReads)
        #result.append([numMistakes, success])
        #os.system("bash ../gepard-1.30/gepardcmd.sh -seq1 " + seq1Name +" -seq2 "+seq2Name+" -matrix matrices/edna.mat -outfile Jane.png")
        #plt.savefig("Resolved_"+str(i)+".png")
        #segmentChopAndTestTest(N,G, L, folderName,30,30)
        transformReads(G,N,L,folderName, 8115)
    print "result: ", result
    

#compareAns()
#logging.savingGenomeSegmentFile("")
#clusteringTestingOfK(G, N, L, folderName, K, threshold, liid)
#branchClearingUnitTest(G, N, L,K, liid, threshold, folderName, 5)
#resolveRepeatsUnitTest(folderName)
#MSAResolverEpeatUnitTest(N, L, folderName)
#numMistakes, success = ECUnitTest(N, L, folderName,12,20) 
#compareAns()
#print numMistakes, success
#listOfTest = [] 
#for index in range(5):
#    listOfTest.append( overallTestRun())
#    print listOfTest
#tuneParamenters()
#logging.fetchResults("largeScaleTest\\", 10, 3)
#compare.transformToFASTA(folderName+"UnitTest_motherGen.txt", folderName+"UnitTest_motherGen.fasta")
#logging.savingLNKFile("")
#cProfile.run('testParallel()')


speedTest()
#debuggingParallel(2)
#print "Time (Sec): ", time.time() - t0


