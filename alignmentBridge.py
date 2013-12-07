import numpy as np 
import bridgeResolve
import math 
from operator import itemgetter
import cleaner
import extender
import common

def MAPDetermination(readList, noisyReads, p, snpRate):
    countArray = [0,0,0,0,0]
    for eachcopy in readList:
        readNum, offset = eachcopy[0], eachcopy[1]
        countArray[noisyReads[readNum][offset] ] = countArray[noisyReads[readNum][offset] ] + 1
    
    totalCount = np.sum(countArray)
    
    possibleChoices = [[1] , [2] , [3] , [4] , [1,2] , [1,3], [1,4], [2,3], [2,4], [3,4]]
    
    #print "countArray", countArray
    
    voteArray = []
    for index in range(len(possibleChoices)):
        setOfTrial = possibleChoices[index] 
        count = 0
        for dummyindex in setOfTrial:
            count = count + countArray[dummyindex]
        #print "count",count
        if len(setOfTrial) == 1:
            score = (1-snpRate) * pow((1-p), count) * pow((p/3), totalCount - count)
            voteArray.append(score)
        elif len(setOfTrial) == 2:
            score = snpRate * pow((1-p)/2 + p/6, count) * pow(p/3, totalCount - count)
            voteArray.append(score)
    
    print "voteArray", voteArray[0]        
    optIndex = np.argmax(voteArray)
    thisResult = possibleChoices[optIndex]
    #print "thisResult", thisResult
    return thisResult

def checkMSABridging(f2, currentNode, noisyReads,p, snpRate, flankinglen):
    
    isMSABridged = True

    matchingList = []
    readList = []
    voteResultList = []
    approxRepeatLen = len(currentNode.nodeIndexList)
    
    # 1 ) Majority at each location to determine SNPs    
    for eachposition in currentNode.nodeIndexList:
        KmerIndex = eachposition
        
        readList = bridgeResolve.obtainReadNum(KmerIndex, f2)
            
        voteResultAtPoistion = MAPDetermination(readList, noisyReads, p, snpRate)

        voteResultList.append(voteResultAtPoistion)
    
    
    #print "voteResultList",voteResultList 
    
    # 2 ) Identify inList , bridgeList , outList
    # a) Anchor Position
    print "voteResultList",voteResultList
    runOfLengthList = []
    runningindex = 0
    tempCount  =0 
    tempStart  =0 
    while (runningindex  < len(voteResultList)):
        if len(voteResultList[runningindex]) == 2:
            runOfLengthList.append([tempStart, tempCount])
            tempCount = 0
            tempStart = runningindex +1 
        elif len(voteResultList[runningindex]) == 1:
            tempCount = tempCount +1 
        
        runningindex = runningindex +1
    
    if tempCount > 0:
        runOfLengthList.append([tempStart, tempCount])
        
    eachindex = 0

    print "approxRepeatLen , flankinglen",approxRepeatLen , flankinglen
    
    while ( eachindex <len(runOfLengthList)):
        if runOfLengthList[eachindex][0] < flankinglen or  runOfLengthList[eachindex][0] > approxRepeatLen-flankinglen: 
            runOfLengthList.pop(eachindex)
        else:
            eachindex = eachindex + 1

    if len(runOfLengthList) == 0:
        isMSABridged, matchingList = False, []
        return isMSABridged, matchingList
        
    runOfLengthList = sorted(runOfLengthList, key = itemgetter(1))
    print "runOfLengthList",runOfLengthList 
    longestRepeat = runOfLengthList[-1]
    
    lrepStart, lrepEnd = longestRepeat[0] - 1 , longestRepeat[0] + longestRepeat[1] 
    
    print "lrepStart, lrepEnd", lrepStart, lrepEnd
    print "longestRepeat[1] , approxRepeatLen" , longestRepeat[1] , approxRepeatLen
    
    ###########################################Special treatment for 1 SNP
    if len(runOfLengthList) == 1:
        print "Special treatment for 1 intermediate SNP"
        inKmerIndex, outKmerIndex = currentNode.nodeIndexList[lrepStart], currentNode.nodeIndexList[lrepStart]
        
        readsAtInList = bridgeResolve.obtainReadNum(inKmerIndex,f2)
        readsAtOutList = bridgeResolve.obtainReadNum(outKmerIndex, f2)
        
        print "len(readsAtInList)", len(readsAtInList)
        inLen ,midLen, outLen, L = lrepStart, 0, approxRepeatLen - lrepStart, len(noisyReads[0])
        inList ,bridgeList , outList = [] , [] ,[]
    
        runningindex  = 0 
        
        while (runningindex < len(readsAtInList)):
            readNum , offset = readsAtInList[runningindex][0], readsAtInList[runningindex][1]
            if offset  - inLen > 0:
                inList.append([readNum, offset])
                readsAtInList.pop(runningindex)
            else:  
                runningindex = runningindex + 1
        
        runningindex  = 0 

        print "readsAtInList",readsAtInList
        
        print "outLen, L", outLen, L
        
        print "len(readsAtOutList)",len(readsAtOutList), readsAtOutList
       
        
        while (runningindex < len(readsAtOutList)):
            readNum , offset = readsAtOutList[runningindex][0], readsAtOutList[runningindex][1]     
            print "offset  ", offset
            if offset + outLen < L :
                outList.append([readNum , offset ])
                readsAtOutList.pop(runningindex)
            else:
                runningindex = runningindex + 1
        
        print "inList,", inList
        print "bridgeList", bridgeList
        print "outList", outList
        
        # 3 ) Establsh H.T. and count vote 
        # a ) Classify into groups
        refinedInList = []
        refinedOutList = []
        for eachinnode in currentNode.listOfPrevNodes:
            tmpgpList = [[],[]]
            edgeWt = eachinnode[1]
            inKmerIndex = eachinnode[0].nodeIndexList[-(edgeWt +1)]
            
            tmpgpList[0]= [inKmerIndex]
            
            associatedReads = bridgeResolve.obtainReadNum(inKmerIndex, f2)
            onlyReadNumList = filterReadNum(associatedReads)
            for eachinread in inList:
                if eachinread[0] in  onlyReadNumList:
                    tmpgpList[1].append(eachinread)
                    
            refinedInList.append(tmpgpList)
        
         
        for eachoutnode in currentNode.listOfNextNodes:
            tmpgpList = [[],[]]
            edgeWt = eachoutnode[1]
            outKmerIndex = eachoutnode[0].nodeIndexList[edgeWt ]
            
            tmpgpList[0]= [outKmerIndex]
            
            associatedReads = bridgeResolve.obtainReadNum(outKmerIndex, f2)
            onlyReadNumList = filterReadNum(associatedReads)
            for eachoutread in outList:
                if eachoutread[0] in  onlyReadNumList:
                    tmpgpList[1].append(eachoutread)
                    
            refinedOutList.append(tmpgpList)    
    
        inSNPResultList = []
        for eachitem in refinedInList:
            kmerIndex = eachitem[0]
            voteChoice = [0,0,0,0,0]
            for eachcopy in eachitem[1]:
                readNum, offset = eachcopy[0], eachcopy[1]
                voteChoice[noisyReads[readNum][offset]] = voteChoice[noisyReads[readNum][offset]]  + 1
            print "voteChoice", voteChoice      
            inSNPResultList.append([kmerIndex, voteChoice])   
            
        outSNPResultList = []
        for eachitem in refinedOutList:
            kmerIndex = eachitem[0]
            voteChoice = [0,0,0,0,0]
            for eachcopy in eachitem[1]:
                readNum, offset = eachcopy[0], eachcopy[1]
                voteChoice[noisyReads[readNum][offset]] = voteChoice[noisyReads[readNum][offset]]  + 1
            print "voteChoice", voteChoice  
            outSNPResultList.append([kmerIndex, voteChoice])   
        
        
        print "inSNPResultList, outSNPResultList",      inSNPResultList , outSNPResultList
        print "voteResultList[lrepStart]", voteResultList[lrepStart]
        
        SNP1, SNP2 = voteResultList[lrepStart][0], voteResultList[lrepStart][1]
        H0CountStart = inSNPResultList[0][1][SNP1] + inSNPResultList[1][1][SNP2]
        H1CountStart = inSNPResultList[0][1][SNP2] + inSNPResultList[1][1][SNP1]
        
        H0CountEnd = outSNPResultList[0][1][SNP1] + outSNPResultList[1][1][SNP2]
        H1CountEnd = outSNPResultList[0][1][SNP2] + outSNPResultList[1][1][SNP1]       
        
        if ( (H0CountStart < H1CountStart) and  (H0CountEnd < H1CountEnd ) )or ( (H0CountStart >H1CountStart) and  (H0CountEnd > H1CountEnd ) ): 
        
            isMSABridged= True
            matchingList = [[inSNPResultList[0][0][0]  , outSNPResultList[0][0][0]],[inSNPResultList[1][0][0], outSNPResultList[1][0][0]]]
        else:
            isMSABridged= True
            matchingList = [[inSNPResultList[0][0][0]  , outSNPResultList[1][0][0]],[inSNPResultList[1][0][0], outSNPResultList[0][0][0]]]
            
            
        print matchingList
        return isMSABridged, matchingList

    
    
    
    ###############################################End special treatment
    
    if 3*longestRepeat[1] < approxRepeatLen or lrepEnd >=len(currentNode.nodeIndexList) or len(noisyReads[0]) <= longestRepeat[1] :
        isMSABridged, matchingList = False, []
        return isMSABridged, matchingList
        
    # b) Find associated Reads
    inKmerIndex, outKmerIndex = currentNode.nodeIndexList[lrepStart], currentNode.nodeIndexList[lrepEnd]
    
    readsAtInList = bridgeResolve.obtainReadNum(inKmerIndex,f2)
    readsAtOutList = bridgeResolve.obtainReadNum(outKmerIndex, f2)
    
    inLen ,midLen, outLen, L = lrepStart, longestRepeat[1], approxRepeatLen - lrepEnd, len(noisyReads[0])
    
    print "inLen ,midLen, outLen, L", inLen ,midLen, outLen, L
    
    inList ,bridgeList , outList = [] , [] ,[]
    
    runningindex  = 0 
    
    while (runningindex < len(readsAtInList)):
        readNum , offset = readsAtInList[runningindex][0], readsAtInList[runningindex][1]
        if offset + midLen < L -1 :
            bridgeList.append([readNum , offset])
            readsAtInList.pop(runningindex)
        elif offset  - inLen > 0:
            inList.append([readNum, offset])
            readsAtInList.pop(runningindex)
        else:  
            runningindex = runningindex + 1
    
    runningindex  = 0 
    
    

    
    while (runningindex < len(readsAtOutList)):
         readNum , offset = readsAtOutList[runningindex][0], readsAtOutList[runningindex][1]       
         if offset + outLen < L :
             outList.append([readNum , offset ])
             readsAtOutList.pop(runningindex)
         else:
             runningindex = runningindex + 1
    
    #print "inList,", inList
    #print "bridgeList", bridgeList
    #print "outList", outList
    
    # 3 ) Establsh H.T. and count vote 
    # a ) Classify into groups
    refinedInList = []
    refinedOutList = []
    for eachinnode in currentNode.listOfPrevNodes:
        tmpgpList = [[],[]]
        edgeWt = eachinnode[1]
        inKmerIndex = eachinnode[0].nodeIndexList[-(edgeWt +1)]
        
        tmpgpList[0]= [inKmerIndex]
        
        associatedReads = bridgeResolve.obtainReadNum(inKmerIndex, f2)
        onlyReadNumList = filterReadNum(associatedReads)
        for eachinread in inList:
            if eachinread[0] in  onlyReadNumList:
                tmpgpList[1].append(eachinread)
                
        refinedInList.append(tmpgpList)
    
     
    for eachoutnode in currentNode.listOfNextNodes:
        tmpgpList = [[],[]]
        edgeWt = eachoutnode[1]
        outKmerIndex = eachoutnode[0].nodeIndexList[edgeWt ]
        
        tmpgpList[0]= [outKmerIndex]
        
        associatedReads = bridgeResolve.obtainReadNum(outKmerIndex, f2)
        onlyReadNumList = filterReadNum(associatedReads)
        for eachoutread in outList:
            if eachoutread[0] in  onlyReadNumList:
                tmpgpList[1].append(eachoutread)
                
        refinedOutList.append(tmpgpList)    

    inSNPResultList = []
    for eachitem in refinedInList:
        kmerIndex = eachitem[0]
        voteChoice = [0,0,0,0,0]
        for eachcopy in eachitem[1]:
            readNum, offset = eachcopy[0], eachcopy[1]
            voteChoice[noisyReads[readNum][offset]] = voteChoice[noisyReads[readNum][offset]]  + 1
        print "voteChoice", voteChoice      
        inSNPResultList.append([kmerIndex, np.argmax(voteChoice)])   
        
    outSNPResultList = []
    for eachitem in refinedOutList:
        kmerIndex = eachitem[0]
        voteChoice = [0,0,0,0,0]
        for eachcopy in eachitem[1]:
            readNum, offset = eachcopy[0], eachcopy[1]
            voteChoice[noisyReads[readNum][offset]] = voteChoice[noisyReads[readNum][offset]]  + 1
        print "voteChoice", voteChoice  
        outSNPResultList.append([kmerIndex, np.argmax(voteChoice)])   
    
    
    print    "inSNPResultList, outSNPResultList",      inSNPResultList , outSNPResultList
    # b) bridge them
    hypothesis =[0, 0]    
    for eachitem in bridgeList:
        startReadNum, startOffset = eachitem[0], eachitem[1]
        endReadNum , endOffset = eachitem[0],  startOffset + midLen + 1 
       # print "startReadNum, startOffset", startReadNum, startOffset
       # print "endReadNum , endOffset ", endReadNum , endOffset 
        
       # print "noisyReads[startReadNum][startOffset], noisyReads[endReadNum][endOffset]" , noisyReads[startReadNum][startOffset], noisyReads[endReadNum][endOffset]
        if noisyReads[startReadNum][startOffset] == inSNPResultList[0][1] and noisyReads[endReadNum][endOffset] == outSNPResultList[0][1]:
            hypothesis[0] = hypothesis[0] + 1
        
        print endReadNum, endOffset
        print noisyReads[endReadNum][endOffset]
       
        
        if noisyReads[startReadNum][startOffset] == inSNPResultList[1][1] and noisyReads[endReadNum][endOffset] == outSNPResultList[1][1]:
            hypothesis[0] = hypothesis[0] + 1
        
        if noisyReads[startReadNum][startOffset] == inSNPResultList[0][1] and noisyReads[endReadNum][endOffset] == outSNPResultList[1][1]:
            hypothesis[1] = hypothesis[1] + 1
        
        
        if noisyReads[startReadNum][startOffset] == inSNPResultList[1][1] and noisyReads[endReadNum][endOffset] == outSNPResultList[0][1]:
            hypothesis[1] = hypothesis[1] + 1
        
    print "hypothesis", hypothesis
    ansIndex = np.argmax(hypothesis)
    
    if ansIndex == 0:
        matchingList = [[inSNPResultList[0][0][0],outSNPResultList[0][0][0]], [inSNPResultList[1][0][0], outSNPResultList[1][0][0]]]
    elif ansIndex == 1:
        matchingList = [[inSNPResultList[0][0][0],outSNPResultList[1][0][0]], [inSNPResultList[1][0][0], outSNPResultList[0][0][0]]]
        
    isMSABridged = True
    
    print matchingList
    return isMSABridged, matchingList


def filterReadNum(listOfReadsDetail):    
    readList = []
    for eachitem in listOfReadsDetail:
        readList.append(eachitem[0])
    return readList
    
def MSAresolve(f2, G3, noisyReads ,snpRate ,parameterRobot):
    G4 = []
    p,flankinglen = parameterRobot.p, parameterRobot.bridgingDepth
    print "---------------"
    print "MSA Resolve : "
    print "p, snpRate ", p , snpRate 
    print "---------------"
    # Flow of the resolving 
    f2 = sorted(f2) 
    
    xNodesList = []
    nodeLabelList = []
    for eachnode in G3:
        if len(eachnode.listOfNextNodes) > 1 and len(eachnode.listOfPrevNodes) > 1:
            xNodesList.append(eachnode)
        nodeLabelList.append(eachnode.nodeIndex)
    
    nodeLabelList = sorted(nodeLabelList)
    
    startingGpNum = nodeLabelList[-1]
    while (len(xNodesList) > 0):

        currentNode = xNodesList.pop(0)
        if parameterRobot.indel == False:
            canResolve, kmerPairsList = checkMSABridging(f2, currentNode,noisyReads,p, snpRate,flankinglen)
        else:
            canResolve, kmerPairsList = indelMSABridging(f2, currentNode,noisyReads,p, snpRate,flankinglen, parameterRobot)
            
        print "canResolve", canResolve 
        newXNodeList,startingGpNum =  bridgeResolve.resolveFramework(currentNode, f2,startingGpNum, G3, canResolve, kmerPairsList )
        xNodesList = xNodesList + newXNodeList    
        

    # Format output 
    G4 = G3
    bridgeResolve.removeEmptyNodes(G4)
    
    return G4 



def indelMSABridging(f2, currentNode,noisyReads,p, snpRate,flankinglen, parameterRobot):
    
    # Need to loop over all the xnodes 
    canResolve, kmerPairsList = False, []
    
    # Using contig Creator 
    indelRobot = common.parameterRobot()
    indelRobot.defaultFolder = parameterRobot.defaultFolder
    indelRobot.setReadStat( Nshort= parameterRobot.N, Nlong=  parameterRobot.N, Lshort= parameterRobot.L, Llong= parameterRobot.L, p= parameterRobot.p , longOnly = True)
    indelRobot.setGenomeStat(G = parameterRobot.G, lrep=500, lsnp=200, lint=50 )
    indelRobot.setThresholdPara(liid = 30, thresForRandom= 0.5,thresForins =0.4, thresFordel=0.4, insMin=4, delMin=4,thresholdForSupport= 0.15, subthreshold= 9, editsub= -10, editins= -1, editdel= -1, editmatch = 1, lookRange =15)
    indelRobot.tunePara()
    indelRobot.snprate = snpRate
    

    # toProcessList : in1IndexList, in2IndexList, out1IndexList, out2IndexList, commonIndexList
    # shortToLongMap : indexlong    indexshort    jstart    jend    istart    iend    
    
    shortToLongMap,toProcessList = [], []
    
    toProcessList = formToProcessList(f2, noisyReads, currentNode, indelRobot, flankinglen)
    if len(toProcessList[4] ) == 0 : 
        return False, []
    
    shortToLongMap = formRelatedMap(f2, noisyReads, currentNode, indelRobot, toProcessList)
    
    
    cleaner.cleaning([noisyReads,noisyReads] ,shortToLongMap, toProcessList,indelRobot, "init")
    in1List, in2List, out1List, out2List, commonList, longReadToUse  = cleaner.cleaning([noisyReads, noisyReads],shortToLongMap, toProcessList,indelRobot, "vote")

    
    extendResult = extender.readExtender(in1List, in2List, out1List, out2List, commonList,indelRobot,longReadToUse, True)
    
    
    
    if extendResult == 0 : 
        canResolve = True
        
        edgeWt = currentNode.listOfPrevNodes[0][1]
        inKmerIndex = currentNode.listOfPrevNodes[0][0].nodeIndexList[-(edgeWt +1)]
        
        edgeWt = currentNode.listOfNextNodes[0][1]
        outKmerIndex = currentNode.listOfNextNodes[0][0].nodeIndexList[edgeWt ]
        
        kmerPairsList.append([inKmerIndex,outKmerIndex ])
        
        edgeWt = currentNode.listOfPrevNodes[1][1]
        inKmerIndex = currentNode.listOfPrevNodes[1][0].nodeIndexList[-(edgeWt +1)]
        
        edgeWt = currentNode.listOfNextNodes[1][1]
        outKmerIndex = currentNode.listOfNextNodes[1][0].nodeIndexList[edgeWt ]
        
        kmerPairsList.append([inKmerIndex,outKmerIndex ])
    elif extendResult == 1 :
        canResolve = True
        edgeWt = currentNode.listOfPrevNodes[0][1]
        inKmerIndex = currentNode.listOfPrevNodes[0][0].nodeIndexList[-(edgeWt +1)]
        
        edgeWt = currentNode.listOfNextNodes[1][1]
        outKmerIndex = currentNode.listOfNextNodes[1][0].nodeIndexList[edgeWt ]
        
        kmerPairsList.append([inKmerIndex,outKmerIndex ])
        
        edgeWt = currentNode.listOfPrevNodes[1][1]
        inKmerIndex = currentNode.listOfPrevNodes[1][0].nodeIndexList[-(edgeWt +1)]
        
        edgeWt = currentNode.listOfNextNodes[0][1]
        outKmerIndex = currentNode.listOfNextNodes[0][0].nodeIndexList[edgeWt ]
        
        kmerPairsList.append([inKmerIndex,outKmerIndex ])
        
    elif extendResult == -1:
        canResolve = False
        kmerPairsList = []
        
    print kmerPairsList
    
    return canResolve, kmerPairsList 



def formToProcessList(f2, noisyReads, currentNode, indelRobot, flankinglen):
    
    print "formToProcessList"
    print "NodeDetail : ", currentNode.nodeIndex, len(currentNode.nodeIndexList) , len(currentNode.listOfPrevNodes), len(currentNode.listOfNextNodes)
    
    searchDepth = 5
    
    in1IndexList, in2IndexList, out1IndexList, out2IndexList, commonIndexList = [] ,[],[],[],[]
    
    for prevNode, i in zip( currentNode.listOfPrevNodes, range(2)):
        #print prevNode
        edgeWt = prevNode[1]
        KmerIndex = prevNode[0].nodeIndexList[-(edgeWt +1)]
        
        if edgeWt+1+flankinglen <= len(prevNode[0].nodeIndexList):
            frankingin = prevNode[0].nodeIndexList[-(edgeWt+1+flankinglen)]
        elif len(prevNode[0].listOfPrevNodes) > 0:
            edgeWt2 = prevNode[0].listOfPrevNodes[0][1]
            maxlen = len(prevNode[0].listOfPrevNodes[0][0].nodeIndexList)
            frankingin = prevNode[0].listOfPrevNodes[0][0].nodeIndexList[-min(edgeWt2 + 1 + edgeWt+flankinglen -len(prevNode[0].nodeIndexList),maxlen)]
        else: 
            frankingin= prevNode[0].nodeIndexList[0]
            
        
        readList = bridgeResolve.obtainReadNum(frankingin, f2)
        if i == 0 : 
            for eachitem in readList:
                in1IndexList.append(eachitem[0]) 
        elif i == 1  : 
            for eachitem in readList:
                in2IndexList.append(eachitem[0]) 

    
    for nextNode, i in zip( currentNode.listOfNextNodes, range(2)):
        #KmerIndex = nextNode[0].nodeIndex
        
        edgeWt = nextNode[1]
        KmerIndex = nextNode[0].nodeIndexList[edgeWt ]
        
        if edgeWt+flankinglen < len(nextNode[0].nodeIndexList) :
            frankingout = nextNode[0].nodeIndexList[edgeWt+flankinglen]
        elif len(nextNode[0].listOfNextNodes) > 0:
            edgeWt2 = nextNode[0].listOfNextNodes[0][1]
            maxlen = len(nextNode[0].listOfNextNodes[0][0].nodeIndexList)
            frankingout = nextNode[0].listOfNextNodes[0][0].nodeIndexList[min(edgeWt2+edgeWt+flankinglen -len(nextNode[0].nodeIndexList), maxlen-1)]
        else: 
            frankingout = nextNode[0].nodeIndexList[-1]
        
        readList = bridgeResolve.obtainReadNum(frankingout, f2)
        if i == 0 : 
            for eachitem in readList:
                out1IndexList.append(eachitem[0]) 
        elif i == 1 :
            for eachitem in readList:
                out2IndexList.append(eachitem[0]) 
    
    lrep= len(currentNode.nodeIndexList) 
    Llong = indelRobot.Llong
    liid = indelRobot.liid
    print "lrep, Llong, liid", lrep, Llong, liid
    
    anchorPoint1 , anchorPoint2 = max ( int(lrep* 0.3 ), liid)  , min ( int(lrep *0.7) , lrep - liid)
    
    print "anchorPoint1, anchorPoint2 ",  anchorPoint1, anchorPoint2 
    
    
    readSet1 = bridgeResolve.findRangeList(currentNode.nodeIndexList[anchorPoint1:anchorPoint1+searchDepth], f2 )
    print "readSet1", readSet1
    
    readSet2 = bridgeResolve.findRangeList(currentNode.nodeIndexList[anchorPoint2: anchorPoint2 + searchDepth], f2)
    print "readSet2", readSet2
    
  
    combinedList = bridgeResolve.distinct(readSet1,"zero") + bridgeResolve.distinct(readSet2, "zero")
    combinedList = sorted(combinedList)
    
    for i in range(len(combinedList) -1 ):
        if combinedList[i][0] == combinedList[i+1][0]:
            commonIndexList.append(combinedList[i][0])
        
    
    commonIndexList = bridgeResolve.distinct(commonIndexList)
    toProcessList =  [in1IndexList, in2IndexList, out1IndexList, out2IndexList, commonIndexList]
    for eachitem in toProcessList:
        print "len(eachitem)" , len(eachitem), eachitem
    
    return toProcessList
    
    
def formRelatedMap(f2, noisyReads, currentNode, indelRobot, toProcessList):
    shortToLongMap = [] 
   
    print "formRelatedMap"
    # shortToLongMap : indexlong    indexshort    jstart    jend    istart    iend   
    
    indexLongList = toProcessList[4]
    indexShortList = toProcessList[0] + toProcessList[1] + toProcessList[2] + toProcessList[3] + toProcessList[4]
    
    for eachlongRead in indexLongList:
        for eachshortRead in indexShortList: 
            if eachlongRead != eachshortRead: 
                
                endNoisy1, endNoisy2 = 0, 0 
                #print noisyReads[eachshortRead]
                while (noisyReads[eachshortRead][endNoisy1] != 0):
                    endNoisy1 += 1
    
                while (noisyReads[eachlongRead][endNoisy2] != 0):
                    endNoisy2 += 1
                
                score, returnalignedSeq1, returnalignedSeq2 , starti, startj , endi, endj = cleaner.SWAlignment(noisyReads[eachshortRead][0:endNoisy1], noisyReads[eachlongRead][0:endNoisy2],indelRobot)
                
                if score > indelRobot.liid :
                    shortToLongMap.append([eachlongRead, eachshortRead,startj,endj,starti,endi])
    
    print "len(shortToLongMap)" , len(shortToLongMap)
    
    return shortToLongMap  


