import graphForm
import logging 

import numpy as np 
import cluster 
from operator import itemgetter

class MBCondensedNode(graphForm.condensedNode):
    xNodeState = False
    naiveForm = True
    def updateXNodeState(self,updateState):
        self.xNodeState = updateState
    
    def checkXNodeState(self):
        return self.xNodeState
        
    def initOverlap(self):
        self.naiveForm = False
        tempList = []

        for eachnextnode, index in zip(self.listOfNextNodes, range(len(self.listOfNextNodes))):
            tempList.append([eachnextnode, 0])
            
        self.listOfNextNodes = tempList
        
        tempList = []
        
        for eachprevnode, index in zip(self.listOfPrevNodes, range(len(self.listOfPrevNodes))):
            tempList.append([eachprevnode, 0])
            
        self.listOfPrevNodes = tempList
        
    def checkPrevList(self):
        print "Check Prev List"
        for eachitem in self.listOfPrevNodes:
            print eachitem[0].nodeIndex
            
    def checkNextList(self):
        print "Check Next List"
        for eachitem in self.listOfNextNodes:
            print eachitem[0].nodeIndex
            
def distinct(oldtmplist, typeOfMode = "all"):
    tmplist = sorted(oldtmplist)
    runningindex = 0 
    thres = 100
    
    while (runningindex < len(tmplist) -1 ):
        if typeOfMode == "all":
            if tmplist[runningindex] == tmplist[runningindex +1]:
                tmplist.pop(runningindex)
            else:
                runningindex = runningindex + 1
                
        elif typeOfMode == "zero" :
            if tmplist[runningindex][0] == tmplist[runningindex +1][0]:
                tmplist.pop(runningindex)
            else:
                runningindex = runningindex + 1
                
        elif typeOfMode == "zeroThres":
            if tmplist[runningindex][0] == tmplist[runningindex +1][0] and abs(tmplist[runningindex][1] - tmplist[runningindex +1][1] ) < thres:
                tmplist.pop(runningindex)
            else:
                runningindex = runningindex + 1
            
            
    return tmplist

def modifyInList(listA, lrepeat):
    newList = []
    
    for index in range(len(listA)):
        newList.append([listA[index][0], listA[index][1] + lrepeat+1])
        
    return newList
    
def isIntersect(listA, listB, lrepeat):
    totalList = modifyInList(listA, lrepeat) + listB
    print "modifyInList(listA, lrepeat) ", modifyInList(listA, lrepeat) 
    print "listB", listB
    totalList = sorted(totalList)
    
    for index in range(len(totalList) -1):
        if totalList[index] == totalList[index +1]:
            return True
            
            
    return False

def combineListAandB(listA, listB):
    combinedList = []
    for eachitem in listA:
        tmp = []
        for eachsubitem in eachitem:
            tmp.append(eachsubitem)
        tmp.append(0)
        combinedList.append(tmp)
        
    for eachitem in listB:
        tmp = []
        for eachsubitem in eachitem:
            tmp.append(eachsubitem)
        tmp.append(1)
        combinedList.append(tmp)
    return combinedList
    

def isIntersectIndel(listA, listB, lrepeat):
   # totalList = combineListAandB(modifyInList(listA, lrepeat), listB) 
    totalList = modifyInList(listA, lrepeat)+ listB 
    
    print "modifyInList(listA, lrepeat) ", modifyInList(listA, lrepeat) 
    print "listB, lrepeat", listB, lrepeat
    totalList = sorted(totalList)
    #threshold = 300
    threshold = 40
    returnList =[] 
    
    for index in range(len(totalList) -1):
        if totalList[index][0] == totalList[index +1][0]  and abs(totalList[index][1]-totalList[index +1][1] ) < threshold:
            returnList.append(totalList[index][0])
    
    if len(returnList) > 0 : 
        returnList = distinct(returnList)
    return returnList

def bisectkk(a, x, lo=0, hi=None):
    """Return the index where to insert item x in list a, assuming a is sorted.

    The return value i is such that all e in a[:i] have e <= x, and all e in
    a[i:] have e > x.  So if x already appears in the list, a.insert(x) will
    insert just after the rightmost x already there.

    Optional args lo (default 0) and hi (default len(a)) bound the
    slice of a to be searched.
    """

    if lo < 0:
        raise ValueError('lo must be non-negative')
    if hi is None:
        hi = len(a)
    while lo < hi:
        mid = (lo+hi)//2
        ### Determine order between x and a[mid]  
        booleanArray = np.equal(x,a[mid])
        #print booleanArray, x, a[mid]
        tempindex = -1
        for index in range(len(booleanArray) -1):
            if booleanArray[index] == False:
                tempindex = index
                break
                
        if x[tempindex] < a[mid][tempindex] and tempindex != -1: hi = mid
        else: lo = mid+1
    return lo

    
def obtainReadNum(KmerIndex, f2):
    readList = []
    startIndex = bisectkk(f2, [KmerIndex, 0 , -1 ,-1, -1]) 
    
    #print "KmerIndex, f2[startIndex-1:startIndex+2]" , KmerIndex, f2[startIndex-1:startIndex+2]
    for index in range(startIndex,len(f2)):
        if f2[index][0] ==  KmerIndex:
            readList.append([f2[index][1], f2[index][2]])
        if f2[index][0] > KmerIndex:
            break
    
    return readList   

def obtainReadDetail(KmerIndex, f2):

    startIndex = bisectkk(f2, [KmerIndex, 0 , -1 ,-1, -1]) 
    endIndex = bisectkk(f2, [KmerIndex+1, 0 , -1 ,-1, -1]) 
    
    return f2[startIndex:endIndex]

    #if KmerIndex == 8002 or KmerIndex == 8001:
    #    print "KmerIndex, f2[startIndex-1:startIndex+2]" , KmerIndex, f2[startIndex-1:startIndex+2]
    #    print len(f2)
    
#    for index in range(startIndex,len(f2)):
#        if f2[index][0] ==  KmerIndex:
#            readList.append(f2[index][1:5])
#        if f2[index][0] > KmerIndex:
#            break
#    
#    return readList 
 
def findRangeList(frankinginList, f2 ):
    inList = []
    frankinginList = sorted(frankinginList)
    for eachfrankin in frankinginList: 
        inList += obtainReadNum(eachfrankin, f2)
    
    inList = sorted(inList)
    inList = distinct(inList, "zero")
    inList = distinct(inList, "zeroThres")
    #inList = distinct(inList, "all")
    #print "inList",inList 
    return inList
 
 
def checkBridgingFromReads(f2, xNode,frankingdepth, noisyReads, parameterRobot): 
    isBridged, matchingList = False, []
    ### Find bridging reads
    bridgingReads = []
    K = 150
    searchdepth = 5
    
    
    
    outList = findRangeList(xNode.nodeIndexList[K: K+searchdepth], f2 )
    inList = findRangeList(xNode.nodeIndexList[frankingdepth:searchdepth + frankingdepth], f2 )
    
    print "my name is Bond", outList
    print inList
    
    bridgingReads = isIntersectIndel(inList, outList, K - frankingdepth)
    
    ### Determine Extension List and isBridged or not
    inputReads = [8 for i in range(len(xNode.listOfPrevNodes))]
    outputReads = [8 for i in range(len(xNode.listOfNextNodes))]
    inMap = []
    outMap = []
    threshold, liid,overhang= parameterRobot.threshold ,parameterRobot.liid ,5 
    speciallen = 40
    
    for i in range(len(inputReads)):
        tempList = obtainReadNum(xNode.listOfPrevNodes[i][0].nodeIndexList[-speciallen], f2)
        tempList = sorted(tempList, key = itemgetter(1))
        inputReads[i] = tempList[0][0]
        print "In ",i , xNode.listOfPrevNodes[i][0].nodeIndexList[-1],inputReads[i], len(xNode.listOfPrevNodes[i][0].nodeIndexList)
    
    for i in range(len(outputReads)):
        tempList = obtainReadNum(xNode.listOfNextNodes[i][0].nodeIndexList[speciallen], f2)
        tempList = sorted(tempList,key = itemgetter(1))
        outputReads[i] = tempList[-1][0]
        print "Out ",i , xNode.listOfNextNodes[i][0].nodeIndexList[0], outputReads[i], len(xNode.listOfNextNodes[i][0].nodeIndexList)
        
    print "hi there", len(inputReads), len(outputReads), len(bridgingReads)
    for prevIndex, dummyindex in zip(inputReads, range(len(inputReads))):
        endNoisy1= 0 
#  print noisyReads[0][0]
#        print prevIndex
#        print endNoisy1
        while (noisyReads[prevIndex][endNoisy1] != 0):
            endNoisy1 += 1

        for bridgeIndex in bridgingReads:     
            endNoisy2 =0 
            while (noisyReads[bridgeIndex][endNoisy2] != 0):
                endNoisy2 += 1        
                
            score, returnalignedSeq1, returnalignedSeq2 , starti, startj , endi, endj = cluster.SWAlignment(noisyReads[prevIndex][0:endNoisy1], noisyReads[bridgeIndex][0:endNoisy2], parameterRobot)

            check = cluster.meetRequirement(score, returnalignedSeq1, returnalignedSeq2 ,starti, startj, endi, endj, threshold, liid,overhang, endNoisy1, endNoisy2 )
            if [prevIndex, bridgeIndex] == [1318, 235]:
                print "strange"
                print " starti, startj , endi, endj",  starti, startj , endi, endj
                cluster.printSeq(returnalignedSeq1)
                cluster.printSeq(returnalignedSeq2) 
                
            if check : 
                inMap.append([prevIndex, bridgeIndex, 0, dummyindex])
           
    
    print "inMap", inMap
    for nextIndex,dummyindex in zip(outputReads,range(len(outputReads))):
        endNoisy1= 0 
        while (noisyReads[nextIndex][endNoisy1] != 0):
            endNoisy1 += 1

        for bridgeIndex in bridgingReads:      
            endNoisy2 =0 
            while (noisyReads[bridgeIndex][endNoisy2] != 0):
                endNoisy2 += 1         
                
            score, returnalignedSeq1, returnalignedSeq2 , starti, startj , endi, endj = cluster.SWAlignment(noisyReads[nextIndex][0:endNoisy1], noisyReads[bridgeIndex][0:endNoisy2], parameterRobot)
            check = cluster.meetRequirement(score, returnalignedSeq1, returnalignedSeq2 ,starti, startj, endi, endj, threshold, liid,overhang, endNoisy1, endNoisy2)
            if [nextIndex, bridgeIndex] == [235, 235]:
                print "strange 2 "
                print " starti, startj , endi, endj",  starti, startj , endi, endj
                cluster.printSeq(returnalignedSeq1)
                cluster.printSeq(returnalignedSeq2)
            if check : 
                outMap.append([nextIndex, bridgeIndex, 1, dummyindex])
                
    print "outMap,", outMap 
                
    print "matchingList", matchingList      
    for eachinmap in inMap:
        for eachoutmap in outMap:
            if eachinmap[1] == eachoutmap[1] :
                dummyIn = eachinmap[3]
                dummyOut = eachoutmap[3]
                print "dummyIn, dummyOut", dummyIn, dummyOut
                #print "eachinmap", eachinmap
                #print "eachoutmap",eachoutmap
                
                edgeWt = xNode.listOfPrevNodes[dummyIn][1]
                inKmerIndex = xNode.listOfPrevNodes[dummyIn][0].nodeIndexList[-(edgeWt+1 )]
                edgeWt = xNode.listOfNextNodes[dummyOut][1]
                outKmerIndex =  xNode.listOfNextNodes[dummyOut][0].nodeIndexList[edgeWt]
                
                matchingList.append([inKmerIndex,outKmerIndex ])
    
    ### Check isBridged                 
    if len(matchingList) > 0 :          
        matchingList = distinct(matchingList)
        reverseList = []
        for eachitem in matchingList:
            reverseList.append([eachitem[1], eachitem[0]])
        
        if len(distinct(matchingList, "zero")) == len(inputReads) and len(distinct(reverseList, "zero")) == len(outputReads):
            isBridged = True
        else:
            isBridged = False
    
    #matchingList = [[28346, 32138], [28346, 4913], [30233, 32138]]
    
    
    return isBridged, matchingList 


def checkBridging(f2, xNode,frankingdepth):
    isBridged = True
    matchingList = []
    searchDepth = 5
    #frankingdepth = 1
    
    lenRepeat = len(xNode.nodeIndexList)
    
    for eachinnode in xNode.listOfPrevNodes:
        edgeWt = eachinnode[1]
        
        print  "len(eachinnode) ",len(eachinnode),  edgeWt, len(eachinnode[0].nodeIndexList) , xNode.nodeIndex
        inKmerIndex = eachinnode[0].nodeIndexList[-(edgeWt+1 )]
        frankinginList = []
        if edgeWt+1+frankingdepth <= len(eachinnode[0].nodeIndexList):
            frankinginList = eachinnode[0].nodeIndexList[-(edgeWt+1+frankingdepth):-(edgeWt+1+frankingdepth) + searchDepth]
            #frankinginList = eachinnode[0].nodeIndexList[-(edgeWt+1+frankingdepth):-1]
            #print "his is run",  len(frankinginList), (edgeWt+1+frankingdepth), len(eachinnode[0].nodeIndexList)
        elif len(eachinnode[0].listOfPrevNodes) > 0:
            edgeWt2 = eachinnode[0].listOfPrevNodes[0][1]
            maxlen = len(eachinnode[0].listOfPrevNodes[0][0].nodeIndexList)
            frankinginList =  eachinnode[0].listOfPrevNodes[0][0].nodeIndexList[-min(edgeWt2 + 1 + edgeWt+frankingdepth -len(eachinnode[0].nodeIndexList),maxlen):-min(edgeWt2 + 1 + edgeWt+frankingdepth -len(eachinnode[0].nodeIndexList),maxlen) +searchDepth ]
        else: 
            frankinginList= [eachinnode[0].nodeIndexList[0] ] 
        
        print "inKmerIndex", inKmerIndex      
        
        print "frankinginList", frankinginList  
        inList = findRangeList(frankinginList, f2 )
        print "inList", inList

 
        testBridged = False
        for eachoutnode in xNode.listOfNextNodes:
            edgeWt = eachoutnode[1]
            outKmerIndex = eachoutnode[0].nodeIndexList[edgeWt]
            frankingoutList = []
            
            if edgeWt+frankingdepth < len(eachoutnode[0].nodeIndexList) :
                frankingoutList = eachoutnode[0].nodeIndexList[edgeWt+frankingdepth- searchDepth: edgeWt+frankingdepth]
                #frankingoutList = eachoutnode[0].nodeIndexList[0: edgeWt+frankingdepth]
            elif len(eachinnode[0].listOfNextNodes) > 0:
                edgeWt2 = eachoutnode[0].listOfNextNodes[0][1]
                maxlen = len(eachoutnode[0].listOfNextNodes[0][0].nodeIndexList)
                frankingoutList = eachoutnode[0].listOfNextNodes[0][0].nodeIndexList[min(edgeWt2+edgeWt+frankingdepth -len(eachoutnode[0].nodeIndexList), maxlen-1)-searchDepth:min(edgeWt2+edgeWt+frankingdepth -len(eachoutnode[0].nodeIndexList), maxlen-1)]
            else: 
                frankingoutList = [eachoutnode[0].nodeIndexList[-1]]
                
                
            print "outKmerIndex", outKmerIndex,len(eachoutnode[0].nodeIndexList)
            
            print "frankingoutList", frankingoutList
            outList = findRangeList(frankingoutList, f2 )
            print "outList",outList
            tmpRdList = isIntersectIndel(inList, outList, lenRepeat+frankingdepth*2)
            if len(tmpRdList) != 0:
                testBridged = True
                #print "[inKmerIndex, outKmerIndex]", [inKmerIndex, outKmerIndex], tmpRdList
                matchingList.append([inKmerIndex, outKmerIndex])
                
        
        if testBridged == False:
            isBridged = False
            #print "Fhere"
            

    for eachoutnode in xNode.listOfNextNodes:
        edgeWt = eachoutnode[1]
        outKmerIndex = eachoutnode[0].nodeIndexList[edgeWt]
        frankingoutList = []
        
        if edgeWt+frankingdepth < len(eachoutnode[0].nodeIndexList) :
            frankingoutList = eachoutnode[0].nodeIndexList[edgeWt+frankingdepth- searchDepth: edgeWt+frankingdepth]
        elif len(eachinnode[0].listOfNextNodes) > 0:
            edgeWt2 = eachoutnode[0].listOfNextNodes[0][1]
            maxlen = len(eachoutnode[0].listOfNextNodes[0][0].nodeIndexList)
            frankingoutList = eachoutnode[0].listOfNextNodes[0][0].nodeIndexList[min(edgeWt2+edgeWt+frankingdepth -len(eachoutnode[0].nodeIndexList), maxlen-1)-searchDepth:min(edgeWt2+edgeWt+frankingdepth -len(eachoutnode[0].nodeIndexList), maxlen-1)]
        else: 
            frankingoutList = [ eachoutnode[0].nodeIndexList[-1] ]
            
            
        #print "outKmerIndex", outKmerIndex
        outList = findRangeList(frankingoutList, f2 )
        #print "outList",outList
                
        testBridged = False
        for eachinnode in xNode.listOfPrevNodes:
            edgeWt = eachinnode[1]
            #print  "len(eachinnode) ",len(eachinnode),  edgeWt, len(eachinnode[0].nodeIndexList) , xNode.nodeIndex
            inKmerIndex = eachinnode[0].nodeIndexList[-(edgeWt+1 )]
            frankinginList = []
            if edgeWt+1+frankingdepth <= len(eachinnode[0].nodeIndexList):
                frankinginList = eachinnode[0].nodeIndexList[-(edgeWt+1+frankingdepth):-(edgeWt+1+frankingdepth) + searchDepth]
                #print len(frankinginList)
            elif len(eachinnode[0].listOfPrevNodes) > 0:
                edgeWt2 = eachinnode[0].listOfPrevNodes[0][1]
                maxlen = len(eachinnode[0].listOfPrevNodes[0][0].nodeIndexList)
                frankinginList =  eachinnode[0].listOfPrevNodes[0][0].nodeIndexList[-min(edgeWt2 + 1 + edgeWt+frankingdepth -len(eachinnode[0].nodeIndexList),maxlen):-min(edgeWt2 + 1 + edgeWt+frankingdepth -len(eachinnode[0].nodeIndexList),maxlen) +searchDepth ]
            else: 
                frankinginList= [ eachinnode[0].nodeIndexList[0] ]

            inList = findRangeList(frankinginList, f2 )
            
            #print "inKmerIndex", inKmerIndex
            
            if len(isIntersectIndel(inList, outList, lenRepeat+frankingdepth*2)) != 0:
                testBridged = True
                #print "Yes ",outKmerIndex
                
        if testBridged == False:
            isBridged = False

    print "isBridgeddd", isBridged
    
#    if not isBridged:
#        isBridged, matchingList = checkBridgingFromReads(f2, xNode,frankingdepth, noisyReads, parameterRobot)
   
    return isBridged, matchingList

def makeDistinctList(tmpList):
    tmpList2 = []
    tmpList2 = sorted(tmpList)
    
    runningindex = 0
    while (runningindex < len(tmpList2)-1):
        if tmpList2[runningindex] == tmpList2[runningindex+1]:
            tmpList2.pop(runningindex)
        else:
            runningindex = runningindex + 1
            
    return tmpList2



def resolveFramework(currentNode, f2,startingGpNum, G3, canResolve, kmerPairsList):
    
    
    print "currentNode.nodeIndex", currentNode.nodeIndex
    print "--------------"
    newXnodes = []            
    if canResolve:
        
        
        print "startingGpNum, len(currentNode.listOfNextNodes), len(currentNode.listOfPrevNodes)", startingGpNum, len(currentNode.listOfNextNodes), len(currentNode.listOfPrevNodes)
        #currentNode.checkPrevList()
        #currentNode.checkNextList()

        #i) create nodes for I/O
        selfLoopinNode = -1 
        selfLoopoutNode = -1 
        selfLoop = False
        selfLoopWt = -1 
        inNodes = []
        outNodes = []   
        
        oldInNodes = []
        oldOutNodes = []
        
        for eachtup in currentNode.listOfPrevNodes:
            oldInNodes.append(eachtup[0])
            
        for eachtup in currentNode.listOfNextNodes:
            oldOutNodes.append(eachtup[0])
        
        lenRepeat = len(currentNode.nodeIndexList)
        
        for eachprevnode in currentNode.listOfPrevNodes:
            # create node and update edge weights
            edgeWt = eachprevnode[1]
            pNode = eachprevnode[0]
            
            newNode = MBCondensedNode(startingGpNum)        
            newNode.nodeIndexList = [pNode.nodeIndexList[-(edgeWt+1)]] + currentNode.nodeIndexList 

            startingGpNum = startingGpNum +1 
            

            
            inNodes.append(newNode)            
            G3.append(newNode)
            
            if pNode == currentNode:
                selfLoop = True
                print "selfLoop = True : pNode", pNode.nodeIndex
                selfLoopinNode = newNode
                selfLoopWt = edgeWt
            else:
                print "Prev node", pNode.nodeIndex            
                pNode.listOfNextNodes.append([newNode, edgeWt +1])
                newNode.listOfPrevNodes.append([pNode, edgeWt +1 ])
        
        print "startingGpNum, len(currentNode.listOfNextNodes), len(currentNode.listOfPrevNodes)", startingGpNum, len(currentNode.listOfNextNodes), len(currentNode.listOfPrevNodes)
            
        for eachnextnode in currentNode.listOfNextNodes:
            edgeWt = eachnextnode[1]
            nNode = eachnextnode[0]
            
            newNode = MBCondensedNode(startingGpNum)      
            newNode.nodeIndexList = currentNode.nodeIndexList + [nNode.nodeIndexList[edgeWt]]

            startingGpNum = startingGpNum +1 
                
            outNodes.append(newNode)
            G3.append(newNode)
            
            if nNode == currentNode:
                selfLoop = True
                print "selfLoop = True : nNode", selfLoop  , nNode.nodeIndex
                
                selfLoopoutNode = newNode
                selfLoopWt = edgeWt
            else:
                print "Next node", nNode.nodeIndex
                nNode.listOfPrevNodes.append([newNode, edgeWt +1])
                newNode.listOfNextNodes.append([nNode, edgeWt +1 ])
        
        
        #ii) handle self-loop
        if selfLoop == True:
            print "self loop: selfLoopWt", selfLoopWt
            selfLoopoutNode.listOfNextNodes.append([selfLoopinNode,selfLoopWt+2] )
            selfLoopinNode.listOfPrevNodes.append([selfLoopoutNode, selfLoopWt +2 ])

        #iii) remove old nodes and edges
        
        for eachprevtup in currentNode.listOfPrevNodes:
            removalIndex = None
            eachprev = eachprevtup[0]
            
            #print "Checking for node : ", eachprev.nodeIndex
            #eachprev.checkNextList()
            #print "--------------"
            
            for eachnode, index in zip(eachprev.listOfNextNodes, range(len(eachprev.listOfNextNodes))):
                if eachnode[0] == currentNode:
                    removalIndex = index
                    
            eachprev.listOfNextNodes.pop(removalIndex)
            
        for eachnexttup in currentNode.listOfNextNodes:
            removalIndex = None
            eachnext = eachnexttup[0]
            
            #print "Checking for node : ", eachnext.nodeIndex
            #eachnext.checkPrevList()
            #print "--------------"
            
            for eachnode, index in zip(eachnext.listOfPrevNodes, range(len(eachnext.listOfPrevNodes))):
                if eachnode[0] == currentNode:
                    removalIndex = index
              
            
            eachnext.listOfPrevNodes.pop(removalIndex)          
            
        
        currentNode.listOfPrevNodes = []
        currentNode.listOfNextNodes = []
        currentNode.nodeIndexList = []
        
        #iv) use read to add extra edges 
        for eachInNode in inNodes:
            inKmerIndex = eachInNode.nodeIndexList[0]

            for eachOutNode in outNodes:
                outKmerIndex = eachOutNode.nodeIndexList[len(eachOutNode.nodeIndexList)-1]
                #print inKmerIndex, outKmerIndex
                temppair = [inKmerIndex, outKmerIndex]
                
                if temppair in kmerPairsList:

                    print "Found"
                    eachInNode.listOfNextNodes.append([eachOutNode, lenRepeat])
                    eachOutNode.listOfPrevNodes.append([eachInNode, lenRepeat])
                
        
        
        #v) locally condense
        print "Condensing "
        nodesToUpdate = inNodes + outNodes + oldInNodes +  oldOutNodes

        nodesToUpdate= distinct(nodesToUpdate)
        #Debug 
        for eachtestnode in nodesToUpdate:
            print "eachtestnode.nodeIndex, len(eachtestnode.nodeIndexList), len(eachtestnode.listOfPrevNodes), len(eachtestnode.listOfNextNodes)", eachtestnode.nodeIndex, len(eachtestnode.nodeIndexList), len(eachtestnode.listOfPrevNodes), len(eachtestnode.listOfNextNodes)

            for eachsubnode in eachtestnode.listOfPrevNodes:
                print "Prev : ", eachsubnode[0].nodeIndex, eachsubnode[1]

            for eachsubnode in eachtestnode.listOfNextNodes:
                print "Next : ", eachsubnode[0].nodeIndex, eachsubnode[1]
          
        # End Debug
                
        for eachnode in nodesToUpdate:
            
            if len(eachnode.listOfNextNodes) == 1:
                targetNextNode = eachnode.listOfNextNodes[0][0]
                
                if len(targetNextNode.listOfPrevNodes) == 1:
                    edgeWt = eachnode.listOfNextNodes[0][1]
                    #print "node index" , eachnode.nodeIndex   , targetNextNode.nodeIndex, len(eachnode.nodeIndexList)
                    targetNextNode.nodeIndexList = eachnode.nodeIndexList + targetNextNode.nodeIndexList[edgeWt:len(targetNextNode.nodeIndexList)]
                    
                    for eachprevnodetup in eachnode.listOfPrevNodes:
                        eachprevnodetup[0].listOfNextNodes.append([ targetNextNode,eachprevnodetup[1] ] )
                        
                        tmpremovalindex = None
                                                 
                        
                        for removalnodeindex in range(len(eachprevnodetup[0].listOfNextNodes)):
                           # print "1" , eachprevnodetup[0].nodeIndex, eachprevnodetup[0].listOfNextNodes[removalnodeindex][0].nodeIndex 
                            if eachnode == eachprevnodetup[0].listOfNextNodes[removalnodeindex][0]:
                                tmpremovalindex =removalnodeindex
                        
                        print "tmpremovalindex",tmpremovalindex
                        if tmpremovalindex != None:
                            eachprevnodetup[0].listOfNextNodes.pop(tmpremovalindex)
                        
                        targetNextNode.listOfPrevNodes.append([eachprevnodetup[0],eachprevnodetup[1]])
                        
                        tmpremovalindex = None
                        for removalnodeindex in range(len(targetNextNode.listOfPrevNodes)):
                           # print "2", targetNextNode.nodeIndex,  targetNextNode.listOfPrevNodes[removalnodeindex][0].nodeIndex
                            if eachnode == targetNextNode.listOfPrevNodes[removalnodeindex][0]:
                                tmpremovalindex =removalnodeindex
                                
                        if tmpremovalindex != None:
                            targetNextNode.listOfPrevNodes.pop(tmpremovalindex)
                        
                    
                    if targetNextNode != eachnode:
                        eachnode.nodeIndexList = []
                        eachnode.listOfPrevNodes =[]
                        eachnode.listOfNextNodes = []
                
        print "Done Condensing "
        
        #vi) update X-node info
        for eachnode in nodesToUpdate:
            if len(eachnode.listOfNextNodes) >1 and len(eachnode.listOfPrevNodes) > 1:
                currentState = True
            else:
                currentState = False
                
            if eachnode.checkXNodeState == False and currentState == True:
                eachnode.updateXNodeState(currentState)
                newXnodes.append(eachnode)
    
    return newXnodes,startingGpNum

def removeEmptyNodes(G3):
    runningindex = 0
    while (runningindex < len(G3)):
        if len(G3[runningindex].nodeIndexList) == 0:
            G3.pop(runningindex)
        else:
            runningindex = runningindex + 1 
            
    
    print "------------------"
    print "index, prevList, nextList, len(nodeindexList)"
    for eachnode in G3:
        print "---------------"
        print eachnode.nodeIndex, len(eachnode.listOfPrevNodes), len(eachnode.listOfNextNodes), len(eachnode.nodeIndexList)
        print "Prev"
        for eachitem in eachnode.listOfPrevNodes:
            print eachitem[0].nodeIndex, eachitem[1]
        print "Next"
        for eachitem in eachnode.listOfNextNodes:
            print eachitem[0].nodeIndex, eachitem[1] 
        print "---------------"
    
    
    

def resolveRepeats(f2, G2,parameterRobot):
    G3 = []
    # f2 = "Group number ","read-id", "offset-id", "fused or not", "prevGroupid"
    # G2 = CondensedNode
    bridgeLimit = parameterRobot.bridgingDepth
    print "bridgeLimit", bridgeLimit
    
    f2 = sorted(f2)
    

    startingGpNum = f2[-1][0] + 1
    print "startingGpNum",startingGpNum
    
    #1 Init Nodes and X-node List
    basicList,seqList, typeOfGraph = [], [], "MB"
    #a Basic List 
    for eachnode in G2:
        nextList, prevList = [], [] 
        
        for eachnext in eachnode.listOfNextNodes:
            nextList.append(eachnext.nodeIndex)
            
        for eachprev in eachnode.listOfPrevNodes:
            prevList.append(eachprev.nodeIndex)
            
        basicList.append([int(eachnode.nodeIndex),prevList,nextList, int(len(eachnode.nodeIndexList))])
        
    #b Seq List 
    for eachnode in G2:
        seqList.append([eachnode.nodeIndex, eachnode.nodeIndexList])
        
    #c Transofrm 
    basicList= makeDistinctList(basicList)
    seqList = makeDistinctList(seqList)
    
    G2 = logging.transformToMBGraph(basicList,seqList, typeOfGraph)

    #d check
    print "-----------------"
    for eachnode in G2:
        
        print "Current node", eachnode.nodeIndex
        print "Next List "
        for eachnextnode in eachnode.listOfNextNodes:
            print eachnextnode[0].nodeIndex, eachnextnode[1]
            
        print "Prev List "
        for eachprevnode in eachnode.listOfPrevNodes:
            print eachprevnode[0].nodeIndex, eachprevnode[1]
        
        print "-----------------"
    
    for eachitem in G2:
        G3.append(eachitem)

        
    #e X-nodeList 
    xNodesList = []   
    for eachnode in G2:
        if len(eachnode.listOfNextNodes) > 1 and len(eachnode.listOfPrevNodes) > 1:
            eachnode.updateXNodeState(True)
            xNodesList.append(eachnode)
        else:
            eachnode.updateXNodeState(False)
        
    
    #2 Resolving Repeats
    print "Resolving Repeats: "
    print "======================"
    round = 0
    while (len(xNodesList) > 0  ):

        currentNode = xNodesList.pop(0)
        
        #print "currentNode.nodeIndex",currentNode.nodeIndex
        
        canResolve, kmerPairsList =  checkBridging(f2, currentNode,bridgeLimit)
        print "kmerPairsList", kmerPairsList
        newXNodeList, startingGpNum =  resolveFramework(currentNode, f2,startingGpNum, G3, canResolve, kmerPairsList )
        xNodesList = xNodesList + newXNodeList
        
        # Hack to get back the selfloop update Xnodes
        #assert(1==2)
        print "Before : len(xNodesList), len(G3)",len(xNodesList), len(G3)
        removeEmptyNodes(G3)
        print "After : len(xNodesList), len(G3)",len(xNodesList), len(G3)
        if len(xNodesList) == 0:
            for eachnode in G3:
                if len(eachnode.listOfPrevNodes) > 0 and len(eachnode.listOfNextNodes) > 0 and len(eachnode.nodeIndexList) > 0:
                    canResolve, kmerPairsList = checkBridging(f2, eachnode,bridgeLimit)
                    if canResolve and not eachnode in xNodesList and  len(eachnode.nodeIndexList)> 0 and len(eachnode.listOfNextNodes) > 1 and len(eachnode.listOfPrevNodes) >1 :
                        xNodesList.append(eachnode) 
                    
        round = round  +1 
    #3 Format Output Graph and output the weighted sequence graph 
    
    removeEmptyNodes(G3)
    print "len(G3)", len(G3)
    return G3
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    