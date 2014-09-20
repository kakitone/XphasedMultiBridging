from operator import itemgetter 
import cluster
import bisect

class condensedNode(object):  
    def __init__(self, index):    
        self.nodeIndex = index
        self.listOfNextNodes = []
        self.listOfPrevNodes = []
        self.visited = False
        self.nodeIndexList = []
        self.numberOfVisit= -1
        self.isAnchored = False
        
        
    def updateNodeList(self):
        self.nodeIndexList.append( self.nodeIndex)
    def visitedNode(self):
        self.visited = True
    def addPrevNodes(self, node):
        self.listOfPrevNodes.append(node)
    def addNextNodes(self, node):
        self.listOfNextNodes.append(node)
    
def  condenseGraph(seqGraphNodes):  
    simGraph, startList  = [], []
    #print "len(seqGraphNodes)",len(seqGraphNodes)
    # Init queue
    queue = []
    for eachitem in seqGraphNodes:
        if len(eachitem.listOfPrevNodes) == 0:
            queue.append(eachitem)
    
    if len(queue) ==0:
        queue.append(seqGraphNodes[0])

    # Init visited state
    for eachnode in seqGraphNodes:
        eachnode.visited = False
    
    # BFS
    while (len(queue) >= 1):

        # a : BFS routine
        currentNode = queue.pop(0)
        currentNode.visited = True
        
        #print "currentNode.nodeIndex",  currentNode.nodeIndex
        
        for eachnode in currentNode.listOfNextNodes:
            
            if eachnode.visited == False and eachnode != currentNode:
                queue.append(eachnode)
                
        # b : Condensing steps     
        #if currentNode.nodeIndex == 85:
        #    print "Found you "   
            
        if len(currentNode.listOfNextNodes) == 1:
            targetnextnode = currentNode.listOfNextNodes[0]

            if len(targetnextnode.listOfPrevNodes) == 1 and targetnextnode != currentNode:

                targetnextnode.nodeIndexList = currentNode.nodeIndexList + targetnextnode.nodeIndexList          
                currentNode.nodeIndexList = []

                
                for eachprevnode in currentNode.listOfPrevNodes:
                    
                    eachprevnode.listOfNextNodes.append(targetnextnode)
                    targetnextnode.listOfPrevNodes.append(eachprevnode)
                    
                    eachprevnode.listOfNextNodes.remove(currentNode)

                targetnextnode.listOfPrevNodes.remove(currentNode)                    
                currentNode.listOfPrevNodes = []
                currentNode.listOfNextNodes = []
                
                
            
             
    # Simplify the graph
    #print "Simplify graph"
    counter = 0
    
    for eachnode in seqGraphNodes:

            
        if len(eachnode.nodeIndexList) != 0:
            counter = counter + len(eachnode.nodeIndexList)
            simGraph.append(eachnode)
           #if len(eachnode.nodeIndexList) == 1 and len(eachnode.listOfNextNodes) > 0  :
           #     print "Strange", len(eachnode.listOfNextNodes), len(eachnode.listOfNextNodes[0].listOfPrevNodes) 
            if len(eachnode.listOfPrevNodes) == 0:
                startList.append(eachnode)
                
                
    if len(startList) == 0:
        #print "circular"
        startList.append(simGraph[0])
    
    #print "counter", counter
    return simGraph, startList  

def newCondensingStep(seqGraph):
    simGraph, startList = [] ,[] 
    toCondenseList = []

    for eachnode in seqGraph:

        if len(eachnode.listOfNextNodes) == 1:
            if len(eachnode.listOfNextNodes[0].listOfPrevNodes) == 1 : 
                toCondenseList.append(eachnode)
                
    
    print "len(toCondenseList)", len(toCondenseList)
    for currentNode in toCondenseList:
        if len(currentNode.listOfNextNodes) > 0 :
            targetnextnode = currentNode.listOfNextNodes[0]
    
            if len(targetnextnode.listOfPrevNodes) == 1 and targetnextnode != currentNode:
    
                targetnextnode.nodeIndexList = currentNode.nodeIndexList + targetnextnode.nodeIndexList          
                currentNode.nodeIndexList = []
    
                
                for eachprevnode in currentNode.listOfPrevNodes:
                    
                    eachprevnode.listOfNextNodes.append(targetnextnode)
                    targetnextnode.listOfPrevNodes.append(eachprevnode)
                    eachprevnode.listOfNextNodes.remove(currentNode)
    
                targetnextnode.listOfPrevNodes.remove(currentNode)                    
                currentNode.listOfPrevNodes = []
                currentNode.listOfNextNodes = []
            
            
    for eachnode in seqGraph:
        if len(eachnode.nodeIndexList) != 0:
            simGraph.append(eachnode)
            if len(eachnode.listOfPrevNodes) == 0:
                startList.append(eachnode)
                
                
    if len(startList) == 0:
        startList.append(simGraph[0])
                

    return simGraph, startList







def linkComponents(seqGraphNodes,fmapping,noisyReads, liid, threshold, K ):
    newSeqGraphNodes = []
    dummyNodesList =[]
    
    dummyNodesMapping =[]
    
    noPrevList =[]
    noNextList =[]
    
    N = len(noisyReads)
    L = len(noisyReads[0])
    

    fmapping = sorted(fmapping)
    
    tempgpid = fmapping[-1][0] + 1
    ### Find no prev/ no next list 
    for eachnode in seqGraphNodes:
        if len(eachnode.listOfPrevNodes) ==0  :
            noPrevList.append(eachnode)
            print "Noprev", eachnode.nodeIndex
        if len(eachnode.listOfNextNodes) ==0:
            noNextList.append(eachnode)
            print "Nonext", eachnode.nodeIndex

    print "len(noPrevList), len(noNextList) ",len(noPrevList), len(noNextList)
    ### Find Best match and add dummy nodes 
    for eachnoprevnode in noPrevList:
        tempNode = []
        tempscore = -1 
        currentKmerIndex = eachnoprevnode.nodeIndex
        findex = bisect.bisect_left(fmapping, [currentKmerIndex, -1,-1])
        groupid, readNum , offset =fmapping[findex]
        #print "No Prev : groupid, readNum , offset ,eachnoprevnode.nodeIndex",groupid, readNum , offset, eachnoprevnode.nodeIndex

        currentKmer = noisyReads[readNum]
        
        done = False
        for score in range(K-1, liid, -1):
            for eachPossibleprev in noNextList:
                prevKmerIndex = eachPossibleprev.nodeIndex
                findex = bisect.bisect_left(fmapping, [prevKmerIndex, -1,-1])
                groupid, readNum , offset =fmapping[findex]
                print "NO next: groupid, readNum , offset ,eachnoprevnode.nodeIndex",groupid, readNum , offset, eachnoprevnode.nodeIndex

                
                prevKmer = noisyReads[readNum]
                
                isMatched = False
                
                str1 = currentKmer[0:score]
                str2 = prevKmer[L-score:L]
                
                isMatched = cluster.matched(str1, str2, threshold, liid)
                if isMatched:
                    ### Link together
                    # @ nodes level and in fmapping level
                    print "Matched", score
                    for tempoffset in range(L-K+1, L-score):
                        dummyNodesMapping.append([tempgpid,readNum,tempoffset ])
                        eachPossibleprev.nodeIndexList = eachPossibleprev.nodeIndexList + [tempgpid]
                        tempgpid = tempgpid + 1
                    
                    print len(eachPossibleprev.nodeIndexList) 
                    eachPossibleprev.listOfNextNodes.append(eachnoprevnode)
                    eachnoprevnode.listOfPrevNodes.append(eachPossibleprev)

                        
                    done = True    
                    break
            if done :
                break

    newSeqGraphNodes =seqGraphNodes 
    return newSeqGraphNodes, dummyNodesMapping


def lowDensityRemove(simGraph, startList, fmapping):
    countList = []
    currentCount = 1 
    runningindex = 0 
    tooLowThres =1
    fmapping = sorted(fmapping)
    
    while (runningindex < len(fmapping) -1 ):
        if fmapping[runningindex][0] != fmapping[runningindex+1][0] :
            countList.append([fmapping[runningindex][0], currentCount])
            currentCount = 1       
        else:
            currentCount += 1 
            
        runningindex += 1  
        
    countList.append([fmapping[runningindex][0], currentCount])
    
    #print "countList",  countList
    removeList = []
    for eachitem in countList:
        if eachitem[1] < tooLowThres: 
            removeList.append(eachitem[0])
    
    for eachitem in simGraph:
        if eachitem.nodeIndex in removeList:
            eachitem.nodeIndexList = []
            for eachsubitem in eachitem.listOfNextNodes:
                eachsubitem.listOfPrevNodes.remove(eachitem)
            for eachsubitem in eachitem.listOfPrevNodes:
                eachsubitem.listOfNextNodes.remove(eachitem)
            eachitem.listOfNextNodes =[] 
            eachitem.listOfPrevNodes = []
    
    dummyindex = 0 

    returnGraph = []
    while ( dummyindex < len(simGraph)):
        if len(simGraph[dummyindex].nodeIndexList) == 0 or (len(simGraph[dummyindex].listOfNextNodes) ==0 and len(simGraph[dummyindex].listOfPrevNodes) == 0):
            dummyindex += 1 
        else:
            returnGraph.append(simGraph[dummyindex])
            dummyindex += 1 
    
    
    newstartlist = [] 
    for eachitem in returnGraph:
        if len(eachitem.listOfPrevNodes) == 0 : 
            newstartlist.append(eachitem)
            
    return returnGraph, newstartlist


def transitiveReduction(simGraph):
    depthOfSearch = 5
    print "Transitive reduction"
    for eachnode in simGraph:
        eachnode.visited = False

    ### Process using forward edge
    for eachnode in simGraph:
        # Find One Hop
        oneHopList = []
        
        
        for eachnextnode in eachnode.listOfNextNodes:
            oneHopList.append(eachnextnode.nodeIndex)
        
        #print "oneHopList", oneHopList
        # Find transitive edges by local search 
        for eachnexthop in oneHopList:
            # Remove transitive edges
            proxyNode = []
            for nnode in eachnode.listOfNextNodes:
                if  nnode.nodeIndex == eachnexthop:
                    proxyNode.append(nnode)
                    nnode.listOfPrevNodes.remove(eachnode)
                    eachnode.listOfNextNodes.remove(nnode)
            
            
            toRemove = False
            tempStack = []
            for eachnicenode in eachnode.listOfNextNodes:
                if eachnicenode.nodeIndex != eachnexthop and len(eachnicenode.nodeIndexList) < depthOfSearch:
                    tempStack.append(eachnicenode)
            
            for i in range(depthOfSearch):
                newStack = []
                
                for eachitem in tempStack:
                    for eachnewnext in eachitem.listOfNextNodes:

                        if eachnewnext.nodeIndex == eachnexthop:
                            toRemove = True
                            
                        if  len(eachnewnext.nodeIndexList) < depthOfSearch:
                            newStack.append(eachnewnext)
                
                tempStack = []
                for eachmyitem in newStack:
                    tempStack.append(eachmyitem)
            
            # Add back if nothing wrong
            if toRemove == False:
                proxyNode[0].listOfPrevNodes.append(eachnode)
                eachnode.listOfNextNodes.append(proxyNode[0])
    
    

    
    newstartList = []
    for eachitem in simGraph:
        if len(eachitem.listOfPrevNodes) == 0:
            newstartList.append(eachitem)
                
    for eachnode in simGraph:
        eachnode.visited = False
        
    returnGraph = []
    for eachnode in simGraph:
        returnGraph.append(eachnode)
    
    return returnGraph, newstartList
 
def endRemoval(simGraph):
    print "End Removal"
    endRemThres = 8
    for eachitem in simGraph:

        if len(eachitem.nodeIndexList) <= endRemThres: 
            if  len(eachitem.listOfPrevNodes) == 0 :
                for eachsubitem in eachitem.listOfNextNodes:
                    eachsubitem.listOfPrevNodes.remove(eachitem)
                eachitem.listOfNextNodes = []
                eachitem.nodeIndexList = []
                
                 
            elif len(eachitem.listOfNextNodes) == 0   :
                for eachsubitem in eachitem.listOfPrevNodes:
                    eachsubitem.listOfNextNodes.remove(eachitem)
                eachitem.listOfPrevNodes = []    
                eachitem.nodeIndexList = []

    newstartList = []
    runningindex = 0
    returnGraph = []
    while ( runningindex < len(simGraph)):
        eachitem = simGraph[runningindex]
        
        if len(eachitem.nodeIndexList) == 0:
            runningindex += 1 
        else:
            returnGraph.append(eachitem)
            runningindex += 1 
    
    
    
    for eachitem in returnGraph:
        if len(eachitem.listOfPrevNodes)== 0 :
            newstartList.append(eachitem)
    
    return    returnGraph, newstartList 

    
    

def getSeqGraph(fmapping,noisyReads, parameterRobot):
    seqGraphNodes = []
    print "getSeqGraph"
    liid, threshold, K = parameterRobot.liid, parameterRobot.threshold, parameterRobot.K
    folderName = parameterRobot.defaultFolder
    
    # initialized nodes
    fmapping = sorted(fmapping, key = itemgetter(0)) 
    totalGroupNumber = fmapping[-1][0] + 1 
    
    for index in range(totalGroupNumber):   
        tempNode = condensedNode(index)
        tempNode.updateNodeList()
        seqGraphNodes.append(tempNode)

    # use reads to connect nodes    
    fmapping = sorted(fmapping, key = itemgetter(slice(1,3)))

    for index in range(len(fmapping)-1):        
        currentgroupindex = fmapping[index][0]
        nextgroupindex = fmapping[index+1][0]
        
        currentreadindex = fmapping[index][1]
        nextreadindex = fmapping[index+1][1]

        currentoffsetindex = fmapping[index][2]
        nextoffsetindex = fmapping[index+1][2]
    
        if currentreadindex == nextreadindex and currentoffsetindex +1 == nextoffsetindex:
            ### Update nodes 
            if not seqGraphNodes[nextgroupindex] in seqGraphNodes[currentgroupindex].listOfNextNodes:
                seqGraphNodes[currentgroupindex].addNextNodes(seqGraphNodes[nextgroupindex])
                                
            if not seqGraphNodes[currentgroupindex] in seqGraphNodes[nextgroupindex].listOfPrevNodes:
                seqGraphNodes[nextgroupindex].addPrevNodes(seqGraphNodes[currentgroupindex])
    
    # Add Dummy nodes
    newSeqGraphNodes = []
    newSeqGraphNodes, dummyNodesMapping = linkComponents(seqGraphNodes,fmapping, noisyReads, liid, threshold, K)
    returnfmapping = fmapping + dummyNodesMapping
    
    returnfmapping = sorted(returnfmapping)    

    # Condensing
    
    simGraph, startList  = condenseGraph(newSeqGraphNodes)   
    print "condenseGraph len(simGraph), len(startList)", len(simGraph), len(startList)   
    


  #  simGraph, startList = lowDensityRemove(simGraph, startList, fmapping)
  #  simGraph, startList  = condenseGraph(simGraph)   
  #  print "lowDensityRemove len(simGraph), len(startList)", len(simGraph), len(startList)    
    
    for trial in range(5):
        simGraph, startList = transitiveReduction(simGraph)
        simGraph, startList  = newCondensingStep(simGraph)   
        print "transitiveReduction len(simGraph), len(startList)", len(simGraph), len(startList)    
    
        simGraph, startList = removeLoopsAndCycles(simGraph)
        simGraph, startList  = newCondensingStep(simGraph) 
        print "removeLoopsAndCycles len(simGraph), len(startList)", len(simGraph), len(startList)    
    
    
        simGraph, startList = combineSelfReferal(simGraph)
        simGraph, startList  = newCondensingStep(simGraph) 
        print "combineSelfReferal len(simGraph), len(startList)", len(simGraph), len(startList)    
        
        
        simGraph, startList = endRemoval(simGraph)
        simGraph, startList  = newCondensingStep(simGraph)   
        print "endRemoval len(simGraph), len(startList)", len(simGraph), len(startList)
  
    
    debugSeqGraph(simGraph)
    
    return simGraph,startList, returnfmapping
    
    
def removeLoopsAndCycles(simGraph):
    # Remove selfLoops
    
    thres = 8
    
    selfLoopList,cycleList = [], []
    for eachitem in simGraph:
        if eachitem in eachitem.listOfNextNodes and len(eachitem.nodeIndexList) < thres :
            selfLoopList.append(eachitem)

    
    for eachnode in selfLoopList:
        if eachnode in eachnode.listOfNextNodes:
            eachnode.listOfNextNodes.remove(eachnode)
            eachnode.listOfPrevNodes.remove(eachnode)
    
    
    for eachitem in simGraph:
        if len(eachitem.listOfPrevNodes) == 1 and len(eachitem.listOfNextNodes ) == 1 and len(eachitem.nodeIndexList) < thres:
            if eachitem.listOfPrevNodes[0] == eachitem.listOfNextNodes[0] :
                cycleList.append(eachitem)  
    
    
    print "len(cycleList)", len(cycleList)
    
    for eachnode in cycleList:
        if len(eachnode.listOfNextNodes) > 0 :
            targetNode = eachnode.listOfNextNodes[0]
            targetNode.listOfPrevNodes.remove(eachnode)
            targetNode.listOfNextNodes.remove(eachnode)
            
            eachnode.listOfNextNodes  = [] 
            eachnode.listOfPrevNodes = []
            eachnode.nodeIndexList = [] 
        
        
    newGraph, startingList  = [],  []
    
    for eachnode in simGraph:
        if len(eachnode.nodeIndexList) > 0 :
            newGraph.append(eachnode)
        if len(eachnode.listOfPrevNodes) == 0 :
            startingList.append(eachnode)
            
    return newGraph, startingList
 
 
def combineSelfReferal(simGraph):
    
   #return simGraph, [simGraph[0]]

    thres = 8
    selfReferalList = [] 
    finishedProcessList = []
    for eachnode in simGraph:
        for eachnextnode in eachnode.listOfNextNodes:
            if eachnextnode in eachnode.listOfPrevNodes and len(eachnode.nodeIndexList) < thres and len(eachnextnode.nodeIndexList) < thres:
                if eachnode.nodeIndex < eachnextnode.nodeIndex:
                    selfReferalList.append([eachnode, eachnextnode])
    
    print "len(selfReferalList)", len(selfReferalList)
    for eachitem in selfReferalList:
        
        currentNode = eachitem[0]
        targetNode = eachitem[1]
        
        if not currentNode.nodeIndex in finishedProcessList and not targetNode.nodeIndex in finishedProcessList:          
            for eachprev in targetNode.listOfPrevNodes :
                if not eachprev in currentNode.listOfPrevNodes and eachprev != currentNode:
                    currentNode.listOfPrevNodes.append(eachprev)
                    eachprev.listOfNextNodes.append(currentNode)
                    
            for eachnext in targetNode.listOfNextNodes :
                if not eachnext in currentNode.listOfNextNodes and eachnext != currentNode:
                    currentNode.listOfNextNodes.append(eachnext)
                    eachnext.listOfPrevNodes.append(currentNode)
    
            for eachresidual in targetNode.listOfNextNodes:
                eachresidual.listOfPrevNodes.remove(targetNode)
            for eachresidual in targetNode.listOfPrevNodes:
                eachresidual.listOfNextNodes.remove(targetNode)
                
            targetNode.listOfNextNodes = []
            targetNode.listOfPrevNodes = [] 
            targetNode.nodeIndexList = []
            
            finishedProcessList.append(targetNode.nodeIndex)
                    
     
    newGraph, startingList  = [],  []
    
    for eachnode in simGraph:
        if len(eachnode.nodeIndexList) > 0 :
            newGraph.append(eachnode)
        if len(eachnode.listOfPrevNodes) == 0 :
            startingList.append(eachnode)
            
    return newGraph, startingList
 


def countNumber(currentKmerIndex, fmapping):
    counter = 0  
    findex = bisect.bisect_left(fmapping, [currentKmerIndex, -1,-1])
    while (findex < len(fmapping) and fmapping[findex][0] == currentKmerIndex):
        findex += 1
        counter+= 1 
        
    return counter  
def flowBalancingTransform(G2, parameterRobot, f1):
    newGraph, startList = [],[]
    liid = parameterRobot.liid
    thres = 8 
    f1 = sorted(f1)
    
    for eachnode in G2:
        if len(eachnode.nodeIndexList) < liid and len(eachnode.listOfPrevNodes) > 1 and len(eachnode.listOfNextNodes) > 1:
            countInSum = 0 
            for eachprevnode in eachnode.listOfPrevNodes:
                countInSum  += countNumber(eachprevnode.nodeIndexList[-1], f1)
            
            
            countMySum = countNumber(eachnode.nodeIndexList[0], f1)
            
            if abs( countMySum - countInSum ) > thres:
                pureIn = []
                pureOut = []
                bothInAndOut  =[]
                for eachprevnode in eachnode.listOfPrevNodes:
                    if not eachprevnode in eachnode.listOfNextNodes:
                        pureIn.append(eachprevnode)
                    else:
                        bothInAndOut.append(eachprevnode)
                for eachnextnode in eachnode.listOfNextNodes:
                    if not eachnextnode in eachnode.listOfPrevNodes:
                        pureOut.append(eachnextnode)
                
                 
                 
                # Process special flow balancing check 
                for eachboth in bothInAndOut:
                    eachboth.listOfPrevNodes.remove(eachnode)
                    eachboth.listOfNextNodes.remove(eachnode)
                    eachnode.listOfPrevNodes.remove(eachboth)
                    eachnode.listOfNextNodes.remove(eachboth)
                    
                    eachboth.listOfPrevNodes= pureIn
                    eachboth.listOfNextNodes = pureOut
                    
                    for eachitem in pureIn:
                        eachitem.listOfNextNodes.append(eachboth)
                        
                    for eachitem in pureOut:
                        eachitem.listOfPrevNodes.append(eachboth)
                
    
    for eachnode in G2:
        if len(eachnode.listOfPrevNodes) == 0:
            startList.append(eachnode)
        
    if len(startList) == 0 :
        startList = [G2[0]]
         
    return G2, startList

            
def debugSeqGraph(simGraph):
    # @kkdebug

    sizeList = []
    noNext = 0 
    noPrev = 0 
    outputNodeList = []
    neighborNodeList = []
    
    
    for eachitem in simGraph:
        
        sizeList.append([eachitem.nodeIndex, len(eachitem.nodeIndexList), len(eachitem.listOfPrevNodes), len(eachitem.listOfNextNodes)] )
        #if len(eachitem.listOfPrevNodes) == 1 : 
        #    print "prev == 1 :", len(eachitem.listOfPrevNodes[0].listOfNextNodes), len(eachitem.nodeIndexList)
        #if len(eachitem.listOfNextNodes) == 1 :
        #    print  "next == 1 :", len(eachitem.listOfNextNodes[0].listOfPrevNodes), len(eachitem.nodeIndexList)
        
        if len(eachitem.listOfNextNodes) == 0:
            noNext +=1 
        if len(eachitem.listOfPrevNodes) == 0 :
            noPrev += 1 
            
        outputNodeList.append(eachitem.nodeIndex)
        
        for dummyNei in eachitem.listOfNextNodes:
            neighborNodeList.append(dummyNei.nodeIndex)
            
        
        for dummyNei in eachitem.listOfPrevNodes:
            neighborNodeList.append(dummyNei.nodeIndex)
            
        
        
    neighborNodeList = filterSameItem(neighborNodeList)
    outputNodeList = filterSameItem(outputNodeList)
    
    #print "neighbNodeList", neighborNodeList
    #print "outputNodeList", outputNodeList
    
    print outputNodeList == neighborNodeList
    print len(neighborNodeList), len(outputNodeList)
    
    for eachitem in outputNodeList:
        if not eachitem in neighborNodeList:
            print "not in neighborNodeList", eachitem 
            
    for eachitem in neighborNodeList:
        if not eachitem in outputNodeList:
            print "not in outputNodeList", eachitem
    
    sizeList = sorted(sizeList)     
    
    counter = 0 
    for eachitem in sizeList:
        counter = counter + eachitem[1]
    print "counter of nodeIndex retained , " , counter
    print "noNext: ", noNext
    print "noPrev: ", noPrev
    
    
    
    # @end kkdebug




def filterSameItem(mylist):
    newList= [] 
    if len(mylist) == 0 :
        return []
    
    mylist = sorted(mylist)
    
    
    for i in range(len(mylist)-1):
        if mylist[i] != mylist[i+1]:
            newList.append(mylist[i])
            
    newList.append(mylist[-1])
    
    return newList





