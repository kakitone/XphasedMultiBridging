import cluster
import graphForm
import logging
import bisect


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
        
def smallNodeEdges(smallNodesList2,levelList, typeOfMode ):
    internalPairsList = []
    print "smallNodesList2", smallNodesList2
    tempList = []
    for eachitem in levelList:
        tempList.append([eachitem[1], eachitem[0]])
    tempList = sorted(tempList)
    
    print "tempList", tempList
    smallNodesList2 = sorted(smallNodesList2)
    internalPairsList = []
    index = 0
    while ( index < len(smallNodesList2) - 1 ):
        if smallNodesList2[index] == smallNodesList2[index + 1]:
            smallNodesList2.pop(index)
        else:
            index = index + 1
    

    for eachitem in smallNodesList2:
        searchkey = [ eachitem[1] , -1 ]
        tempindex = bisect.bisect_left(tempList, searchkey)
        
        if tempindex < len(tempList)  :
            print "tempList[tempindex][0] ,  eachitem[1]", tempList[tempindex][0] ,  eachitem[1]
            if tempList[tempindex][0] == eachitem[1]:
                #if eachitem[0] != tempList[tempindex][1]:
                if typeOfMode == "next":
                    internalPairsList.append([eachitem[0], tempList[tempindex][1]])
                elif typeOfMode == "prev":
                    internalPairsList.append([ tempList[tempindex][1],eachitem[0]])
        
    return internalPairsList



def clearResidual(returnfmapping, G1, parameterRobot):
    
    logging.fmappingSave(returnfmapping, parameterRobot.defaultFolder)
    
    UBd = parameterRobot.brachingDepth               
    for index in range(UBd-1,UBd,1):    
        print "index", index 
        parameterRobot.brachingDepth = index
        returnfmapping, G1= clearResidualMain(returnfmapping,G1,parameterRobot)
    
    
    f2, G2 = returnfmapping, G1
    
    G2 = sorted(G2)
    eachnodeindex =0 
    while ( eachnodeindex < len( G2) -1 ):
        if G2[eachnodeindex] == G2[eachnodeindex +1]:
            G2.pop(eachnodeindex)
        else:
            eachnodeindex = eachnodeindex +1 
            
    f2 = sorted(f2) 
    logging.fmapfusedSave(f2, parameterRobot.defaultFolder)
    logging.storeGraph(G2, parameterRobot.defaultFolder )
    return f2, G2 
    
    
def clearResidualMain(f1, G1,parameterRobot):

    f2, G2 = [], []
    branchLimit = parameterRobot.brachingDepth    
    queue = []
    
    # Classify Big or Small sequences
    print "Classify Big or Small sequences"
    startList = []
    for v in G1:
        v.visited = False
        if len(v.listOfPrevNodes) == 0:
            startList.append(v)
      
    if len(startList) != 0 :
        queue = startList
    else: 
        runningindex = 0

        while len(G1[runningindex].listOfNextNodes) == 0:
            runningindex = runningindex+1
        queue  = [G1[runningindex]]
        #print len(G1[runningindex].listOfNextNodes)
        
    bigList = []
    smallList = []
    for eachelem in queue:
        eachelem.visited = True
        
    while (len(queue) > 0):

        currentNode = queue.pop(0)
        print "currentNode.nodeIndex, len(currentNode.listOfNextNodes), len(currentNode.nodeIndexList)  "    ,currentNode.nodeIndex, len(currentNode.listOfNextNodes), len(currentNode.nodeIndexList)
        
        for eachnextnode in currentNode.listOfNextNodes:
            if eachnextnode.visited == False:
                queue.append(eachnextnode)
                eachnextnode.visited = True
        
        print len(currentNode.nodeIndexList) 
        if len(currentNode.nodeIndexList) > branchLimit:
            bigList.append(currentNode)
        else:
            smallList.append(currentNode)
        
    # Clear Branches originated from big branches
    print "Clear Branches originated from big branches"
    #countGroup = 0
    for v in G1:
        v.visited = False
    
    f1 = sorted(f1)
    countGroup = int(f1[-1][0] ) +1
    #countGroup = countGroup + len(v.nodeIndexList)
    #print "countGroup", countGroup
    
    clusterList = []
    for index in range(countGroup):
        clusterList.append(cluster.clusterElem(index))
    
    print "len(bigList)", len(bigList)
    
    # remove small dead-end in 
    for v in bigList:
        eachitemindex = 0
        while (eachitemindex < len( v.listOfPrevNodes)):
            if len(v.listOfPrevNodes[eachitemindex].nodeIndexList) < branchLimit and len(v.listOfPrevNodes[eachitemindex].listOfPrevNodes) == 0:
                print "eachitemindex", eachitemindex
                v.listOfPrevNodes[eachitemindex].listOfNextNodes = []
                v.listOfPrevNodes[eachitemindex].nodeIndexList = []
                v.listOfPrevNodes.remove(v.listOfPrevNodes[eachitemindex])
                
                
            else:
                eachitemindex = eachitemindex + 1
    
    for v in bigList:
        
        
        levelList = []
        inAtLevelList = []
        outAtLevelList = []
        
        smallNodesListnext = []
        smallNodesListprev = []
        # Collect Associated clusters and put them into levels 
        queue = []
        
        print "v.nodeIndex, len(v.listOfNextNodes)", v.nodeIndex, len(v.listOfNextNodes)
        
        runningindex = 0
        while (runningindex < len(v.listOfNextNodes)):
            eachnode = v.listOfNextNodes[runningindex]
            #print " eachnode.nodeIndex", eachnode.nodeIndex
            if not eachnode in bigList :
                # Cut deadend
                if  len(eachnode.listOfNextNodes) > 0:
                    queue.append([eachnode,0])
                v.listOfNextNodes.remove(eachnode)
                eachnode.listOfPrevNodes.remove(v)
            else:
                runningindex  = runningindex +1 

        inAtLevelList.append([0, v])
        print "len(v.listOfNextNodes),len(v.listOfPrevNodes)", len( v.listOfNextNodes), len(v.listOfPrevNodes)
        #print "len(queue)", len(queue)
        
        while ( len(queue)  > 0):
            currentNode, cumLvl = queue.pop(0)
            currentNode.visited = True
            
            for eachprevnode in currentNode.listOfPrevNodes:
                if eachprevnode in bigList and not [cumLvl, eachprevnode] in inAtLevelList :
                    inAtLevelList.append([cumLvl, eachprevnode])
                if not eachprevnode in bigList :
                    smallNodesListprev.append([ cumLvl, eachprevnode.nodeIndexList[-1],eachprevnode])
                #if eachprevnode.nodeIndex == 3161:
                #    inAtLevelList.append([cumLvl, eachprevnode])
                #    print "???", v.nodeIndex

            for eachnextnode in currentNode.listOfNextNodes:
                if eachnextnode in bigList and not [cumLvl, eachnextnode] in outAtLevelList:
                    outAtLevelList.append([cumLvl+ len(currentNode.nodeIndexList) -1, eachnextnode])
                    
                if not eachnextnode in bigList :
                    smallNodesListnext.append([ cumLvl+len(currentNode.nodeIndexList) -1, eachnextnode.nodeIndexList[0], eachnextnode])
                
            
            for eachindex, runningindex in zip(currentNode.nodeIndexList, range(len(currentNode.nodeIndexList))):
                levelList.append([runningindex + cumLvl, eachindex])
                
            for eachnode in currentNode.listOfNextNodes :
                
                if ( not eachnode in bigList) and ( eachnode.visited == False):
                    queue.append([eachnode, cumLvl+ len(currentNode.nodeIndexList)])
                
                
            # remove edges
            runningindex = 0
            while (runningindex < len(currentNode.listOfNextNodes)):
                eachnode = currentNode.listOfNextNodes[runningindex]
                
                if currentNode in eachnode.listOfPrevNodes:
                    eachnode.listOfPrevNodes.remove(currentNode)
                    currentNode.listOfNextNodes.remove(eachnode)
                else:
                    runningindex = runningindex + 1
            
            runningindex = 0
            while (runningindex < len(currentNode.listOfPrevNodes)):    
                eachnode = currentNode.listOfPrevNodes[runningindex]   
                if currentNode in eachnode.listOfNextNodes:
                    eachnode.listOfNextNodes.remove(currentNode)
                    currentNode.listOfPrevNodes.remove(eachnode)
                else:
                    runningindex = runningindex  +1 
                    
            currentNode.nodeIndexList = []
            
            
        
        # Find backward edges : 
        # internalPairsList : Formats: { (inLvl, outLvl)  }  e.g. {(0,1), (1,2), (2,3), (3,4)... }
        # smallNodesList.append([ cumLvl+len(currentNode.nodeIndexList) -1, eachnextnode.nodeIndex])]
        
        #Filtering of small node missing
        runningindex = 0
        while (runningindex < len(smallNodesListnext) ):
            
            eachitem = smallNodesListnext[runningindex]
            nodeIndex = eachitem[1]
            
            found = False
            for dummy in levelList:
                if nodeIndex == dummy[1]:
                    found = True
            
            if not found:
                outAtLevelList.append([eachitem[0], eachitem[2]])
                smallNodesListnext.pop(runningindex)
            else:
                runningindex = runningindex +1 


        runningindex = 0
        while (runningindex < len(smallNodesListprev) ):
            eachitem = smallNodesListprev[runningindex]
            #print eachitem[2]
            nodeIndex = eachitem[1]
            found = False
            for dummy in levelList:
                if nodeIndex == dummy[1]:
                    found = True
            
            if not found:
                inAtLevelList.append([eachitem[0], eachitem[2]])
                smallNodesListprev.pop(runningindex)
            else:
                runningindex = runningindex +1                 
        #End filtering
        
        ### Special treatment for indel : no backedge added for small nodes
        internalPairsList = []
        #internalPairsList = smallNodeEdges(smallNodesListnext,levelList, "next" ) 
        #internalPairsList = filterSameItem(internalPairsList)
        #print "internalPairsList" , internalPairsList
        
        
        
        #internalPairsList = internalPairsList + smallNodeEdges(smallNodesListprev,levelList, "prev" ) 
        #print "internalPairsList", internalPairsList


        ### Special treatment for indel : no backedge added for small nodes End
        # End Find backward edges 
                

        # Merge Nodes
        levelList= sorted(levelList)
        print "levelList",levelList
        print "inAtLevelList",inAtLevelList
        print "outAtLevelList",outAtLevelList
        inAtLevelList = sorted(inAtLevelList)
        outAtLevelList = sorted(outAtLevelList)
        
        ### list hacks :
        if len(outAtLevelList) > 0 :
            finalOut = outAtLevelList[-1][0]
        else: 
            finalOut = -1
            
        finalinAtLevellist = []
        for eachinlvl in inAtLevelList:
            if eachinlvl[0] > finalOut:
                finalinAtLevellist.append([finalOut, eachinlvl[1]])
            else:
                finalinAtLevellist.append(eachinlvl)
                
        inAtLevelList = sorted(finalinAtLevellist)        
        # End list hacks 
        
        if len(levelList) > 0:
            numberOfLevels = levelList[-1][0]
            toMergeList = [[[],[],[]] for i in range(numberOfLevels+1)]
            
            for item in levelList:
                index = item[0]
                content = item[1]
                #print index
                toMergeList[index][0].append(content)
                
            for item in inAtLevelList:
                index = item[0]
                content = item[1]
                toMergeList[index][1].append(content)
                
            for item in outAtLevelList:
                index = item[0]
                content = item[1]
                toMergeList[index][2].append(content)        
                
                
            print "toMergeList",toMergeList
            
            # init nodes array
            vArray = []
            for i in range(len(toMergeList)):
                idOfNode = toMergeList[i][0][0]
                v = []
                v = graphForm.condensedNode(idOfNode)
                v.updateNodeList()
                vArray.append(v)
                
            for i in range(len(toMergeList)):
                idOfNode = toMergeList[i][0][0]
                mylistOfcluster = toMergeList[i][0]
                mylistOfPrevNodes = toMergeList[i][1]
                mylistOfNextNodes = toMergeList[i][2]
                
                if i > 0:
                    mylistOfPrevNodes = mylistOfPrevNodes + [vArray[i-1]]
                    
                if i < len(toMergeList) -1:  
                    mylistOfNextNodes = mylistOfNextNodes + [vArray[i+1]]
                
                v = vArray[i]
                
                for eachnode in mylistOfPrevNodes:
                    if not eachnode in v.listOfPrevNodes:
                        v.addPrevNodes(eachnode)
                    if not v in eachnode.listOfNextNodes:
                        eachnode.addNextNodes(v)
                    
                    
                for eachnode in mylistOfNextNodes:
                    if not eachnode in v.listOfNextNodes:
                        v.addNextNodes(eachnode)
                    if not v in eachnode.listOfPrevNodes:
                        eachnode.addPrevNodes(v)
                
                
                #print mylistOfcluster
                for eachindex in mylistOfcluster:
                    cluster.union(clusterList[idOfNode], clusterList[eachindex])
            
            
            print "internalPairsList", internalPairsList
            for eachitem in internalPairsList:
                tmpprevnode = eachitem[0]
                tmpnextnode = eachitem[1]

                if not vArray[tmpprevnode] in vArray[tmpnextnode].listOfPrevNodes:
                    vArray[tmpnextnode].listOfPrevNodes.append(vArray[tmpprevnode])
                if not vArray[tmpnextnode] in vArray[tmpprevnode].listOfNextNodes:
                    vArray[tmpprevnode].listOfNextNodes.append(vArray[tmpnextnode])
                        
            for i in range(len(toMergeList)):
                print "vArray[i].nodeIndex, len(vArray[i].listOfPrevNodes), len(vArray[i].listOfNextNodes)",vArray[i].nodeIndex, len(vArray[i].listOfPrevNodes), len(vArray[i].listOfNextNodes)
            
        
    # Formatting Return 
    print "Formatting Return "
    seqGraphNodes = []
    print len(bigList)
    #for eachitem in bigList:
    #    print len(eachitem.nodeIndexList), len(eachitem.listOfPrevNodes), len(eachitem.listOfNextNodes), eachitem.listOfPrevNodes[0].nodeIndex
    queue = [bigList[0]]
    
    while len(queue) > 0:
        currentNode = queue.pop(0)
        currentNode.visited = True
        #print "currentNode.nodeIndex ", currentNode.nodeIndex
        if len(currentNode.nodeIndexList) > 0:
            seqGraphNodes.append(currentNode)
        
        for eachnode in currentNode.listOfNextNodes:
            if eachnode.visited == False:
                queue.append(eachnode)

     
    #G2 =            seqGraphNodes
    G2,startList2 = graphForm.condenseGraph(seqGraphNodes)  
    sizeOfGraph = len(G2) 
    for index in range(sizeOfGraph):
        G2,startList2 = graphForm.condenseGraph(G2)  

    # Hacking the deadends :- ??? 
    for eachnode in G2:
        if len(eachnode.listOfPrevNodes) == 0 :
            print "no prev"
            runningindex =0 
            while (runningindex < len(eachnode.listOfNextNodes)):
                if eachnode in eachnode.listOfNextNodes[runningindex].listOfPrevNodes:
                    eachnode.listOfNextNodes[runningindex].listOfPrevNodes.remove(eachnode)
                else:
                    runningindex = runningindex+1
            eachnode.listOfNextNodes = []
            eachnode.nodeIndexList = []
            
        if len(eachnode.listOfNextNodes) == 0:
            print "nonext", eachnode.nodeIndex, len(eachnode.nodeIndexList)
            runningindex =0 
            while (runningindex < len(eachnode.listOfPrevNodes)):
                if eachnode in eachnode.listOfPrevNodes[runningindex].listOfNextNodes:
                    eachnode.listOfPrevNodes[runningindex].listOfNextNodes.remove(eachnode)
                    print "eachnode.nodIndex, len(eachnode.nodeIndexList)", eachnode.nodeIndex, len(eachnode.nodeIndexList)
                    print "eachnode.listOfPrevNodes[runningindex].nodeIndex",eachnode.listOfPrevNodes[runningindex].nodeIndex
                else:
                    runningindex = runningindex +1 
            eachnode.listOfPrevNodes = []
            eachnode.nodeIndexList = []    
    
      
    #print   "len(G2)", len(G2)

    runningindex =0 
    while (runningindex < len(G2)):
        if len(G2[runningindex].nodeIndexList) == 0:
            G2.pop(runningindex)
        else:
            runningindex = runningindex + 1 
    
    print "len(G2)", len(G2)
     
    #G2, startList2 = graphForm.newCondensingStep(G2)
    
    for trial in range(5):
        G2, startList = graphForm.transitiveReduction(G2)
        G2, startList  = graphForm.newCondensingStep(G2)   
 
    
        G2, startList = graphForm.removeLoopsAndCycles(G2)
        G2, startList  = graphForm.newCondensingStep(G2) 
   
    
        G2, startList = graphForm.combineSelfReferal(G2)
        G2, startList  = graphForm.newCondensingStep(G2) 
        
        G2, startList = graphForm.endRemoval(G2)
        G2, startList  = graphForm.newCondensingStep(G2)   

    
    G2, startList = graphForm.flowBalancingTransform(G2, parameterRobot, f1)
    ### Finish hack
    # OutputFormat :  Gp id , read #, offset #, fusedOrNot, prevGroup id 
    f2 = []
    for eachitem in f1:
        oldGpid = eachitem[0]

        newGp = cluster.find(clusterList[oldGpid])
        newGpid = newGp.id

        readNum = eachitem[1]
        offset = eachitem[2]
        
        if len(cluster.familyList(clusterList[newGpid])) == 1:
            fused = False
        else:
            fused = True
            
        rowRecord = [newGpid, readNum, offset, fused, oldGpid]
        
        f2.append(rowRecord)
    
    f2 = sorted(f2) 

    return f2, G2