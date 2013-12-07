def explore(vtup):
    tempSeq = []
    returnNodes = []
    v = vtup[0]
    
    returnNodes.append(vtup)
    
    ttuple = v.listOfNextNodes.pop(0)
    t = ttuple[0]
    
    returnNodes.append(ttuple)
    
    while t.nodeIndex != v.nodeIndex and len(t.listOfNextNodes) > 0 :
        tempnext = t.listOfNextNodes.pop(0)
        returnNodes.append(tempnext) 
        t = tempnext[0]
        

    return returnNodes

def findEC(G4):
    ### Tranverse Graph
    nodeList = [ [G4[0], 0 ] ]
    tempSeq = []
    finishedNodes = []
    
    while len(nodeList) > 0: 
        currentNodeTup = nodeList.pop(0)
        currentNode = currentNodeTup[0]
        
        if len(currentNode.listOfNextNodes) > 0:
            tempSeq = explore(currentNodeTup)
            
            nodeList = tempSeq + nodeList
            
            while (len(nodeList) > 0 and len(nodeList[0][0].listOfNextNodes ) == 0 ):
                tempnode = nodeList.pop(0)
                finishedNodes.append(tempnode)


    ### List of all the Kmers needed
    listOfKmerIndexToBeAssembled = []
    
    print "---------------------"
    print "Check EC Finding"
    
    for eachitem in finishedNodes:
        print "eachitem[0].nodeIndex, len(eachitem[0].nodeIndexList), eachitem[1]",     eachitem[0].nodeIndex, len(eachitem[0].nodeIndexList), eachitem[1]
    
    for eachnode in finishedNodes:
        listOfKmerIndexToBeAssembled = listOfKmerIndexToBeAssembled + eachnode[0].nodeIndexList[eachnode[1]:len( eachnode[0].nodeIndexList)]
        print "len(listOfKmerIndexToBeAssembled)", len(listOfKmerIndexToBeAssembled)

    recovSeq = listOfKmerIndexToBeAssembled
    return recovSeq





















