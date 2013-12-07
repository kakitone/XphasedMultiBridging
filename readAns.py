import numpy as np
import bridgeResolve
import alignmentBridge
from operator import itemgetter
import cleaner

def majorityVote(readsList, noisyReads):
    countArray = [0,0,0,0,0]
    consensusBase = -1
    
    for eachitem in readsList:
        readNum , offset, fused, prevgpNum = eachitem[1:5]
        #print readNum,offset
        chosenOne = noisyReads[readNum][offset]
        countArray[chosenOne] = countArray[chosenOne] +1
    
    #print "countArray",countArray
    optimalIndex = np.argmax(countArray)
    totalCount = np.sum(countArray)
    
    if countArray[optimalIndex] > 0.7*totalCount:
        consensusBase =  optimalIndex
    else:
        consensusBase =  0
        
    return consensusBase


def reportRecovSeqOriginal(recovSeq, f2, noisyReads, parameterRobot):
    
    thresTrust = 50 
    G = parameterRobot.G
    
    f2 = sorted(f2)
    recoveredGenome =  np.zeros(2*G,dtype = np.int8)
    correctedGen = np.zeros(G, dtype = np.int8)
    # -1 : fused, 0 : Snps

    print "len(recovSeq)", len(recovSeq)
    ### MajorityVote

    runningIndex = 0 
    for eachindex in range(len(recovSeq)):

        KmerIndex = recovSeq[eachindex]
        readsList =bridgeResolve.obtainReadDetail(KmerIndex, f2)
        if eachindex % thresTrust == 0 :
            readNum, offset = readsList[0][1:3]
            recoveredGenome[runningIndex : runningIndex+ thresTrust] = noisyReads[readNum][offset: offset + thresTrust]
            runningIndex += thresTrust
          
    correctedGen= recoveredGenome[0:runningIndex] 
    print correctedGen[-10:-1], len(correctedGen)
    countfused = 0
    countsnps = 0
    for eachitem in correctedGen:
        if eachitem == -1 :
            countfused = countfused  +1
        if eachitem == 0: 
            countsnps = countsnps +1 
    
    print "correctedGen[0:30]" ,correctedGen[0:30]  , len(correctedGen)
    print "Not corrected yet", countfused, countsnps
        
    return correctedGen[0:int(G*1.1)]


def countDiff(currentList, nextList,thres = 5):
    countDiffMax  = 0 
    
    readListOfInterest = []
    countList = [] 
    for eachcurrent in currentList:
        currentRead, currentOffset = eachcurrent[1:3]
        for eachnext in nextList: 
            nextRead, nextOffset = eachnext[1:3]
            
            if currentRead == nextRead and nextOffset > currentOffset:
                countList.append(nextOffset- currentOffset)
                
    countList = sorted(countList)
    countDetail = np.bincount(countList)
    
    if len(countDetail) > 0 and np.max(countDetail) > thres:
        countDiffMax = np.argmax(countDetail)
        #print "np.max(countDetail)", np.max(countDetail)
        for eachcurrent in currentList:
            currentRead, currentOffset = eachcurrent[1:3]
            for eachnext in nextList: 
                nextRead, nextOffset = eachnext[1:3]
                
                chk = nextOffset - currentOffset
                if currentRead == nextRead and (chk == countDiffMax):
                    readListOfInterest.append([eachcurrent[0], currentRead, currentOffset])
                    #print "here"
                    
    elif  len(countDetail) > 0 and  np.max(countDetail) <= thres:
        countDiffMax = 1            
    else:
        countDiffMax = 0
        
    return countDiffMax, readListOfInterest

def searchFromRO(readSortedf2,searchitem):

    searchIndex = bridgeResolve.bisectkk(readSortedf2, searchitem) 
    #print "searchitem, searchresult", searchitem, readSortedf2[searchIndex-1]
    newKmerIndex = readSortedf2[searchIndex-1][2]
    return newKmerIndex

def  findIndex(currentList, f2,readSortedf2, dummy):
    newKmerIndex = 0 
    mycurrentList = sorted(currentList, key = itemgetter(2))
    #print "mycurrentList", mycurrentList

    readNum, offset = mycurrentList[0][1:3]
    
    offset = offset + dummy
    #print "dummy",dummy
    newKmerIndex= searchFromRO(readSortedf2, [readNum, offset, -1])
    #print "-------"
    #print newKmerIndex
    return newKmerIndex

def addBackDeleted(recovSeq, f2,readSortedf2, noisyReads, parameterRobot):
    newRecovSeq = np.zeros(len(recovSeq)*2, dtype= np.int64)

    counter = 0 
    for eachindex in range(len(recovSeq)-1):
        KmerIndex = recovSeq[eachindex]
        currentList =bridgeResolve.obtainReadDetail(KmerIndex, f2) 
    
        KmerIndex = recovSeq[eachindex+1 ]
        nextList =bridgeResolve.obtainReadDetail(KmerIndex, f2) 
        
        countDiffMax,readListOfInterest = countDiff(currentList, nextList)
        
        if countDiffMax >1 :
            #print "countDiffMax", countDiffMax
            
            for dummy in range(countDiffMax):
              
                newKmerIndex = findIndex(readListOfInterest, f2,readSortedf2, dummy)
                newRecovSeq[counter+ dummy] = newKmerIndex 
            counter = counter + countDiffMax
        elif countDiffMax == 1 :
            newRecovSeq[counter] = currentList[0][0]
            counter = counter + 1 

    
    return newRecovSeq[0:counter]

def deleteExtra(recovSeq, f2, noisyReads, parameterRobot):
    newRecovSeq = np.zeros(len(recovSeq), dtype = np.int64)
    thresTrue = 3
    counter = 0 
    for eachindex in range(len(recovSeq)):
        KmerIndex = recovSeq[eachindex]
        readsList =bridgeResolve.obtainReadDetail(KmerIndex, f2)
        if len(readsList) >= thresTrue: 
            newRecovSeq[counter] = KmerIndex
            counter = counter + 1 
    
    return newRecovSeq[0:counter]

def convertSeqtoBase(recovSeq, f2, noisyReads, parameterRobot):
    correctedGen = np.zeros(len(recovSeq), dtype = np.int8)  
    
    for eachindex in range(len(recovSeq)): 
        KmerIndex = recovSeq[eachindex]
        readsList =bridgeResolve.obtainReadDetail(KmerIndex, f2)
        if len(readsList) == 0:
            print len(readsList), eachindex
        readNum, offset = readsList[0][1:3] 
        correctedGen[eachindex] = noisyReads[readNum][offset]
    
    return  correctedGen


def reorder(f2):
    newList = []
    for eachitem in f2:
        newList.append([eachitem[1], eachitem[2], eachitem[0]])
        
    newList = sorted(newList)
    
    return newList


def findCommon(currentList, nextList):
    commonRead = 0
    for eachcur in currentList:
        for eachnext in nextList:
            if eachcur[1] == eachnext[1]:
                commonRead,offset1, offset2 = eachcur[1], eachcur[2], eachnext[2]
    
    return commonRead, offset1, offset2


def checkSkip(currentList):
    skip = False
    myCurrentList = sorted(currentList, key = itemgetter(1))
    thres = 1
    
    prev = -1
    count = 0
    
    for i in range(len(myCurrentList)-1):
        if myCurrentList[i][1]!= prev:
            if count > thres:
                return True
            else:
                count = 1
                prev = myCurrentList[i][1]
        else:
            count = count +1 
    
    return skip 
def segmentedAdd(recovSeq, f2 ,readSortedf2,  noisyReads, parameterRobot):
    recovSeqNew = np.zeros(2*len(recovSeq), dtype = np.int64)
    thres = 2
    W = 100

    counter = 0 
    runningSum = 0
    
    while counter< len(recovSeq)-W:
        KmerIndex = recovSeq[counter]
        currentList =bridgeResolve.obtainReadDetail(KmerIndex, f2) 
    
        KmerIndex = recovSeq[counter+W ]
        nextList =bridgeResolve.obtainReadDetail(KmerIndex, f2) 
        
        countDiffMax, readListOfInterest = countDiff(currentList, nextList, 5)
        
        skip = checkSkip(currentList)
        
        if (not skip and -thres <=countDiffMax - W <= thres) or countDiffMax ==0 or len(readListOfInterest) == 0:
            recovSeqNew[runningSum] = recovSeq[counter]
            
            counter += 1
            runningSum += 1
        else:
            for dummy in range(countDiffMax):
                newKmerIndex = findIndex(readListOfInterest, f2,readSortedf2, dummy)
                recovSeqNew[runningSum+ dummy] = newKmerIndex 

            counter = counter + W
            runningSum += countDiffMax
            
    for counter in range(len(recovSeq)-W, len(recovSeq)):
        recovSeqNew[runningSum] = recovSeq[counter] 
        runningSum += 1
        
    return recovSeqNew[0:runningSum]


def reportRecovSeq(recovSeq, f2, noisyReads, parameterRobot):
    readSortedf2 = reorder(f2)
    
    # Delete the extra bases for confident Kmers 
    for i in range(1): 
        recovSeq = deleteExtra(recovSeq, f2, noisyReads, parameterRobot)
        print "len(recovSeq)", len(recovSeq)
        # Add back the deleted items 
        ### Todo : when adding, use only the read that fit the criterion
        ### Bug in adding back deleted still
        ### Add a filter to make sure thing added back is not too few frequency
        
        ## Remaining are related to the run of long alphabets 
        
        #recovSeq = addBackDeleted(recovSeq, f2, readSortedf2, noisyReads, parameterRobot)
        print "len(recovSeq)", len(recovSeq)
        
        recovSeq = segmentedAdd(recovSeq, f2 ,readSortedf2, noisyReads, parameterRobot)
        print "len(recovSeq)",len(recovSeq) 
        
        #recovSeq = segmentedAdd(recovSeq, f2 ,readSortedf2, noisyReads, parameterRobot)
        #print "len(recovSeq)",len(recovSeq) 
    # Form the corrected Sequences 
    correctedGen = convertSeqtoBase(recovSeq, f2, noisyReads, parameterRobot)
    
    return correctedGen
    
    

    




