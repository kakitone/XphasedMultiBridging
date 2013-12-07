import numpy as np 
import logging
import time
import csv 
#from blist import sortedlist, blist
import bisect
import random
from scipy import weave
from operator import itemgetter

### Disjoint Union Data Structure
class clusterElem(object):
    def __init__(self,index):
        self.rank = 0
        self.parent = self
        self.id = index
        self.childList =[]
        
        #self.size = 1

def find(x):
    #if x != x.parent:
    #    x.parent = find(x.parent)
    #return x.parent
    if x.parent == x:
        return x
    else:
        return find(x.parent)
       
def union(x,y):
    xRoot = find(x)
    yRoot = find(y)
    
    if xRoot == yRoot:
        return 0
    
    if xRoot.rank < yRoot.rank:
        xRoot.parent = yRoot
        yRoot.childList.append(xRoot)
        
        #yRoot.size = yRoot.size + xRoot.size
        
    elif xRoot.rank > yRoot.rank:
        yRoot.parent = xRoot
        xRoot.childList.append(yRoot)
        
        #xRoot.size = yRoot.size + xRoot.size
        
    else:
        yRoot.parent = xRoot
        xRoot.childList.append(yRoot)
        
        xRoot.rank = xRoot.rank + 1 
        
        #xRoot.size = yRoot.size + xRoot.size

    return 1

def familyList(x):
    root = find(x)
    stack = []
    familyKmers = [root]
    
    stack.append(root)
    while (len(stack) >0 ):
        item = stack.pop(0)
        for eachsubitem in item.childList:
            familyKmers.append(eachsubitem)
            stack.append(eachsubitem)
            
    return familyKmers
                    
    
### Taylor made sorting 
def itemgetterkk(items):
    if len(items) == 1:
        item = items[0]
        def g(obj):
            return obj[0][item]
    else:
        def g(obj):
            return tuple(obj[0][item] for item in items)
    return g



### Linear Tranversal 
def matched(str1, str2, threshold, liid):
    K = len(str1)
    wholeThreshold = int(threshold* K/liid) 
    finalMatch= np.zeros(1,dtype = np.int8)
    lenPass = min(len(str1) , len(str2))
    
    code = """
    int numberOfMismatch ;
    int i ; 
    int startCounter ;
    int endCounter ;
    
    startCounter = 0; 
    endCounter = 0 ;
    
    numberOfMismatch = 0 ;
    for (i = 0 ; i < lenPass ;i++){
        if (str1[i] != str2[i]){
            numberOfMismatch = numberOfMismatch+ 1 ;
            if (i < liid){
                startCounter = startCounter + 1;
            }
            
            if (i >= lenPass - liid){
                endCounter = endCounter + 1 ;
            }
        } 
        
        
    }
    
    if (numberOfMismatch <= wholeThreshold && startCounter <= threshold && endCounter <= threshold ) {
        finalMatch[0] = 1;
    }else{
        finalMatch[0] = 0;
    }
    
    """
    weave.inline(code, ['str1', 'str2', 'threshold','finalMatch', 'lenPass', 'wholeThreshold', 'liid'])

    if finalMatch[0] == 1:
        return True
    elif finalMatch[0] == 0:
        return False
        
        

    
 #   isMatched = True
 #   K = len(str1)
 #   wholeThreshold = int(threshold* K/liid) 
    
 #   numberOfMismatch = 0
    
    
 #   for index in range(min(len(str1), len(str2))):
 #       if str1[index] != str2[index]:
 #           numberOfMismatch = numberOfMismatch + 1

#    if numberOfMismatch > wholeThreshold:
#        isMatched = False 
    
    
#    numberOfMismatch = 0
 #   for index in range(liid):
 #       if str1[index] != str2[index]:
 #           numberOfMismatch = numberOfMismatch + 1

    
#    if numberOfMismatch > threshold:
#        isMatched = False 
#        
#    numberOfMismatch = 0
#    for index in range(liid):
#        if str1[K-1-index] != str2[K-1-index]:
#            numberOfMismatch = numberOfMismatch + 1
    
#    if numberOfMismatch > threshold:
#        isMatched = False     
        
        
#    return isMatched    
    
def fastMatch(str1, str2, threshold, startIndex):
    numberOfMismatch = 0
    for index in range(min(len(str1), len(str2))):
        if str1[index] != str2[index]:
            numberOfMismatch = numberOfMismatch + 1
    #print numberOfMismatch
    if numberOfMismatch <= threshold:
        return True
    else:
        return False      

def linearTranverseList(kmerList, clusterList, threshold, liid):
    for  index in range(len(kmerList) -1 ):

        matchKmer1 = kmerList[index][1][2]
        matchKmer2 = kmerList[index+1][1][2]      
        if  find(clusterList[matchKmer1]) != find(clusterList[matchKmer2]):
            if matched(kmerList[index][0], kmerList[index+1][0], threshold, liid):
                #print matchKmer1, matchKmer2
                union(clusterList[matchKmer1], clusterList[matchKmer2])
                #print clusterList[matchKmer1]
          
    return clusterList


### Use reads to assist clustering                 
def checkAndMerge(clusterList, read,readid , offsetid, targetread,targetreadid, targetoffsetid, threshold, K,clusterTreeSize, liid):              
    lengthreq = -1    
    L = len(read)    

    if offsetid <= targetoffsetid :

        lengthreq = L -K +1 - targetoffsetid + offsetid
        
        #print "lengthreq",lengthreq
        for index in range(lengthreq):
            kmer1 = read[index:index+K]
            kmer2 = targetread[index + targetoffsetid - offsetid:index + targetoffsetid - offsetid+K]
            #selfRoot = find(clusterList[readid*(L-K+1)+ index ])


            if  find(clusterList[readid*(L-K+1)+ index]) != find(clusterList[targetreadid*(L-K+1) + index + targetoffsetid - offsetid]):
            #if (selfRoot.size <= clusterTreeSize):
                if matched(kmer1, kmer2, threshold, liid):
                
                    union(clusterList[readid*(L-K+1)+ index ] ,clusterList[targetreadid*(L-K+1) + index + targetoffsetid - offsetid] )
        

    else: 
        lengthreq = L -K + 1 - offsetid  + targetoffsetid
                  
        for index in range(lengthreq):
            kmer1 = read[index + offsetid - targetoffsetid:index + offsetid - targetoffsetid+K]            
            kmer2 = targetread[index:index+K]
            #selfRoot = find(clusterList[readid*(L-K+1) + index + offsetid - targetoffsetid])


            if find(clusterList[readid*(L-K+1) + index + offsetid - targetoffsetid]) != find(clusterList[targetreadid*(L-K+1)  + index]):
            #if (selfRoot.size <= clusterTreeSize):
                if matched(kmer1, kmer2, threshold, liid):
                
                    union(clusterList[readid*(L-K+1) + index + offsetid - targetoffsetid ], clusterList[targetreadid*(L-K+1)  + index ])
        
    
        
        
        
### Formatting 
def formatClusteringMap(clusterList):
    fmapping = []
    
    headOfListStack = []
    for eachitem in clusterList:
        if eachitem.parent == eachitem :
            headOfListStack.append(eachitem)

    for eachitem, index in zip(headOfListStack,range(len(headOfListStack))):
        fmapping.append([index, familyList(eachitem)])

    return fmapping 

def multiplier(indexStart,indexEnd):
    def cmpFn(kmer1,kmer2):
        i =indexStart
        threshold = indexEnd

        while i< threshold and kmer1[0][i]==kmer2[0][i]:
                i= i+1
        
        if i == threshold :
            return 0
        else:
            return cmp(kmer1[0][i] ,kmer2[0][i])

    return cmpFn


### Fast Clustering 
def fastClusteringAlgo(N, L, K, kmerList, clusterList, noisyReads, threshold,clusterRounds, fingerprint,clusterTreeSize, liid):
    # first order clustering by liid match
    #fingerprint = 6
    
    ## KK Debug
    fout = open("clusterDetail.csv", 'wb')
    mywriter = csv.writer(fout)
    mywriter.writerow(["Tree size","nodeIndex"])
    ## end KK Debug
    
    
    numberOfRounds = min ( int ( K/fingerprint ), clusterRounds)
  
    t0 = time.time()
    
    print "Before linear Trans"  
    for index in range(numberOfRounds):
        #kmerList = sorted(kmerList, key = itemgetterkk(range(fingerprint*index, fingerprint*(index+1))))
        kmerList = sorted(kmerList, cmp=multiplier(fingerprint*index,fingerprint*(index+1)) )
        linearTranverseList(kmerList, clusterList,threshold, liid)
    
    for index in range(numberOfRounds):
        #kmerList = sorted(kmerList, key = itemgetterkk(range(K -fingerprint*(index+1), K -fingerprint*index)))
        #print K -fingerprint*(index+1), K -fingerprint*index, len(kmerList[0][0]), K
        kmerList = sorted(kmerList, cmp=multiplier(K -fingerprint*(index+1), K -fingerprint*index) )
        linearTranverseList(kmerList, clusterList,threshold, liid)
        
    
    print "After linear Trans", time.time() - t0
    mywriter.writerow([ time.time() - t0])
    t0 = time.time()
    
    # second order clustering by  by reads 
        # Read# = div (L-K+1) ; Offset# = mod (L-K+1)
    
    #finishedPairsList = blist([])
    
    for indexN in range(N):

        if (random.random() < 0.8) :
            associatedReads = []
            for eachkmerindex in range(L-K+1):  
                tmpList= []
                tmpList = familyList(clusterList[eachkmerindex + indexN*(L-K+1)])
                for eachitem in tmpList:
                    associatedReads.append(int(eachitem.id / (L-K+1) ))
                    
            associatedReads = sorted(associatedReads)
            
            index = 0
            while (index < len(associatedReads) - 1 ):
                if associatedReads[index] == associatedReads[index +1]:
                    associatedReads.pop(index)
                else:
                    index = index +1
            
            
            for loc in range(L-K+1):
    
                
                #if (random.random() < 2 ):
                tmpList= []
                tmpList = familyList(clusterList[loc + indexN*(L-K+1)])
    
                for eachitem in tmpList:
                    tmpid = int(eachitem.id / (L-K+1) )    
                    targetreadid, targetoffsetid, targetkmerid = int (eachitem.id/ (L-K+1) ), np.mod(eachitem.id,L-K+1), eachitem.id
                    readid , offsetid, kmerid = indexN, loc, indexN * (L-K+1) + loc 
                   # finishIndex = bisect.bisect_left(finishedPairsList, [min(readid,targetreadid), max(readid,targetreadid)] )
    
                    if tmpid in associatedReads:
    
                        #if  (finishIndex >= len(finishedPairsList)) or  ( finishIndex < len(finishedPairsList) and  [min(readid,targetreadid), max(readid,targetreadid)] != finishedPairsList[finishIndex])  :   
                        checkAndMerge(clusterList, noisyReads[readid],readid , offsetid, noisyReads[targetreadid],targetreadid, targetoffsetid, threshold, K,clusterTreeSize, liid)
                        associatedReads.remove(tmpid)
                        #finishedPairsList.append([min(readid,targetreadid), max(readid,targetreadid)])
    
                        #print "finishedPairsList", finishedPairsList
                
    print "After using reads", time.time() - t0    
    mywriter.writerow([ time.time() - t0])
    fout.close()  
                   
### Main Clustering Flow    

### Indel Treatment
def filterDuplicate(mylist):
    
    ### Important parameter to tune
    thres = 10
    
    if len(mylist) == 0:
        return []
    
    mylist = sorted(mylist)
    newlist = [mylist[0]] 
    
    counter = 1
    
    for i in range(len(mylist)-1):
        if mylist[i] != mylist[i+1]:
            if counter > thres:
                newlist.append(mylist[i])
            counter = 1 
            
        else:
            counter += 1 
            
    if counter > thres:
        newlist.append(mylist[-1])
        
    return newlist
            
            
def filterFingerPrint(mylist):
    
    ### Important parameter to tune
    thres = 3
    
    if len(mylist) == 0:
        return []
    
    mylist = sorted(mylist)
    #print "mylist[0:30]", mylist[0:30]
    
    newlist = [] 
    
    counter = 1
    startIndex = 0
    
    for i in range(len(mylist)-1):
        if mylist[i][0] != mylist[i+1][0]:
            if counter > thres:
                targetReadNumSmallest, targetOffsetSmallest, myOffsetSmallest = mylist[startIndex]
                targetReadNumLargest, targetOffsetLargest, myOffsetLargest = mylist[startIndex+counter -1]
                
                assert(targetReadNumSmallest== targetReadNumLargest)
                newlist.append([targetReadNumSmallest, [targetOffsetSmallest, myOffsetSmallest],[targetOffsetLargest, myOffsetLargest]])
                
            counter = 1 
            startIndex = i+1 
            
        else:
            counter += 1 
            
    if counter > thres:
        targetReadNumSmallest, targetOffsetSmallest, myOffsetSmallest = mylist[startIndex]
        targetReadNumLargest, targetOffsetLargest, myOffsetLargest = mylist[startIndex+counter -1]
        
        assert(targetReadNumSmallest== targetReadNumLargest)
        newlist.append([targetReadNumSmallest, [targetOffsetSmallest, myOffsetSmallest],[targetOffsetLargest, myOffsetLargest]])

    return newlist
def fromCompareList(toCompareList, activeKmerList, liid, threshold):
    
    runningList = []
    
    readNum ,offset = activeKmerList[0][1][0], activeKmerList[0][1][1]
    tempgroup = [[readNum ,offset]]

    
    for i in range(len(activeKmerList) -1 ):
        
        if np.array_equal(activeKmerList[i][0], activeKmerList[i+1][0]):
            readNum ,offset = activeKmerList[i+1][1][0], activeKmerList[i+1][1][1]
            tempgroup.append([readNum,offset])
        
        else:
            runningList.append(tempgroup)
            readNum ,offset =activeKmerList[i+1][1][0], activeKmerList[i+1][1][1] 
            tempgroup = [[readNum ,offset ]]
        
        
    runningList.append(tempgroup)
    
    for i in range(len(runningList)):
        runningList[i] = sorted(runningList[i])
    
    for i in range(len(runningList)):
        for j in range(0,len(runningList[i])-1):
            myReadNum, myOffset = runningList[i][j][0],runningList[i][j][1] 
            for k in range(j+1, len(runningList[i])):
                targetReadNum, targetOffset = runningList[i][k][0],runningList[i][k][1]
                toCompareList[myReadNum] += [[targetReadNum, targetOffset, myOffset]] 
                      
    returnList    = []
    for eachitem in  toCompareList:
        returnList.append(filterFingerPrint(eachitem))
        
        
    #print "returnList[0:3]: ", returnList[0:3]
    
            
    return returnList


def reverseString(str):
    newstr = []
    for index in range(len(str)):
        newstr.append(str[len(str)-1-index])
    
    return newstr

def meetRequirement(score, returnalignedSeq1, returnalignedSeq2 ,starti, startj, endi, endj, threshold, liid,overhang, L1, L2 ):
    check = True
    if ( endi - starti ) < liid :
        check = False 
        
    numEdit = ( len(returnalignedSeq1) - score ) /2
    if numEdit > ( threshold/float(liid) ) * len(returnalignedSeq1):
        check = False
    
    if min(L1 - endi, L2 - endj)  > overhang or min(starti, startj) > overhang :
        check  =False 
    

    return check 


def canDoFast(startFingerPrint, endFingerPrint, parameterRobot):
    check = False
    myOffsetStart, targetOffsetStart = startFingerPrint[0], startFingerPrint[1]
    myOffsetEnd, targetOffsetEnd = endFingerPrint[0], endFingerPrint[1]
    
    mylength = myOffsetEnd - myOffsetStart
    targetlength = targetOffsetEnd - targetOffsetStart
    
    #if  0 < mylength < 60 and 0 < targetlength < 60 and abs(mylength - targetlength)< parameterRobot.fingerPrint and not startFingerPrint== endFingerPrint:
    if  0 < mylength  and 0 < targetlength  and abs(mylength - targetlength)< parameterRobot.fingerPrint and not startFingerPrint== endFingerPrint:
     
        check = True
        
    return check


### This function "SWAlignmentBanded" contains bugs... please resolve it ASAP 
def SWAlignmentBanded(startFingerPrint, endFingerPrint, read1, read2, parameterRobot):
    score, returnalignedSeq1, returnalignedSeq2 , starti, startj , endi, endj = 0,[],[],0,0,0,0
    
    myOffsetStart, targetOffsetStart = startFingerPrint[0], startFingerPrint[1]
    myOffsetEnd, targetOffsetEnd = endFingerPrint[0], endFingerPrint[1]
    bandSearch = 20
    #print " myOffsetStart, targetOffsetStart, myOffsetEnd, targetOffsetEnd", myOffsetStart, targetOffsetStart, myOffsetEnd, targetOffsetEnd
    if myOffsetStart < targetOffsetStart:
        suffixRead, prefixRead =  read1 ,read2
        startSkip, endSkip = myOffsetStart+bandSearch , len(read2) - targetOffsetEnd+ bandSearch
        
        
        suffixStart, suffixEnd  = 0 , min( myOffsetEnd + endSkip, len(read1)-1)
        prefixStart, prefixEnd  = max( targetOffsetStart-startSkip, 0), len(read2) -1
         
        
    else:
        suffixRead, prefixRead = read2, read1
        startSkip, endSkip = targetOffsetStart +bandSearch, len(read1) - myOffsetEnd +bandSearch
        
        suffixStart, suffixEnd = 0 , min(targetOffsetEnd + endSkip, len(read2)-1)
        prefixStart, prefixEnd = max(myOffsetStart- startSkip, 0) , len(read1) -1 
        
        

    wts = parameterRobot.editsub
    wti = parameterRobot.editins 
    wtd = parameterRobot.editdel
    wtm=  parameterRobot.editmatch
    

    
    #print "prefixStart, suffixStart: ", prefixStart, suffixStart
    seq1 ,seq2 =  prefixRead[prefixStart: prefixEnd+1], suffixRead[suffixStart:suffixEnd +1 ]
    
    #print "len(seq1),len(seq2)", len(seq1),len(seq2)
    #print prefixRead, suffixRead
    
    m = len(seq1) + 1
    n = len(seq2) +1
    

    H = np.zeros([m,n], dtype = np.float64)
    B = np.zeros([m,n], dtype = np.float64)

    # Assign weights 
    for i in range(m):
        H[i][0] = 0
        B[i][0] = 4
    for j in range(n):
        H[0][j]  =0
        B[0][j] = 4 
        
        
        
    seq1NP = np.zeros(m-1, dtype = np.float64)
    seq2NP = np.zeros(n-1, dtype = np.float64)
    
    for i in range(m-1):
        seq1NP[i] = seq1[i]
    for j in range(n-1) :
        seq2NP[j] = seq2[j]
        

    
    code =\
        """
        
        int i; 
        int j ;
        double w;

        int startIndex ;
        int endIndex ; 
        
        for (i =1 ;i <m ; i++){
            startIndex = i - 2*bandSearch ;
            endIndex = i + 2*bandSearch ;
             
            if (startIndex <= 0 ){
                startIndex = 1 ; 
            }
            
            if (endIndex >=n ){
                endIndex = n-1; 
            }
            
            
            
            for (j=startIndex; j<= endIndex; j++){
                if (seq1NP[i-1] == seq2NP[j-1]){
                    w = wtm ;
                }
                else{
                    w= wts ;
                }
                
                    H2(i,j) = 0 ;
                    
                    if (H2(i,j) < H2(i-1,j-1) + w) {
                        H2(i,j) = H2(i-1,j-1) + w ; 
                    }

                    if (  (j> startIndex    || startIndex == 1 ) && H2(i,j) < H2(i-1,j)+wtd ) {
                        H2(i,j) = H2(i-1,j)+wtd ;
                    }

                    if  (   (j<  endIndex  || endIndex == n-1) &&  H2(i,j) < H2(i,j-1) + wti){
                        H2(i,j) = H2(i,j-1) + wti;
                    }


                    if (H2(i-1,j-1) + w == H2(i,j)){
                        B2(i,j) = 1 ; 
                    }
                    else if ( (j> startIndex    || startIndex == 1 ) && H2(i-1,j)+wtd == H2(i,j)){
                        B2(i,j) = 2;
                    }
                    
                    else if  ((j<  endIndex  || endIndex == n-1) &&  H2(i,j-1) + wti == H2(i,j)){
                        B2(i,j) = 3;
                    }
                    else if (0 == H2(i,j)) {
                        B2(i,j) = 4 ;
                    }
            }
        }

        """
    
    weave.inline(code, ['H','B','m', 'n', 'wtm','wts','wti','wtd', 'seq1NP', 'seq2NP', 'bandSearch'])
    

    # Backtrack 
    alignedSeq1 = []
    alignedSeq2 = []
    
    bestindex = np.argmax(H)

 
    
    endi = bestindex / n
    endj = bestindex%n
    
    score = H[endi][endj]
    scoremax = np.max(H)
    
    assert(score == scoremax)
    tempi, tempj = endi , endj
    while (B[tempi][tempj] != 4):
        if B[tempi][tempj] == 1:
            alignedSeq1.append(seq1[tempi -1 ])
            alignedSeq2.append(seq2[tempj-1 ]) 
            tempi = tempi -1 
            tempj = tempj -1           

        elif B[tempi][tempj] == 2 :
            alignedSeq1.append(seq1[tempi-1 ])
            alignedSeq2.append(0) 
            tempi = tempi - 1

        elif B[tempi][tempj] == 3 :
            alignedSeq1.append(0)
            alignedSeq2.append(seq2[tempj-1 ])
            tempj = tempj -1

    
    starti, startj = tempi , tempj 
    
    returnalignedSeq1 = reverseString(alignedSeq1)
    returnalignedSeq2 = reverseString(alignedSeq2)
    
    scoretmp , returnalignedSeq1tmp, returnalignedSeq2tmp , startitmp, startjtmp , enditmp, endjtmp  = SWAlignment(seq1NP,seq2NP, parameterRobot)
    #print "starti, startj , endi, endj , m , n", starti, startj , endi, endj , m , n
    if [starti, startj , endi, endj ]!= [startitmp, startjtmp , enditmp, endjtmp]:
        print "--------------"
        print "This fcn : score, starti, startj , endi, endj , m , n", score, starti, startj , endi, endj , m , n
        print "This fcn: [prefixStart,suffixStart, prefixEnd+1. suffixEnd +1 ]", prefixStart,  suffixStart,prefixEnd, suffixEnd

        print "SWAlignment : scoretmp, startitmp, startjtmp , enditmp, endjtmp , m , n",scoretmp, startitmp, startjtmp , enditmp, endjtmp , m , n
        scoretmp , returnalignedSeq1tmp, returnalignedSeq2tmp , startitmp, startjtmp , enditmp, endjtmp  = SWAlignment(prefixRead,suffixRead, parameterRobot)
        print "SWAlignment: ", startitmp, startjtmp , enditmp, endjtmp
        
        
        print "startSkip, endSkip",startSkip, endSkip
        print "myOffsetStart, targetOffsetStart",myOffsetStart, targetOffsetStart
        print "myOffsetEnd, targetOffsetEnd ",myOffsetEnd, targetOffsetEnd 
        print "--------------"

    
    starti, endi = starti + prefixStart, endi + prefixStart
   # startj , endj = startj  ,endj 
    #starti , endi , startj , endj = prefixStart , len(prefixRead), 0, suffixEnd
    
    
    if myOffsetStart < targetOffsetStart:
        reverseOrder = score , returnalignedSeq1, returnalignedSeq2 , starti, startj , endi, endj 
        returnOrder = reverseOrder[0], reverseOrder[2], reverseOrder[1], reverseOrder[4], reverseOrder[3], reverseOrder[6], reverseOrder[5] 
    else:
        returnOrder = score , returnalignedSeq1, returnalignedSeq2 , starti, startj , endi, endj 
        
    return returnOrder

def SWAlignment(seq1 , seq2, parameterRobot):
    score  = 0 

    wts = parameterRobot.editsub
    wti = parameterRobot.editins 
    wtd = parameterRobot.editdel
    wtm=  parameterRobot.editmatch
    
    
    
    m = len(seq1) + 1
    n = len(seq2) +1
    

    H = np.zeros([m,n], dtype = np.float64)
    B = np.zeros([m,n], dtype = np.float64)

    # Assign weights 
    for i in range(m):
        H[i][0] = 0
        B[i][0] = 4
    for j in range(n):
        H[0][j]  =0
        B[0][j] = 4 
        
        
        
    seq1NP = np.zeros(m-1, dtype = np.float64)
    seq2NP = np.zeros(n-1, dtype = np.float64)
    
    for i in range(m-1):
        seq1NP[i] = seq1[i]
    for j in range(n-1) :
        seq2NP[j] = seq2[j]
        

    
    code =\
        """
        int i; 
        int j ;
        double w;

        for (i =1 ;i <m ; i++){
            for (j=1; j<n; j++){
                if (seq1NP[i-1] == seq2NP[j-1]){
                    w = wtm ;
                }
                else{
                    w= wts ;
                }
                
                    H2(i,j) = 0 ;
                    
                    if (H2(i,j) < H2(i-1,j-1) + w) {
                        H2(i,j) = H2(i-1,j-1) + w ; 
                    }

                    if (H2(i,j) < H2(i-1,j)+wtd ) {
                        H2(i,j) = H2(i-1,j)+wtd ;
                    }

                    if  (H2(i,j) < H2(i,j-1) + wti){
                        H2(i,j) = H2(i,j-1) + wti;
                    }

                    if (H2(i-1,j-1) + w == H2(i,j)){
                        B2(i,j) = 1 ; 
                    }
                    else if (H2(i-1,j)+wtd == H2(i,j)){
                        B2(i,j) = 2;
                    }
                    
                    else if  (H2(i,j-1) + wti == H2(i,j)){
                        B2(i,j) = 3;
                    }
                    else if (0 == H2(i,j)) {
                        B2(i,j) = 4 ;
                    }
            }
        }

        """
    
    weave.inline(code, ['H','B','m', 'n', 'wtm','wts','wti','wtd', 'seq1NP', 'seq2NP'])
    

    # Backtrack 
    alignedSeq1 = []
    alignedSeq2 = []
    
    bestindex = np.argmax(H)

 
    
    endi = bestindex / n
    endj = bestindex%n
    
    score = H[endi][endj]
    

    tempi, tempj = endi , endj
    while (B[tempi][tempj] != 4):
        if B[tempi][tempj] == 1:
            alignedSeq1.append(seq1[tempi -1 ])
            alignedSeq2.append(seq2[tempj-1 ]) 
            tempi = tempi -1 
            tempj = tempj -1           

        elif B[tempi][tempj] == 2 :
            alignedSeq1.append(seq1[tempi-1 ])
            alignedSeq2.append(0) 
            tempi = tempi - 1

        elif B[tempi][tempj] == 3 :
            alignedSeq1.append(0)
            alignedSeq2.append(seq2[tempj-1 ])
            tempj = tempj -1

    
    starti, startj = tempi , tempj 
    
    returnalignedSeq1 = reverseString(alignedSeq1)
    returnalignedSeq2 = reverseString(alignedSeq2)
    

    return score, returnalignedSeq1, returnalignedSeq2 , starti, startj , endi, endj

def printSeq(str2):
    printStr = ""

    for eachchar in str2:
        if eachchar != 0:
            printStr = printStr + str(eachchar)
        elif eachchar == 0 :
            printStr = printStr + " "
    print printStr + "\n"
    
def groupIndelNoisyKmers(noisyReads, parameterRobot, typeOfOpt = "fast"):

    ### Setup 
    returnfmapping = [] 
    N = parameterRobot.N
    L = parameterRobot.L

    G=  parameterRobot.G
    threshold = parameterRobot.threshold
    liid = parameterRobot.liid
    

    folderName = parameterRobot.defaultFolder
    
    #clusterRounds, fingerPrint, clusterTreeSize = 3 , 6, int(N*L/G)
    clusterRounds, fingerPrint, clusterTreeSize = parameterRobot.clusterRounds , parameterRobot.fingerPrint, int(parameterRobot.clusterRatio*N*L/G)
    overhang = 5

    # Find fingerPrint
    K = 10

    print "clusterRounds, fingerPrint, clusterTreeSize",clusterRounds, fingerPrint ,clusterTreeSize
    
    specification = str(K)+'int8, 3int64'
    kmerList = np.zeros(N*(L-K+1)*2, dtype=specification)
    tempKmer = np.zeros(1, dtype=specification)
    
    endIndexArray = []
    print "len(noisyReads)", len(noisyReads)
    
    runningindex =0 
    for indexN in range(N):
        indexL = 0
        #print indexN, indexL, K 
        #print noisyReads
        while (noisyReads[indexN][indexL+K-1] !=0 ):
            tempKmer[0][0][:] = noisyReads[indexN][indexL:indexL+K]                
            tempKmer[0][1][0]  = indexN
            tempKmer[0][1][1] =  indexL
            tempKmer[0][1][2] =  runningindex 
            kmerList[runningindex] = tempKmer[0]
            
            runningindex += 1
            indexL = indexL +1 
            
        endIndexArray.append(runningindex - 1)

    activeKmerUbdd = runningindex -1
    activeKmerList = kmerList[0:activeKmerUbdd+1]
    
    ### Computing overlap using fingerprint
    toCompareList = [[] for i in range(N)]
    print "len(activeKmerList): ",activeKmerList[0]
    #activeKmerList = sorted(activeKmerList, cmp=multiplier(0,K))
    #activeKmerList = sorted(activeKmerList, key = itemgetterkk(range(0,K)))
    #activeKmerList = sorted(activeKmerList, key = itemgetter(slice(0)))
    activeKmerList.sort()
    print "len(activeKmerList): ", activeKmerList[0]
    
    
    toCompareList = fromCompareList(toCompareList, activeKmerList, liid, threshold)
    
    print "activeKmerList[0] : ", activeKmerList[0]
    print "toCompareList[0] : " , toCompareList[0]     
    
    ### Alignment and establish K-mer link graphs
    K = parameterRobot.K
    runningindex =0 
    endIndexArray = []
    endOfEachReadArr = []
    for indexN in range(N):
        indexL = 0
        while (noisyReads[indexN][indexL+K-1] !=0 ):
            runningindex += 1
            indexL = indexL +1 
            
        endOfEachReadArr.append(indexL+K-1)  
        endIndexArray.append(runningindex - 1)
        
    activeKmerUbdd = runningindex -1
    kmerlinkGraph =  [[] for i in range(activeKmerUbdd +1 ) ]
    
    
    
    
    print "Check the slow checking:"
    tkk = time.time()
    for i in range(N):
        for eachtargetRead in toCompareList[i]:
            
            j = eachtargetRead[0]
            startFingerPrint = eachtargetRead[1]
            endFingerPrint = eachtargetRead[2]
            
            
            #endNoisy1, endNoisy2 = 0, 0 
            #while (noisyReads[i][endNoisy1] != 0):
            #    endNoisy1 += 1

            #while (noisyReads[j][endNoisy2] != 0):
            #    endNoisy2 += 1
            endNoisy1, endNoisy2 = endOfEachReadArr[i],endOfEachReadArr[j]
                
            if canDoFast(startFingerPrint, endFingerPrint, parameterRobot):
                #scoretmp, returnalignedSeq2tmp ,returnalignedSeq1tmp  , startjtmp,  startitmp , endjtmp ,enditmp  = SWAlignmentBanded(startFingerPrint, endFingerPrint,  noisyReads[j][0:endNoisy2],noisyReads[i][0:endNoisy1], parameterRobot)
                score, returnalignedSeq2 ,returnalignedSeq1  , startj,  starti , endj ,endi  = SWAlignmentBanded(startFingerPrint, endFingerPrint,  noisyReads[j][0:endNoisy2],noisyReads[i][0:endNoisy1], parameterRobot)
            
            else:
                score, returnalignedSeq1, returnalignedSeq2 , starti, startj , endi, endj = SWAlignment(noisyReads[i][0:endNoisy1], noisyReads[j][0:endNoisy2], parameterRobot)
    
            '''
            if [scoretmp   , startjtmp,  startitmp , endjtmp ,enditmp ]!= [score  , startj,  starti , endj ,endi ]: 
                print "Important:::"
                print "[scoretmp, returnalignedSeq2tmp ,returnalignedSeq1tmp  , startjtmp,  startitmp , endjtmp ,enditmp ]"
                print [scoretmp  , startjtmp,  startitmp , endjtmp ,enditmp ]
                printSeq(returnalignedSeq1tmp)
                printSeq(returnalignedSeq2tmp)
                
                print "[score, returnalignedSeq2 ,returnalignedSeq1  , startj,  starti , endj ,endi ] "
                print [score  , startj,  starti , endj ,endi ]
                printSeq(returnalignedSeq1)
                printSeq(returnalignedSeq2)
            '''
                        
            if i != j and meetRequirement(score, returnalignedSeq1, returnalignedSeq2,starti, startj, endi, endj, threshold, liid, overhang , endNoisy1, endNoisy2):
                #establish Links 
                
                read1Start , read2Start =  endIndexArray[i-1] +1 + starti, endIndexArray[j-1]+1 + startj
                if i == 0 : 
                    read1Start =   starti
                if j == 0:
                    read2Start =  startj
                    
                runningindex1, runningindex2, dummyIndex = 0,0, 0 
                

                     
                while ( runningindex1+ starti + K < endNoisy1+1 and runningindex2+ startj +K < endNoisy2+1):
                    #print "dummyIndex, runningindex1 , K, runningindex2 ,endNoisy1, endNoisy2", dummyIndex, runningindex1 , K, runningindex2 ,endNoisy1, endNoisy2
                    if returnalignedSeq1[dummyIndex] ==returnalignedSeq2[dummyIndex]  :
                    
                        startOffset1 = starti + runningindex1
                        startOffset2 = startj + runningindex2

                        kmerlinkGraph[read1Start+runningindex1 ].append(read2Start+runningindex2)
                        kmerlinkGraph[read2Start+runningindex2].append(read1Start+runningindex1 )
                    
                    
                    if returnalignedSeq1[dummyIndex] != 0 :
                        runningindex1 += 1
                    if returnalignedSeq2[dummyIndex] != 0:
                        runningindex2 += 1

                    dummyIndex += 1
    
    print "taken: (sec)", time.time() - tkk                
    ### Form clusters and form fmapping
    tempfmapping = formClusterMapping(kmerlinkGraph)
    
    ### Format Return  
    returnfmapping = []
    for i in range(N):
        if i == 0 :
            startIndex = 0
        else:
            startIndex = endIndexArray[i-1] + 1
            
        for j in range(endIndexArray[i] - startIndex + 1 ):
            readNum = i 
            offset = j 
            clusterIndex = tempfmapping[j+startIndex][1]
            returnfmapping.append([clusterIndex, readNum, offset])
    
    
    logging.fmappingSave(returnfmapping, folderName)


    return returnfmapping

def formRoughMapping(kmerlinkGraph):
    tempfmapping = [] 
    for i in range(len(kmerlinkGraph)):
        tempfmapping.append([i, -1])
    
    runningClusterIndex =0 

    for i in range(len(kmerlinkGraph)):
        kmerindex = i
        clusterExist = False
        targetclusterIndex = -1 
        for j in kmerlinkGraph[i]:
            if tempfmapping[j][1] != -1:
                clusterExist = True
                targetclusterIndex = tempfmapping[j][1]
                break
        
        if clusterExist == True:
            tempfmapping[kmerindex][1] = targetclusterIndex
        else:
            tempfmapping[kmerindex][1] = runningClusterIndex
            runningClusterIndex += 1
    return tempfmapping


def formClusterMapping(kmerlinkGraph):
    returnMapping = []
    tmpList = []
    for i in range(len(kmerlinkGraph)):
        tmpList.append(clusterElem(i))
        
    for i in range(len(kmerlinkGraph)):
        for j in kmerlinkGraph[i]:
            union(tmpList[i], tmpList[j])
    
    fmapping = formatClusteringMap(tmpList)
    for eachitem in fmapping:
        clusterindex = eachitem[0]
        for eachsubitem in eachitem[1]:
            returnMapping.append([eachsubitem.id, clusterindex])
            
    returnMapping = sorted(returnMapping)
    
    return returnMapping  
def groupNoisyKmers(noisyReads,parameterRobot, typeOfOpt= 'fast'):
    fmapping = []
    N = parameterRobot.N
    L = parameterRobot.L
    K = parameterRobot.K
    G=  parameterRobot.G
    threshold = parameterRobot.threshold
    liid = parameterRobot.liid
    

    #clusterRounds, fingerPrint, clusterTreeSize = 3 , 6, int(N*L/G)
    clusterRounds, fingerPrint, clusterTreeSize = parameterRobot.clusterRounds , parameterRobot.fingerPrint, int(parameterRobot.clusterRatio*N*L/G)

    
    print "clusterRounds, fingerPrint, clusterTreeSize",clusterRounds, fingerPrint ,clusterTreeSize
    
    specification = str(K)+'int8, 3int64'
    kmerList = np.empty(N*(L-K+1), dtype=specification)
    tempKmer = np.zeros(1, dtype=specification)
    
    folderName = parameterRobot.defaultFolder
    
   # list of K mers
    for indexN in range(N):
        for indexL in range(L- K+1):
            tempKmer[0][0][:] = noisyReads[indexN][indexL:indexL+K]                
            tempKmer[0][1][0]  = indexN
            tempKmer[0][1][1] =  indexL
            tempKmer[0][1][2] =  indexN*(L-K+1)+indexL

            kmerList[indexN*(L-K+1)+indexL] = tempKmer[0]

    
    clusterList = []
    for index in range(N*(L-K+1)):
        clusterList.append(clusterElem(index))
   
    if typeOfOpt == 'fast':
        print "Fast Clustering"
        fastClusteringAlgo(N, L, K, kmerList, clusterList, noisyReads, threshold, clusterRounds, fingerPrint,clusterTreeSize, liid)
    elif typeOfOpt == 'pair':
        print "Pairwise Clustering"
        pairwiseCompareKmers(kmerList,clusterList, N, L, K, threshold, liid)

    fmapping = formatClusteringMap(clusterList)       
    
    returnfmapping = [] 
    
    checksum = 0
    for eachitem in fmapping : 
        #print "Group ", eachitem[0] , len(eachitem[1])
        checksum = checksum + len(eachitem[1])
        for eachsubitem in eachitem[1]:
            #print "kmer id, read num", eachsubitem.id, int( eachsubitem.id / (L-K+1) )
            returnfmapping.append([eachitem[0], int(eachsubitem.id/(L-K+1)), int(np.mod(eachsubitem.id, L-K+1))])
            
    logging.fmappingSave(returnfmapping, folderName)
    
    print "correctNum", N*(L-K+1)
    print "checksum" , checksum
    return returnfmapping              

# as a check
def pairwiseCompareKmers(kmerList,clusterList, N, L, K, threshold, liid):
    total = len(kmerList)

    for index1 in range(total -1 ):
        for index2 in range(index1+1, total):
            if matched(kmerList[index1][0], kmerList[index2][0], threshold, liid):
                union(clusterList[index1], clusterList[index2])
                

                



        
 
