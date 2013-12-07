import numpy as np
import os
import cluster
import bridgeResolve
import common
import cleaner

def transformToFASTA(infile, outfile):

    f2 = open(infile, 'r')

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
    fout= open(outfile, 'w')
    
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

    f2.close()




def outputToFastaFiles(reconstructedGenome, motherGenome, parameterRobot):
    frecov = open(parameterRobot.defaultFolder+"rec.txt", 'w')
    for eachbase in reconstructedGenome:
        frecov.write(str(eachbase))
    frecov.close()
    
    
    transformToFASTA(parameterRobot.defaultFolder+"rec.txt", parameterRobot.defaultFolder+"rec.fasta")
    transformToFASTA(parameterRobot.defaultFolder+"UnitTest_motherGen.txt", parameterRobot.defaultFolder+"UnitTest_motherGen.fasta")
    
    os.chdir("../gepard-1.30")
    os.system("./gepardcmd.sh -seq1 ../indelAssembler/"+parameterRobot.defaultFolder+"UnitTest_motherGen.fasta -seq2 ../indelAssembler/"+parameterRobot.defaultFolder+"rec.fasta -matrix matrices/edna.mat -outfile "+"../indelAssembler/"+parameterRobot.defaultFolder+"jane.png")
    os.chdir("../indelAssembler")
    
def itemgetterkk(items):
    if len(items) == 1:
        item = items
        def g(obj):
            return obj[item]
    else:
        def g(obj):
            return tuple(obj[item] for item in items)
    return g

def arrangeSeqBasedOnRef(motherGenome, reconstructedGenome, parameterRobot):
    
    # 1) Extract components and sort
    G = parameterRobot.G 
    fingerPrint = 20
    searchNumber = 1000

    allFingerPrint = np.zeros(( G-fingerPrint+1) *(fingerPrint+1) , dtype= np.int32).reshape(G-fingerPrint+1, fingerPrint+1)
    
    for i in range(G-fingerPrint +1 ):
        allFingerPrint[i][0:-1] = motherGenome[i:i+fingerPrint]
        allFingerPrint[i][-1] = i
        
    allFingerPrint = sorted(allFingerPrint,key = itemgetterkk(range(0, fingerPrint+1)))

    # 2) Search for appropriate hash items 
    searchItemList = []
    starterIndex = -1
    oneItemList = [] 
    hashSearch = np.zeros(fingerPrint+1 , dtype = np.int32)
    
    for i in range(searchNumber):
        foundIndex = int( min(i*len(reconstructedGenome)/float(searchNumber), len(reconstructedGenome)- fingerPrint -1  ) )
        hashSearch[0:fingerPrint+1] = reconstructedGenome[foundIndex: foundIndex +fingerPrint+1]
        hashSearch[-1] = -10

        indexInFPList = bridgeResolve.bisectkk(allFingerPrint, hashSearch) -1

        tmpItem = [] 
        #print indexInFPList
        while (0<=indexInFPList < len(allFingerPrint) and np.array_equal(allFingerPrint[indexInFPList][0:-1],hashSearch[0:-1])):
            tmpItem.append(allFingerPrint[indexInFPList][-1])
            indexInFPList = indexInFPList - 1 
        print "tmpItem", tmpItem
        searchItemList.append(tmpItem)
        if len(tmpItem) == 1 : 
            oneItemList.append(i)
        #else:
        #    print len(tmpItem)
    
    starterIndexList = []      
    print "oneItemList",oneItemList
    print "searchItemList", searchItemList
    
    for dummyIndex in oneItemList:
        tmp = int(searchItemList[dummyIndex][0] - dummyIndex*len(reconstructedGenome)/float(searchNumber))
        #print tmp
        if tmp < 0 : 
            tmp = tmp + G 
        starterIndexList.append(tmp)
    
    print "starterIndexList", starterIndexList
    starterIndex = int(np.median(starterIndexList))
    starterIndex = G - starterIndex
    print starterIndex
    
    # 3) Arrange the starter of reconstructed guy
    newArranged = np.zeros(len(reconstructedGenome))
    
    newArranged[0:len(reconstructedGenome)-starterIndex] = reconstructedGenome[starterIndex:len(reconstructedGenome)]
    newArranged[len(reconstructedGenome)-starterIndex:len(reconstructedGenome)]= reconstructedGenome[0:starterIndex]
    
    return newArranged 
    
def findMismatchNumber(motherGenome, reconstructedGenome):
    ## Fill in the algorithm here 
    # 1)  Sliding window matching and count
    # 2)  Make sure to not include the loop around thing. 
    
    W= 50 
    
    runningI = 0 
    runningJ = 0 
    G = len(motherGenome)
    

    totalScore = 0 
    parameterRobot = common.parameterRobot()
    counter = 0
    tmprunningI = -1
    print "----------------------"
    overallSum = 0
    while runningI < len(motherGenome) -1 and tmprunningI != runningI:
        score, returnalignedSeq1, returnalignedSeq2 , starti, startj , endi, endj = cleaner.SWAlignmentFixRef(motherGenome[runningI:runningI+W] , reconstructedGenome[runningJ:runningJ+W], parameterRobot)
   #     print "Mother: "
   #     cluster.printSeq(returnalignedSeq1)
   #     print "Reconstructed: "
   #     cluster.printSeq(returnalignedSeq2)
   #     print "----------"
        
        tmprunningI = runningI
        
        runningI = endi  + runningI 
        runningJ = endj +  runningJ  
        
        counter  = counter +1 
        totalScore = totalScore + score 
        #print "score", score
        overallSum += score
    print "----------------------"
    numberOfMismatch = max((G - totalScore)/2, 0 )
    print overallSum / counter

    return numberOfMismatch 

def arrangeSeqBasedOnRefEasy(motherGenome, reconstructedGenome, parameterRobot):
    runningI = 0 
    runningJ = 0
    G = len(motherGenome)
    
    W = 200
    totalScore = 0 
    parameterRobot = common.parameterRobot()
    counter = 0
    scoreList = []
    
    while runningI < len(motherGenome) -1 and counter < 10000:
        score, returnalignedSeq1, returnalignedSeq2 , starti, startj , endi, endj = cleaner.SWAlignment(motherGenome[runningI:runningI+W] , reconstructedGenome[runningJ:runningJ+W], parameterRobot)
        scoreList.append([score,starti, startj , endi, endj,runningI ])
        
        runningI = runningI + W
        counter  = counter +1 
        
    scoreList = sorted(scoreList)

    testCases = [] 
    for i in range(2):
        score, starti, startj , endi, endj,runningI =  scoreList[-(i+1)]
    
        newGenome = np.zeros(G, dtype = np.int32)
        myindex = runningI +starti - startj-runningJ
                
        if myindex > 0 :
            newGenome[0:G-myindex] = motherGenome[myindex:G] 
            newGenome[G-myindex:G] = motherGenome[0:myindex]
        else:
            myposindex = G+ myindex
            newGenome[0:G-myposindex] = motherGenome[myposindex:G]
            newGenome[G-myposindex:G] = motherGenome[0:myposindex]

        testCases.append(newGenome)
        
    return testCases 

def subAlignCompare(reconstructedGenome, motherGenome,parameterRobot):
    # Output and print to .fasta files 
    outputToFastaFiles(reconstructedGenome, motherGenome, parameterRobot)    

    
    print "len(motherGenome)", len(motherGenome)
    # Find the hashing starter 
    testCases = arrangeSeqBasedOnRefEasy(motherGenome, reconstructedGenome, parameterRobot)
    
    # Do the iterative alignment 
    finalAns= -1
    tmpMin = len(motherGenome)
    for motherGenome, i in zip(testCases,range(len(testCases))):
        numberOfMismatch = findMismatchNumber(motherGenome, reconstructedGenome)
        if numberOfMismatch < tmpMin:
            tmpMin = numberOfMismatch 
            finalAns = i 
            
   
    p, G = parameterRobot.p, parameterRobot.G
    
    if tmpMin < 2*p*G:
        check = True 
    else:
        check = False
    
    print "tmpMin: ",tmpMin 
    return check , tmpMin 


