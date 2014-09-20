
import numpy as np 
import csv 
import graphForm
import bridgeResolve
import numericalCompute
import math

### Logging raw Data
def rawDataSave(filename, motherGen, reads, noisyReads, typeOfSave = 'a'): 
    if typeOfSave == 'a':
        f = open(filename + "_motherGen.txt",'w')
    
        for eachindex in range(len(motherGen)):
            f.write(str(motherGen[eachindex]))
        
        f.close()

    if typeOfSave == 'a':
        f = open(filename + "_reads.txt" , 'w')
    
        for eachread in reads:
            for eachcharacter in eachread: 
                f.write(str(eachcharacter))            
            f.write("\n")
    
        f.close()

    if typeOfSave == 'a' or typeOfSave == 'n':    
        f = open(filename + "_noisyReads.txt" , 'w')
        
        
        for eachread in noisyReads:
            for eachcharacter in eachread: 
                f.write(str(eachcharacter))            
            f.write("\n")
    
        f.close()    


def rawDataLoad(filename,G,N,L, typeOfSave = 'a'):
    motherGen, reads, noisyReads = "", "", ""
    print G, N, L
    
    if typeOfSave == 'a'or typeOfSave == "dn" or typeOfSave == 'd':
        # Logging mother genome    
        f = open(filename + "_motherGen.txt",'r')
        
        motherGenStr = f.read()
        motherGen = np.zeros(G,dtype = np.int8)
        
        
        for eachindex in range(len(motherGen)):
            motherGen[eachindex] = int( motherGenStr[eachindex] )    
        f.close()
    
    if typeOfSave == 'a' or typeOfSave == 'd' or typeOfSave == "dn" :
        # Logging noiseless reads
        f = open(filename + "_reads.txt" , 'r')
        
        reads = np.zeros(N*L*2, dtype = np.int8).reshape(N,2*L)
        
        eachread = f.readline()
        readindex = 0
        while (len(eachread) > 0 ):
            for characterindex in range(len(eachread)): 
                if eachread[characterindex] != '\n':
                    #print characterindex, readindex, eachread
                    reads[readindex][characterindex] = eachread[characterindex]
                    
            eachread = f.readline()
            readindex = readindex +1 
                
        f.close()
    
    if typeOfSave == 'a' or typeOfSave == 'n':
        # Logging noisy reads
        f = open(filename + "_noisyReads.txt" , 'r')
            
        noisyReads = np.zeros(N*L, dtype = np.int8).reshape(N,L)
        
        eachread = f.readline()
        readindex = 0
        while (len(eachread) > 0 ):
            for characterindex in range(len(eachread)): 
                if eachread[characterindex] != '\n':
                    noisyReads[readindex][characterindex] = eachread[characterindex]
                    
            eachread = f.readline()
            readindex = readindex +1 
                
        f.close()
        
    if typeOfSave == 'd' or typeOfSave == "dn":
        # Logging noisy indel reads
        f = open(filename + "_noisyReads.txt" , 'r')
            
        noisyReads = np.zeros(N*2*L, dtype = np.int8).reshape(N,2*L)
        
        eachread = f.readline()
        readindex = 0
        while (len(eachread) > 0 ):
            for characterindex in range(len(eachread)): 
                if eachread[characterindex] != '\n':
                    noisyReads[readindex][characterindex] = eachread[characterindex]
                    
            eachread = f.readline()
            readindex = readindex +1 
                
        f.close()
        
    
    return motherGen, reads, noisyReads

### Logging fmapping before branch clearing 
def fmappingSave(returnfmapping, folderName):
    ofile  = open(folderName+'clusteredGroup.csv', "wb")
    writer = csv.writer(ofile)   
    writer.writerow(["Group number ","read-id", "offset-id"])
    
    for eachitem in returnfmapping:
        writer.writerow([eachitem[0],eachitem[1], eachitem[2]])
        
        
def fmappingLoad(filename):
    ifile = open(filename, 'r')
    myreader = csv.reader(ifile)
    
    returnfmapping = []
    for row in myreader:
        if row[0] != 'Group number ':
            returnfmapping.append([int(row[0]), int(row[1]), int(row[2])])
    
    return returnfmapping
    
### Logging fmapping after branch clearing
def fmapfusedSave(returnfmapping,folderName= ""):
    ofile  = open(folderName + 'clusteredGroup2.csv', "wb")
    writer = csv.writer(ofile)   
    writer.writerow(["Group number ","read-id", "offset-id", "fused or not", "prevGroupid"])
        # OutputFormat :     Gp id , read #, offset #, fusedOrNot, prevGroup id 
    for eachitem in returnfmapping:
        writer.writerow([eachitem[0],eachitem[1], eachitem[2], eachitem[3], eachitem[4]])
        
def fmapfusedLoad(filename):
    ifile = open(filename, 'r')
    myreader = csv.reader(ifile)
    
    returnfmapping = []
    for row in myreader:
        if row[0] != 'Group number ':
            returnfmapping.append([int(row[0]), int(row[1]), int(row[2]), bool(row[3] == "True"), int(row[4])])

    return returnfmapping  

### Logging Graph Structure for MB   
def storeGraph(G,folderName ):    
    ofile  = open(folderName+'basicMapping.csv', "wb")
    mywriter = csv.writer(ofile)  
    mywriter.writerow(["nodeIndex",  "prev List","next List", "next List length"])

    for eachnode in G:
        nextList = []
        prevList = []
        for eachnextNode in eachnode.listOfNextNodes:
            nextList.append(eachnextNode.nodeIndex)
        
        for eachprevNode in eachnode.listOfPrevNodes:
            prevList.append(eachprevNode.nodeIndex)
        mywriter.writerow([eachnode.nodeIndex, prevList,nextList , len(eachnode.nodeIndexList)])
        
    
    ofile.close()
    
    ofile = open(folderName+'seqMapping.txt', "w")
    for eachnode in G:
        ofile.write(str(eachnode.nodeIndex))
        ofile.write(str(eachnode.nodeIndexList))
        ofile.write("\n")
        
    ofile.close()
    

def transformToMBGraph(basicList,seqList, typeOfGraph):
    nodeIndexing = []
    G = []
    for eachitem in seqList:
        temp = []
        if typeOfGraph == 'simple':
            temp = graphForm.condensedNode(eachitem[0])
        elif typeOfGraph == 'MB':
            temp = bridgeResolve.MBCondensedNode(eachitem[0])
            
        temp.nodeIndexList = eachitem[1]
        G.append(temp)
        nodeIndexing.append(eachitem[0])
        
    
    
    for index in range(len(basicList)):
        currentNode = G[nodeIndexing.index(basicList[index][0])]
        
        for eachprevnodeindex in basicList[index][1]:
            prevNode = G[nodeIndexing.index(eachprevnodeindex)]
            
            if not currentNode in prevNode.listOfNextNodes:
                prevNode.listOfNextNodes.append(currentNode)
                
            if not prevNode in currentNode.listOfPrevNodes:
                currentNode.listOfPrevNodes.append(prevNode)
            
        for eachnextnodeindex in basicList[index][2]:
            nextNode = G[nodeIndexing.index(eachnextnodeindex)]
            
            if not currentNode in  nextNode.listOfPrevNodes:
                nextNode.listOfPrevNodes.append(currentNode)
            if not nextNode in currentNode.listOfNextNodes:
                currentNode.listOfNextNodes.append(nextNode)
                
    if typeOfGraph == "MB":
        if G[0].naiveForm == True:
            for eachnode in G:
                eachnode.initOverlap()
                
    return G
    
    
def loadGraph(basicmapping, seqmapping,typeOfGraph):
    infile = open(basicmapping, "r")
    
    myreader = csv.reader(infile)
    
    basicList = []
    for eachrow in myreader:
        if eachrow[0] != "nodeIndex":
            prevList , nextList = [], []
            prevList = eachrow[1][1:-1].split(',')
            nextList = eachrow[2][1:-1].split(',')
            for index in range(len(prevList)):
                prevList[index] = int(prevList[index])
            for index in range(len(nextList)):
                print nextList[index]
                nextList[index] = int(nextList[index])
                
            basicList.append([int(eachrow[0]),prevList,nextList, int(eachrow[3])])
    
    infile.close()
    infile = open(seqmapping, 'r')
    temp = infile.readline()
    seqList = []
    
    
    while (len(temp) > 0):
        
        list1 = temp.split('[')
        nodeIndex = int(list1[0])
        
        #print list1
        list2 = list1[1].split(',')
        
        #print list2
        if len(list2) > 1:

            nodeIndexList = []
            for index in range(len(list2) - 1):
                nodeIndexList.append(int(list2[index]))
            nodeIndexList.append(int(list2[index+1][0:-2]))
                
        elif len(list2) <= 1:
            nodeIndexList = []
            nodeIndexList.append(int(list2[0][0:-2]))
        
        seqList.append([nodeIndex, nodeIndexList])

        temp = infile.readline()
            
        
    infile.close()   
    
    G = transformToMBGraph(basicList,seqList, typeOfGraph )

    
    return G

        
        
### Batch Processing
def savingLNKFile(folderName = ""):
    fout = open(folderName+ "dataPoints.csv", 'wb')
    
    mywriter = csv.writer(fout)
    mywriter.writerow(["G", "N", "L", "p", "epsilon", "K", "liid", "threshold", "NKcov", "Nbridge", "Ncov", "Nratio" ,"numberOfClusterRounds","brachingDepth", "bridgingDepth", "msaWidth" , "Nbridgenoiseless", "ratioNoiseless" , "clusterRounds", "fingerPrint", "clusterRatio"])
    
    
    
    for index in range(5):
        
        G, L, p= 50000, 200, 0.015
        linter , ltriple = 100, 10 
        L = L - index*20
        
        epsilon = 0.05
    
    
        ### Noisy Compute
        calculator =numericalCompute.thresholdCompute(p, G)
        liid , threshold = calculator.findRoot()
        threshold = threshold - 2
        
        K =  int(liid*1.3) 

        calculator = numericalCompute.Ncompute(G,L,epsilon)
        Ncov = int(calculator.findRoot()) 
                 
        calculator = numericalCompute.Ncompute(G,L- K,epsilon/3)
        NKcov = int(calculator.findRoot()) 

        Nbridge = int ( G*math.log(3/epsilon)/float(L-max(linter, ltriple) - 10) )
        
        #N = int ( max(NKcov, Nbridge)*1.5) 
        N = int ( max(NKcov, Nbridge)*1.5) 
        
        numberOfClusterRounds,brachingDepth,bridgingDepth,msaWidth = 6 , liid*2/3, liid*2/3 , liid*2/3
        
        
        ### Noiseless Compute
        calculator = numericalCompute.Ncompute(G,L,epsilon)
        Nbridgenoiseless = int ( G*math.log(3/epsilon)/float(L-max(linter, ltriple)  ) )
        
        Noiseless = max(Nbridgenoiseless, Ncov)
        ratioNoiseless = Noiseless/ float(Ncov)
        
        clusterRounds, fingerPrint, clusterRatio = 2 , 6 , 1 
        
        mywriter.writerow([G, N, L, p, epsilon, K, liid, threshold,NKcov, Nbridge, Ncov, N/float(Ncov),numberOfClusterRounds,brachingDepth,bridgingDepth,msaWidth, Nbridgenoiseless,ratioNoiseless,clusterRounds, fingerPrint, clusterRatio])
    
    fout.close()
    

def loadingLNKFile(folderName = ""):    
    fin = open(folderName+"dataPoints.csv", 'r')
    myreader = csv.reader(fin)
    
    listOfNLKDataPts = []
    for eachrow in myreader:
        if eachrow[0] != "G":
            G, N, L, p, epsilon, K, liid, threshold, NKcov, Nbridge, Ncov, Nratio,numberOfClusterRounds,brachingDepth,bridgingDepth,msaWidth  = int(eachrow[0]), int(eachrow[1]), int(eachrow[2]), float(eachrow[3]), float(eachrow[4]), int(eachrow[5]), int(eachrow[6]), int(eachrow[7]),  int(eachrow[8]), int(eachrow[9]), int(eachrow[10]), float(eachrow[11]), int(eachrow[12]),  int(eachrow[13]),  int(eachrow[14]),  int(eachrow[15])
            clusterRounds, fingerPrint, clusterRatio = int(eachrow[18]), int(eachrow[19]), float(eachrow[20])
            listOfNLKDataPts.append([G, N, L, p, epsilon, K, liid, threshold, NKcov, Nbridge, Ncov, Nratio,numberOfClusterRounds,brachingDepth,bridgingDepth,msaWidth,clusterRounds, fingerPrint, clusterRatio  ])
            
    fin.close()
    
    return listOfNLKDataPts
        
class parameterObj(object):
    def __init__(self,G=0, N=0, L=0, p=0, epsilon=0, K=0, liid=0, threshold=0,NKcov=0, Nbridge=0, Ncov=0, ratio=0, numberOfClusterRounds=3,brachingDepth=20,bridgingDepth=20,msaWidth=20,defaultFolder = "", clusterRounds = 2, fingerPrint = 6, clusterRatio = 1, indel= False, editsub=-10, editins=-1, editdel = -1, editmatch =1  ):
        self.G, self.N, self.L, self.p, self.epsilon, self.K, self.liid, self.threshold,self.NKcov, self.Nbridge, self.Ncov, self.ratio, self.numberOfClusterRounds,self.brachingDepth,self.bridgingDepth,self.msaWidth = G, N, L, p, epsilon, K, liid, threshold,NKcov, Nbridge, Ncov, ratio, numberOfClusterRounds,brachingDepth,bridgingDepth,msaWidth 
        self.defaultFolder = defaultFolder
        self.clusterRounds , self.fingerPrint, self.clusterRatio = clusterRounds, fingerPrint, clusterRatio
        self.indel = indel
        self.editsub, self.editins, self.editdel , self.editmatch = editsub, editins, editdel , editmatch

def logBatch(resultList):
    fout = open("batchResults.csv", 'r')
    mywriter = csv.writer(fout)
    
    mywriter.writerow(["N", "L", "roundNum" , "numMistakes", "success"])
    for eachitem in resultList:
        mywriter.writerow(eachitem)
    
    fout.close()
    
    
def transformReadsToFasta(filename, outputFilename):
    f = open(filename, 'r')
    fout = open(outputFilename, 'w')
    
    
    temp = f.readline()
    
    runningindex = 0
    while (len(temp) > 0 ):
        fout.write(">Seq "+ str(runningindex)+ "\n")
        for eachcharacter in temp:
            if eachcharacter == '1' :
                fout.write('A')
            elif eachcharacter == '2' :
                fout.write('C')
            elif eachcharacter == '3' :
                fout.write('G')
            elif eachcharacter == '4' :
                fout.write('T')
                
        runningindex = runningindex + 1
        
        fout.write("\n")
        temp = f.readline()
    
    
    f.close()
    fout.close()
        
    
    
    

def savingGenomeSegmentFile(folderName):
    fin = open(folderName+ "genomeStat.csv", 'r')
    myreader = csv.reader(fin)
   
    dataList = []
    firstTime = True
    for eachrow in myreader:
        if firstTime:
            firstTime = False   
        else:
            dataList.append([int(eachrow[0]), int(eachrow[1]), int(eachrow[2]), int(eachrow[3]), int(eachrow[4])  , int(eachrow[5]), int(eachrow[6])]) 
    
    fin.close()
    
    
    fout = open(folderName+ "dataPoints.csv", 'wb')
    
    
    
    mywriter = csv.writer(fout)
    #mywriter.writerow(["G", "N", "L", "p", "epsilon", "K", "liid", "threshold", "NKcov", "Nbridge", "Ncov", "Nratio" ,"numberOfClusterRounds","brachingDepth", "bridgingDepth", "msaWidth" , "Nbridgenoiseless", "ratioNoiseless" , "clusterRounds", "fingerPrint", "clusterRatio", "approx repeat", "Lcrit", "approxinter"])
    mywriter.writerow(["G", "N", "L", "p", "epsilon", "K", "liid", "threshold", "NKcov", "Nbridge", "Ncov", "Nratio" ,"numberOfClusterRounds","brachingDepth", "bridgingDepth", "msaWidth" , "Nbridgenoiseless", "ratioNoiseless" , "clusterRounds", "fingerPrint", "clusterRatio", "startIndex", "endIndex", "approxinter"])
    
    

    for index in range(len(dataList)):
        

        linter , ltriple = dataList[index][3], dataList[index][4]
        G, L, p= dataList[index][1] - dataList[index][0], dataList[index][5], 0.015

        
        epsilon = 0.05
    
    
        ### Noisy Compute
        calculator =numericalCompute.thresholdCompute(p, G)
        liid , threshold = calculator.findRoot()

        
        #K =  liid*2
        K = 600

        calculator = numericalCompute.Ncompute(G,L,epsilon)
        Ncov = int(calculator.findRoot()) 
                 
        calculator = numericalCompute.Ncompute(G,L- liid,epsilon)
        NKcov = int(calculator.findRoot()) 

        Nbridge = int ( G*math.log(9/epsilon)/float(L-max(linter, ltriple) - liid) )
        
        N = max(NKcov, Nbridge)
        
        #numberOfClusterRounds,brachingDepth,bridgingDepth,msaWidth = 6 , liid*1/3, liid*1/3 , liid*2/3
        numberOfClusterRounds,brachingDepth,bridgingDepth,msaWidth = 6 , liid*2/3, liid*2/3 , liid*2/3
        
        
        ### Noiseless Compute
        calculator = numericalCompute.Ncompute(G,L,epsilon)
        Nbridgenoiseless = int ( G*math.log(3/epsilon)/float(L-max(linter, ltriple)  ) )
        
        Noiseless = max(Nbridgenoiseless, Ncov)
        ratioNoiseless = Noiseless/ float(Ncov)
        
        clusterRounds, fingerPrint, clusterRatio = 2 , 6 , 1 
        
        #mywriter.writerow([G, N, L, p, epsilon, K, liid, threshold,NKcov, Nbridge, Ncov, N/float(Ncov),numberOfClusterRounds,brachingDepth,bridgingDepth,msaWidth, Nbridgenoiseless,ratioNoiseless,clusterRounds, fingerPrint, clusterRatio,dataList[index][2], dataList[index][3] ,dataList[index][6]])
        mywriter.writerow([G, N, L, p, epsilon, K, liid, threshold,NKcov, Nbridge, Ncov, N/float(Ncov),numberOfClusterRounds,brachingDepth,bridgingDepth,msaWidth, Nbridgenoiseless,ratioNoiseless,clusterRounds, fingerPrint, clusterRatio,dataList[index][0], dataList[index][1] ,dataList[index][6]])
    
    fout.close()
    
    
    
def loadingGenomeSegmentFile(folderName = ""):    
    fin = open(folderName+"dataPoints.csv", 'r')
    myreader = csv.reader(fin)
    
    listOfNLKDataPts = []
    for eachrow in myreader:
        if eachrow[0] != "G":
            G, N, L, p, epsilon, K, liid, threshold, NKcov, Nbridge, Ncov, Nratio,numberOfClusterRounds,brachingDepth,bridgingDepth,msaWidth  = int(eachrow[0]), int(eachrow[1]), int(eachrow[2]), float(eachrow[3]), float(eachrow[4]), int(eachrow[5]), int(eachrow[6]), int(eachrow[7]),  int(eachrow[8]), int(eachrow[9]), int(eachrow[10]), float(eachrow[11]), int(eachrow[12]),  int(eachrow[13]),  int(eachrow[14]),  int(eachrow[15])
            clusterRounds, fingerPrint, clusterRatio, startIndex , endIndex = int(eachrow[18]), int(eachrow[19]), float(eachrow[20]), int(eachrow[21]), int(eachrow[22])
            listOfNLKDataPts.append([G, N, L, p, epsilon, K, liid, threshold, NKcov, Nbridge, Ncov, Nratio,numberOfClusterRounds,brachingDepth,bridgingDepth,msaWidth,clusterRounds, fingerPrint, clusterRatio,startIndex , endIndex   ])
            
    fin.close()
    
    return listOfNLKDataPts


def generateGenomeStatFile():
    fout  = open("genomeStat.csv", 'wb')
    mywriter = csv.writer(fout)
    mywriter.writerow(["Segment in",    "Segment in",     "lrepeat",    "linter",    "ltriple"])
    
    mywriter.writerow([1210000   , 1220000  ,  1784 ,   20,    770])    
    mywriter.writerow([1205000 ,   1215000,   1784  ,  20  ,  770])

    
    for eachindex in range(0, 1400000, 5000):
        mywriter.writerow([eachindex,eachindex + 10000,   1784  ,  20  ,  770])
    fout.close()



def fetchResults(header, numFile, roundNum):
    for index1 in range(numFile):
        filename = header + "sample_point_"+str(index1) + "\\result.txt"
        f = open(filename, 'r')
        temp = f.read()
        print "sample_point_"+str(index1) 
        print temp
        print "---------------------------------"
        f.close()
#savingLNKFile()
#transformReadsToFasta("sample_point_0//round_0//UnitTest_noisyReads.txt", "abc.fasta")    


#savingGenomeSegmentFile("")   
    
    
    
    