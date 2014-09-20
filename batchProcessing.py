import assemblerMain
import logging
import time
import os
import cProfile


def batchProcessingUnitTest(numberOfRounds,parameterRobot,snpRate, typeOfGen , detail):
    print "Batch Processing Unit Test : "   
    #snpRate, typeOfGen, detail = 0.001 ,'d', "genome.fasta-10000-20000" 
    
    resultList = []

    oldName = parameterRobot.defaultFolder
    
    
    for index in range(numberOfRounds):
        tempfile = open(oldName+"/result.txt", 'a')
        os.system("mkdir "+ oldName + "/round_" +str(index) )
        
        parameterRobot.defaultFolder =  oldName + "/round_" +str(index) + "/"
        #try:
        numMistakes, success = assemblerMain.runAssembler(snpRate, typeOfGen, detail,parameterRobot)
       # except:
       #     numMistakes, success= parameterRobot.G, False
            
        print [numMistakes, success]
        resultList.append([numMistakes, success])
            
        print "index, resultList",index, resultList

        tempfile.write(str(parameterRobot.N) +" " + str(parameterRobot.L) +" " +  str(index) +" " + str(numMistakes) + ", " + str(success) + "\n")
        tempfile.close()
    # Print result
    print "results( N, L )",parameterRobot.N, parameterRobot.L 
    
    for eachitem in resultList:
        print "eachitem",  eachitem 
        
    
    
def batchProcessingLNKTest():
    print "Batch Processing LNK Test"
    headerName = "synthetic_reads/"
    os.system("mkdir "+headerName)
    
    logging.savingLNKFile(headerName)
    listOfNLKDataPts = logging.loadingLNKFile(headerName)
    numberOfRounds = 1
    listOfNLKDataPts = [listOfNLKDataPts[0]]

    
    
    for testPoint,roundNum in  zip(listOfNLKDataPts, range(len(listOfNLKDataPts))):
        folderName = headerName+"sample_point_" + str(roundNum)
        os.system("mkdir " + folderName )
        

        snpRate, typeOfGen , detail = 0.001,'m', "500-200-50" 
        [G, N, L, p, epsilon, K, liid, threshold,NKcov, Nbridge, Ncov, ratio, numberOfClusterRounds,brachingDepth,bridgingDepth,msaWidth,clusterRounds, fingerPrint, clusterRatio ] = testPoint
        parameterRobot = logging.parameterObj(G, N, L, p, epsilon, K, liid, threshold,NKcov, Nbridge, Ncov, ratio, numberOfClusterRounds,brachingDepth,bridgingDepth,msaWidth,folderName,clusterRounds, fingerPrint, clusterRatio )
        parameterRobot.indel = True
 
        batchProcessingUnitTest(numberOfRounds,parameterRobot,snpRate, typeOfGen , detail)
    
    
    
def batchProcessingGenomeSegTest():
    print "Batch Processing LNK Test"
    headerName = "largeScaleTest/"
    os.system("mkdir "+headerName)
    
    logging.savingGenomeSegmentFile(headerName)
    listOfNLKDataPts = logging.loadingGenomeSegmentFile(headerName)
    numberOfRounds = 1
    
    for testPoint,roundNum in  zip(listOfNLKDataPts, range(len(listOfNLKDataPts))):
        folderName = headerName+"sample_point_" + str(roundNum)
        os.system("mkdir " + folderName )
                
        [G, N, L, p, epsilon, K, liid, threshold,NKcov, Nbridge, Ncov, ratio, numberOfClusterRounds,brachingDepth,bridgingDepth,msaWidth,clusterRounds, fingerPrint, clusterRatio, startindex, endindex ] = testPoint
        print "G", G
        parameterRobot = logging.parameterObj(G, N, L, p, epsilon, K, liid, threshold,NKcov, Nbridge, Ncov, ratio, numberOfClusterRounds,brachingDepth,bridgingDepth,msaWidth,folderName,clusterRounds, fingerPrint, clusterRatio )
        snpRate, typeOfGen , detail = 0.001,'d', "genome.fasta-"+str(startindex)+"-"+str(endindex) 
        parameterRobot.indel = True
        temptime = time.time()
        batchProcessingUnitTest(numberOfRounds,parameterRobot,snpRate, typeOfGen , detail)
        print "time per sample point ",  time.time() - temptime
        

t0 = time.time()
#batchProcessingGenomeSegTest()
#logging.generateGenomeStatFile()    
#batchProcessingLNKTest()

cProfile.run("batchProcessingLNKTest()")
print "Time (Sec)", time.time() - t0






