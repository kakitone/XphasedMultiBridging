import numpy as np
import random 
import logging

def generateData( typeOfGen,detail, parameterRobot): 
    N, G, L,p = parameterRobot.N, parameterRobot.G, parameterRobot.L,parameterRobot.p 
    
    motherGen = np.zeros(G, dtype = np.int8 )
    
    reads = np.zeros(N*L, dtype = np.int8).reshape(N,L) 
    noisyReads = np.zeros(N*L, dtype = np.int8).reshape(N,L)
    detailArr = detail.split('-')
    
    for index in range(G): 
        motherGen[index] = random.randint(1,4) 

    if typeOfGen == 'r' :         
        
        longestRepeatlength = int(detailArr[0])
        
        repeatLoc1 = random.randint(0, G-2*longestRepeatlength)
        repeatLoc2 = random.randint(repeatLoc1 +longestRepeatlength ,G-longestRepeatlength )
        
        print "repeatLoc1, repeatLoc2, detail", repeatLoc1, repeatLoc2, detail
        
        
        for eachindex in range(repeatLoc2, repeatLoc2 + longestRepeatlength):
            motherGen[eachindex] = motherGen[eachindex - repeatLoc2 + repeatLoc1]           
            
    elif typeOfGen == 'i': 
        
        lrepeat = int(detailArr[0])
        linter = int(detailArr[1])
        repeatLoc1 = random.randint(0  , G- 2*lrepeat -2*linter  )
        intrepeatLoc1 = random.randint(repeatLoc1 + lrepeat ,G- lrepeat -2*linter )
        repeatLoc2 = random.randint(intrepeatLoc1 + linter , G- lrepeat - linter )
        intrepeatLoc2 = random.randint(repeatLoc2 + lrepeat, G-linter)
        
        for eachindex in range(repeatLoc2, repeatLoc2+ lrepeat):
            motherGen[eachindex] = motherGen[eachindex - repeatLoc2 + repeatLoc1]
        
        for eachindex in range(intrepeatLoc2, intrepeatLoc2+ linter):
            motherGen[eachindex] = motherGen[eachindex -intrepeatLoc2 + intrepeatLoc1 ]
        
        print "repeatLoc1,repeatLoc2,intrepeatLoc1,intrepeatLoc2, lrepeat,linter",repeatLoc1,repeatLoc2,intrepeatLoc1,intrepeatLoc2, lrepeat,linter
    
    elif typeOfGen == 't' :
        ltriple = int(detailArr[0])
        
        repeatLoc1 = random.randint(0, G-3*ltriple)
        repeatLoc2 = random.randint(repeatLoc1+ ltriple, G- 2*ltriple)
        repeatLoc3 = random.randint(repeatLoc2+ ltriple, G- ltriple)
        
        for eachindex in range(repeatLoc2, repeatLoc2 + ltriple):
            motherGen[eachindex] = motherGen[eachindex -repeatLoc2 + repeatLoc1 ]
            
        for eachindex in range(repeatLoc3, repeatLoc3 + ltriple):
            motherGen[eachindex] = motherGen[eachindex - repeatLoc3 + repeatLoc1]
            
            
        print "repeatLoc1, repeatLoc2, repeatLoc3, ltriple", repeatLoc1, repeatLoc2, repeatLoc3, ltriple

    elif typeOfGen == 'd' :
        filename = detailArr[0] 
        truncatestart = int(detailArr[1])
        truncateend = int(detailArr[2])
        
        
        f = open(filename , 'r')
        
        line = f.readline()
        
        motherStr = ""
        print "Genome Detail ", line
        
        while(len(line) > 0):   
            line= f.readline()
            motherStr = motherStr + line[0:-1]
        
        print "genomeLen", len(motherStr)
        
        G = truncateend - truncatestart 
        print "G:", G
        motherGen = np.zeros(G, dtype = np.int8)
        
        for eachindex in range(truncatestart, truncateend):
            if motherStr[eachindex] == 'A':
                motherGen[eachindex- truncatestart] = 1
            elif motherStr[eachindex] == 'C':
                motherGen[eachindex- truncatestart] = 2
            elif motherStr[eachindex] == 'G':
                motherGen[eachindex- truncatestart] = 3
            elif motherStr[eachindex] == 'T':
                motherGen[eachindex- truncatestart] = 4

        f.close()
    
    
    elif typeOfGen == 'm':
        lrepeat = int(detailArr[0])
        lsnp = int(detailArr[1])
        lint = int(detailArr[2])

        randLoc1 = random.randint(0, G - 2*lrepeat - 2*lsnp - 2*lint)
        randLocSnp1 = random.randint(randLoc1+ lrepeat, G - lrepeat - 2*lsnp - 2*lint)
        randLocint1 = random.randint(randLocSnp1 + lsnp , G - lrepeat - lsnp - 2*lint)
        
        randLoc2 = random.randint(randLocint1 + lint, G - lrepeat - lsnp - lint)
        randLocSnp2 = random.randint(randLoc2 + lrepeat, G - lsnp - lint)
        randLocint2 = random.randint(randLocSnp2 + lsnp, G - lint)
        
        
        for eachindex in range(randLoc2, randLoc2+ lrepeat):
            motherGen[eachindex] = motherGen[eachindex - randLoc2 + randLoc1]
            
        for eachindex in range(randLocSnp2, randLocSnp2+ lsnp):
            motherGen[eachindex] = motherGen[eachindex - randLocSnp2 + randLocSnp1]
            
        for eachindex in range(randLocint2, randLocint2 + lint):
            motherGen[eachindex] = motherGen[eachindex - randLocint2 + randLocint1]
            #print eachindex, eachindex - randLocint2 + randLocint1

        #print motherGen[randLocint1:randLocint1+lint] == motherGen[randLocint2:randLocint2+lint]
        # introduce 1 SNP on both sides
        
        motherGen[randLocSnp2 + lsnp/4] = motherGen[randLocSnp2 + lsnp/4] + 1
        if motherGen[randLocSnp2 + lsnp/4] == 5 :
            motherGen[randLocSnp2 + lsnp/4] = 1
            
        motherGen[randLocSnp2 + lsnp*3/4] = motherGen[randLocSnp2 + lsnp*3/4] + 1
        if motherGen[randLocSnp2 + lsnp*3/4] == 5:
            motherGen[randLocSnp2 + lsnp*3/4] = 1
        
        print "randLocSnp1 + lsnp/4, randLocSnp2 + lsnp/4",randLocSnp1 + lsnp/4, randLocSnp2 + lsnp/4
        print motherGen[randLocSnp1 + lsnp/4], motherGen[randLocSnp2 + lsnp/4]

        print "randLocSnp1 + 3*lsnp/4, randLocSnp2 + 3*lsnp/4",randLocSnp1 + 3*lsnp/4, randLocSnp2 + 3*lsnp/4
        print motherGen[randLocSnp1 + 3*lsnp/4], motherGen[randLocSnp2 + 3*lsnp/4]
        
 #       motherGen[randLocSnp2 + lsnp/2] = motherGen[randLocSnp2 + lsnp/2] + 1
#        if motherGen[randLocSnp2 + lsnp/2] == 5 :
#            motherGen[randLocSnp2 + lsnp/2] = 1

#        print motherGen[randLocSnp1 + lsnp/2], motherGen[randLocSnp2 + lsnp/2]


    elif typeOfGen == 'a':
        print "Tandem Repeat"
        ltandem1 = int(detailArr[0])
        lcopyNum1 = int(detailArr[1])
        ltandem2 = int(detailArr[2])
        lcopyNum2 = int(detailArr[3])
        
        lrepeat1 = ltandem1*lcopyNum1
        lrepeat2 = ltandem2*lcopyNum2
        
        randLoc1 = random.randint(0, G - lrepeat1 - lrepeat2)
        randLoc2 = random.randint(randLoc1 + lrepeat1, G  - lrepeat2)
        
        for copyindex in range(lcopyNum1):
            for eachindex in range(randLoc1+copyindex*ltandem1, randLoc1 + (copyindex+1) *ltandem1 ):
                motherGen[eachindex]= motherGen[eachindex- copyindex*ltandem1]
        
        for copyindex in range(lcopyNum2):
            for eachindex in range(randLoc2+copyindex*ltandem2, randLoc2 + (copyindex+1) *ltandem2 ):
                motherGen[eachindex]= motherGen[eachindex- copyindex*ltandem2]
                
        print "randLoc1, randLoc2", randLoc1, randLoc2
        
    else:
        
        print "Error in Type "
        
    indel = parameterRobot.indel
    reads, noisyReads = readGen(N,L,p, motherGen, indel)
    
    logging.rawDataSave(parameterRobot.defaultFolder+"UnitTest", motherGen, reads, noisyReads)
    
    return motherGen, reads, noisyReads
    
def readGen(N,L, p,motherGen, indel):

    reads = np.zeros(N*L, dtype = np.int8).reshape(N,L)
    noisyReads = np.zeros(N*L, dtype = np.int8).reshape(N,L)
    
    G = len(motherGen)
    
    cicularGen = np.zeros(G + L, dtype = np.int8)
    cicularGen[0:G] = motherGen
    cicularGen[G:G+L] = motherGen[0:L]
    
    for indexN in range(N):
        randomLoc = random.randint(0, G-1)

        reads[indexN][:] = cicularGen[randomLoc:randomLoc+L]
 
    if indel == False:
        noisyReads = noisyReading(reads, p)
    else:
        noisyReads = addIndelNoise(reads, p)
        
    return reads, noisyReads
    
def noisyReading(reads, p):
    
    N = len(reads)
    L = len(reads[0])
    
    noisyReads = np.zeros(N*L, dtype = np.int8).reshape(N,L)
    
    for indexN in range(N):
        
        for indexL in range(L):
            randombit = random.random()
            corruptOrNot = 4
            
            if p < randombit :
                corruptOrNot = 0            
            elif 2*p/3 < randombit <= p :
                corruptOrNot = 1 
            elif p/3<  randombit <= 2*p/3:
                corruptOrNot = 2
            elif randombit <= p/3:
                corruptOrNot = 3
   
            if corruptOrNot >3 :
                print "Base transform error"
            
            noisyReads[indexN][indexL] = np.mod( reads[indexN][indexL] + corruptOrNot , 4) 
            
            if noisyReads[indexN][indexL]  == 0:
                noisyReads[indexN][indexL]  = 4
            
    return noisyReads



            
        
def addIndelNoise(readList ,p):
    noisyReadList = []
    #L = 1.3*len(readList[0])
    
    for eachread in readList:
        tempread = []
        for eachbase in eachread:
           
            randDel = random.random()
            randIns = random.random()
            
            if randDel > p :
                tempread.append(eachbase)
            
            if randIns< p/4 : 
                corruptOrNot = 1
            elif randIns < 2*p/4 and randIns >= p/4:
                corruptOrNot = 2
            elif randIns < 3*p/4 and randIns >= 2*p/4:
                corruptOrNot = 3
            elif  randIns < 4*p/4 and randIns >= 3*p/4:
                corruptOrNot = 4
            else:
                corruptOrNot = 0 
     
            if corruptOrNot != 0 :
                tempread.append(corruptOrNot)
        
        #while (len(tempread) < L):
        #    tempread.append(0)
            
        noisyReadList.append(tempread)
        
    return noisyReadList
             
    
    
    
    
    