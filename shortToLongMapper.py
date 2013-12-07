import loggingIndel

def findPerfectMap(readsMapping, i, parameterRobot, longRepeatStat):
    
    # Output Format : "indexlong", "indexshort",  "jstart", "jend", "istart", "iend"
    tempMappingi = []
    Nshort , Nlong, Lshort, Llong = parameterRobot.Nshort , parameterRobot.Nlong, parameterRobot.Lshort, parameterRobot.Llong
    
    shortreadMap = readsMapping[0]
    longreadMap = readsMapping[1]
    
    [ startLoc1, startLoc2, lsnp ]= longRepeatStat[1]
    
    
    overlapthreshold = parameterRobot.liid 
    
    longnum, longstart, longend, longtype = longreadMap[i][0], longreadMap[i][1] , longreadMap[i][1] + Llong, longreadMap[i][2]
    #print "longnum, longstart, longend, longtype", longnum, longstart, longend, longtype
    copyNum = 0
    if min( abs(longstart-startLoc1), abs(longstart-startLoc1- lsnp) ) < min(abs(longstart-startLoc2), abs(longstart-startLoc2- lsnp) ):
        copyNum = 1 
    else:
        copyNum = 2
    
    
    c1, c2 = 0 ,0 
    
    for j in range(len(shortreadMap)):
        shortnum, shortstart, shortend, shorttype = shortreadMap[j][0], shortreadMap[j][1] , shortreadMap[j][1] + Lshort, shortreadMap[j][2]
        
        
        if copyNum == 2 :
            correction = startLoc2 - startLoc1 
            
        elif copyNum == 1:
            correction = startLoc1 - startLoc2
        

        if (shortstart < longend - overlapthreshold and longstart < shortend - overlapthreshold ) or (shortstart + correction< longend - overlapthreshold and longstart  < shortend+ correction - overlapthreshold ):

            if   (shortstart < longend - overlapthreshold and longstart < shortend - overlapthreshold ) :
                
                if shortstart < longstart:
                    jstart, jend, istart, iend = 0,Lshort - longstart + shortstart,longstart - shortstart ,Lshort-1
                elif shortend > longend:
                    jstart, jend, istart, iend = Llong - Lshort + shortend - longend,Llong,0,Lshort - shortend + longend
                else:
                    jstart, jend, istart, iend = shortstart - longstart,shortstart - longstart + Lshort,0,Lshort

            elif (shortstart + correction< longend - overlapthreshold and longstart  < shortend+ correction - overlapthreshold ):
                
                shortstart, shortend = shortstart + correction,shortend+ correction 
                
                if shortstart < longstart:
                    jstart, jend, istart, iend = 0,Lshort - longstart + shortstart,longstart - shortstart ,Lshort-1
                elif shortend > longend:
                    jstart, jend, istart, iend = Llong - Lshort + shortend - longend,Llong,0,Lshort - shortend + longend
                else:
                    jstart, jend, istart, iend = shortstart - longstart,shortstart - longstart + Lshort,0,Lshort
                
                shortstart, shortend = shortstart - correction,shortend- correction 

            if longtype == shorttype : 
                tempMappingi.append([longnum, shortnum,jstart, jend, istart, iend ])
            else:
                tempMappingi.append([longnum, shortnum + Nshort,jstart, jend, istart, iend])
                
                
            if (shortstart < longend - overlapthreshold and longstart < shortend - overlapthreshold ) :
                c1 += 1 
            elif  (shortstart + correction< longend - overlapthreshold and longstart < shortend + correction - overlapthreshold ):
                c2 += 1

    return tempMappingi

def createIdealMapping(readsMapping, longRepeatStat,listOfXnodesInfo, parameterRobot, debug = "map"):
    
    print "createIdealMapping" 
    # Parameters and setting 
    Nlong, Nshort, Llong, Lshort = parameterRobot.Nlong, parameterRobot.Nshort, parameterRobot.Llong, parameterRobot.Lshort
    liid = parameterRobot.liid
    thresForRandom = parameterRobot.thresForRandom
    overhang = parameterRobot.liid 
    editParameters = [parameterRobot.editsub, parameterRobot.editins, parameterRobot.editdel, parameterRobot.editmatch]
    
    
    toProcessList = []
    for item in listOfXnodesInfo[0]:
        toProcessList += item
    
    tempMapping = [[] for i in range(Nlong)]
    
    # Format :   "indexlong", "indexshort",  "jstart", "jend", "istart", "iend"
    for i in toProcessList:
        tempMapping[i] = findPerfectMap(readsMapping, i, parameterRobot, longRepeatStat)

    if debug == "map":
        shortToLongMap = []
        for i in range(Nlong):
            shortToLongMap += tempMapping[i]
            
        loggingIndel.logshortToLongMap(parameterRobot.defaultFolder,shortToLongMap)
        
    elif debug == "vote": 
        shortToLongMap = loggingIndel.loadShortToLongMap(parameterRobot.defaultFolder)
    
    return shortToLongMap








