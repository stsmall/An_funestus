import sys, os

cDir, outBgFileName = sys.argv[1:]
strictThreshold=0.05
lenientThreshold=0.1

def overlap(x1,y1,x2,y2):
    r = 0
    if (y2 < x2 or y1 < x1):
        raise Exception
    elif (y1 <= y2 and y1 >= x2):
        if (x1 > x2):
            r = (x1,y1);
        else:
            r = (x2,y1)
    elif (x1 <= y2 and x1 >= x2):
        if (y1 < y2):
            r = (x1,y1)
        else:
            r = (x1,y2)
    elif (y2 <= y1 and y2 >= x1):
        if (x2 > x1):
            r = (x2,y2)
        else:
            r = (x1,y2)
    elif (x2 <= y1 and x2 >= x1):
        if (y2 < y1):
            r = (x2,y2)
        else:
            r = (x2,y1)
    return r

def writeBed(wins, bedFileName):
    with open(bedFileName, "w") as bedFile:
        for c, s, e in wins:
            bedFile.write("%s\t%d\t%d\n" %(c, s, e))

runs = []
bgWins = []
for cFileName in os.listdir(cDir):
    with open(cDir + "/" + cFileName) as cFile:
        inRun = False
        for line in cFile:
            c, s, e, numSites, classNum, noMigProb, mig12Prob, mig21Prob = line.strip().split()
            s, e = int(s), int(e)
            bgWins.append((c, s, e))
            noMigProb = float(noMigProb)
            if inRun:
                if s == runE:
                    if noMigProb < lenientThreshold:
                        runE = e
                        if noMigProb < minNoMigProb:
                            minNoMigProb = noMigProb
                    else:
                        inRun = False
                        if minNoMigProb < strictThreshold:
                            runs.append((c, runS, runE))
                else:
                    if minNoMigProb < strictThreshold:
                        runs.append((c, runS, runE))
                    if noMigProb < lenientThreshold:
                        inRun = True
                        runS, runE = s, e
                        minNoMigProb = noMigProb
                    else:
                        inRun = False
            else:
                if noMigProb < lenientThreshold:
                    inRun = True
                    runS, runE = s, e
                    minNoMigProb = noMigProb
        if inRun:
            if minNoMigProb < strictThreshold:
                runs.append((c, runS, runE))

writeBed(bgWins, outBgFileName)
