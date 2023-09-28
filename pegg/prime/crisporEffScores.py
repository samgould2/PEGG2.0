# main functions: calcAllScores and calcMutSeqs 

# this library re-implements the efficiency scoring functions of these articles in calcAllScores():

# - WangSvm: Wang et al, Science 2014, PMID 24336569, no website
# - Doench: Doench et al, Nat Biotech 2014, PMID 25184501, http://www.broadinstitute.org/rnai/public/analysis-tools/sgrna-design
# - CrisprScan: Moreno-Mateos, Nat Meth 2015, PMID 26322839, http://crisprscan.org
# - ssc: Xu et al, Gen Res 2015, PMID 26063738, http://crispr.dfci.harvard.edu/SSC/
# - Chari: Chari et al, PMID 26167643 http://crispr.med.harvard.edu/sgRNAScorer
# - Fusi: Fusi et al, prepublication manuscript on bioarxiv, http://dx.doi.org/10.1101/021568 http://research.microsoft.com/en-us/projects/azimuth/, only available as a web API
# - Housden: Housden et al, PMID 26350902, http://www.flyrnai.org/evaluateCrispr/
# - Wu-Crispr: Wong et al, http://www.genomebiology.com/2015/16/1/218
# - DeepCpf1, Kim et al, PMID 29431740, https://www.ncbi.nlm.nih.gov/pubmed/29431740
# - SaCas9 efficiency score (no name), Najm et al, https://www.ncbi.nlm.nih.gov/pubmed/29251726

# Also includes the prediction of DSB-repair outcome in calcMutSeqs:
# - OOF: Microhomology and out-of-frame score from Bae et al, Nat Biotech 2014 , PMID24972169 http://www.rgenome.net/mich-calculator/
# - Wei Chen et al: 

# the input are 100bp sequences that flank the basepair just 5' of the PAM +/-50bp.
# so 50bp 5' of the PAM, and 47bp 3' of the PAM -> 100bp

# this module uses pipes to feed data into some programs
# If you run too many sequences at once, it may hang. Increase the BUFSIZE variable in this case.

from subprocess import Popen, PIPE, STDOUT, check_output, CalledProcessError, call
import platform, math, tempfile, bisect, sys, os, logging, types, optparse, shutil
from os.path import dirname, join, basename, isfile, expanduser, isdir, abspath
from math import log10

import urllib.request, urllib.error, urllib.parse, pickle
import json

myDir = dirname(__file__)

aziDir = join(myDir, "bin/Azimuth-2.0/")
sys.path.append(aziDir)

global binDir
binDir = None

# the name of a directory to use for caching some efficiency values that are slow to calculate
# deactivated by default
cacheDir = None

# by default bindir is relative to the location of this library
if binDir is None:
    binDir = join(dirname(__file__), "bin")

BUFSIZE = 10000000

def setBinDir(path):
    global binDir
    binDir = path

def setCacheDir(path):
    global cacheDir
    cacheDir = path

def getBinPath(name, isDir=False):
    """
    get the full pathname of a platform-specific binary, in the bin/ directory relative to this directory
    """
    currPlatform = platform.system()
    binPath = join(binDir, currPlatform, name)
    if isDir and not isdir(binPath):
        raise Exception("Could not find directory %s" % binPath)
    if not isDir and not isfile(binPath):
        raise Exception("Could not find file %s" % binPath)
    return binPath

def seqToVec(seq, offsets={"A":0,"C":1,"G":2,"T":3}):
    """ convert a x bp sequence to a 4 * x 0/1 vector
    >>> seqToVec("AAAAATTTTTGGGGGCCCCC")
    [1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0]
    """
    assert(len(seq)==20)
    row = [0]*len(seq)*4
    pseudoOffset = offsets["A"]
    for pos, nucl in enumerate(seq):
        nucl = nucl.upper()
        # treat N, Y, etc like "A". Happens very rarely.
        nuclOffset = offsets.get(nucl, pseudoOffset)
        vecPos = (pos*len(offsets))+nuclOffset
        #if vecPos not in range(len(row)):
            #ofh = open("temp.txt", "a")
            #ofh.write(str(vecPos)+" "+seq+" "+str(row)+"pos %d, nucl %s" % (pos, nucl)+"\n")
            #assert(False)
        row[vecPos] = 1
    return row

def vecToSeqDicts(coefs):
    " convert a list of 80 floats to 20 dictionaries with A/C/T/G -> float "
    freqs = []
    for i in range(0,20):
        charFreqs = {}
        for nucl, x in zip("ACGT", list(range(0,4))):
            freq = coefs[i*4+x]
            if freq==0.0:
                continue
            charFreqs[nucl] = freq
        freqs.append(charFreqs)
    return freqs


def writeDict(d, fname):
    " write dict as a tab file "
    if not isdir(dirname(fname)):
        logging.debug("Cannot write %s, no caching of efficiency scores" % fname)
        return

    ofh = open(fname, "w")
    for k, v in d.items():
        if type(v)==tuple:
            ofh.write("%s\t%s\n" % (k, "\t".join([str(x) for x in v])))
        else:
            ofh.write("%s\t%s\n" % (k, str(v)))
    ofh.close()

def readDict(fname, isFloat=True):
    " read dict from a tab sep file "
    if not isfile(fname):
        logging.debug("%s does not exist. Returning empty dict" % fname)
        return {}

    logging.info("Reading %s" %fname)
    data = {}
    for line in open(fname):
        fs = line.rstrip("\n").split("\t")
        if len(fs)==2:
            k, v = fs
            if isFloat:
                v = float(v)
        else:
            k = fs[0]
            v = tuple(fs[1:])
            if isFloat:
                v = tuple([float(x) for x in v])
        data[k] = v
    return data

class ScoreCache:
    """
    a cache of eff scores, kept on disk. Can avoid slow calculations by keeping
    the value of the score in a tab-sep file.
    """

    def __init__(self, scoreName):
        self.cacheFname = join(cacheDir, "%s.tab" % scoreName)
        scoreCache = readDict(self.cacheFname, isFloat=True)
        self.scoreCache = scoreCache

    def findNewSeqs(self, seqs):
        """ get seqs that are not in cache. If all are, return the list of scores.
        Otherwise return None for the scores.
        Returns tuple (seqs, scores)
        """
        self.allSeqs = seqs
        newSeqs = set()
        for s in seqs:
            if not s in self.scoreCache:
                newSeqs.add(s)

        scoreList = None
        if len(newSeqs)==0:
            scoreList = [self.scoreCache[s] for s in seqs]
        self.newSeqs = newSeqs
        return newSeqs, scoreList

    def mergeIntoCache(self, newScores):
        # create final result merging cache and newly obtained scores
        scoreList = []
        assert(len(newScores)==len(self.newSeqs))
        newScoreDict = dict(list(zip(self.newSeqs, newScores)))

        for s in self.allSeqs:
            if s in newScoreDict:
                scoreList.append(newScoreDict[s])
            else:
                scoreList.append(self.scoreCache[s])

        self.scoreCache.update(newScoreDict)
        writeDict(self.scoreCache, self.cacheFname)
        return scoreList

def calcAziScore(seqs):
    " the official implementation of the Doench2016 (aka Fusi) score from Microsoft "
    import numpy
    import azimuth.model_comparison
    res = []
    for seq in seqs:
        if "N" in seq:
            res.append(-1) # can't do Ns
            continue

        pam = seq[25:27]
        # pam_audit = do not check for NGG PAM
        seq = seq.upper()
        score = azimuth.model_comparison.predict(numpy.array([seq]), None, None, pam_audit=False)[0]
        res.append(int(round(100*score)))
    return res

# ----------- MAIN --------------
if __name__=="__main__":
    import cProfile
    #cProfile.runctx('print runLindel(["test"], ["CCCTGGCGGCCTAAGGACTCGGCGCGCCGGAAGTGGCCAGGGCGGGGGCGACCTCGGCTCACAG"])', globals(), locals())
    #pr = cProfile.Profile()
    #pr.enable()
    #runLindel(["test"], ["CCCTGGCGGCCTAAGGACTCGGCGCGCCGGAAGTGGCCAGGGCGGGGGCGACCTCGGCTCACAG"])
    #pr.disable()
    #pr.print_stats(sort='tottime')
    #sys.exit(0)


    args, options = parseArgs()
    if options.test:
        test()
        sys.exit(0)

    #setBinDir("../crispor/bin")
    setBinDir("./bin")
    inFname = sys.argv[1]
    seqs = readSeqs(inFname)
    if len(seqs)==0:
        logging.error("No sequences in input left")
    else:
        printScoreTabSep(seqs, options.all, options.enzyme)
