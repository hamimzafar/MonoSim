from collections import defaultdict
import sys
import random
import bisect
import cPickle
import gzip
import os
from trees import *
from generator_helper import *

def bisect_choiceTUP(items):
    """Returns a function that makes a weighted random choice from a list of tuples."""
    added_weights = []
    last_sum = 0.0
    for item, weight in items:
        weight = float(weight)
        last_sum += weight
        added_weights.append(last_sum)

    def choice(rnd=random.random, bis=bisect.bisect):

        return items[bis(added_weights, rnd() * last_sum)][0]

    return choice

def parseModel(gzipFile, readlen):
    """prepares error models for input to mkErrors."""
    try:
        file = gzip.open(gzipFile, 'rb')
    except IOError:
        sys.exit()
    modReadLen = cPickle.load(file)
    if readlen != 'd' and readlen > modReadLen:
        file.close()
        sys.exit()
    mx = cPickle.load(file)
    insD = cPickle.load(file)
    delD = cPickle.load(file)
    gQualL = cPickle.load(file)
    bQualL = cPickle.load(file)
    iQualL = cPickle.load(file)
    readCount = cPickle.load(file)
    rdLenD = cPickle.load(file)
    file.close()
    return mx, insD, delD, gQualL, bQualL, iQualL, readCount, rdLenD

def mkErrors(read, readLen, mx, gQ, bQ, qual):
    """Adds random errors to read."""
    inds = {'A': 0, 'T': 1, 'G': 2, 'C': 3, 'N': 4, 'a': 0, 't': 1, 'g': 2, 'c': 3, 'n': 4}
    pos = 0
    quals = ''
    pos += 1
    read2 = ""
    for s in read:
        if s in inds:
            read2 += s
    read = read2
    while pos <= readLen and pos < len(read) - 4:
        #print "read",read
        #print "type",type(read)
        prev = read[pos:pos + 4]
        after = read[pos + 4]
        d0 = pos
        d1 = inds[prev[3]]
        d2 = inds[prev[2]]
        d3 = inds[prev[1]]
        d4 = inds[prev[0]]
        d5 = inds[after]
        tot = float(mx[d0][d1][d2][d3][d4][d5][5])
        # print tot
        Mprobs = mx[d0][d1][d2][d3][d4][d5] / tot
        val = random.random()
        a = Mprobs[0]
        t = Mprobs[1] + a
        g = Mprobs[2] + t
        c = Mprobs[3] + g
        n = Mprobs[4] + c
        success = False
        if val > n or tot == 0:
            gPos = pos - 1
            while gPos >= 0:
                try:
                    quals += gQ[gPos]()
                    success = True
                    break
                except:
                    gPos -= 1
            if success == False:
                quals += chr(30 + qual)
        elif val > c:
            read = read[:pos + 3] + 'N' + read[pos + 4:]
            bPos = pos - 1
            while bPos >= 0:
                try:
                    quals += bQ[bPos]()
                    success = True
                    break
                except:
                    bPos - 1
                if success == False:
                    quals += chr(2 + qual)
        elif val > g:
            read = read[:pos + 3] + 'C' + read[pos + 4:]
            bPos = pos - 1
            while bPos >= 0:
                try:
                    quals += bQ[bPos]()
                    success = True
                    break
                except:
                    bPos - 1
                if success == False:
                    quals += chr(2 + qual)
        elif val > t:
            read = read[:pos + 3] + 'G' + read[pos + 4:]
            bPos = pos - 1
            while bPos >= 0:
                try:
                    quals += bQ[bPos]()
                    success = True
                    break
                except:
                    bPos - 1
                if success == False:
                    quals += chr(2 + qual)
        elif val > a:
            read = read[:pos + 3] + 'T' + read[pos + 4:]
            bPos = pos - 1
            while bPos >= 0:
                try:
                    quals += bQ[bPos]()
                    success = True
                    break
                except:
                    bPos - 1
                if success == False:
                    quals += chr(2 + qual)
        else:
            read = read[:pos + 3] + 'A' + read[pos + 4:]
            bPos = pos - 1
            while bPos >= 0:
                try:
                    quals += bQ[bPos]()
                    success = True
                    break
                except:
                    bPos - 1
                if success == False:
                    quals += chr(2 + qual)
        pos += 1
    if quals == "":
        #sys.exit()
        return read, quals
    #print "quals",quals,"quals"
    quals += quals[-1]
    read = read[4:readLen + 4]
    quals = quals[:readLen]
    if len(quals) != len(read):
        sys.exit()
    return read, quals

def Comp(sequence):
    """ complements a sequence, preserving case."""
    d = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'a': 't', 't': 'a', 'c': 'g', 'g': 'c', 'N': 'N', 'n': 'n'}
    cSeq = ''
    for s in sequence:
        if s in d:
            cSeq += d[s]
    return cSeq

def GenReads(fasta = "", length = '', models = '', qual = '', out = ''):
    out1=open(out+'.fastq','a')
    mx1,insD1,delD1,gQualL,bQualL,iQualL,readCount,rdLenD=parseModel(models,length)
    #choose good quality bases
    gQList=[]
    for i in (gQualL):
        gL=[]
        keys=i.keys()
        keys.sort()
        for k in keys:
            gL.append((chr(k+qual),i[k]))
        gQList.append(bisect_choiceTUP(gL))
    #choose bad quality bases
    bQList=[]
    for i in (bQualL):
        bL=[]
        keys=i.keys()
        keys.sort()
        for k in keys:
            bL.append((chr(k+qual),i[k]))
        bQList.append(bisect_choiceTUP(bL))
    end_pos = 100
    cur_pos = 0
    while end_pos < len(fasta):
        #print "end_pos",end_pos
        read = fasta[cur_pos:end_pos]
        read1, quals1 = mkErrors(read, length, mx1, gQList, bQList, qual)
        #print "read1",read1
        #print "quals1",quals1
        head1='@'+'r'+'_from_'+'_#0/1\n'
        cur_pos = end_pos
        end_pos += 100
        out1.write(head1)
        out1.write(read1+'\n')
        out1.write('+\n')
        out1.write(quals1+'\n')
    return out1
