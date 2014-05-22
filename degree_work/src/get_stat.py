import math
from decimal import Decimal

import numpy as np

# paths to input files
assembly_path = "/home/sergey/GAGE/ecoli/Ecoli_first350k1.fna"
sam_file_path = "Ecoli_first350k.map.sam"
# paths to output files
gc_stat_file = "depth-stat"
depth_stat_file = "gc-stat"

# dictionaries/maps for gc-stat and depth-stat
full_depth = {}
full_gc = {}
# dictionary, that matches contig name and line in the file
# ctg_line = {}

def get_data():
    # initialize all the dictionaries
    ass = open(assembly_path, "r")
    lines = ass.readlines()
    for i in range(len(lines)):
        line = lines[i].split(' ')
        contig_name = line[0]
        if (contig_name[0] == '>'):
            name = contig_name[1:]
            # ctg_line[name] = i + 1
            ll = len(lines[i + 1])
            full_depth[name] = np.zeros(ll)
            full_gc[name] = np.zeros(ll)
            full_gc[name] = get_gc(lines[i+1])
        else:
            continue
    ass.close()
    #count depth stats
    with open(sam_file_path, 'r') as aln:
        for l in aln:
            line = l.split('\t')
            if (line[0][0] == '@'):
                ind = 0
            else:
                # if (ind % 2 == 0):
                #     st = line[0] + '/1'
                # else:
                #     st = line[0] + '/2'
                ctg = line[2].lstrip()
                if (ctg != '*'):
                    readlen = len(line[9])
                    pos = int(line[3].lstrip())
                    l = min(pos + readlen, len(full_depth[ctg]))
                    for j in range(pos, l):
                        full_depth[ctg][j] += 1
                ind += 1


def get_gc(contig):
    gc = calculateContigGCcont(contig, 101)
    return gc

def isGC(seq):
    return (seq == 'G' or seq == 'C' or seq == 'g' or seq == 'c')

def getGCtotal(seq1, seq1len):
    GCtot = 0
    for i in range(0, seq1len):
        if (isGC(seq1[i])):
            GCtot += 1
    return (GCtot)

def calculateContigGCcont(contig, windowSize):
    j = 0
    baseGC = 0
    ctglen = len(contig)
    GCpast = np.zeros(windowSize)
    GCcont = np.zeros(ctglen)

    if (ctglen < 2 * windowSize):
        # contig is too small to estimate per-base windowing
        baseGC = getGCtotal(contig, ctglen)
        baseGC = math.floor(100.0 * (float)(baseGC / ctglen))

        for j in range(0, ctglen):
            GCcont[j] = baseGC
    else:
        baseGC = getGCtotal(contig, windowSize)
        GCpast[0] = baseGC

        for j in range(0, windowSize - 1):
            GCcont[j] = math.floor(100.0 * baseGC / (float)((j + 1) * windowSize))
            if (GCcont[j] > 100):
                GCcont[j] = 100
            #   printf("GC correction out of range %d %s %d len %d '%c'\n", baseGC, contig->name, j, contig->seqLen, contig->seq[j]);
            GCpast[(j + 1) % windowSize] = GCpast[j % windowSize]

            if (isGC(contig)):
                GCpast[(j + 1) % windowSize] -= 1

            if (isGC(contig[j + windowSize])):
                GCpast[(j + 1) % windowSize] += 1;

            baseGC += GCpast[(j + 1) % windowSize]

        for j in range(windowSize - 1, ctglen - windowSize):
            GCcont[j] = math.floor(100.0 * baseGC / (float)(windowSize * windowSize))

            if (GCcont[j] > 100):
                GCcont[j] = 100
            #   printf("GC correction out of range %d %s %d len %d '%c'\n", baseGC, contig->name, j, contig->seqLen, contig->seq[j]);

            baseGC -= GCpast[(j + 1) % windowSize]
            GCpast[(j + 1) % windowSize] = GCpast[j % windowSize]

            if (isGC(contig[j])):
                GCpast[(j + 1) % windowSize] -= 1

            if (isGC(contig[j + windowSize])):
                GCpast[(j + 1) % windowSize] += 1

            baseGC += GCpast[(j + 1) % windowSize]

        for j in range(ctglen - windowSize, ctglen):
            GCcont[j] = math.floor(100.0 * baseGC / (float)((ctglen - j) * windowSize));
            if (GCcont[j] > 100):
                GCcont[j] = 100
            #   printf("GC correction out of range %d %s %d len %d '%c'\n", baseGC, contig->name, j, contig->seqLen, contig->seq[j]);
            baseGC -= GCpast[(j + 1) % windowSize]
    return GCcont

def main():
    get_data()

    gcstat = open(gc_stat_file, 'w')
    depthstat = open(depth_stat_file, 'w')

    for k, v in sorted(full_depth.items()):
        depthstat.write(">" + k + '\n')
        gcstat.write(">" + k + '\n')

        depth_cur = v
        gc_cur = full_gc[k]

        for i, x in enumerate(v):
            gcstat.write(str(int(gc_cur[i])) + ' ')
            depthstat.write(str(int(depth_cur[i])) + ' ')

        depthstat.write('\n')
        gcstat.write('\n')

    gcstat.close()
    depthstat.close()

if __name__ == '__main__':
    main()