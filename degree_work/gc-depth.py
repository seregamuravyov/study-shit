import numpy as np
import math
import pylab as p
from decimal import *
from matplotlib.pyplot import plot, draw, show
from decimal import Decimal
from os.path import exists
import scipy.stats as st

assembly = {}
reads = {}
alig = []
depth = {}
gc = {}

depthAvgSum = 0.0
depthAvgNorm = 0.0

contigs = []
GCcont = np.zeros(101)
depthNormalizer = np.zeros(102) 
negBinomParam_r = np.zeros(102)
GC = 0

class Alignment:
	def __init__(self, r, ctg, p):
		self.read = r
		self.contig = ctg
		self.position = p

def get_data():
	ass = open('genome.ctg.fasta')
	r1 = open('frag_1.fastq')
	r2 = open('frag_2.fastq')
	aln = open('alig.sam')

	#parse assembly
	lines = ass.readlines()
	for i in range(len(lines)):
		line = lines[i].split(' ')
		contig_name = line[0]
		if (contig_name[0] == '>'):
			assembly[contig_name[1:]] = lines[i + 1]
			contigs.append(contig_name[1:])
		else:
			continue

	#parse mate 1
	lines = r1.readlines()
	parse_reads(lines)

	#parse mate 2
	lines = r2.readlines()
	parse_reads(lines)

	#parse alignment
	lines = aln.readlines()
	ind = 0
	for i in range(len(lines)):
		line = lines[i].split('\t');
		if (line[0][0] == '@'):
			ind = 0
		else:
			if (ind % 2 == 0):
				st = line[0]+'/1'
			else:
				st = line[0]+'/2'
			x = Alignment(st, line[2].lstrip(), int(line[3].lstrip()))
			alig.append(x)
			ind += 1

def get_depth():
	for k, v in assembly.items():
		depth[k] = np.zeros(len(assembly[k]))

	for i in range(len(alig)):
		rkey = alig[i].read
		ckey = alig[i].contig
		pos = alig[i].position
		rlen = reads[rkey]

		if (ckey != '*'):
			l = min(pos+rlen, len(depth[ckey]))
			for j in range(pos, l):
				depth[ckey][j] += 1

def get_gc():
	for k, v in assembly.items():
		gc[k] = calculateContigGCcont(assembly[k], 101)
		#print(GCcont)

def parse_reads(lines):
	for i in range(len(lines)):
		line = lines[i]
		if (i % 4 == 0):
			reads[line[1:-1]] = len(lines[i+1]) - 1
		else:
			continue

def isGC(seq):
    return (seq == 'G'or seq == 'C' or seq == 'g' or seq == 'c')

def getGCtotal(seq1, seq1len):
    GCtot = 0
    for i in range(0, seq1len):
        if(isGC(seq1[i])):
            GCtot+=1
    return (GCtot)

def calculateContigGCcont(contig, windowSize):
	j = 0
	baseGC = 0
	ctglen = len(contig)
	#print(contig)
	GCpast = np.zeros(windowSize)
	GCcont = np.zeros(ctglen)
	#print(GCcont)
	#print(ctglen)

	if (ctglen < 2 * windowSize):
		# contig is too small to estimate per-base windowing
		baseGC = getGCtotal(contig, ctglen)
		baseGC = math.floor(100.0 * (float)(baseGC / ctglen))
		
		for j in range(0, ctglen):
			GCcont[j] = baseGC;
	else:
		baseGC = getGCtotal(contig, windowSize);
		GCpast[0] = baseGC;

		for j in range(0, windowSize - 1):
			GCcont[j] = math.floor(100.0 * baseGC / (float)((j+1) * windowSize));
			if (GCcont[j] > 100):
				GCcont[j] = 100
			#	printf("GC correction out of range %d %s %d len %d '%c'\n", baseGC, contig->name, j, contig->seqLen, contig->seq[j]);
			GCpast[(j+1) % windowSize] = GCpast[j % windowSize];

			if(isGC(contig)):
				GCpast[(j+1) % windowSize]-=1;
			
			if(isGC(contig[j+windowSize])):
				GCpast[(j+1) % windowSize]+=1;
			
			baseGC += GCpast[(j+1) % windowSize]

		for j in range(windowSize - 1, ctglen - windowSize):
			GCcont[j] = math.floor(100.0*baseGC / (float)(windowSize * windowSize));

			if (GCcont[j] > 100):
				GCcont[j] = 100
			#	printf("GC correction out of range %d %s %d len %d '%c'\n", baseGC, contig->name, j, contig->seqLen, contig->seq[j]);
			
			baseGC -= GCpast[(j+1) % windowSize];
			GCpast[(j+1) % windowSize] = GCpast[j % windowSize];

			if(isGC(contig[j])):
				GCpast[(j+1) % windowSize]-=1;
			
			if(isGC(contig[j+windowSize])):
				GCpast[(j+1) % windowSize]+=1;

			baseGC += GCpast[(j+1) % windowSize];

		for j in range(ctglen - windowSize, ctglen):
			GCcont[j] = math.floor(100.0*baseGC / (float)((ctglen - j) * windowSize));
			if (GCcont[j] > 100):
				GCcont[j] = 100
			# 	printf("GC correction out of range %d %s %d len %d '%c'\n", baseGC, contig->name, j, contig->seqLen, contig->seq[j]);
			baseGC -= GCpast[(j+1) % windowSize];
	return GCcont

###not used
def plotDistribution(gc_const, maxdp, gc, dp):
	xasix = np.arange(maxdp)
	depth_pos = np.zeros(maxdp)
	for i in range(len(gc)):
		if (gc[i] == gc_const):
			depth_pos[dp[i]]+=1
	plt.plot(xasix, depth_pos, 'ro')
	plt.show()
##########################################


def negBinomPMF(k, r, p):
	#ans = Decimal(0.0)
	ans = 0.0
	for i in range(1, k+1):
		#ans *= Decimal((r - 1 + i)/i)
		ans += math.log(r - 1 + i) - math.log(i)
	#ans *= Decimal((1.0 - p) ** r)
	#ans *= Decimal(p ** k)
	ans += r * math.log(1.0 - p)
	ans += (float)(k) * math.log(p)
	return ans

def computeNormaliziedDepthGCParameters(depthNormalizer, depthNormalizerCount, negBinomParam_r, avgDepth):
	#int j;

	for j in range(0, 101):
		if(depthNormalizerCount[j] > 0):
			depthNormalizer[j] = depthNormalizer[j] / depthNormalizerCount[j]
		else:
			depthNormalizer[j] = avgDepth

		if(depthNormalizer[j] < 10.0):
			depthNormalizer[j] = 10.0

		negBinomParam_r[j] = depthNormalizer[j]


def main():
	get_data()
	get_depth()
	get_gc()

	out = open('depth_per_contig.out', 'w')
	out_gc = open('gc_per_contig.out', 'w')
	out_pl = open('plot.dat', 'w')
	out_pl_nb = open('plot_nb.dat', 'w')

	pl1 = open('pl1.txt', 'w')
	pl2 = open('pl2.txt', 'w')
	pl11 = open('pl11.txt', 'w')
	

	#out1 = open('pairs_of_depth.out', 'w')
	#out2 = open('probabilities.out', 'w')
	#wg = open('res.my.depth.wig', 'w')
	

	maxdepth = 0
	GC = 0
	gcval = 0
	kgc = ''
	gccount = np.zeros(102)
	depthNormalizer = np.zeros(102)
	depthAvgSum = 0
	depthAvgNorm = 0

	for k, v in assembly.items():
		gccount = np.zeros(101)

		d = depth[k]
		g = gc[k]
		out.write('\n>>>' + k + ' ' + str(int(np.amax(d))) + ' \n')
		out_gc.write('\n>>>' + k + ' ' + str(int(len(g))) +' \n')
		for i in range(len(d)):
			out.write(str(int(d[i])) + ' ')
			if (d[i] > maxdepth):
				maxdepth = d[i]
			depthAvgSum += d[i]
			depthAvgNorm += 1.0

		for i in range(len(g)):
			out_gc.write(str(int(g[i])) + ' ')
			gccount[g[i]] += 1
			depthNormalizer[g[i]] += d[i]

			if (gccount[g[i]] > GC):
				GC = gccount[g[i]]
				gcval = g[i]
				kgc = k
		
	getcontext().prec = 7
	##Эмпирическое распределение##
	##out_pl.write(str(GC) + "\n" + str(maxdepth) + "\n")
	out_pl_nb.write("# This file is called plot.dat\n")
	out_pl_nb.write("# Force-Deflection data for a beam and a bar\n")
	out_pl_nb.write("# Depth  Position\n")

	xasix = np.arange(maxdepth + 1)
	depth_pos = np.zeros(maxdepth + 1)
	f_obs_n = np.zeros(maxdepth + 1)
	f_obs = np.zeros(maxdepth + 1)
	f_exp = np.zeros(maxdepth + 1)
	d = depth[kgc]
	g = gc[kgc]
	dpsum = 0
	out_pl.write(kgc + "\n")
	for i in range(len(g)):
	 	if (g[i] == gcval):
	 		depth_pos[d[i]]+=1

	for i in range(len(xasix)):
		dpsum += depth_pos[i]

	for i in range(len(xasix)):
		out_pl.write(str(int(xasix[i])) + " " +  str(depth_pos[i]) + '\n')
		pl11.write(str(depth_pos[i]) + '\n')
		pl1.write(str(depth_pos[i] / dpsum) + '\n')

	f_obs_n = [depth_pos[i] / dpsum for i in range(len(xasix))]
	f_obs = depth_pos
	##Neg_Binom из статьи##
	#out_pl_nb.write(str(GC) + "\n" + str(maxdepth) + "\n")
	#out_pl_nb.write(kgc + '\n')

	out_pl.write("# This file is called plot_nb.dat\n")
	out_pl.write("# Force-Deflection data for a beam and a bar\n")
	out_pl.write("# Depth  Position\n")

	depth_pos = np.zeros(maxdepth + 1)
	computeNormaliziedDepthGCParameters(depthNormalizer, gccount, negBinomParam_r, depthAvgSum/depthAvgNorm)

	for i in range(len(g)):
		if (g[i] == gcval):
			depth_pos[d[i]]=negBinomPMF((int)(d[i]), negBinomParam_r[g[i]], 0.5)


	for i in range(len(xasix)):
		out_pl_nb.write(str(int(xasix[i])) + " " +  str(Decimal(depth_pos[i])) + '\n')
		pl2.write((str(math.exp(depth_pos[i])) if depth_pos[i] != 0 else '0.0') + '\n')


	f_exp = [math.exp(depth_pos[i]) if depth_pos[i] != 0 else '0.0' for i in range(len(xasix))]

	obs = np.array([f_obs_n[8:49], f_exp[8:49]])
	#print(obs)

	print(st.chisquare(f_obs_n[8:49], f_exp[8:49], ddof=[0, 0.00000000001, 0.1, 0.02]))
	#print(st.chisquare(f_obs[8:49], f_exp[8:49], ddof=[0, 0.00000000001, 0.1, 0.02]))
	print(st.chi2_contingency(obs))

	#for i in range(len(depth_pos)):
	#	out_pl.write(str(int(depth_pos[i])) + "\n")
	#p.plot(xasix, depth_pos, 'ro')
	#draw()
	#p.show()
	#plotDistribution(GC, maxdepth, g, d)




if __name__ == '__main__':
    main()