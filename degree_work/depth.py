import numpy as np
import math
from decimal import Decimal
from os.path import exists

assembly = {}
reads = {}
alig = []
depth = {}
contigs = []

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
			assembly[contig_name[1:]] = len(lines[i + 1])
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
		depth[k] = np.zeros(assembly[k])

	for i in range(len(alig)):
		rkey = alig[i].read
		ckey = alig[i].contig
		pos = alig[i].position
		rlen = reads[rkey]

		if (ckey != '*'):
			l = min(pos+rlen, len(depth[ckey]))
			for j in range(pos, l):
				depth[ckey][j] += 1

def parse_reads(lines):
	for i in range(len(lines)):
		line = lines[i]
		if (i % 4 == 0):
			reads[line[1:-1]] = len(lines[i+1]) - 1
		else:
			continue

def main():
	get_data()
	get_depth()

	out = open('depth_per_contig.out', 'w')
	out1 = open('pairs_of_depth.out', 'w')
	out2 = open('probabilities.out', 'w')
	wg = open('res.my.depth.wig', 'w')
	

	maxdepth = 0
	for k, v in assembly.items():
		d = depth[k]
		out.write('\n>>>' + k + ' ' + str(int(np.amax(d))) + ' \n')
		for i in range(len(d)):
			out.write(str(int(d[i])) + ' ')
			if (d[i] > maxdepth):
				maxdepth = d[i]

	n = maxdepth + 1
	pairs = np.zeros((n, n))
	probs = np.zeros((n, n))

	for k, v in assembly.items():
		d = depth[k]
		for j in range(len(d) - 1):
			pairs[d[j], d[j+1]] += 1

	for i in range(int(n)):
		sums = 0
		for j in range(int(n)):
			sums += pairs[i, j]
		for j in range(int(n)):
			#out1.write(str(i) + ' ' + str(j) + ' | ' + str(int(pairs[i, j])) + ' | ' + '\n')
			if sums != 0:
				probs[i, j] = float(pairs[i, j]) / sums
			else:
				probs[i, j] = 0
			if (pairs[i, j] > 0):
				out1.write(str(i) + ' ' + str(j) + ' | ' + str(pairs[i, j]) + ' 	| ' + str("%.4f" % probs[i, j]) + '\n')

	#подсчёт исходной вероятности
	dp = 0
	#for k, v in assembly.items():
	for k in contigs:
		d = depth[k]
		res = 0
		wg.write('track name=My-depth color=255,0,0 group=ALE priority=2\nfixedStep chrom=' + k + ' start=1 step=1\n')
		for i in range(len(d)-1):
			t = math.log(probs[d[i], d[i+1]])
			res = res + t
			wg.write(str(100 * t) + '\n')
		out2.write('>>>' + k + ' ' + str(int(np.amax(d))) + ' \n')
		out2.write(str(Decimal(res)) + '\n')
		dp =  Decimal(dp) + Decimal(res)
	print(str(Decimal(dp)) + '\n')
if __name__ == '__main__':
    main()