import math
from decimal import Decimal
import znorm
import numpy as np

# paths to input files
gc_stat_file = "depth-stat"
depth_stat_file = "gc-stat"
# paths to output files
distr_file = "res.my100.depth.wig"

# dictionaries/maps for all stuff
full_depth = {}
full_gc = {}
full_ans = {}

# some constants
lnfactconst2 = 0.918938533204672741780329

# functions for normalization
def lnfact2(input):
    input2 = input * input
    input3 = input * input2
    input5 = input3 * input2
    input7 = input5 * input2
    return (input - 0.5) * math.log(input) - input + lnfactconst2 + 1.0 / (12.0 * input) - 1.0 / (
        360.0 * input3) + 1.0 / (
               1260.0 * input5) - 1.0 / (1680.0 * input7)


def getNegBinomZnorm(r):
    ans = 0.0
    n = 0
    diff = 1.0
    while ((diff > 1e-8 or n < r) and n < 10000000):
        diff = math.exp(2.0 * lnfact2(r + n) - 2.0 * lnfact2(r) - 2.0 * lnfact2(n + 1.0) + (r + n) * math.log(4.0));
        ans += diff
        n += 1
    return ans


def computeNormaliziedDepthGCParameters(depthNormalizer, depthNormalizerCount, negBinomParam_r,
                                        negBinomParamZnorm_r, avgDepth):
    for j in range(0, 101):
        if (depthNormalizerCount[j] > 0):
            depthNormalizer[j] = depthNormalizer[j] / depthNormalizerCount[j]
        else:
            depthNormalizer[j] = avgDepth

        if (depthNormalizer[j] < 10.0):
            depthNormalizer[j] = 10.0

        negBinomParam_r[j] = depthNormalizer[j]

        if (math.floor(negBinomParam_r[j]) < 2047):
            negBinomParamZnorm_r[j] = znorm.negBinomZ[math.floor(negBinomParam_r[j])]
        else:
            negBinomParamZnorm_r[j] = getNegBinomZnorm(negBinomParam_r[j]);

    return negBinomParam_r, negBinomParamZnorm_r


def main():
    gcstat = open(gc_stat_file, 'w')
    depthstat = open(depth_stat_file, 'w')

    dp_lines = depthstat.readlines()
    gc_lines = gcstat.readlines()

    # parse data and initialize variables
    for i, x in enumerate(dp_lines):
        if(x[0] == '>'):
            name = x[1:]
            full_depth[name] = map(int, dp_lines[i + 1][:-1].split(' '))
            full_gc[name] = map(int, gc_lines[i + 1][:-1].split(' '))
            full_ans[name] = len(full_gc[name])
        else:
            continue

    maxdepth = 0
    depthAvgSum = 0.0
    depthAvgNorm = 0
    gc_count = np.zeros(102)
    depthNormalizer = np.zeros(102)
    negBinomParam_r = np.zeros(102)
    negBinomParamZnorm_r = np.zeros(102)

    # count some additional parameters, incl. Z-normalization
    for k, v in full_depth.items():
        d = np.max(v)
        gc_cur = full_gc[k]
        if d > maxdepth:
            maxdepth = d
        # if (np.sum(v) != 0):
        for i, x in enumerate(v):
            g = gc_cur[i]
            depthAvgSum += d
            depthAvgNorm += 1
            gc_count[g] += 1
            depthNormalizer[g] += d

    negBinomParam_r, negBinomParamZnorm_r = computeNormaliziedDepthGCParameters(depthNormalizer, gc_count,
                                                                                negBinomParam_r, negBinomParamZnorm_r,
                                                                                depthAvgSum / depthAvgNorm)
    # count distribution for depth
    depth_pos = np.zeros((101, maxdepth + 1))
    for k, v in full_depth.items():
        depth_cur = v
        gc_cur = full_gc[k]
        for i, x in enumerate(depth_cur):
            if (np.sum(depth_cur) != 0):
                d = x
                g = gc_cur[i]
                for j in range(0, 101):
                    if (g == j):
                        depth_pos[j, d] += 1

    # count log-likelihood and z-norm, write distribution table to file (optional)
    pl_emp = open('pl1.txt', 'w')
    for i in range(0, 101):
        sm = np.sum(depth_pos[i])
        for j in range(len(depth_pos[i])):
            depth_pos[i, j] = math.log(depth_pos[i, j] / sm) - negBinomParamZnorm_r[i] if sm != 0 and depth_pos[
                i, j] != 0 else -120.0
            pl_emp.write(str(depth_pos[i, j] / sm if sm != 0 else 0.0) + ' ');
        pl_emp.write('\n')
    pl_emp.close()

    #write all the stuff to *.wig file (to run it in IGV afterwards)
    final_stat = open(distr_file, 'w')
    for k, v in sorted(full_ans.items()):
        final_stat.write("track name=100-ALE-depth color=255,0,0 group=ALE priority=2 \n" +
                     "fixedStep chrom=" + k + " start=1  step=1\n")
        d = full_depth[k]
        g = full_gc[k]
        for i, x in enumerate(v):
            full_ans[k][i] = depth_pos[g[i], d[i]]
            final_stat.write(str(full_ans[k][i]) + '\n')
    final_stat.close()



if __name__ == '__main__':
    main()

