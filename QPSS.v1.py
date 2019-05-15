#!/usr/bin/env python

## QPSS v1.0
## updated 4/3/2019

def create_rarebfile(filename):
    ## create plink files for specified chromosome
    cmd_bfile = plink + " --bfile " + args.bfile + " --allow-no-sex --prune --chr " + chrn + " --make-bed --out " + cwd + "/tmp/" + filename + ".rare" + " --max-maf " + args.max_maf + " --mac " + args.min_ac
   
    if args.G_position:
        cmd_bfile = cmd_bfile + " --from-bp " + args.G_position[0] + " --to-bp " + args.G_position[1]
    if args.pheno:
        cmd_bfile = cmd_bfile + " --pheno " + args.pheno
    if args.pheno_name:
        cmd_bfile = cmd_bfile + " --pheno-name " + args.pheno_name
    if args.extract:
        cmd_bfile = cmd_bfile + " --extract " + args.extract
    if args.exclude:
        cmd_bfile = cmd_bfile + " --exclude " + args.exclude
    if args.filter:
        cmd_bfile = cmd_bfile + " --filter " + args.filter[0] + " " + args.filter[1]
    if args.keep:
        cmd_bfile = cmd_bfile + " --keep " + args.keep
    if args.remove:
        cmd_bfile = cmd_bfile + " --remove " + args.remove

    subprocess.call(cmd_bfile, shell = True, stdout=devnull)
    

    if args.G_position and (int(args.G_position[1]) - int(args.G_position[0])) < 1000000:
        #print(args.Gposition[0])
        cmd_rareG = plink + " --bfile " + cwd + "/tmp/" + filename + ".rare --allow-no-sex --recode A --out " + cwd + "/tmp/" + filename + ".rareG"
        subprocess.call(cmd_rareG, shell = True, stdout=devnull)
        rawG = np.genfromtxt(cwd + "/tmp/" + filename + ".rareG.raw", skip_header = 1, missing_values = 'NA', unpack = True)
        col_sum = np.nansum(rawG[6:], axis=0)
        rawG_index = np.where(col_sum >= 1)
        idG = np.genfromtxt(cwd + "/tmp/" + filename + ".rareG.raw", dtype='str', usecols=[0,1], skip_header = 1, unpack = True)
        ids = idG.transpose()[rawG_index]
        #np.savetxt(cwd + "/tmp/" + fname + ".rare.rawG.ids", ids, fmt='%i')
        np.savetxt(cwd + "/tmp/" + fname + ".rare.rawG.ids", ids, delimiter=" ", fmt="%s")
        cmd_rareG2 = plink + " --bfile " + cwd + "/tmp/" + filename + ".rare --allow-no-sex --keep " + cwd + "/tmp/" + filename + ".rare.rawG.ids --make-bed --out " + cwd + "/tmp/" + filename + ".rare"
        subprocess.call(cmd_rareG2, shell = True, stdout=devnull)


def create_raw(W_start, W_end, filename):
    cmd_raw = plink + " --bfile " + cwd + "/tmp/" + filename + ".rare --allow-no-sex --chr " + chrn + " --from-bp " + str(W_start) + " --to-bp " + str(W_end) + " --recode A --out " + cwd + "/tmp/" + filename + ".rare"
    #cmd_raw = ["plink19b3x", "--bfile", "tmp/" + filename + ".rare", "--allow-no-sex", "--chr", chrn, "--from-bp", str(W_start), "--to-bp", str(W_end), "--recode", "A", "--out", "tmp/" + filename + ".rare"]
    
    subprocess.call(cmd_raw, shell = True, stdout=devnull)
    #subprocess.Popen(cmd_raw, stdout=FNULL)

def scan_stat(chrn, W_start, W_end, col_sum, phe):
    nG = len(phe)
    #phe = (raw[5] - np.mean(raw[5])) / np.std(raw[5], ddof=0)
    rawin = phe[np.where(col_sum >= 1)]
    rawout = np.delete(phe, np.where(col_sum >= 1))
    re =[int(chrn), W_start, W_end]
    meanG = np.mean(phe)
    sigmaGsq = np.sum((phe - meanG)**2)/nG
    meanin = np.mean(rawin)
    meanout = np.mean(rawout)
    sigmaWsq = (np.sum((rawin - meanin)**2) + np.sum((rawout - meanout)**2))/nG
    
    #indicator = int(meanin > meanout)
    if meanin > meanout:
        indicator = 1
    else:
        indicator = -1
        
    lnLRW = nG/2*math.log(float(sigmaGsq)/(float(sigmaWsq))) ## calculate part of LRW (1/sigmaW^2)
    re = re + [round(lnLRW,4), indicator, round(meanin,4), round(meanout,4)]

    return(re)

## simulation in each thread
def MonteCarlo(i, j, l, ssl, phe):
    count = 0
    sim_ss_list = []
    #[fi for fi in map(fn, lst) if fi is not None]
    #counta = [int(scan_stat(args.chrn, w_start, w_end, row_sum, np.random.permutation(phe))[3] >= ss[3]) for num in xrange(i)]
   
    for _ in xrange(i):
        phe = np.random.permutation(phe)
        sim_ss = scan_stat(chrn, w_start, w_end, row_sum, np.random.permutation(phe))
        sim_ss_list = sim_ss_list + [sim_ss[3]]
        count += int(sim_ss[3] >= ss[3])

    l[j] = count
    #sim_ss_list = np.array(sorted(sim_ss_list, reverse=True))
    ssl[j] = np.trim_zeros(sim_ss_list)
    #print(l)
    #print(ssl)
    #l[j] = np.sum(counta)
    

## multiple threads	
def mp(k, sim, phe):
    processes = []
    manager = multiprocessing.Manager() # to share memories between processes
    l = manager.list([0]*k)
    ssl = manager.list([0]*k)
    
    # replication number
    rep = [(sim + x)//k  for x in range(k)]

    for j in xrange(k):
        i = rep[j]
            
        if __name__ == '__main__':
            phe = np.random.permutation(phe)
            p = multiprocessing.Process(target = MonteCarlo, args = (i, j, l, ssl, np.random.permutation(phe)))
            processes.append(p)
            p.start()
    
    [p.join() for p in processes]

    gc.collect()
    return(sum(l[:]), ssl[:])

## GPD approximation
def params(mexc, y): # y = np.array
    t = (y[(mexc-1)] + y[mexc])/2
    z = y[:mexc] - t
    z2 = z**2
    m = np.mean(z)
    m2 = np.mean(z2)
    a_hat = m*m2/(m2-m**2)/2
    k_hat = (m**2/(m2-m**2)-1)/2
    
    return([a_hat, k_hat, t])
    
def goodness(mexc, y, filename):
    #print(y)
    t = (y[(mexc-1)] + y[mexc])/2
    #print(t)
    z = y[:mexc] - t
    #print(z)
    np.savetxt(cwd + "/tmp/" + filename + ".y.gpd", z, fmt='%f')
    #cmd_gpd = "Rscript --vanilla --slave " + abs_path + "/gpd.r " + cwd + "/tmp/" + fname + ".y.gpd"
    cmd_gpd = ["Rscript", "--vanilla", "--slave", abs_path + "/gpd.r", cwd + "/tmp/" + fname + ".y.gpd"]
    gpd = subprocess.Popen(cmd_gpd, stdout=subprocess.PIPE, stderr=devnull)
    gpd.wait()
    mexc_re = gpd.communicate()[0].split()
    return(mexc_re) ##[p, mexc]
    

def gpd(a_hat, k_hat, z0, nsim, mexc): 
    F = 1-(1-k_hat*z0/a_hat)**(1/k_hat)
    p = float(mexc)/nsim*(1-F)
    
    return(p)
        
if __name__ == '__main__':
    from fractions import Fraction
    import argparse
    import os
    import sys
    import time
    import subprocess
    import multiprocessing
    import numpy as np
    import math
    import ntpath
    import shutil
    import random
    import gc
    import datetime
    
    #import scipy.stats as sp

    parser = argparse.ArgumentParser()
    parser.add_argument("--bfile", required=True, help = "input plink bfile")
    parser.add_argument("--chr", required=True, help = "chromosome #")
    parser.add_argument("--G-position", "--G", help = "start and end positions of G (e.g., 12345 56789)", nargs = 2)
    parser.add_argument("--pheno", help = "phenotype data with FID, IID, phenotype")
    parser.add_argument("--pheno-name", help = "column name for phenotype")
    parser.add_argument("--max-maf", default = '0.05', type = str, help = "max threshold of MAF to define rare variant")
    parser.add_argument("--min-ac", default = '5', type = str, help = "min threshold of MAC to define rare variant")
    parser.add_argument("--extract", help = "to extract variants")
    parser.add_argument("--exclude", help = "to exclude variants")
    parser.add_argument("--filter", help = "filename value", nargs = 2)
    parser.add_argument("--keep", help = "to keep individuals")
    parser.add_argument("--remove", help = "to remove individuals")
    parser.add_argument("--W-position", "--W", type = int, help = "start and end positions of W (e.g., 12345 56789)", nargs = 2)
    parser.add_argument("--W-file", help = "set ss.out file")
    parser.add_argument("--W-fixed", "--wf", default = 2000, type = int, help = "fixed window size")
    parser.add_argument("--W-slide", "--ws", default = 1000, type = int, help = "slide window size")
    parser.add_argument("--out", help = "output file name")
    parser.add_argument("--perm", "--p", const = 'True', help = "whether p is computed. If GPD approximation is used, add 'gpd' as --perm gpd", nargs = '?')
    parser.add_argument("--threads", default = 10, type = int, help = "# of threads")
    parser.add_argument("--max-sim", default = 10000, type = int, help = "max# of simulation")
    args = parser.parse_args()

    ##
    plink = "plink19b67"

    ## output path and set the file name
    try:        
        cwd = os.path.dirname(os.path.abspath(args.out))
        fname = ntpath.basename(args.out)
    except:
        cwd = os.getcwd()
        fname = ntpath.basename(args.bfile)
    #print(cwd)

    ## absolutive path
    abs_path = os.path.dirname(os.path.abspath(__file__))
    
    # make tmp directory for working under the current directory
    if not os.path.exists(cwd + "/tmp/"):
        os.mkdir(cwd + "/tmp/")

    if os.path.isfile(cwd + "/tmp/" + fname + ".rare.raw"):
        os.remove(cwd + "/tmp/" + fname + ".rare.raw")

    if os.path.exists(cwd + '/' + fname + '.ss.out'):
        os.remove(cwd + '/' + fname + '.ss.out')

    log = open(cwd + '/' + fname + '.log', 'w')    
    devnull = open(os.devnull, 'wb')

    argsstr = ' '.join(sys.argv)
    argsstr = argsstr.replace ('--', '\n --')

    sys.stdout = log
    print("QPSS v1.0 \(@_@)/")
    print("Start time: " + str(datetime.datetime.today().strftime("%Y/%m/%d %H:%M:%S")))
    print("")
    print(argsstr)
    sys.stdout.flush()
    
    chrn = args.chr
    
    create_rarebfile(fname)

    ## create chr, start, end list
    if args.W_file:
        ssout = np.loadtxt(args.W_file, dtype="int", usecols=(range(0, 3)))
        if ssout.ndim == 1:
            ssout = [ssout]
    elif args.W_position:
        ssout = [[int(chrn), int(args.W_position[0]), int(args.W_position[1])]]
    else:
        bim = np.loadtxt(cwd + "/tmp/" + fname + ".rare.bim", usecols=[3])
        G_start = int(bim[0])
        G_end = int(bim[-1])
        w_start = G_start
        ssout = []

        while w_start < G_end:
            w_end = w_start + args.W_fixed - 1
            ssout.append([int(chrn), w_start, w_end])
            w_start = w_start + args.W_slide
    sys.stdout = sys.__stdout__
    print("G = " + str(ssout[0][1]) + " - " + str(ssout[-1][1]))
    for term in ssout:
        w_start = int(term[1])
        w_end = int(term[2])
        if os.path.exists(cwd + "/tmp/" + fname + ".rare.raw"):
            os.remove(cwd + "/tmp/" + fname + ".rare.raw")
        create_raw(w_start, w_end, fname)
        sys.stdout = log
        print("")
        print("---")
        print("Chromosome " + chrn)
        print("Start position = " + str(w_start))
        print("End position = " + str(w_end))
        print("---")
        sys.stdout.flush()
        
        try:
            raw = np.genfromtxt(cwd + "/tmp/" + fname + ".rare.raw", skip_header = 1, missing_values = 'NA', unpack = True)
            row_sum = np.nansum(raw[6:], axis=0)
            ncol = np.size(raw[6:], axis=0)
            phe = raw[5]
            ss = scan_stat(chrn, w_start, w_end, row_sum, phe)

            print("# of variants within the window = " + str(ncol))
            print("Mean within the window = " + str(ss[5]))
            print("Mean outside of the window = " + str(ss[6]))
            sys.stdout.flush()
            
            ## permutation test
            if args.perm:
                nsim = 1000
                nsimlist = [nsim]
                list_sum = [0]
                list_ss = []
                sys.stdout = sys.__stdout__
                print('--- CHR '+ chrn + ': ' + str(w_start) + ' - ' + str(w_end))
                print("Computing permutation p-value "),
                while(list_sum[0] < 100):
                    if(sum(nsimlist) > args.max_sim):
                        break
                    sys.stdout.write('\033[33m=>\033[1D\033[0m')
                    #print('\033[33m' + str(sum(nsimlist)) + '\033[0m'),
                    mpre = mp(args.threads, nsim, phe)
                    #list_sum.append(sum(mpre[0]))
                    list_sum.append(mpre[0])
                    list_ss.extend(mpre[1])
                    list_sum = [sum(list_sum)]
                    nsimlist = [sum(nsimlist)]
                    nsim = nsimlist[0]*10 - nsimlist[0]
                    nsimlist.append(nsim)
                nsimp = sum(nsimlist[:-1])
                sys.stdout.write('\033[33m> ' + str(sum(list_sum)) + '/' + str(nsimp) + '\033[0m' + '\n')    
                #print("")
                if args.perm == "gpd" and list_sum[0] < 100:
                    i = 0
                    if list_sum[0] == 0:
                        sys.stdout = sys.__stdout__
                        #print("")
                        sys.stdout.write("Computing permutation p-value with GPD approximation ")
                        nsim = args.max_sim
                        while(list_sum[0] < 1):
                            i = i + 1
                            mpre = mp(args.threads, nsim, phe)
                            #list_sum.append(sum(mpre[0]))
                            list_sum.append(mpre[0])
                            list_ss.extend(mpre[1])
                            list_sum = [sum(list_sum)]
                            sys.stdout.write('\033[33m' + str(sum(list_sum)) + '>\033[1D\033[0m')
                        nsimp = nsimp + args.max_sim*i
                        sys.stdout.write('\033[33m> ' + str(sum(list_sum)) + '/' + str(nsimp) + '\033[0m' + '\n')    
                    
                    sys.stdout = sys.__stdout__
                    #print("")
                    sys.stdout.write("Performing GPD goodness-of-fit test" + '\n')
                    y = np.array(sorted([item for sublist in list_ss for item in sublist], reverse=True))
                    mexc = 250 ## start mexc
                    a_hat = params(mexc, y)[0]
                    k_hat = params(mexc, y)[1]
                    t = params(mexc, y)[2]
                    z0 = ss[3] - t
                    m = len(y)
                    mexc_re = goodness(mexc, y, fname)
                    sys.stdout = log
                    print("---")
                    print("GPD goodness-of-fit test p-value = " + mexc_re[0])
                    print("Nexc = " + mexc_re[1])
                    sys.stdout.flush()
                    a_hat = params(int(mexc_re[1]), y)[0]
                    k_hat = params(int(mexc_re[1]), y)[1]
                    t = params(int(mexc_re[1]), y)[2]
                    z0 = ss[3] - t
                    p = gpd(a_hat, k_hat, z0, nsimp, int(mexc_re[1]))
                    if(mexc_re[0] < 0.05):
                        gpdtxt = "GPD*"
                    else:
                        gpdtxt = "GPD"
                    ss = ss + [list_sum[0], str(nsimp), p, gpdtxt]
                    
                else:
                    p = float(list_sum[0])/nsimp
                    ss = ss + [list_sum[0], nsimp, p, "Permutation"]
              
                sys.stdout = log
                print("Number of permutation values = " + str(ss[7]))
                print("Number of permutation replicates = " + str(ss[8]))
                print("p-value = " + str(ss[9]))
                print("Method = " + str(ss[10]))
                sys.stdout.flush()
            
                            
            #with open(cwd + '/' + fname + '.ss.out', mode = "a") as fout:
                #np.savetxt(fout, [ss], delimiter=" ")
            #print("---")
            #np.savetxt(cwd + "/tmp/" + fname + ".rare.rawG.ids", ids, delimiter=" ", fmt="%s")
            
            with open(cwd + '/' + fname + '.ss.out', mode = "a") as fout:
                for term in ss:
                    fout.write(str(term) + ' ')
                fout.write('\n')
            #sys.stdout = sys.__stdout__
            #print("")
        
        except IOError:
            print("# of variants within the window = 0")
            
            
    log.close()            
            
            
        

        
