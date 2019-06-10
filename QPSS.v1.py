#!/usr/bin/env python

## QPSS v1.0
## updated 4/3/2019

def create_rarebfile(filename):
    ## create plink files for specified chromosome
    cmd_bfile = [plink, "--bfile", args.bfile, "--allow-no-sex", "--prune", "--chr", chrn, "--make-bed", "--out", cwd + "/tmp/" + filename + ".rare", "--max-maf", args.max_maf, "--mac", args.min_ac]
   
    if args.G_position:
        cmd_bfile = cmd_bfile + ["--from-bp", args.G_position[0], "--to-bp", args.G_position[1]]
    if args.pheno:
        cmd_bfile = cmd_bfile + ["--pheno", args.pheno]
    if args.pheno_name:
        cmd_bfile = cmd_bfile + ["--pheno-name", args.pheno_name]
    if args.extract:
        cmd_bfile = cmd_bfile + ["--extract", args.extract]
    if args.exclude:
        cmd_bfile = cmd_bfile + ["--exclude", args.exclude]
    if args.filter:
        cmd_bfile = cmd_bfile + ["--filter", args.filter[0], args.filter[1]]
    if args.keep:
        cmd_bfile = cmd_bfile + ["--keep", args.keep]
    if args.remove:
        cmd_bfile = cmd_bfile + ["--remove", args.remove]
    
    try:
        subprocess.check_call(cmd_bfile, stdout=devnull)
    except subprocess.CalledProcessError as e:
        print e.output


def create_raw(W_start, W_end, filename):
    cmd_raw = [plink, "--bfile", cwd + "/tmp/" + filename + ".rare", "--allow-no-sex", "--chr", chrn, "--from-bp", str(W_start), "--to-bp",  str(W_end), "--recode",  "A",  "--out", cwd + "/tmp/" + filename + ".rare"]
   
    try:
        subprocess.call(cmd_raw, stdout=devnull, stderr=devnull)
    except subprocess.CalledProcessError as e:
        print e.output

def scan_stat(chrn, W_start, W_end, col_sum, phe):
    nG = len(phe)
    rawin = phe[np.where(col_sum >= 1)]
    rawout = np.delete(phe, np.where(col_sum >= 1))
    re =[int(chrn), W_start, W_end]
    meanG = np.mean(phe)
    sigmaGsq = np.sum((phe - meanG)**2)/nG
    meanin = np.mean(rawin)
    meanout = np.mean(rawout)
    sigmaWsq = (np.sum((rawin - meanin)**2) + np.sum((rawout - meanout)**2))/nG
    
    #sign function
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
  
    for _ in xrange(i):
        phe = np.random.permutation(phe)
        sim_ss = scan_stat(chrn, w_start, w_end, row_sum, np.random.permutation(phe))
        sim_ss_list = sim_ss_list + [sim_ss[3]]
        count += int(sim_ss[3] >= ss[3])

    l[j] = count
    ssl[j] = sim_ss_list
    
## multiple threads
def mp(k, sim, phe): ## k = args.threads
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

# run R for goodness of fit test -- > requires "goft" package
def goodness(mexc, y, filename):
    t = (y[(mexc-1)] + y[mexc])/2
    z = y[:mexc] - t
    np.savetxt(cwd + "/tmp/" + filename + ".y.gpd", z, fmt='%f')
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
    import progressbar
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
    parser.add_argument("--perm", "--p", const = 'True', help = "whether p is computed. Add 'gpd' as --perm gpd for GPD approximation", nargs = '?')
    parser.add_argument("--threads", default = 10, type = int, help = "# of threads")
    parser.add_argument("--max-sim", default = 10000, type = int, help = "max# of simulation")
    args = parser.parse_args()

    ##
    plink = "plink19b67" ## <-- change to plink command on your server

    ## output path and set the file name
    try:        
        cwd = os.path.dirname(os.path.abspath(args.out))
        fname = ntpath.basename(args.out)
    except:
        cwd = os.getcwd()
        fname = ntpath.basename(args.bfile)

    ## absolutive path to this python file
    abs_path = os.path.dirname(os.path.abspath(__file__))
    
    # make tmp directory for working under the current directory
    if not os.path.exists(cwd + "/tmp/"):
        os.mkdir(cwd + "/tmp/")
    ## delete .raw file if it exists
    if os.path.exists(cwd + "/tmp/" + fname + ".rare.raw"):
        os.remove(cwd + "/tmp/" + fname + ".rare.raw")
    
    ## delete .ss.out file if it exists
    if os.path.exists(cwd + '/' + fname + '.ss.out'):
        os.remove(cwd + '/' + fname + '.ss.out')
    
    ## open log file
    log = open(cwd + '/' + fname + '.log', 'w')    
    devnull = open(os.devnull, 'wb')
    
    ## get arguments
    argsstr = ' '.join(sys.argv)
    argsstr = argsstr.replace ('--', '\n --')
    
    chrn = args.chr
    
    ## output STDOUT
    sys.stdout = sys.__stdout__
    message = '''
    Welcome to   ___  ____  ____ ____  
                / _ \|  _ \/ ___/ ___| 
               | | | | |_) \___ \___ \    
               | |_| |  __/ ___) |__) | 
                \__\_\_|   |____/____/   
    '''
    print(message)
    
    ## output start time and arguments in the log file
    sys.stdout = log
    print("QPSS v1.0 \(@_@)/")
    print("Start time: " + str(datetime.datetime.today().strftime("%Y/%m/%d %H:%M:%S")))
    print("Working directory:" + str(os.getcwd()))
    print("")
    print(argsstr)
    sys.stdout.flush()

    ## create plink files for rare varaints
    create_rarebfile(fname)
    
    ## get G_start and G_end
    bim = np.loadtxt(cwd + "/tmp/" + fname + ".rare.bim", usecols=[3])
    G_start = int(bim[0])
    G_end = int(bim[-1])
        
    ## create list with [chr, start, end] list for window(s)
    if args.W_file:
        ssout = np.loadtxt(args.W_file, dtype="int", usecols=(range(0, 3)))
        if ssout.ndim == 1:
            ssout = [ssout]
    elif args.W_position:
        ssout = [[int(chrn), int(args.W_position[0]), int(args.W_position[1])]]
    else:
        w_start = G_start
        ssout = []
        while w_start < G_end:
            w_end = w_start + args.W_fixed - 1
            ssout.append([int(chrn), w_start, w_end])
            w_start = w_start + args.W_slide
    
    ## run scan_stat for each window
    sys.stdout = log
    print("")
    print("Large genetic region G: " + str(G_start) + "-" + str(G_end))
    print("# of windows within G: " + str(len(ssout)))
    sys.stdout.flush()
    
    ## progressbar
    widgets = [
        'QPSS progress: ', progressbar.Percentage(),
        ' ', progressbar.Bar(),
        ' ', progressbar.ETA(),
        ' ', progressbar.FileTransferSpeed(),
    ]
    pb = progressbar.ProgressBar(widgets=widgets, maxval=len(ssout)).start()
    
    i = 0 ## for progressbar
    k = 0 ## count # of windows contining no variants
    fout = open(cwd + '/' + fname + '.ss.out', mode = "a")
    fout.write("chr start end loglr sign meanin meanout n_variants p method p_goodness_fit\n")
    for term in ssout:
        w_start = int(term[1])
        w_end = int(term[2])
        if os.path.exists(cwd + "/tmp/" + fname + ".rare.raw"):
            os.remove(cwd + "/tmp/" + fname + ".rare.raw")
        create_raw(w_start, w_end, fname)

        try:
            raw = np.genfromtxt(cwd + "/tmp/" + fname + ".rare.raw", skip_header = 1, missing_values = 'NA', unpack = True)
            row_sum = np.nansum(raw[6:], axis=0)
            ncol = np.size(raw[6:], axis=0)
            phe = raw[5]
            ss = scan_stat(chrn, w_start, w_end, row_sum, phe)
            ss = ss + [ncol]
            
            ## permutation test
            if args.perm:
                totalnsim = 0
                nsim = 1000
                sumcount = 0
                list_ss = []
                while(sumcount < 100):
                    totalnsim = totalnsim + nsim
                    if(totalnsim > args.max_sim):
                        break
                    mpre = mp(args.threads, nsim, phe)
                    sumcount = sumcount + mpre[0]
                    list_ss = sorted(list_ss + sum(mpre[1], []), reverse = True)[:500]
                nsimp = totalnsim

                ## GPD approximation
                if args.perm == "gpd" and sumcount < 100:
                    if sumcount == 0:
                        while(sumcount < 1):
                            totalnsim = totalnsim + nsim
                            mpre = mp(args.threads, nsim, phe)
                            sumcount = sumcount + mpre[0]
                            list_ss = sorted(list_ss + sum(mpre[1], []), reverse = True)[:500]
                        nsimp = totalnsim
                    
                    y = np.array(sorted(list_ss, reverse=True))
                    mexc = 250 ## start mexc
                    a_hat = params(mexc, y)[0]
                    k_hat = params(mexc, y)[1]
                    t = params(mexc, y)[2]
                    z0 = ss[3] - t
                    m = len(y)
                    ## goodness-fit test
                    mexc_re = goodness(mexc, y, fname)
                    sys.stdout = log
                    a_hat = params(int(mexc_re[1]), y)[0]
                    k_hat = params(int(mexc_re[1]), y)[1]
                    t = params(int(mexc_re[1]), y)[2]
                    z0 = ss[3] - t
                    p = gpd(a_hat, k_hat, z0, nsimp, int(mexc_re[1]))
                    #ss = ss + [sumcount, nsimp, round(p, 3), "GPD", round(float(mexc_re[0]), 3)]
                    ss = ss + ["{:.3g}".format(p), "GPD", "{:.3g}".format(float(mexc_re[0]))]
                    
                else:
                    p = float(sumcount)/nsimp
                    ss = ss + ["{:.3g}".format(p), "Permutation", "."]

            
            with open(cwd + '/' + fname + '.ss.out', mode = "a") as fout:
                for term in ss:
                    fout.write(str(term) + ' ')
                fout.write('\n')
        except IOError:
            k = k + 1
        i = i + 1
        pb.update(i)
        time.sleep(0.001)
    pb.finish()
    print(str(k) + " windows contain 0 variants")
    print("")
    print("End time: " + str(datetime.datetime.today().strftime("%Y/%m/%d %H:%M:%S")))
            
    log.close()            
