#!/usr/bin/env python

## boundaries for MLL (KMT2A) exon 3 and 27
## in human genome reference GRCh37p13

LABEL_EX3 = "MLL_EXON3"
LABEL_EX27 = "MLL_EXON27"

BEDFIL_EXON_LABELS = {
    LABEL_EX3:
    ("11", 118342353, 118345053), # chromosome, start, end
    LABEL_EX27:
    ("11", 118373091, 118377381)
    }

## samtools executable (version 1.3)
CMD_SAMTOOLS = "samtools"

FILEXT_BAM = ".bam"
FILEXT_BAI = ".bai"
FILPRFX_BED = "MFITB"
FILPRFX_TMP = "MFITT"
FILPRFX_SRT = "MFITS"

import sys
import os

def createBEDfile():
    fnam_bed = os.tempnam(os.getcwd(), FILPRFX_BED)
    infh = open(fnam_bed, 'w')
    for label in BEDFIL_EXON_LABELS.keys():
        (chrom, start, end) = BEDFIL_EXON_LABELS[label]
        infh.write("%s\t%i\t%i\t%s\n" % (chrom, start, end, label))
    infh.close()
    return fnam_bed

def loadFileNames(fnam_list, logfh = sys.stderr):
    bamfnamd = {}
    bamfnams = []
    infh = open(fnam_list, 'r')
    while 1:
        lin = infh.readline()
        if not lin:
            break
        fld = lin.split()
        nf = len(fld)
        if nf > 0 and fld[0][0] != '#':
            label = fld[0]
            fnam = fld[-1]
            fnprfx, fnext = os.path.splitext(fnam)
            if fnext and fnext != FILEXT_BAM:
                sys.exit("ERROR: unexpected file name extension '%s' for BAM files.\n" %
                         (fnext))
            if fnprfx in bamfnams:
                logfh.write("WARNING: BAM file '%s'.bam occurs multiple times in list '%s'\n" %
                            (fnprfx, fnam_list))
            else:
                bamfnams.append(fnprfx)
            if bamfnamd.has_key(label):
                sys.exit("ERROR: sample label '%s' occurs multiple times in list '%s'\n" %
                         (label, fnam_list))
            bamfnamd[label] = fnprfx
    infh.close()
    logfh.write("Read %i file names from '%s'." % (len(bamfnams), fnam_list))
    
    return bamfnamd
        
def bedcov(fnprfx, bedfnam, oufh = sys.stdout, logfh = sys.stderr):
    from subprocess import call

    oufnprfx = os.path.basename(fnprfx)
    bamfnam_sorted = " %s.srt%s" % (oufnprfx, FILEXT_BAM)

    cmd1str = "%s sort -m 2G -T %s -o %s %s%s" % \
              (CMD_SAMTOOLS, os.tempnam(os.getcwd(), FILPRFX_SRT), bamfnam_sorted, fnprfx, FILEXT_BAM)
    cmd2str = "%s index %s" % (CMD_SAMTOOLS, bamfnam_sorted)
    cmd3str = "%s bedcov %s %s" % (CMD_SAMTOOLS, bedfnam, bamfnam_sorted)

    is_sorted = os.access(bamfnam_sorted, os.F_OK)
    if is_sorted:
        logfh.write("File %s already present ...\n" % (bamfnam_sorted))
    else:
        logfh.write(cmd1str + '\n')
        rc = call(cmd1str.split())
        if rc:
            sys.exit("ERROR code %i when executing: %s\n" % (rc, cmd1str))

    bainam = bamfnam_sorted + FILEXT_BAI
    if is_sorted and os.access(bainam, os.F_OK):
        logfh.write("File %s already present ...\n" % (bainam))
    else:
        logfh.write(cmd2str + '\n')
        rc = call(cmd2str.split())
        if rc:
            sys.exit("ERROR code %i when executing: %s\n" % (rc, cmd2str))
    
    logfh.write(cmd3str + '\n')
    rc = call(cmd3str.split(), stdout=oufh, stderr=logfh)
    if rc:
        sys.exit("ERROR code %i when executing: %s\n" % (rc, cmd3str))

    return

    

def processFiles(bamfnamd, datadir=None):

    fnam_bed = createBEDfile()

    for label in bamfnamd.keys():
        fnprfx = bamfnamd[label]
        if datadir is None:
            fpathprfx = fnprfx
        else:
            fpathprfx = path.join(datadir, fnprfx)

        oufnam = label + '.out'
        oufh = open(oufnam, 'w')
        bedcov(fpathprfx, fnam_bed, oufh=oufh)
        oufh.close()
    

    # os.remove(fnam_bed)
    return

if __name__ == '__main__':
    nargs = len(sys.argv)

    if nargs < 2 or nargs > 3:
        exit("usage %s <list of file names> [<data dir>]\n" % (sys.argv[0]))

    fnam_list = sys.argv[1]
    datadir = None
    if nargs > 2:
        datadir = sys.argv[2]
    
    bamfnamd = loadFileNames(fnam_list)

    processFiles(bamfnamd, datadir)
    sys.exit()
