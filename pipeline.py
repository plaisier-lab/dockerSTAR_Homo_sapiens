##########################################################
## OncoMerge:  alignment.py                             ##
##  ______     ______     __  __                        ##
## /\  __ \   /\  ___\   /\ \/\ \                       ##
## \ \  __ \  \ \___  \  \ \ \_\ \                      ##
##  \ \_\ \_\  \/\_____\  \ \_____\                     ##
##   \/_/\/_/   \/_____/   \/_____/                     ##
## @Developed by: Plaisier Lab                          ##
##   (https://plaisierlab.engineering.asu.edu/)         ##
##   Arizona State University                           ##
##   242 ISTB1, 550 E Orange St                         ##
##   Tempe, AZ  85281                                   ##
## @github: plaisier-lab/dockerSTAR_Homo_sapeins        ##
## @Author:  Chris Plaisier                             ##
## @License:  GNU GPLv3                                 ##
##                                                      ##
## If this program is used in your analysis please      ##
## mention who built it. Thanks. :-)                    ##
##########################################################

from subprocess import *
from shutil import move
import os

## Read in manifest of samples to process
samples = []
with open('manifest.csv','r') as inFile:
    while 1:
        inLine = inFile.readline()
        if not inLine:
            break
        splitUp = inLine.strip().split(',')
        samples.append(splitUp[0])

## Make output folder if doesn't exist
if not os.path.exists('/fastq/output'):
    os.mkdir('/fastq/output')

### RNA-seq pipeline overview:
##   1. Run cutadapt to remove adaptor sequences
##   2. Align with STAR to create SAM file
##   3. Tabulate counts for each transcript using htseq-count
gexpMatrix = {} # gexpMatrix[transcript][sample] = {'coverage':<>, 'FPKM':<>}
output = {}
outOrder = []
for s1 in samples:
    if not os.path.exists('/fastq/raw_data/'+s1+'_1_trimmed.fq.gz'):
        # Second remove adaptor sequences
        cmd = ['cutadapt',
                '-a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC',
                '-A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT',
                '-j 12',
                '-m 15',
                '--max-n=5',
                '-o /fastq/raw_data/'+s1+'_1_trimmed.fq.gz',
                '-p /fastq/raw_data/'+s1+'_2_trimmed.fq.gz',
                '/fastq/raw_data/'+s1+'_1.fq.gz',
                '/fastq/raw_data/'+s1+'_2.fq.gz']
        print('Running = '+' '.join(cmd))
        cutadaptProc = Popen(' '.join(cmd), shell=True,stdout=PIPE,stderr=PIPE)
        cutadaptProc.communicate()

    if not os.path.exists('/fastq/output/'+s1+'Aligned.out.counts'):
        # Run STAR on each experiment
        cmd = ['STAR',
               '--genomeDir /index',
               '--runMode alignReads --runThreadN 12',
               '--outSAMstrandField intronMotif',
               '--outFilterIntronMotifs RemoveNoncanonicalUnannotated',
               '--readFilesIn /fastq/raw_data/'+s1+'_1_trimmed.fq.gz',
                             '/fastq/raw_data/'+s1+'_2_trimmed.fq.gz',
               '--readFilesCommand zcat',
               '--outFileNamePrefix /fastq/output/'+s1]
        print('Running = '+' '.join(cmd))
        starProc = Popen(' '.join(cmd), shell=True,stdout=PIPE,stderr=PIPE)
        starProc.communicate()

        # Run htseq on each experiment
        cmd = ['htseq-count',
               '--quiet',
               '--format=sam',
               '--mode=intersection-strict',
               '--order=pos',
               '--stranded=no',
               '/fastq/output/'+s1+'Aligned.out.sam',
               '/GRCh38.p12/gencode.v31.primary_assembly.annotation.gtf',
               '>',
               '/fastq/output/'+s1+'Aligned.out.counts']
        print('Running = '+' '.join(cmd))
        htseqProc = Popen(' '.join(cmd), shell=True,stdout=PIPE,stderr=PIPE)
        htseqProc.communicate()

        # gzip Aligned.out.sam file
        print('Gzipping '+s1)
        gzipProc = Popen('pigz /fastq/output/'+s1+'Aligned.out.sam', shell=True,stdout=PIPE) #,stderr=errOut)
        gzipProc.communicate()

        # Done
        print('Done with '+s1+'.\n')

    print('Reading'+s1+'...')
    with open('/fastq/output/'+s1+'Aligned.out.counts','r') as inFile:
        while 1:
            line = inFile.readline()
            if not line:
                break
            splitUp = line.strip().split('\t')
            if not splitUp[0] in gexpMatrix:
                gexpMatrix[splitUp[0]] = {}
            gexpMatrix[splitUp[0]][s1] = splitUp[1]

    with open('/fastq/output/'+s1+'Log.final.out','r') as inFile:
        inLines = [i.strip() for i in inFile.readlines() if i]
        for inLine in inLines:
            if not inLine.find('\t')==-1:
                splitUp = [i.strip() for i in inLine.split('\t')]
                key1 = splitUp[0].rstrip(' |').replace(',',' ')
                if not key1 in output:
                    output[key1] = {}
                    outOrder.append(key1)
                output[key1][s1] = splitUp[1]

print('Writing full count matrix...')
# Dump out matrix of coverages
with open('/fastq/output/gexp_counts.csv','w') as outFile:
    outFile.write('transcript_id,'+','.join(samples)+'\n')
    outFile.write('\n'.join([transcript+','+','.join([gexpMatrix[transcript][sample] for sample in samples]) for transcript in gexpMatrix]))

# Write out qc file
print('Writing QC file...')
writeMe = ['Sample,'+','.join(outOrder)]
for sample in samples:
    writeMe += [sample+','+','.join([output[i][sample] for i in outOrder])]

with open('/fastq/output/qc_alignment.csv','w') as outFile:
    outFile.write('\n'.join(writeMe))

print('Done.')
