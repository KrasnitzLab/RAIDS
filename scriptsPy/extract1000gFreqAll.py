import vcf
from datetime import date
import os
import sys
import numpy as np

def main():
    fileVCF = sys.argv[1]
    fileBASE = sys.argv[2]
    freqCutOff = -1
    today = date.today()
    displayLOG(fileBASE, today, fileVCF, freqCutOff)
    extractFreq(fileBASE, fileVCF, freqCutOff, today)

def displayLOG(fileBASE, today, fileVCF, freqCutOff):
    try:
        fileLOG = "log/" + fileBASE + "_" + str(today) + ".log"
        FLOG = open(fileLOG, "w")
        FLOG.write("fileVCF " + fileVCF + "\n")
        FLOG.write("fileBASE " + fileBASE + "\n")
        FLOG.write("FreqCutOff " + str(freqCutOff) + "\n")
        FLOG.close()
    except (OSError, IOError) as e:
        sys.stderr.write( 'Problem with the file: '+ fileLOG + '\n')
        sys.stderr.write(e + "\n")
        sys.exit(1)

def extractFreq(fileBASE, fileVCF, freqCutOff, today):
    try:
        fileName = fileBASE + "_" + str(today) + "_f_All" + ".txt"

        fileOUT = os.path.join("mat1000gFAll", fileName)
        sep = ";"
        vcf_reader = vcf.Reader(filename=fileVCF)
        OUT = open(fileOUT, "w")

        for record in vcf_reader:
            # At least the frequency in one super population
            # higher or equal to freqCutOff
            if (record.INFO['EAS_AF'][0] >= freqCutOff or
                record.INFO['EUR_AF'][0] >= freqCutOff or
                record.INFO['AFR_AF'][0] >= freqCutOff or
                record.INFO['AMR_AF'][0] >= freqCutOff or
                record.INFO['SAS_AF'][0] >= freqCutOff):

                OUT.write( record.CHROM +"\t" + str( record.start ) + "\t" + record.REF)

                flag = 0
                for g in record.ALT:
                    if flag == 1:
                        OUT.write(sep)
                    OUT.write("\t" + str(g))
                    flag = 1
                OUT.write("\t" + str(record.INFO['AF'][0]))
                OUT.write("\t" + str(record.INFO['EAS_AF'][0]))
                OUT.write("\t" + str(record.INFO['EUR_AF'][0]))
                OUT.write("\t" + str(record.INFO['AFR_AF'][0]))
                OUT.write("\t" + str(record.INFO['AMR_AF'][0]))
                OUT.write("\t" + str(record.INFO['SAS_AF'][0]))
                OUT.write("\n")

        OUT.close


    except (OSError, IOError) as e:
        sys.stderr.write( 'Problem with the file: '+ fileBASE + '_' + fileVCF + '\n')
        sys.stderr.write(e + "\n")
        sys.exit(1)

if __name__ == "__main__":
    main()

