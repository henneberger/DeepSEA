#Example usage: python 1_fasta2input.py infile.fasta 
#Convert sequences to hdf5 format. Used for processing bed or fasta format input.

from Bio import SeqIO
import numpy as np
import sys
import h5py

    
inputwindow=1000
mutpos=inputwindow/2-1
fasta_sequences = SeqIO.parse(open(sys.argv[1]),'fasta')
np.random.seed(1)




#write mat / h5 file

def writeh5(seqs,filename,offset_values=None):
    seqsnp=np.zeros((len(seqs),4,1000),np.bool_)

    mydict={'A':np.asarray([1,0,0,0]),'G':np.asarray([0,1,0,0]),'C':np.asarray([0,0,1,0]),'T':np.asarray([0,0,0,1]),'N':np.asarray([0,0,0,0]),'H':np.asarray([0,0,0,0]),'a':np.asarray([1,0,0,0]),'g':np.asarray([0,1,0,0]),'c':np.asarray([0,0,1,0]),'t':np.asarray([0,0,0,1]),'n':np.asarray([0,0,0,0])}
    n=0
    offset_values = np.zeros(len(seqs))
    for line,o in zip(seqs,offset_values):
        if len(line)<1000:
            print len(line)
            continue
            raise Exception("Each fasta sequence has to be at least 1000bp.")
        #if the user specified region/sequence is longer than 1000bp, use the center 1000bp
        cline = line[((int(o) + (len(line)-1000)/2)):(int(o)+(len(line)+1000)/2)]
        for c,i in zip(cline,range(len(cline))):
            seqsnp[n,:,i]=mydict[c]
        n=n+1
    seqsnp=seqsnp[:n,:,:]
    dataflip=seqsnp[:,::-1,::-1];
    seqsnp=np.concatenate([seqsnp, dataflip],axis=0)
    seqsnp = seqsnp.astype(np.uint8)
    f=h5py.File(filename,'w')
    f.create_dataset('testxdata', data= seqsnp,compression="gzip")
    f.close()

    
    
seqs=[str(fasta.seq) for fasta in fasta_sequences]


if True:
    writeh5(seqs,sys.argv[1]+'.ref.h5')
    
        
    
