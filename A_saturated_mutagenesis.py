# in silico mutagenesis
from Bio import SeqIO
import numpy as np
import sys
import h5py



BASENUM = 4
mutwindow = 1000

fasta_sequences = SeqIO.parse(open(sys.argv[1]), 'fasta')
np.random.seed(1)

seqs = [str(fasta.seq) for fasta in fasta_sequences]

useAlternativeAllele = False



mutwindow = 1000
BASENUM = 4

if len(sys.argv)>2:
    inputwindow = int(sys.argv[2])
    mutpos = inputwindow / 2 - 1
    oris = []
    muts = []
    chrs = []
    poss = []
    annos = []
    names = []
    fasta_sequences = SeqIO.parse(open(sys.argv[1]), 'fasta')
    for fasta in fasta_sequences:
        anno = fasta.name.split('_')
        annos.append(fasta.name)
        oris.append(anno[0])
        muts.append(anno[1])
        chrs.append(anno[2])
        poss.append(int(anno[3]) + mutpos)
        if len(anno) > 5:
            names.append('_'.join(anno[5:]))

    oris = np.asarray(oris)
    muts = np.asarray(muts)
    chrs = np.asarray(chrs)
    seqs = np.asarray(seqs)
    poss = np.asarray(poss)
    annos = np.asarray(annos)
    names = np.asarray(names)

    seqsmut = []
    inds = []
    for i in range(len(seqs)) or seqs[i][mutpos:(mutpos+len(oris[i]))]!=oris[i]:
        if type(oris[i]) is not np.string_:
            continue
        else:
            inds.append(i)
            seqsmut.append(
                seqs[i][:mutpos] + muts[i] + seqs[i][(mutpos + len(oris[i])):])
            seqs[i] = seqs[i][:mutpos] + oris[i] + seqs[i][(mutpos + 1):]
    seqs=seqsmut
    #write mut sequence
    wfile=open(sys.argv[1].replace('wt','mut'),'w')
    wfile.write('>\n')
    wfile.write(seqs[0][(( (len(seqs[0]) - 1000) / 2)):
                 ( len(seqs[0]) - (len(seqs[0]) - 1000) / 2)])
    wfile.close()



# write mat / h5 file

filename = sys.argv[1]
offset_values = None

seqmut = np.zeros((len(seqs) * (BASENUM - 1) * mutwindow, 4, 1000), np.bool_)
seqori = np.zeros((len(seqs), 4, 1000), np.bool_)

mydict = {'A': np.asarray([1, 0, 0, 0]), 'G': np.asarray([0, 1, 0, 0]),
          'C': np.asarray([0, 0, 1, 0]), 'T': np.asarray([0, 0, 0, 1]),
          'N': np.asarray([0, 0, 0, 0]), 'H': np.asarray([0, 0, 0, 0]),
          'a': np.asarray([1, 0, 0, 0]), 'g': np.asarray([0, 1, 0, 0]),
          'c': np.asarray([0, 0, 1, 0]), 't': np.asarray([0, 0, 0, 1]),
          'n': np.asarray([0, 0, 0, 0])}

n = 0
for line in seqs:
    if len(line) < 1000:
        raise Exception("Each fasta sequence has to be at least 1000bp.")
    cline = line[( (len(line) - 1000) / 2):  ( ( len(line) + 1000) / 2)]
    for c, i in zip(cline, range(len(cline))):
        seqori[n, :, i] = mydict[c]

    # saturated mutagenesis
    seqmut[n * (BASENUM - 1) * mutwindow : (n + 1) * (BASENUM - 1) * mutwindow, :, :] = \
        np.tile(seqori[n, :, :],
                ((BASENUM - 1) * mutwindow, 1, 1))
    for x in range(mutwindow):
        for y in range(BASENUM - 1):
            mutagenesis_offset = int((1000 - mutwindow) / 2)
            seqmut[n * (BASENUM - 1) * mutwindow + y + x * (BASENUM - 1), :, mutagenesis_offset + x] = \
                np.roll(seqmut[n * (BASENUM - 1) * mutwindow + y + x * (BASENUM - 1), :, mutagenesis_offset + x], y + 1)

    n = n + 1

dataflip = seqmut[:, ::-1, ::-1]
seqmut = np.concatenate([seqmut, dataflip], axis=0)
seqmut = seqmut.astype(np.uint8)
f = h5py.File(filename + '.saturatedmut.h5', 'w')
f.create_dataset('testxdata', data=seqmut, compression="gzip")
f.close()
