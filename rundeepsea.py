#python rundeepsea.py infilename outdir
#The output files will be in outdir.

from subprocess import *
from tempfile import mkdtemp
import sys
import os 


infilename=sys.argv[1]
outdir=sys.argv[2]+'/'
cpoutdir=True
check_call(['mkdir','-p',outdir])

#for fasta input 
if infilename.endswith('fasta'):
    try:
        tempdir = mkdtemp()
        check_call(['cp',infilename,tempdir+'/infile.fasta'])
        print "Successfully copied input to working directory."
        check_call(["grep '>'  "+tempdir+'/infile.fasta '+">"+ tempdir+'/infile.fasta.name'],shell=True)
        try:
            check_call(["python 1_fasta2input.nomut.py  "+tempdir+"/infile.fasta"],shell=True)
        except:
            raise Exception("Fasta format error.")
        print "Successfully converted to input format"
        check_call(["luajit 2_DeepSEA.lua -test_file_h5 "+tempdir+"/infile.fasta.ref.h5"],shell=True)
        print "Finished running DeepSEA. Now prepare output files..."
        check_call(["python 3_h5ToOutput.py "+tempdir+"/infile.fasta "+tempdir+"/infile.fasta.ref.h5.pred.h5"],shell=True) 
        if cpoutdir:
            check_call(["cp", tempdir+"/infile.fasta.out" ,outdir ])
    finally:
        print "Finished creating output file. Now clean up..."
        call(['rm',tempdir,'-r'])
        print "Everything done."
#for bed input 
elif infilename.endswith('bed'):
    try:
        tempdir = mkdtemp()
        check_call(['cp',infilename,tempdir+'/infile.bed'])
        print "Successfully copied input to working directory "+ tempdir
        try:
            check_call(["Rscript  0_coor2fasta.R "+tempdir+'/infile.bed '],shell=True)
        except:
            raise Exception('Bed format error')
        check_call(["python 1_fasta2input.nomut.py  "+tempdir+"/infile.bed.wt1000.fasta"],shell=True) 
        check_call(["grep '>' "+tempdir+"/infile.bed.wt1000.fasta |sed 's/>__//g'|sed 's/_/\t/g' > "+tempdir+"/infile.bed.wt1000.fasta.ref.bed"],shell=True)
        print "Successfully converted to input format"
        check_call(["luajit 2_DeepSEA.lua -test_file_h5 "+tempdir+"/infile.bed.wt1000.fasta.ref.h5"],shell=True)
        print "Finished running DeepSEA. Now prepare output files..."
        check_call(["python 3_h5ToOutput.py "+tempdir+"/infile.bed "+tempdir+"/infile.bed.wt1000.fasta.ref.h5.pred.h5"],shell=True)  #(output infile.bed.out)
        if cpoutdir:
            check_call(["cp", tempdir+"/infile.bed.out" ,outdir ])
    finally:
        print "Finished creating output file. Now clean up..."
        call(['rm',tempdir,'-r'])
        print "Everything done."
#for vcf input 
elif infilename.endswith('vcf'):
    try:
        tempdir = mkdtemp()
        check_call(['cp',infilename,tempdir+'/infile.vcf'])
        print "Successfully copied input to working directory " + tempdir 
        try:
            check_call(["Rscript  0_coor2fasta.R "+tempdir+'/infile.vcf '],shell=True)
        except:
            raise Exception('Vcf format error.')
        #retrieve 1100bp instead of 1000bp for supporting deletion variants (<100bp) 
        check_call(["python 1_fasta2input.py  "+tempdir+"/infile.vcf.wt1100.fasta 1100"],shell=True) 
        print "Successfully converted to input format"
        check_call(["luajit 2_DeepSEA.lua -test_file_h5 "+tempdir+"/infile.vcf.wt1100.fasta.ref.h5"],shell=True)
        check_call(["luajit 2_DeepSEA.lua -test_file_h5 "+tempdir+"/infile.vcf.mut1100.fasta.ref.h5"],shell=True)
        print "Finished running DeepSEA. Now prepare output files..."
        check_call(["python 3_h5ToOutput.py "+tempdir+"/infile.vcf "+tempdir+"/infile.vcf.wt1100.fasta.ref.h5.pred.h5  "+tempdir+"/infile.vcf.mut1100.fasta.ref.h5.pred.h5"],shell=True) 
        if cpoutdir:
            evoinfo_exist=os.path.isfile('./resources/Gerp/gerp_elements.tsv.gz') and os.path.isfile('./resources/Gerp/gerp_scores.tsv.gz') and  os.path.isfile('./resources/phastCons/primates_nohuman.tsv.gz') and os.path.isfile('./resources/phyloP/primates_nohuman.tsv.gz')
            check_call(["cp", tempdir+"/infile.vcf.out.alt" ,outdir ])
            check_call(["cp", tempdir+"/infile.vcf.out.ref" ,outdir ])
            check_call(["cp", tempdir+"/infile.vcf.out.diff" ,outdir ])
            check_call(["cp", tempdir+"/infile.vcf.out.logfoldchange" ,outdir ])
            check_call(["cp", tempdir+"/infile.vcf.out.evalue" ,outdir ])
            if evoinfo_exist:
                check_call(["cp", tempdir+"/infile.vcf.wt1100.fasta.ref.vcf.evoall" ,outdir ])
                check_call(["cp", tempdir+"/infile.vcf.wt1100.fasta.ref.vcf.evo.evalues" ,outdir ])
                check_call(["cp", tempdir+"/infile.vcf.out.snpclass" ,outdir ])
                check_call(["cp", tempdir+"/infile.vcf.out.funsig" ,outdir ])
                check_call(["cp", tempdir+"/infile.vcf.out.summary" ,outdir ])
    finally:
        print "Finished creating output file. Now clean up..."
        call(['rm',tempdir,'-r'])
        print "Everything done."
