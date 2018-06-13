#python rundeepsea.py infilename outdir

from subprocess import *
from tempfile import mkdtemp
import sys
import os 


infilename=sys.argv[1]
outdir=sys.argv[2]+'/'
cpoutdir=True
check_call(['mkdir','-p',outdir])

feature_ind=sys.argv[3]

if infilename.endswith('fasta'):
    try:
        tempdir = mkdtemp()
        check_call(['cp',infilename,tempdir+'/infile.fasta'])
        print "Successfully copied input to working directory."
        check_call(["grep '>'  "+tempdir+'/infile.fasta '+">"+ tempdir+'/infile.fasta.name'],shell=True)
        try:
            check_call(["python 1_fasta2input.nomut.py  "+tempdir+"/infile.fasta"],shell=True) #(output tempdir/infile.fasta.ref.h5)
        except:
            raise Exception("Fasta format error.")
        #
        check_call(["python A_saturated_mutagenesis.py  "+tempdir+"/infile.fasta"],shell=True) #(output tempdir/infile.fasta.saturatedmut.h5)

        print "Successfully converted to input format"
        check_call(["luajit 2_DeepSEA.lua -test_file_h5 "+tempdir+"/infile.fasta.ref.h5"],shell=True)
        #
        check_call(["luajit 2_DeepSEA.lua -test_file_h5 "+tempdir+"/infile.fasta.saturatedmut.h5"],shell=True)
        print "Finished running DeepSEA. Now prepare output files..."
        #
        check_call(["python B_MutaVis.py "+tempdir+"/infile.fasta "+tempdir+"/infile.fasta.ref.h5.pred.h5 "+tempdir+"/infile.fasta.saturatedmut.h5.pred.h5 "+feature_ind],shell=True)  #(output tempdir/infile.fasta.out)
        if cpoutdir:
            check_call(["cp", tempdir+"/infile.fasta.preview.png" ,outdir+'/preview.png' ])
            check_call(["cp", tempdir+"/infile.fasta.vis.png" ,outdir+'/vis.png' ])
            check_call(["cp", tempdir+"/infile.fasta.colorbar.png" ,outdir+'/colorbar.png' ])
            check_call(["cp", tempdir+"/infile.fasta.csv" ,outdir+'/log2foldchange_profile.csv' ])
        print "Finished creating output file. Now clean up..."
        print "Everything done."
    finally:
        call(['rm',tempdir,'-r'])
elif infilename.endswith('bed'):
    try:
        tempdir = mkdtemp()
        print tempdir
        check_call(['cp',infilename,tempdir+'/infile.bed'])
        print "Successfully copied input to working directory "+ tempdir
        try:
            check_call(["Rscript  0_coor2fasta.R "+tempdir+'/infile.bed '],shell=True)
        except:
            raise Exception('Bed format error')
        check_call(["python 1_fasta2input.nomut.py  "+tempdir+"/infile.bed.wt1000.fasta"],shell=True) #(output tempdir/infile.bed.wt1000.fasta.ref.h5)
        check_call(["grep '>' "+tempdir+"/infile.bed.wt1000.fasta |sed 's/>__//g'|sed 's/_/\t/g' > "+tempdir+"/infile.bed.wt1000.fasta.ref.bed"],shell=True)
        #
        check_call(["python A_saturated_mutagenesis.py  "+tempdir+"/infile.bed.wt1000.fasta"],shell=True) #(output tempdir/infile.fasta.saturatedmut.h5)

        print "Successfully converted to input format"
        check_call(["luajit 2_DeepSEA.lua -test_file_h5 "+tempdir+"/infile.bed.wt1000.fasta.ref.h5"],shell=True)
        #
        check_call(["luajit 2_DeepSEA.lua -test_file_h5 "+tempdir+"/infile.bed.wt1000.fasta.saturatedmut.h5"],shell=True)
        print "Finished running DeepSEA. Now prepare output files..."
        #
        check_call(["python B_MutaVis.py "+tempdir+"/infile.bed.wt1000.fasta "+tempdir+"/infile.bed.wt1000.fasta.ref.h5.pred.h5 "+tempdir+"/infile.bed.wt1000.fasta.saturatedmut.h5.pred.h5 "+feature_ind],shell=True)  #(output tempdir/infile.fasta.out)
        if cpoutdir:
            check_call(["cp", tempdir+"/infile.bed.wt1000.fasta.preview.png" ,outdir+'/preview.png' ])
            check_call(["cp", tempdir+"/infile.bed.wt1000.fasta.vis.png" ,outdir+'/vis.png' ])
            check_call(["cp", tempdir+"/infile.bed.wt1000.fasta.colorbar.png" ,outdir+'/colorbar.png' ])
            check_call(["cp", tempdir+"/infile.bed.wt1000.fasta.csv" ,outdir+'/log2foldchange_profile.csv' ])

        print "Finished creating output file. Now clean up..."
        print "Everything done."
    finally:
        #pass
        call(['rm',tempdir,'-r'])
elif infilename.endswith('vcf'):
    try:
        tempdir = mkdtemp()
        check_call(['cp',infilename,tempdir+'/infile.vcf'])
        print "Successfully copied input to working directory " + tempdir
        try:
            check_call(["Rscript  0_coor2fasta.R "+tempdir+'/infile.vcf '],shell=True)
        except:
            raise Exception('Vcf format error.')
        #
        check_call(["python A_saturated_mutagenesis.py  "+tempdir+"/infile.vcf.wt1100.fasta 1100"],shell=True) #(output tempdir/infile.fasta.saturatedmut.h5)

        check_call(["python 1_fasta2input.py  "+tempdir+"/infile.vcf.wt1100.fasta 1100"],shell=True) #(output tempdir/infile.vcf.wt1100.fasta.ref.h5)
        print "Successfully converted to input format"
        check_call(["luajit 2_DeepSEA.lua -test_file_h5 "+tempdir+"/infile.vcf.mut1100.fasta.ref.h5"],shell=True)
        #
        check_call(["luajit 2_DeepSEA.lua -test_file_h5 "+tempdir+"/infile.vcf.wt1100.fasta.saturatedmut.h5"],shell=True)

        print "Finished running DeepSEA. Now prepare output files..."
        #
        check_call(["python B_MutaVis.py "+tempdir+"/infile.vcf.mut1100.fasta "+tempdir+"/infile.vcf.mut1100.fasta.ref.h5.pred.h5 "+tempdir+"/infile.vcf.wt1100.fasta.saturatedmut.h5.pred.h5 "+feature_ind],shell=True)  #(output tempdir/infile.fasta.out)
        if cpoutdir:
            check_call(["cp", tempdir+"/infile.vcf.mut1100.fasta.preview.png" ,outdir+'/preview.png' ])
            check_call(["cp", tempdir+"/infile.vcf.mut1100.fasta.vis.png" ,outdir+'/vis.png' ])
            check_call(["cp", tempdir+"/infile.vcf.mut1100.fasta.colorbar.png" ,outdir+'/colorbar.png' ])
            check_call(["cp", tempdir+"/infile.vcf.mut1100.fasta.csv" ,outdir+'/log2foldchange_profile.csv' ])

        print "Finished creating output file. Now clean up..."
        print "Everything done."
    finally:
        pass
        #call(['rm',tempdir,'-r'])




