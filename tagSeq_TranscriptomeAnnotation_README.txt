# Transcriptome Annotation
# Created by Misha Matz (matz@utexas.edu), modified by Michael Studivan (mstudiva@fau.edu) 
# for use in generating transcriptome annotation files for Montastraea cavernosa

#------------------------------
# BEFORE STARTING, replace, in this whole file:
#	- mstudiva@fau.edu by your actual email;
#	- mstudiva with your KoKo user name.

# The idea is to copy the chunks separated by empty lines below and paste them into your cluster 
# terminal window consecutively. 

# The lines beginning with hash marks (#) are explanations and additional instructions - 
# please make sure to read them before copy-pasting. 

#------------------------------
# To install Bioperl in your bin directory, please follow these instructions:
cd bin
module load python/anaconda
conda create -y -n bioperl perl-bioperl

# getting scripts

cd ~/bin
git clone https://github.com/z0on/annotatingTranscriptomes.git
mv annotatingTranscriptomes/* .
rm -rf annotatingTranscriptomes
rm launcher_creator.py

# creating annotation dir
mkdir annotate
cd annotate

# getting transcritpome to play with
wget http://meyerlab:coral@files.cgrb.oregonstate.edu/Meyer_Lab/transcriptomes/Mcav/Mcav_transcriptome_v1.fasta.gz
gunzip Mcav_transcriptome_v1.fasta.gz
mv Mcav_transcriptome_v1.fasta mcav.fasta

# statistics: (takes a long time to run for some reason but works):
source activate bioperl
seq_stats.pl mcav.fasta
source deactivate bioperl

# getting uniprot_swissprot KB database
wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz

# getting annotations (this file is over 3G, will take a while)
echo 'wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/idmapping_selected.tab.gz '> getz
launcher_creator.py -j getz -n getz -t 2:00:00 -q shortq7 -e mstudiva@fau.edu
sbatch getz.slurm

# unzipping
gunzip uniprot_sprot.fasta.gz &
gunzip idmapping_selected.tab.gz &

# indexing the fasta database
module load blast
echo "makeblastdb -in uniprot_sprot.fasta -dbtype prot" >mdb
launcher_creator.py -j mdb -n mdb -q shortq7 -t 2:00:00 -e mstudiva@fau.edu
sbatch mdb.slurm

# splitting the transcriptome into 190 chunks
splitFasta.pl mcav.fasta 190

# blasting all 190 chunks to uniprot in parallel, 4 cores per chunk
module load blast
ls subset* | perl -pe 's/^(\S+)$/blastx -query $1 -db uniprot_sprot\.fasta -evalue 0\.0001 -num_threads 4 -num_descriptions 5 -num_alignments 5 -out $1.br/'>bl
launcher_creator.py -j bl -n blast -t 2:00:00 -q shortq7 -e mstudiva@fau.edu
sbatch blast.slurm

# watching progress:
grep "Query= " subset*.br | wc -l

# combining all blast results
cat subset*br > myblast.br
rm subset*

# for trinity-assembled transcriptomes: annotating with "isogroup" (=component)
grep ">" mcav.fasta | perl -pe 's/>comp(\d+)(\S+)\s.+/comp$1$2\tisogroup$1/' >mcav_seq2iso.tab
cat mcav.fasta | perl -pe 's/>comp(\d+)(\S+).+/>comp$1$2 gene=isogroup$1/' >mcav_iso.fasta

#-------------------------

# extracting gene names (per isogroup):
getGeneNameFromUniProtKB.pl blast=myblast.br prefix=mcav fastaQuery=mcav_iso.fasta

# extracting GO annotations (per isogroup)
echo "getGOfromUniProtKB.pl blast=myblast.br prefix=mcav fastaQuery=mcav_iso.fasta" >getgo
launcher_creator.py -j getgo -n getgo -l gg -q shortq7 -t 2:00:00 -e mstudiva@fau.edu
sbatch gg

# extracting coding sequences and corresponding protein translations:
source activate bioperl
echo "perl ~/bin/CDS_extractor_v2.pl mcav_iso.fasta myblast.br allhits bridgegaps" >cds
launcher_creator.py -j cds -n cds -l cddd -t 2:00:00 -q shortq7 -e mstudiva@fau.edu
sbatch cddd
source deactivate bioperl

# calculating contiguity:
contiguity.pl hits=mcav_iso_hits.tab threshold=0.75
# 0.23

# core gene set form korflab: to characterize representation of genes:
wget http://korflab.ucdavis.edu/Datasets/genome_completeness/core/248.prots.fa.gz
gunzip 248.prots.fa.gz

makeblastdb -in mcav.fasta -dbtype nucl
echo 'tblastn -query 248.prots.fa -db mcav.fasta -evalue 1e-10 -outfmt "6 qseqid sseqid evalue bitscore qcovs" -max_target_seqs 1 -num_threads 12 >mcav_248.brtab' >bl248
launcher_creator.py -j bl248 -n bl -l blj -q shortq7 -t 02:00:00 -e mstudiva@fau.edu
sbatch blj
# calculating fraction of represented KOGs:
cat mcav_248.brtab | perl -pe 's/.+(KOG\d+)\s.+/$1/' | uniq | wc -l | awk '{print $1/248}'
# 0.991935

#------------------------------
# KOG annotation
# scp your *_PRO.fas file to laptop, submit it to
http://weizhong-lab.ucsd.edu/metagenomic-analysis/server/kog/
cd /Users/Mike/Documents/Grad_School/Dissertation/Data/RNA_Seq/Transcriptome
scp mstudiva@koko-login.fau.edu:/scratch/02475/mstudiva/annotate/*_PRO.fas .

# copy link to job ID status and output file, paste it below instead of current link:
# check status: go on web to http://weizhong-lab.ucsd.edu/metagenomic-analysis/result/?jobid=95827620170616072315007524
# once it is done, download results:
wget http://weizhong-lab.ucsd.edu/metagenomic-analysis/output/95827620170616072315007524/output.zip

unzip output.zip
mv output.2 mcav.kog.tab

# generates iso2kogClass and iso2kogDef (another kind of gene names)
getKOGs.pl fastaQuery=mcav_iso.fasta prefix=mcav kogMatch=mcav.kog.tab 

# removing "multiple classes" annotation, renaming comp to isogroup
grep -v "Multiple classes" mcav_iso2kogClass.tab | perl -pe 's/^comp/isogroup/' > mcav_iso2kogClassNR.tab

#------------------------------
# KEGG annotations:

# selecting the longest contig per isogroup:
fasta2SBH.pl mcav_iso.fasta >mcav_4kegg.fasta

# scp mcav_4kegg.fasta to your laptop
cd /Users/Mike/Documents/Grad_School/Dissertation/Data/RNA_Seq/Transcriptome
scp mstudiva@ls5.tacc.utexas.edu:/scratch/02475/mstudiva/annotate/mcav_4kegg.fasta .
# use web browser to submit mcav_4kegg.fasta file to KEGG's KAAS server ( http://www.genome.jp/kegg/kaas/ )
# select SBH algorithm, upload nucleotide query
# Once it is done, download the 'text' output from KAAS, name it query.ko (default)

wget http://www.genome.jp/tools/kaas/files/dl/1497646795/query.ko

# selecting only the lines with non-missing annotation:
cat query.ko | awk '{if ($2!="") print }' > mcav_iso2kegg.tab

# the KEGG mapping result can be explored for completeness of transcriptome in terms of genes found,
# use 'html' output link from KAAS result page, see how many proteins you have for conserved complexes and pathways,
# such as ribosome, spliceosome, proteasome etc

#------------------------------
# copy all files to laptop
cd /Users/Mike/Documents/Grad_School/Dissertation/Data/RNA_Seq/Transcriptome
scp mstudiva@ls5.tacc.utexas.edu:/scratch/02475/mstudiva/annotate/* .

