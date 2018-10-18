# Transcriptome Annotation, version Dec 31, 2017
# Created by Misha Matz (matz@utexas.edu), modified by Michael Studivan (mstudiva@fau.edu) 
# for use in generating transcriptome annotation files for Montastraea cavernosa
# also includes the concatention of M. cav and Symbiodinium Clade C transcriptomes

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

# M. cavernosa transcriptome, v1 (2014)
# wget http://meyerlab:coral@files.cgrb.oregonstate.edu/Meyer_Lab/transcriptomes/Mcav/Mcav_transcriptome_v1.fasta.gz
# gunzip Mcav_transcriptome_v1.fasta.gz
# mv Mcav_transcriptome_v1.fasta mcav.fasta

# updated (July 2018) M. cavernosa genome with transcriptome
wget https://www.dropbox.com/s/0inwmljv6ti643o/Mcavernosa_genome.tgz
tar -xzf Mcavernosa_genome.tgz
cd Mcav_genome/Mcavernosa_annotation/
cp Mcavernosa.maker.transcripts.fasta ../..
mv Mcavernosa.maker.transcripts.fasta mcav.fasta

# Symbiodinium Clade C transcriptome (mid 2017)
# wget ftp://marchanon:anon@rc-ns-ftp.its.unc.edu/CladeC_Symbiodinium_transcriptome.zip
# unzip CladeC_Symbiodinium_transcriptome.zip
# mv CladeC_Symbiodinium_transcriptome/* .
# rm -rf CladeC_Symbiodinium_transcriptome/ CladeC_Symbiodinium_transcriptome.zip _MACOSX/

# Cladocopium sp. (formerly Symbiodinium Clade C) transcriptome as of November 2017
wget http://sites.bu.edu/davieslab/files/2017/11/CladeC_Symbiodinium_transcriptome.zip
unzip CladeC_Symbiodinium_transcriptome.zip
mv CladeC_Symbiodinium_transcriptome/* .
rm -rf CladeC_Symbiodinium_transcriptome/ CladeC_Symbiodinium_transcriptome.zip __MACOSX/
mv davies_cladeC_feb.fasta symC.fasta
mv davies_cladeC_iso2go.tab symC_iso2go.tab
mv davies_cladeC_iso2gene.tab symC_iso2gene.tab
mv davies_cladeC_seq2iso.tab symC_seq2iso.tab

# use the stream editor to find and replace all instances of "comp" and "isogroup" with "sym" in the symbiont transcriptome files
sed -i 's/comp/sym/g' symC_seq2iso.tab
sed -i 's/comp/sym/g' symC.fasta
sed -i 's/isogroup/sym/g' symC_iso2gene.tab
sed -i 's/isogroup/sym/g' symC_iso2go.tab
sed -i 's/isogroup/sym/g' symC_seq2iso.tab

# concatenate the host and symbiont transcriptomes into a holobiont transcriptome
cat symC.fasta mcav.fasta > mcav_holobiont.fasta

# transcriptome statistics
source activate bioperl
seq_stats.pl mcav_holobiont.fasta
source deactivate bioperl

# mcav_holobiont.fasta
# -------------------------
# 90980 sequences.
# 1465 average length.
# 43960 maximum length.
# 75 minimum length.
# N50 = 1798
# 133.3 Mb altogether (133273996 bp).
# 0 ambiguous Mb. (1637 bp, 0%)
# 0 Mb of Ns. (1637 bp, 0%)
# -------------------------

# getting uniprot_swissprot KB database
wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz

# getting annotations (this file is large, may take a while)
echo 'wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/idmapping_selected.tab.gz '> getz
launcher_creator.py -j getz -n getz -t 2:00:00 -q shortq7 -e mstudiva@fau.edu
sbatch getz.slurm

# if the US mirror is down, use the line below, then run the getz script as normal
echo 'wget ftp://ftp.ebi.ac.uk/pub/databases/uniprot/current_release/knowledgebase/idmapping/idmapping_selected.tab.gz '> getz

# unzipping
gunzip uniprot_sprot.fasta.gz &
gunzip idmapping_selected.tab.gz &

# indexing the fasta database
module load blast
echo "makeblastdb -in uniprot_sprot.fasta -dbtype prot" >mdb
launcher_creator.py -j mdb -n mdb -q shortq7 -t 2:00:00 -e mstudiva@fau.edu
sbatch mdb.slurm

# splitting the transcriptome into 190 chunks
splitFasta.pl mcav_holobiont.fasta 190

# blasting all 190 chunks to uniprot in parallel, 4 cores per chunk
module load blast
ls subset* | perl -pe 's/^(\S+)$/blastx -query $1 -db uniprot_sprot\.fasta -evalue 0\.0001 -num_threads 4 -num_descriptions 5 -num_alignments 5 -out $1.br/'>bl
launcher_creator.py -j bl -n blast -t 2:00:00 -q shortq7 -N 8 -e mstudiva@fau.edu
sbatch blast.slurm

# watching progress:
grep "Query= " subset*.br | wc -l

# combining all blast results
cat subset*br > myblast.br
rm subset*

# for trinity-assembled transcriptomes: annotating with "sym" or "Mcavernosa" depending on if component is from symbiont or host (=component)
grep ">" mcav_holobiont.fasta | perl -pe 's/>sym(\d+)(\S+)\s.+/sym$1$2\tsym$1/' | perl -pe 's/>Mcavernosa(\d+)(\S+)\s.+/Mcavernosa$1$2\tMcavernosa$1/'>mcav_holobiont_seq2iso.tab
cat mcav_holobiont.fasta | perl -pe 's/>sym(\d+)(\S+).+/>sym$1$2 gene=sym$1/' | perl -pe 's/>Mcavernosa(\d+)(\S+).+/>Mcavernosa$1$2 gene=Mcavernosa$1/'>mcav_holobiont_iso.fasta
#-------------------------

# extracting gene names (per isogroup):
getGeneNameFromUniProtKB.pl blast=myblast.br prefix=mcav_holobiont fastaQuery=mcav_holobiont_iso.fasta

# extracting GO annotations (per isogroup)
echo "getGOfromUniProtKB.pl blast=myblast.br prefix=mcav_holobiont fastaQuery=mcav_holobiont_iso.fasta" >getgo
launcher_creator.py -j getgo -n getgo -l gg -q shortq7 -t 2:00:00 -e mstudiva@fau.edu
sbatch gg

# extracting coding sequences and corresponding protein translations:
source activate bioperl
echo "perl ~/bin/CDS_extractor_v2.pl mcav_holobiont_iso.fasta myblast.br allhits bridgegaps" >cds
launcher_creator.py -j cds -n cds -l cddd -t 2:00:00 -q shortq7 -e mstudiva@fau.edu
sbatch cddd
source deactivate bioperl

# calculating contiguity:
contiguity.pl hits=mcav_holobiont_iso_hits.tab threshold=0.75
# 0.38

# core gene set form korflab: to characterize representation of genes:
wget http://korflab.ucdavis.edu/Datasets/genome_completeness/core/248.prots.fa.gz
gunzip 248.prots.fa.gz

makeblastdb -in mcav_holobiont.fasta -dbtype nucl
echo 'tblastn -query 248.prots.fa -db mcav_holobiont.fasta -evalue 1e-10 -outfmt "6 qseqid sseqid evalue bitscore qcovs" -max_target_seqs 1 -num_threads 12 >mcav_holobiont_248.brtab' >bl248
launcher_creator.py -j bl248 -n bl -l blj -q shortq7 -t 02:00:00 -e mstudiva@fau.edu
sbatch blj
# calculating fraction of represented KOGs:
cat mcav_holobiont_248.brtab | perl -pe 's/.+(KOG\d+)\s.+/$1/' | uniq | wc -l | awk '{print $1/248}'
# 0.959677

#------------------------------
# KOG annotation
# scp your *_PRO.fas file to laptop, submit it to
http://weizhong-lab.ucsd.edu/webMGA/server/kog/
cd /path/to/local/directory
scp mstudiva@koko-login.fau.edu:/path/to/HPC/directory/*_PRO.fas .

# copy link to job ID status and output file, paste it below instead of current link:
# check status: go on web to http://weizhong-lab.ucsd.edu/webMGA/result/?jobid=20180912075910950229025904
# once it is done, download results to HPC:
wget http://weizhong-lab.ucsd.edu/webMGA/result/output/20180912075910950229025904.zip

unzip 20180912075910950229025904.zip
mv kog/* .
rm -rf kog/
# As of 9/13/18, there is a problem with the webMGA KOG annotation where the correct output files are not being produced. I was informed that they are working on it.
# But, as a result, the following lines of code in this section do not currently work.

# generates iso2kogClass and iso2kogDef (another kind of gene names)
getKOGs.pl fastaQuery=mcav_holobiont_iso.fasta prefix=mcav_holobiont kogMatch=mcav_holobiont.kog.tab 

# removing "multiple classes" annotation, renaming comp to isogroup
grep -v "Multiple classes" mcav_holobiont_iso2kogClass.tab | perl -pe 's/^comp/isogroup/' > mcav_holobiont_iso2kogClassNR.tab

#------------------------------
# KEGG annotations:

# selecting the longest contig per isogroup:
fasta2SBH.pl mcav_holobiont_iso.fasta >mcav_holobiont_4kegg.fasta

# scp mcav_holobiont_4kegg.fasta to your laptop
cd /path/to/local/directory
scp mstudiva@koko-login.fau.edu:/path/to/HPC/directory/mcav_holobiont_4kegg.fasta .
# use web browser to submit mcav_holobiont_4kegg.fasta file to KEGG's KAAS server ( http://www.genome.jp/kegg/kaas/ )
# select SBH method, upload nucleotide query
# Once it is done, download the 'text' output from KAAS, name it query.ko (default)
https://www.genome.jp/kaas-bin/kaas_main?mode=user&id=1536766154&key=v6N9KqeR

wget https://www.genome.jp/tools/kaas/files/dl/1536766154/query.ko

# selecting only the lines with non-missing annotation:
cat query.ko | awk '{if ($2!="") print }' > mcav_holobiont_iso2kegg.tab

# the KEGG mapping result can be explored for completeness of transcriptome in terms of genes found,
# use 'html' output link from KAAS result page, see how many proteins you have for conserved complexes and pathways,
# such as ribosome, spliceosome, proteasome etc

#------------------------------
# copy all files to laptop
cd /path/to/local/directory
scp mstudiva@koko-login.fau.edu:/path/to/HPC/directory/* .