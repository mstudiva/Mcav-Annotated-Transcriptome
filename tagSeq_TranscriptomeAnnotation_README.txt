# Transcriptome Annotation, version Feb 7, 2019
# Created by Misha Matz (matz@utexas.edu), modified by Michael Studivan (mstudiva@fau.edu)
# for use in generating transcriptome annotation files for Montastraea cavernosa
# also includes the concatention of M. cavernosa and Cladocopium sp. (formerly Symbiodinium Clade C) transcriptomes

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

git clone https://github.com/z0on/emapper_to_GOMWU_KOGMWU.git
mv emapper_to_GOMWU_KOGMWU/* .
rm -rf emapper_to_GOMWU_KOGMWU

# creating annotation dir
cd
mkdir annotate
cd annotate

# old version, commented out for now
# M. cavernosa transcriptome, v1 (2014)
# wget http://meyerlab:coral@files.cgrb.oregonstate.edu/Meyer_Lab/transcriptomes/Mcav/Mcav_transcriptome_v1.fasta.gz
# gunzip Mcav_transcriptome_v1.fasta.gz
# mv Mcav_transcriptome_v1.fasta Mcavernosa.fasta

# updated (July 2018) M. cavernosa genome with transcriptome
wget https://www.dropbox.com/s/0inwmljv6ti643o/Mcavernosa_genome.tgz
tar -xzf Mcavernosa_genome.tgz
cd Mcav_genome/Mcavernosa_annotation/
cp Mcavernosa.maker.transcripts.fasta ../..
cd ../..
mv Mcavernosa.maker.transcripts.fasta Mcavernosa.fasta

# old version, commented out for now
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
mv davies_cladeC_feb.fasta Cladocopium.fasta
mv davies_cladeC_iso2go.tab Cladocopium_iso2go.tab
mv davies_cladeC_iso2gene.tab Cladocopium_iso2gene.tab
mv davies_cladeC_seq2iso.tab Cladocopium_seq2iso.tab

# use the stream editor to find and replace all instances of "comp" and "isogroup" with "sym" in the symbiont transcriptome files
sed -i 's/comp/Cladocopium/g' Cladocopium_seq2iso.tab
sed -i 's/comp/Cladocopium/g' Cladocopium.fasta
sed -i 's/isogroup/Cladocopium/g' Cladocopium_iso2gene.tab
sed -i 's/isogroup/Cladocopium/g' Cladocopium_iso2go.tab
sed -i 's/isogroup/Cladocopium/g' Cladocopium_seq2iso.tab

# concatenate the host and symbiont transcriptomes into a holobiont transcriptome
cat Cladocopium.fasta Mcavernosa.fasta > Mcavernosa_Cladocopium.fasta

# transcriptome statistics
source activate bioperl
seq_stats.pl Mcavernosa_Cladocopium.fasta
source deactivate bioperl

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
splitFasta.pl Mcavernosa_Cladocopium.fasta 190

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
grep ">" Mcavernosa_Cladocopium.fasta | perl -pe 's/>sym(\d+)(\S+)\s.+/sym$1$2\tsym$1/' | perl -pe 's/>Mcavernosa(\d+)(\S+)\s.+/Mcavernosa$1$2\tMcavernosa$1/'>Mcavernosa_Cladocopium_seq2iso.tab
cat Mcavernosa_Cladocopium.fasta | perl -pe 's/>sym(\d+)(\S+).+/>sym$1$2 gene=sym$1/' | perl -pe 's/>Mcavernosa(\d+)(\S+).+/>Mcavernosa$1$2 gene=Mcavernosa$1/'>Mcavernosa_Cladocopium_iso.fasta

#-------------------------
# old code, commenting out for now
# extracting gene names (per isogroup):
# getGeneNameFromUniProtKB.pl blast=myblast.br prefix=Mcavernosa_Cladocopium fastaQuery=Mcavernosa_Cladocopium_iso.fasta

# extracting GO annotations (per isogroup)
# echo "getGOfromUniProtKB.pl blast=myblast.br prefix=Mcavernosa_Cladocopium fastaQuery=Mcavernosa_Cladocopium_iso.fasta" >getgo
# launcher_creator.py -j getgo -n getgo -l gg -q shortq7 -t 2:00:00 -e mstudiva@fau.edu
# sbatch gg

# extracting coding sequences and corresponding protein translations:
source activate bioperl
echo "perl ~/bin/CDS_extractor_v2.pl Mcavernosa_Cladocopium_iso.fasta myblast.br allhits bridgegaps" >cds
launcher_creator.py -j cds -n cds -l cddd -t 2:00:00 -q shortq7 -e mstudiva@fau.edu
sbatch cddd
source deactivate bioperl

# calculating contiguity:
contiguity.pl hits=Mcavernosa_Cladocopium_iso_hits.tab threshold=0.75
# 0.38

# core gene set from korflab: to characterize representation of genes:
wget http://korflab.ucdavis.edu/Datasets/genome_completeness/core/248.prots.fa.gz
gunzip 248.prots.fa.gz

makeblastdb -in Mcavernosa_Cladocopium.fasta -dbtype nucl
echo 'tblastn -query 248.prots.fa -db Mcavernosa_Cladocopium.fasta -evalue 1e-10 -outfmt "6 qseqid sseqid evalue bitscore qcovs" -max_target_seqs 1 -num_threads 12 >Mcavernosa_Cladocopium_248.brtab' >bl248
launcher_creator.py -j bl248 -n bl -l blj -q shortq7 -t 02:00:00 -e mstudiva@fau.edu
sbatch blj
# calculating fraction of represented KOGs:
cat Mcavernosa_Cladocopium_248.brtab | perl -pe 's/.+(KOG\d+)\s.+/$1/' | uniq | wc -l | awk '{print $1/248}'
# 0.959677

#------------------------------
# GO annotation
# updated based on Misha Matz's new GO and KOG annotation steps on github: https://github.com/z0on/emapper_to_GOMWU_KOGMWU

# scp your *_PRO.fas file to laptop, submit it to
http://eggnogdb.embl.de/#/app/emapper
cd /path/to/local/directory
scp mstudiva@koko-login.fau.edu:~/path/to/HPC/directory/*_PRO.fas .

# copy link to job ID status and output file, paste it below instead of current link:
# check status: go on web to http://eggnogdb.embl.de/#/app/emapper?jobname=MM_7_Z9m9
# once it is done, download results to HPC:
wget http://eggnogdb.embl.de/MM_7_Z9m9/Mcavernosa_Cladocopium_iso_PRO.fas.emapper.annotations

# GO:
awk -F "\t" 'BEGIN {OFS="\t" }{print $1,$6 }' Mcavernosa_Cladocopium_iso_PRO.fas.emapper.annotations | grep GO | perl -pe 's/,/;/g' >Mcavernosa_Cladocopium_iso2go.tab
# gene names:
awk -F "\t" 'BEGIN {OFS="\t" }{print $1,$13 }' Mcavernosa_Cladocopium_iso_PRO.fas.emapper.annotations | grep -Ev "\tNA" >Mcavernosa_Cladocopium_iso2geneName.tab

#------------------------------
# KOG annotation
# updated based on Misha Matz's new GO and KOG annotation steps on github: https://github.com/z0on/emapper_to_GOMWU_KOGMWU

cp ~/bin/kog_classes.txt .

#  KOG classes (single-letter):
awk -F "\t" 'BEGIN {OFS="\t" }{print $1,$12 }' Mcavernosa_Cladocopium_iso_PRO.fas.emapper.annotations | grep -Ev "[,#S]" >Mcavernosa_Cladocopium_iso2kogClass1.tab
# converting single-letter KOG classes to text understood by KOGMWU package (must have kog_classes.txt file in the same dir):
awk 'BEGIN {FS=OFS="\t"} NR==FNR {a[$1] = $2;next} {print $1,a[$2]}' kog_classes.txt Mcavernosa_Cladocopium_iso2kogClass1.tab > Mcavernosa_Cladocopium_iso2kogClass.tab

# old code, commenting out for now
# generates iso2kogClass and iso2kogDef (another kind of gene names)
# getKOGs.pl fastaQuery=Mcavernosa_Cladocopium_iso.fasta prefix=Mcavernosa_Cladocopium kogMatch=Mcavernosa_Cladocopium.kog.tab

# removing "multiple classes" annotation, renaming comp to isogroup
# grep -v "Multiple classes" Mcavernosa_Cladocopium_iso2kogClass.tab | perl -pe 's/^comp/isogroup/' > Mcavernosa_Cladocopium_iso2kogClassNR.tab

#------------------------------
# KEGG annotations:

# selecting the longest contig per isogroup:
fasta2SBH.pl Mcavernosa_Cladocopium_iso.fasta >Mcavernosa_Cladocopium_4kegg.fasta

# scp Mcavernosa_Cladocopium_4kegg.fasta to your laptop
cd /path/to/local/directory
scp mstudiva@koko-login.fau.edu:~/path/to/HPC/directory/Mcavernosa_Cladocopium_4kegg.fasta .
# use web browser to submit Mcavernosa_Cladocopium_4kegg.fasta file to KEGG's KAAS server ( http://www.genome.jp/kegg/kaas/ )
# select SBH method, upload nucleotide query
# Once it is done, download the 'text' output from KAAS, name it query.ko (default)
https://www.genome.jp/kaas-bin/kaas_main?mode=user&id=1550845868&key=h7cIpELd

wget https://www.genome.jp/tools/kaas/files/dl/1550845868/query.ko

# selecting only the lines with non-missing annotation:
cat query.ko | awk '{if ($2!="") print }' > Mcavernosa_Cladocopium_iso2kegg.tab

# the KEGG mapping result can be explored for completeness of transcriptome in terms of genes found,
# use 'html' output link from KAAS result page, see how many proteins you have for conserved complexes and pathways,
# such as ribosome, spliceosome, proteasome etc

#------------------------------
# copy all files to laptop
cd /path/to/local/directory
scp mstudiva@koko-login.fau.edu:~/path/to/HPC/directory/* .
