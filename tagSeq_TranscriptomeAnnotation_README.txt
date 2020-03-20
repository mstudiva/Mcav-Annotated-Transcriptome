# Transcriptome Annotation, version March 20, 2020
# Created by Misha Matz (matz@utexas.edu), modified by Michael Studivan (studivanms@gmail.com)
# for use in generating transcriptome annotation files for Montastraea cavernosa
# also includes the concatention of M. cavernosa and Cladocopium sp. (formerly Symbiodinium Clade C) transcriptomes

#------------------------------
# BEFORE STARTING, replace, in this whole file:
#	- studivanms@gmail.com by your actual email;
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

git clone https://github.com/mstudiva/Mcav-Annotated-Transcriptome
mv Mcav-Annotated-Transcriptome/* .
rm -rf Mcav-Annotated-Transcriptome

# creating backup directory
mkdir backup

# creating annotation directory
cd
mkdir annotate
cd annotate

# original transcriptome from Kitchen et al (2015)
# M. cavernosa transcriptome, v1 (2014)
# wget http://meyerlab:coral@files.cgrb.oregonstate.edu/Meyer_Lab/transcriptomes/Mcav/Mcav_transcriptome_v1.fasta.gz
# gunzip Mcav_transcriptome_v1.fasta.gz
# mv Mcav_transcriptome_v1.fasta Mcavernosa.fasta

# updated (July 2018) M. cavernosa genome with transcriptome
wget https://www.dropbox.com/s/0inwmljv6ti643o/Mcavernosa_genome.tgz
tar -xzf Mcavernosa_genome.tgz
cp Mcav_genome/Mcavernosa_annotation/Mcavernosa.maker.transcripts.fasta .
mv Mcavernosa.maker.transcripts.fasta Mcavernosa.fasta

# Cladocopium sp. (formerly Symbiodinium Clade C) transcriptome as of November 2017
wget http://sites.bu.edu/davieslab/files/2017/11/CladeC_Symbiodinium_transcriptome.zip
unzip CladeC_Symbiodinium_transcriptome.zip
cp CladeC_Symbiodinium_transcriptome/davies_cladeC_feb.fasta .
rm -rf __MACOSX/
mv davies_cladeC_feb.fasta Cladocopium.fasta

# use the stream editor to find and replace all instances of "comp" with "Cladocopium" in the symbiont transcriptome
sed -i 's/comp/Cladocopium/g' Cladocopium.fasta

# concatenate the host and symbiont transcriptomes into a holobiont transcriptome
cat Cladocopium.fasta Mcavernosa.fasta > Mcavernosa_Cladocopium.fasta

# transcriptome statistics
seq_stats.pl Mcavernosa_Cladocopium.fasta > seqstats_Mcavernosa_Cladocopium.txt
nano seqstats_Mcavernosa_Cladocopium.txt

# Mcavernosa_Cladocopium.fasta
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
launcher_creator.py -j getz -n getz -t 6:00:00 -q shortq7 -e studivanms@gmail.com
sbatch getz.slurm

# if the US mirror is down, uncomment the line below, then run the getz script as normal
# echo 'wget ftp://ftp.ebi.ac.uk/pub/databases/uniprot/current_release/knowledgebase/idmapping/idmapping_selected.tab.gz '> getz

# unzipping
gunzip uniprot_sprot.fasta.gz &
gunzip idmapping_selected.tab.gz &

# indexing the fasta database
module load blast
echo "makeblastdb -in uniprot_sprot.fasta -dbtype prot" >mdb
launcher_creator.py -j mdb -n mdb -q shortq7 -t 6:00:00 -e studivanms@gmail.com
sbatch mdb.slurm

# splitting the transcriptome into 190 chunks
splitFasta.pl Mcavernosa_Cladocopium.fasta 190

# blasting all 190 chunks to uniprot in parallel, 4 cores per chunk
ls subset* | perl -pe 's/^(\S+)$/blastx -query $1 -db uniprot_sprot\.fasta -evalue 0\.0001 -num_threads 4 -num_descriptions 5 -num_alignments 5 -out $1.br/'>bl
launcher_creator.py -j bl -n blast -t 6:00:00 -q shortq7 -e studivanms@gmail.com
sbatch blast.slurm

# watching progress:
grep "Query= " subset*.br | wc -l
# you should end up with the same number of queries as sequences from the seq_stats script (90980 sequences)

# combining all blast results
cat subset*br > myblast.br
mv subset* ~/backup/

# for trinity-assembled transcriptomes: annotating with "Cladocopium" or "Mcavernosa" depending on if component is from symbiont or host (=component)
grep ">" Mcavernosa_Cladocopium.fasta | perl -pe 's/>Cladocopium(\d+)(\S+)\s.+/Cladocopium$1$2\tCladocopium$1/' | perl -pe 's/>Mcavernosa(\d+)(\S+)\s.+/Mcavernosa$1$2\tMcavernosa$1/'>Mcavernosa_Cladocopium_seq2iso.tab
cat Mcavernosa_Cladocopium.fasta | perl -pe 's/>Cladocopium(\d+)(\S+).+/>Cladocopium$1$2 gene=Cladocopium$1/' | perl -pe 's/>Mcavernosa(\d+)(\S+).+/>Mcavernosa$1$2 gene=Mcavernosa$1/'>Mcavernosa_Cladocopium_iso.fasta

#-------------------------
# extracting coding sequences and corresponding protein translations:
echo "perl ~/bin/CDS_extractor_v2.pl Mcavernosa_Cladocopium_iso.fasta myblast.br allhits bridgegaps" >cds
launcher_creator.py -j cds -n cds -l cddd -t 6:00:00 -q shortq7 -e studivanms@gmail.com
sbatch cddd

# use the stream editor to remove all instances of "gene=" from the query IDs in the CDS_extractor outputs
# not needed unless running CDS_extractor_MS on the transcriptome_iso file
# sed -i 's/gene=//g' Mcavernosa_Cladocopium_iso_CDSends.fas
# sed -i 's/gene=//g' Mcavernosa_Cladocopium_iso_CDS.fas
# sed -i 's/gene=//g' Mcavernosa_Cladocopium_iso_hits.tab
# sed -i 's/gene=//g' Mcavernosa_Cladocopium_iso_PROends.fas
# sed -i 's/gene=//g' Mcavernosa_Cladocopium_iso_PRO.fas

# calculating contiguity:
contiguity.pl hits=Mcavernosa_Cladocopium_iso_hits.tab threshold=0.75
# contiguity at 0.75 threshold: 0.38

# core gene set from korflab: to characterize representation of genes:
wget http://korflab.ucdavis.edu/Datasets/genome_completeness/core/248.prots.fa.gz
gunzip 248.prots.fa.gz

module load blast
makeblastdb -in Mcavernosa_Cladocopium_iso.fasta -dbtype nucl
echo 'tblastn -query 248.prots.fa -db Mcavernosa_Cladocopium_iso.fasta -evalue 1e-10 -outfmt "6 qseqid sseqid evalue bitscore qcovs" -max_target_seqs 1 -num_threads 12 >Mcavernosa_Cladocopium_248.brtab' >bl248
launcher_creator.py -j bl248 -n bl -l blj -q shortq7 -t 06:00:00 -e studivanms@gmail.com
sbatch blj

# calculating fraction of represented KOGs:
cat Mcavernosa_Cladocopium_248.brtab | perl -pe 's/.+(KOG\d+)\s.+/$1/' | uniq | wc -l | awk '{print $1/248}'
# 0.959677

#------------------------------
# GO annotation
# updated based on Misha Matz's new GO and KOG annotation steps on github: https://github.com/z0on/emapper_to_GOMWU_KOGMWU

# selecting the longest contig per isogroup (also renames using isogroups based on Mcavernosa and Cladocopium annotations):
fasta2SBH_MS.pl Mcavernosa_Cladocopium_iso_PRO.fas >Mcavernosa_Cladocopium_out_PRO.fas

# scp your *_out_PRO.fas file to laptop, submit it to
http://eggnogdb.embl.de/#/app/emapper
cd /path/to/local/directory
scp mstudiva@koko-login.fau.edu:~/path/to/HPC/directory/*_out_PRO.fas .

# copy link to job ID status and output file, paste it below instead of current link:
# check status: go on web to http://eggnogdb.embl.de/#/app/emapper?jobname=MM_w4_mOZ
# once it is done, download results to HPC:
wget http://eggnogdb.embl.de/MM_w4_mOZ/Mcavernosa_Cladocopium_out_PRO.fas.emapper.annotations

# GO:
awk -F "\t" 'BEGIN {OFS="\t" }{print $1,$6 }' Mcavernosa_Cladocopium_out_PRO.fas.emapper.annotations | grep GO | perl -pe 's/,/;/g' >Mcavernosa_Cladocopium_iso2go.tab
# gene names:
awk -F "\t" 'BEGIN {OFS="\t" }{print $1,$13 }' Mcavernosa_Cladocopium_out_PRO.fas.emapper.annotations | grep -Ev "\tNA" >Mcavernosa_Cladocopium_iso2geneName.tab

#------------------------------
# KOG annotation
# updated based on Misha Matz's new GO and KOG annotation steps on github: https://github.com/z0on/emapper_to_GOMWU_KOGMWU

cp ~/bin/kog_classes.txt .

#  KOG classes (single-letter):
awk -F "\t" 'BEGIN {OFS="\t" }{print $1,$12 }' Mcavernosa_Cladocopium_out_PRO.fas.emapper.annotations | grep -Ev "[,#S]" >Mcavernosa_Cladocopium_iso2kogClass1.tab
# converting single-letter KOG classes to text understood by KOGMWU package (must have kog_classes.txt file in the same dir):
awk 'BEGIN {FS=OFS="\t"} NR==FNR {a[$1] = $2;next} {print $1,a[$2]}' kog_classes.txt Mcavernosa_Cladocopium_iso2kogClass1.tab > Mcavernosa_Cladocopium_iso2kogClass.tab

#------------------------------
# KEGG annotations:

# selecting the longest contig per isogroup:
fasta2SBH_MS.pl Mcavernosa_Cladocopium_iso.fasta >Mcavernosa_Cladocopium_4kegg.fasta

# scp Mcavernosa_Cladocopium_4kegg.fasta to your laptop
cd /path/to/local/directory
scp mstudiva@koko-login.fau.edu:~/path/to/HPC/directory/Mcavernosa_Cladocopium_4kegg.fasta .
# use web browser to submit Mcavernosa_Cladocopium_4kegg.fasta file to KEGG's KAAS server ( http://www.genome.jp/kegg/kaas/ )
# select SBH method, upload nucleotide query
# Once it is done, download the 'text' output from KAAS, name it query.ko (default)
https://www.genome.jp/kaas-bin/kaas_main?mode=user&id=1556130193&key=KCa76L0E

wget https://www.genome.jp/tools/kaas/files/dl/1556130193/query.ko

# selecting only the lines with non-missing annotation:
cat query.ko | awk '{if ($2!="") print }' > Mcavernosa_Cladocopium_iso2kegg.tab

# the KEGG mapping result can be explored for completeness of transcriptome in terms of genes found,
# use 'html' output link from KAAS result page, see how many proteins you have for conserved complexes and pathways, such as ribosome, spliceosome, proteasome etc

#------------------------------
# move the very large idmapping_selected.tab to backup
mv idmapping_selected.tab ~/backup/

# copy all files to laptop
cd /path/to/local/directory
scp mstudiva@koko-login.fau.edu:~/path/to/HPC/directory/* .
