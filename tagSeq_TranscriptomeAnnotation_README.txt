# Transcriptome Annotation, version November 4, 2021
# Created by Misha Matz (matz@utexas.edu), modified by Michael Studivan (studivanms@gmail.com) for use on FAU's HPC (KoKo)
# for use in generating transcriptome annotation files for Montastraea cavernosa
# also includes the separation of reads associated with M. cavernosa and Cladocopium transcriptomes


#------------------------------
# BEFORE STARTING, replace, in this whole file:
#	- email@gmail.com by your actual email;
#	- username with your KoKo user name.

# The idea is to copy the chunks separated by empty lines below and paste them into your cluster
# terminal window consecutively.

# The lines beginning with hash marks (#) are explanations and additional instructions -
# please make sure to read them before copy-pasting.


#------------------------------
# setup

# To install Bioperl in your bin directory, please follow these instructions:
cd bin
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

git clone https://github.com/mstudiva/Mcav-Cladocopium-Annotated-Transcriptome.git
mv Mcav-Cladocopium-Annotated-Transcriptome/* .
rm -rf Mcav-Cladocopium-Annotated-Transcriptome

# creating backup directory
mkdir backup

# creating annotation directory
cd
mkdir annotate
cd annotate


#------------------------------
# getting transcriptomes

# updated (July 2018) M. cavernosa genome with transcriptome
wget https://www.dropbox.com/s/zj0rbpc35b6x707/Mcavernosa_genome.tgz
tar -xzf Mcavernosa_genome.tgz
cp Mcav_genome/Mcavernosa_annotation/Mcavernosa.maker.transcripts.fasta .
mv Mcavernosa.maker.transcripts.fasta Mcavernosa.fasta

# Cladocopium spp. (formerly Symbiodinium Clade C) transcriptome as of November 2017
wget http://sites.bu.edu/davieslab/files/2017/11/CladeC_Symbiodinium_transcriptome.zip
unzip CladeC_Symbiodinium_transcriptome.zip
cp CladeC_Symbiodinium_transcriptome/davies_cladeC_feb.fasta .
rm -rf __MACOSX/
mv davies_cladeC_feb.fasta Cladocopium.fasta

# use the stream editor to find and replace all instances of "comp" with "Cladocopium" in the symbiont transcriptome
sed -i 's/comp/Cladocopium/g' Cladocopium.fasta

# transcriptome statistics
echo "seq_stats.pl Mcavernosa.fasta > seqstats_Mcavernosa.txt" > seq_stats
echo "seq_stats.pl Cladocopium.fasta > seqstats_Cladocopium.txt" >> seq_stats
launcher_creator.py -j seq_stats -n seq_stats -q shortq7 -t 6:00:00 -e email@gmail.com
sbatch seq_stats.slurm

nano seqstats_Mcavernosa.txt

Mcavernosa.fasta
-------------------------
25142 sequences.
1420 average length.
43960 maximum length.
75 minimum length.
N50 = 2020
35.7 Mb altogether (35692498 bp).
0 ambiguous Mb. (1637 bp, 0%)
0 Mb of Ns. (1637 bp, 0%)
-------------------------

nano seqstats_Cladocopium.txt

Cladocopium.fasta
-------------------------
65838 sequences.
1482 average length.
18168 maximum length.
500 minimum length.
N50 = 1746
97.6 Mb altogether (97581498 bp).
0 ambiguous Mb. (0 bp, 0%)
0 Mb of Ns. (0 bp, 0%)
-------------------------


#------------------------------
# uniprot annotations with blast

# getting uniprot_swissprot KB database
wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz

# unzipping
gunzip uniprot_sprot.fasta.gz &

# indexing the fasta database
echo "makeblastdb -in uniprot_sprot.fasta -dbtype prot" >mdb
launcher_creator.py -j mdb -n mdb -q shortq7 -t 6:00:00 -e email@gmail.com
sbatch mdb.slurm

# splitting the transcriptome into 200 chunks
splitFasta.pl Mcavernosa.fasta 200
splitFasta.pl Cladocopium.fasta 200

# blasting all 200 chunks to uniprot in parallel, 4 cores per chunk
ls subset* | perl -pe 's/^(\S+)$/blastx -query $1 -db uniprot_sprot\.fasta -evalue 0\.0001 -num_threads 4 -num_descriptions 5 -num_alignments 5 -out $1.br/'>bl
launcher_creator.py -j bl -n blast -t 6:00:00 -q shortq7 -e email@gmail.com
sbatch blast.slurm

# watching progress:
grep "Query= " subset*.br | wc -l
# you should end up with the same number of queries as sequences from the seq_stats script

# combining all blast results
cat subset*br > myblast.br
mv subset* ~/annotate/backup/

# for trinity-assembled transcriptomes: annotating with "Cladocopium" or "Mcavernosa" depending on if component is from symbiont or host (=component)
grep ">" Mcavernosa.fasta | perl -pe 's/>Mcavernosa(\d+)(\S+)\s.+/Mcavernosa$1$2\tMcavernosa$1/'>Mcavernosa_seq2iso.tab
cat Mcavernosa.fasta | perl -pe 's/>Mcavernosa(\d+)(\S+).+/>Mcavernosa$1$2 gene=Mcavernosa$1/'>Mcavernosa_iso.fasta

grep ">" Cladocopium.fasta | perl -pe 's/>Cladocopium(\d+)(\S+)\s.+/Cladocopium$1$2\tCladocopium$1/' > Cladocopium_seq2iso.tab
cat Cladocopium.fasta | perl -pe 's/>Cladocopium(\d+)(\S+).+/>Cladocopium$1$2 gene=Cladocopium$1/' > Cladocopium_iso.fasta


#-------------------------
# extracting coding sequences and corresponding protein translations:
echo "perl ~/bin/CDS_extractor_v2.pl Mcavernosa_iso.fasta myblast.br allhits bridgegaps" >cds
launcher_creator.py -j cds -n cds -l cddd -t 6:00:00 -q shortq7 -e email@gmail.com
sbatch cddd

echo "perl ~/bin/CDS_extractor_v2.pl Cladocopium_iso.fasta myblast.br allhits bridgegaps" >cds
launcher_creator.py -j cds -n cds -l cddd -t 6:00:00 -q shortq7 -e email@gmail.com
sbatch cddd


#------------------------------
# GO annotation
# updated based on Misha Matz's new GO and KOG annotation steps on github: https://github.com/z0on/emapper_to_GOMWU_KOGMWU

# selecting the longest contig per isogroup (also renames using isogroups based on Mcavernosa and Cladocopium annotations):
fasta2SBH_Mcav.pl Mcavernosa_iso_PRO.fas >Mcavernosa_out_PRO.fas

fasta2SBH_Mcav.pl Cladocopium_iso_PRO.fas >Cladocopium_out_PRO.fas

# scp your *_out_PRO.fas file to laptop, submit it to
http://eggnog-mapper.embl.de
cd /path/to/local/directory
scp username@koko-login.hpc.fau.edu:~/path/to/HPC/directory/*_out_PRO.fas .

# copy link to job ID status and output file, paste it below instead of current link:
# Mcav status: go on web to http://eggnog-mapper.embl.de/job_status?jobname=MM_s3chyymz
# symC status: go on web to http://eggnog-mapper.embl.de/job_status?jobname=MM_0wawlio7

# once it is done, download results to HPC:
wget http://eggnog-mapper.embl.de/MM_s3chyymz/out.emapper.annotations
wget http://eggnog-mapper.embl.de/MM_0wawlio7/out.emapper.annotations

# GO:
awk -F "\t" 'BEGIN {OFS="\t" }{print $1,$10 }' out.emapper.annotations | grep GO | perl -pe 's/,/;/g' >Mcavernosa_iso2go.tab
# gene names:
awk -F "\t" 'BEGIN {OFS="\t" }{print $1,$8 }' out.emapper.annotations | grep -Ev "\tNA" >Mcavernosa_iso2geneName.tab

awk -F "\t" 'BEGIN {OFS="\t" }{print $1,$10 }' out.emapper.annotations | grep GO | perl -pe 's/,/;/g' >Cladocopium_iso2go.tab
# gene names:
awk -F "\t" 'BEGIN {OFS="\t" }{print $1,$8 }' out.emapper.annotations | grep -Ev "\tNA" >Cladocopium_iso2geneName.tab


#------------------------------
# KOG annotation
# updated based on Misha Matz's new GO and KOG annotation steps on github: https://github.com/z0on/emapper_to_GOMWU_KOGMWU

cp ~/bin/kog_classes.txt .

#  KOG classes (single-letter):
awk -F "\t" 'BEGIN {OFS="\t" }{print $1,$7 }' out.emapper.annotations | grep -Ev "[,#S]" >Mcavernosa_iso2kogClass1.tab
# converting single-letter KOG classes to text understood by KOGMWU package (must have kog_classes.txt file in the same dir):
awk 'BEGIN {FS=OFS="\t"} NR==FNR {a[$1] = $2;next} {print $1,a[$2]}' kog_classes.txt Mcavernosa_iso2kogClass1.tab > Mcavernosa_iso2kogClass.tab

awk -F "\t" 'BEGIN {OFS="\t" }{print $1,$7 }' out.emapper.annotations | grep -Ev "[,#S]" >Cladocopium_iso2kogClass1.tab
awk 'BEGIN {FS=OFS="\t"} NR==FNR {a[$1] = $2;next} {print $1,a[$2]}' kog_classes.txt Cladocopium_iso2kogClass1.tab > Cladocopium_iso2kogClass.tab


#------------------------------
# KEGG annotations:

# selecting the longest contig per isogroup:
srun fasta2SBH_Mcav.pl Mcavernosa_iso.fasta >Mcavernosa_4kegg.fasta

srun fasta2SBH_Mcav.pl Cladocopium_iso.fasta >Cladocopium_4kegg.fasta

# scp *4kegg.fasta to your laptop
cd /path/to/local/directory
scp mstudiva@koko-login.hpc.fau.edu:~/path/to/HPC/directory/*4kegg.fasta .
# use web browser to submit 4kegg.fasta file to KEGG's KAAS server (http://www.genome.jp/kegg/kaas/)
# select SBH method, upload nucleotide query
https://www.genome.jp/kaas-bin/kaas_main?mode=user&id=1636136731&key=fzxEy2PS
https://www.genome.jp/kaas-bin/kaas_main?mode=user&id=1636139487&key=VECJ5NEi

# Once it is done, download to HPC - it is named query.ko by default
wget https://www.genome.jp/tools/kaas/files/dl/1636136731/query.ko
wget https://www.genome.jp/tools/kaas/files/dl/1636139487/query.ko

# selecting only the lines with non-missing annotation:
cat query.ko | awk '{if ($2!="") print }' > Mcavernosa_iso2kegg.tab

cat query.ko | awk '{if ($2!="") print }' > Cladocopium_iso2kegg.tab

# the KEGG mapping result can be explored for completeness of transcriptome in terms of genes found,
# use 'html' output link from KAAS result page, see how many proteins you have for conserved complexes and pathways, such as ribosome, spliceosome, proteasome etc


#------------------------------
# file transfer

# copy all files to local machine
cd /path/to/local/directory
scp mstudiva@koko-login.fau.edu:~/path/to/HPC/directory/* .
