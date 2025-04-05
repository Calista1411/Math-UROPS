# Running MAFFT for multiple sequence alignment
mafft --auto --thread -1 --keeplength --quiet --mapout --preservecase --addfragments asia_seq.fasta refseq.fasta > asia_seq_aligned.fasta    

# Getting SNP counts using SNP-sites
docker run -v $(pwd):/data --rm biocontainers/snp-sites:v2.4.1-1-deb_cv1 snp-sites -v -o /data/snp_asia..vcf  /data/asia_seq_aligned.fasta
Rscript positions.R snp_asia.vcf snp_asia.txt    
wc -w < snp_asia.txt

# Running PHI Test
./Phi -f asia_seq_aligned.fasta

# Running GARD
git clone https://github.com/veg/hyphy.git
cd hyphy
docker build -t hyphy:latest .
docker run --rm -v ~/Downloads/snp-sites:/hyphy/data -it hyphy:latest
./HYPHYMP ENV=TOLERATE_NUMERICAL_ERRORS=1 LIBPATH=/hyphy/res /hyphy/res/TemplateBatchFiles/GARD.bf --alignment /hyphy/data/asia_seq_aligned.fasta --type nucleotide

# Running KwARG
./kwarg -T50,30 -Q5 \
-S0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1,0.01,1 \
-M0.91,0.81,0.71,0.61,0.51,0.41,0.31,0.21,0.11,0.02,1.1 \
-R1,1,1,1,1,1,1,1,1,1,-1 \
-C2,2,2,2,2,2,2,2,2,2,-1 \
-n -f < asia_seq_aligned.fasta > asia_seq_kwarg.fasta
