mkdir -p $1
cd $1

wget https://sam.s3.climb.ac.uk/dehumanizer/GCA_000786075.2_hs38d1_genomic.fna.mmi
wget https://sam.s3.climb.ac.uk/dehumanizer/GCA_000001405.27_GRCh38.p12_genomic.fna.mmi
wget https://sam.s3.climb.ac.uk/dehumanizer/hla_gen.fasta.mmi

echo "asm10" > manifest.txt
echo "hs38d1	$(pwd)/GCA_000786075.2_hs38d1_genomic.fna.mmi" >> manifest.txt
echo "GRCh38	$(pwd)/GCA_000001405.27_GRCh38.p12_genomic.fna.mmi" >> manifest.txt
echo "HLA	$(pwd)/hla_gen.fasta.mmi" >> manifest.txt
