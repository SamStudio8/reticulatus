mkdir -p $1
cd $1

wget https://sam.s3.climb.ac.uk/dehumanizer/20200421/GCA_000786075.2_hs38d1_genomic.mmi
wget https://sam.s3.climb.ac.uk/dehumanizer/20200421/GCF_000001405.39_GRCh38.p13_genomic.mmi
wget https://sam.s3.climb.ac.uk/dehumanizer/20200421/ipd-imgt-3_39_0.hla_gen.mmi

echo "hs38d1 $(pwd)/GCA_000786075.2_hs38d1_genomic.mmi" >> manifest.txt
echo "GRCh38 $(pwd)/GCF_000001405.39_GRCh38.p13_genomic.mmi" >> manifest.txt
echo "HLA $(pwd)/ipd-imgt-3_39_0.hla_gen.mmi" >> manifest.txt
