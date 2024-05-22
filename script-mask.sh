IN=/home/mdelima/share/accessibility/

for CHR in 2L 2R 3L 3R X; do

	bcftools view accessibility.$CHR.vcf.gz -f PASS | bcftools query -f '%CHROM\t%POS\n'>  accessibility.$CHR.mask

done
