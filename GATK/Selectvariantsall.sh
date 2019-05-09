#
# Original script by Marni Tausen https://github.com/MarniTausen/CloverAnalysisPipeline
#

#!/bin/bash

source /com/extra/java/8/load.sh
source /com/extra/GATK/3.6/load.sh

java -Xmx64g -jar /com/extra/GATK/3.6/jar-bin/GenomeAnalysisTK.jar \
     -R /home/runewind/NChain/faststorage/WHITE_CLOVER/RNASEQ/REFERENCE/TrR.v5.fasta \
     -T SelectVariants \
     -o /home/runewind/NChain/faststorage/WHITE_CLOVER/RNASEQ/VARIANT_CALLING/GATK/WC_SNP_Qual.vcf \
     --variant /home/runewind/NChain/faststorage/WHITE_CLOVER/RNASEQ/VARIANT_CALLING/GATK/WC_RNA_SNP.vcf \
     --restrictAllelesTo ALL -select "MQ>30.00 && DP>140 && QUAL>20.00"
