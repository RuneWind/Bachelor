#
# Original script by Marni Tausen https://github.com/MarniTausen/CloverAnalysisPipeline
#

#!/bin/bash

source /com/extra/java/8/load.sh
source /com/extra/GATK/3.6/load.sh

java -Xmx64g -jar /com/extra/GATK/3.6/jar-bin/GenomeAnalysisTK.jar \
     -T UnifiedGenotyper \
     -I $1 \
     -o $2 \
     -stand_call_conf 50 \
     -stand_emit_conf 20.0 \
     --sample_ploidy 2 \
     -nct 12 --genotyping_mode DISCOVERY \
     -R /faststorage/project/NChain/WHITE_CLOVER/RNASEQ/REFERENCE/TrR.v5.fasta \
     --output_mode EMIT_VARIANTS_ONLY \
     -rf BadCigar 2>&1 > log
