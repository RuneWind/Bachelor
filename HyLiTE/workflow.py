#
# Original script by Marni Tausen https://github.com/MarniTausen/CloverAnalysisPipeline
#

from gwf import Workflow

gwf = Workflow()


def bowtie2_index(reference, reference_name):
    inputs = [reference]
    outputs = [reference_name+".1.bt2"]
    options = {
        'cores': 1,
        'memory': '16g',
        'account': 'NChain',
        'walltime': '12:00:00'
    }

    spec = '''
    source activate HyLiTE

    bowtie2-build {ref} {ref_n}
    '''.format(ref=reference, ref_n=reference_name)

    return inputs, outputs, options, spec

def star_index(reference, output):
    inputs = [reference]
    outputs = [output]
    options = {"cores":8, "memory":"64g", "account":"NChain", "walltime": "12:00:00"}

    directory = "/".join(output.split("/")[:-1])

    spec = """
    source /com/extra/STAR/2.5.2b/load.sh
    STAR --runMode genomeGenerate --runThreadN 8 --genomeDir {dir} --genomeFastaFiles {ref} --limitGenomeGenerateRAM=64000000000 --genomeSAindexNbases 3
    """.format(ref=reference, dir=directory)

    return inputs, outputs, options, spec


reference_file = "./references/To/To.v5.gDNA.fasta"
index_ref_file = "./references/To/To.ref"

raw_gDNA = {"To": ["/faststorage/project/NChain/20181120_clover_180bp_gDNA/Trifolium_occidentale_180bp_1.fastq.gz",
                   "/faststorage/project/NChain/20181120_clover_180bp_gDNA/Trifolium_occidentale_180bp_2.fastq.gz"],
            "Tp": ["/faststorage/project/NChain/20181120_clover_180bp_gDNA/Trifolium_pallescens_180_1.fastq.gz",
                   "/faststorage/project/NChain/20181120_clover_180bp_gDNA/Trifolium_pallescens_180_2.fastq.gz"],
            "TrR": ["/faststorage/project/NChain/20181120_clover_180bp_gDNA/Trifolium_repens_180bp_1.fastq.gz",
                    "/faststorage/project/NChain/20181120_clover_180bp_gDNA/Trifolium_repens_180bp_2.fastq.gz"]}

raw_RNA_old = {"To": {"floral": ["../BRAKER_ANNOTATION/pipeline/RNAdata/occidentale/To_F2_pooled_1.fastq",
                             "../BRAKER_ANNOTATION/pipeline/RNAdata/occidentale/To_F2_pooled_2.fastq"],
                  "leaf": ["/faststorage/project/NChain/WHITE_CLOVER/20181211_RNASEQ/old_rnaseq/To/To_1_1.fastq",
                           "/faststorage/project/NChain/WHITE_CLOVER/20181211_RNASEQ/old_rnaseq/To/To_1_2.fastq"],
                  "stolon": ["/faststorage/project/NChain/WHITE_CLOVER/20181211_RNASEQ/old_rnaseq/To/To_3_1.fastq",
                             "/faststorage/project/NChain/WHITE_CLOVER/20181211_RNASEQ/old_rnaseq/To/To_3_2.fastq"],
                  "root": ["/faststorage/project/NChain/WHITE_CLOVER/20181211_RNASEQ/old_rnaseq/To/To_6_1.fastq",
                           "/faststorage/project/NChain/WHITE_CLOVER/20181211_RNASEQ/old_rnaseq/To/To_6_2.fastq"]},
           "Tp": {"floral": ["/faststorage/project/NChain/WHITE_CLOVER/20181211_RNASEQ/old_rnaseq/Tp/Tp_F1_pooled_1.fastq",
                             "/faststorage/project/NChain/WHITE_CLOVER/20181211_RNASEQ/old_rnaseq/Tp/Tp_F1_pooled_2.fastq"],
                  "leaf": ["/faststorage/project/NChain/WHITE_CLOVER/20181211_RNASEQ/old_rnaseq/Tp/Tp_L4_1_1.fastq",
                           "/faststorage/project/NChain/WHITE_CLOVER/20181211_RNASEQ/old_rnaseq/Tp/Tp_L4_1_2.fastq"],
                  "stolon": ["/faststorage/project/NChain/WHITE_CLOVER/20181211_RNASEQ/old_rnaseq/Tp/T_pal_stolon_1_1_9_1.fastq",
                             "/faststorage/project/NChain/WHITE_CLOVER/20181211_RNASEQ/old_rnaseq/Tp/T_pal_stolon_1_1_9_2.fastq"],
                  "root": ["/faststorage/project/NChain/WHITE_CLOVER/20181211_RNASEQ/old_rnaseq/Tp/Tp_root_3e1_47_1.fastq",
                           "/faststorage/project/NChain/WHITE_CLOVER/20181211_RNASEQ/old_rnaseq/Tp/Tp_root_3e1_47_2.fastq"]},
           "TrR": {"floral": ["../BRAKER_ANNOTATION/pipeline/RNAdata/repens/floral/TRfloral_1.fastq",
                              "../BRAKER_ANNOTATION/pipeline/RNAdata/repens/floral/TRfloral_2.fastq"],
                   "leaf": ["../BRAKER_ANNOTATION/pipeline/RNAdata/repens/leaf/TRleaf_1.fastq",
                            "../BRAKER_ANNOTATION/pipeline/RNAdata/repens/leaf/TRleaf_2.fastq"],
                   "stolon": ["../BRAKER_ANNOTATION/pipeline/RNAdata/repens/stolon/TRstolon_1.fastq",
                              "../BRAKER_ANNOTATION/pipeline/RNAdata/repens/stolon/TRstolon_2.fastq"],
                   "root": ["../BRAKER_ANNOTATION/pipeline/RNAdata/repens/root/TRroot_1.fastq",
                            "../BRAKER_ANNOTATION/pipeline/RNAdata/repens/root/TRroot_2.fastq"]}}



raw_RNA = {"To": {"YL1": ["/faststorage/project/NChain/WHITE_CLOVER/20181211_RNASEQ/new_rnaseq/raw_data/To_1_YL_1.fq.gz",
                          "/faststorage/project/NChain/WHITE_CLOVER/20181211_RNASEQ/new_rnaseq/raw_data/To_1_YL_2.fq.gz"],
                  "YL2": ["/faststorage/project/NChain/WHITE_CLOVER/20181211_RNASEQ/new_rnaseq/raw_data/To_2_YL_1.fq.gz",
                          "/faststorage/project/NChain/WHITE_CLOVER/20181211_RNASEQ/new_rnaseq/raw_data/To_2_YL_2.fq.gz"],
                  "YL3": ["/faststorage/project/NChain/WHITE_CLOVER/20181211_RNASEQ/new_rnaseq/raw_data/To_3_YL_1.fq.gz",
                          "/faststorage/project/NChain/WHITE_CLOVER/20181211_RNASEQ/new_rnaseq/raw_data/To_3_YL_2.fq.gz"],
                  "OL1": ["/faststorage/project/NChain/WHITE_CLOVER/20181211_RNASEQ/new_rnaseq/raw_data/To_1_OL_1.fq.gz",
                          "/faststorage/project/NChain/WHITE_CLOVER/20181211_RNASEQ/new_rnaseq/raw_data/To_1_OL_2.fq.gz"],
                  "OL2": ["/faststorage/project/NChain/WHITE_CLOVER/20181211_RNASEQ/new_rnaseq/raw_data/To_2_OL_1.fq.gz",
                          "/faststorage/project/NChain/WHITE_CLOVER/20181211_RNASEQ/new_rnaseq/raw_data/To_2_OL_2.fq.gz"],
                  "OL3": ["/faststorage/project/NChain/WHITE_CLOVER/20181211_RNASEQ/new_rnaseq/raw_data/To_3_OL_1.fq.gz",
                          "/faststorage/project/NChain/WHITE_CLOVER/20181211_RNASEQ/new_rnaseq/raw_data/To_3_OL_2.fq.gz"],
                  "St1": ["/faststorage/project/NChain/WHITE_CLOVER/20181211_RNASEQ/new_rnaseq/raw_data/To_1_St_1.fq.gz",
                          "/faststorage/project/NChain/WHITE_CLOVER/20181211_RNASEQ/new_rnaseq/raw_data/To_1_St_2.fq.gz"],
                  "St2": ["/faststorage/project/NChain/WHITE_CLOVER/20181211_RNASEQ/new_rnaseq/raw_data/To_2_St_1.fq.gz",
                          "/faststorage/project/NChain/WHITE_CLOVER/20181211_RNASEQ/new_rnaseq/raw_data/To_2_St_2.fq.gz"],
                  "St3": ["/faststorage/project/NChain/WHITE_CLOVER/20181211_RNASEQ/new_rnaseq/raw_data/To_3_St_1.fq.gz",
                          "/faststorage/project/NChain/WHITE_CLOVER/20181211_RNASEQ/new_rnaseq/raw_data/To_3_St_2.fq.gz"],
                  "R1": ["/faststorage/project/NChain/WHITE_CLOVER/20181211_RNASEQ/new_rnaseq/raw_data/To_1_R_1.fq.gz",
                         "/faststorage/project/NChain/WHITE_CLOVER/20181211_RNASEQ/new_rnaseq/raw_data/To_1_R_2.fq.gz"],
                  "R2": ["/faststorage/project/NChain/WHITE_CLOVER/20181211_RNASEQ/new_rnaseq/raw_data/To_2_R_1.fq.gz",
                         "/faststorage/project/NChain/WHITE_CLOVER/20181211_RNASEQ/new_rnaseq/raw_data/To_2_R_2.fq.gz"],
                  "R3": ["/faststorage/project/NChain/WHITE_CLOVER/20181211_RNASEQ/new_rnaseq/raw_data/To_3_R_1.fq.gz",
                         "/faststorage/project/NChain/WHITE_CLOVER/20181211_RNASEQ/new_rnaseq/raw_data/To_3_R_2.fq.gz"]},
           "Tp": {"YL1": ["/faststorage/project/NChain/WHITE_CLOVER/20181211_RNASEQ/new_rnaseq/raw_data/Tp_1_YL_1.fq.gz",
                          "/faststorage/project/NChain/WHITE_CLOVER/20181211_RNASEQ/new_rnaseq/raw_data/Tp_1_YL_2.fq.gz"],
                  "YL2": ["/faststorage/project/NChain/WHITE_CLOVER/20181211_RNASEQ/new_rnaseq/raw_data/Tp_2_YL_1.fq.gz",
                          "/faststorage/project/NChain/WHITE_CLOVER/20181211_RNASEQ/new_rnaseq/raw_data/Tp_2_YL_2.fq.gz"],
                  "YL3": ["/faststorage/project/NChain/WHITE_CLOVER/20181211_RNASEQ/new_rnaseq/raw_data/Tp_3_YL_1.fq.gz",
                          "/faststorage/project/NChain/WHITE_CLOVER/20181211_RNASEQ/new_rnaseq/raw_data/Tp_3_YL_2.fq.gz"],
                  "OL1": ["/faststorage/project/NChain/WHITE_CLOVER/20181211_RNASEQ/new_rnaseq/raw_data/Tp_1_OL_1.fq.gz",
                          "/faststorage/project/NChain/WHITE_CLOVER/20181211_RNASEQ/new_rnaseq/raw_data/Tp_1_OL_2.fq.gz"],
                  "OL2": ["/faststorage/project/NChain/WHITE_CLOVER/20181211_RNASEQ/new_rnaseq/raw_data/Tp_2_OL_1.fq.gz",
                          "/faststorage/project/NChain/WHITE_CLOVER/20181211_RNASEQ/new_rnaseq/raw_data/Tp_2_OL_2.fq.gz"],
                  "OL3": ["/faststorage/project/NChain/WHITE_CLOVER/20181211_RNASEQ/new_rnaseq/raw_data/Tp_3_OL_1.fq.gz",
                          "/faststorage/project/NChain/WHITE_CLOVER/20181211_RNASEQ/new_rnaseq/raw_data/Tp_3_OL_2.fq.gz"],
                  "St1": ["/faststorage/project/NChain/WHITE_CLOVER/20181211_RNASEQ/new_rnaseq/raw_data/Tp_1_St_1.fq.gz",
                          "/faststorage/project/NChain/WHITE_CLOVER/20181211_RNASEQ/new_rnaseq/raw_data/Tp_1_St_2.fq.gz"],
                  "St2": ["/faststorage/project/NChain/WHITE_CLOVER/20181211_RNASEQ/new_rnaseq/raw_data/Tp_2_St_1.fq.gz",
                          "/faststorage/project/NChain/WHITE_CLOVER/20181211_RNASEQ/new_rnaseq/raw_data/Tp_2_St_2.fq.gz"],
                  "St3": ["/faststorage/project/NChain/WHITE_CLOVER/20181211_RNASEQ/new_rnaseq/raw_data/Tp_3_St_1.fq.gz",
                          "/faststorage/project/NChain/WHITE_CLOVER/20181211_RNASEQ/new_rnaseq/raw_data/Tp_3_St_2.fq.gz"],
                  "R1": ["/faststorage/project/NChain/WHITE_CLOVER/20181211_RNASEQ/new_rnaseq/raw_data/Tp_1_R_1.fq.gz",
                         "/faststorage/project/NChain/WHITE_CLOVER/20181211_RNASEQ/new_rnaseq/raw_data/Tp_1_R_2.fq.gz"],
                  "R2": ["/faststorage/project/NChain/WHITE_CLOVER/20181211_RNASEQ/new_rnaseq/raw_data/Tp_2_R_1.fq.gz",
                         "/faststorage/project/NChain/WHITE_CLOVER/20181211_RNASEQ/new_rnaseq/raw_data/Tp_2_R_2.fq.gz"],
                  "R3": ["/faststorage/project/NChain/WHITE_CLOVER/20181211_RNASEQ/new_rnaseq/raw_data/Tp_3_R_1.fq.gz",
                         "/faststorage/project/NChain/WHITE_CLOVER/20181211_RNASEQ/new_rnaseq/raw_data/Tp_3_R_2.fq.gz"]},
          
          'TrR': {'001': ['/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/5samples/raw_data/WC_001_1.fq.gz', '/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/5samples/raw_data/WC_001_2.fq.gz'],
                  '002': ['/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_002_1.fq.gz', '/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_002_2.fq.gz'],
                  '003': ['/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_003_1.fq.gz', '/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_003_2.fq.gz'],
                  '004': ['/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_004_1.fq.gz', '/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_004_2.fq.gz'],
                  '005': ['/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_005_1.fq.gz', '/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_005_2.fq.gz'],
                  '006': ['/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/5samples/raw_data/WC_006_1.fq.gz', '/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/5samples/raw_data/WC_006_2.fq.gz'],
                  '007': ['/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/5samples/raw_data/WC_007_1.fq.gz', '/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/5samples/raw_data/WC_007_2.fq.gz'],
                  '008': ['/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_008_1.fq.gz', '/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_008_2.fq.gz'],
                  '009': ['/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_009_1.fq.gz', '/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_009_2.fq.gz'],
                  '010': ['/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_010_1.fq.gz', '/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_010_2.fq.gz'],
                  '011': ['/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_011_1.fq.gz', '/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_011_2.fq.gz'],
                  '012': ['/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_012_1.fq.gz', '/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_012_2.fq.gz'],
                  '013': ['/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_013_1.fq.gz', '/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_013_2.fq.gz'],
                  '014': ['/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_014_1.fq.gz', '/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_014_2.fq.gz'],
                  '015': ['/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_015_1.fq.gz', '/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_015_2.fq.gz'],
                  '016': ['/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_016_1.fq.gz', '/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_016_2.fq.gz'],
                  '017': ['/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_017_1.fq.gz', '/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_017_2.fq.gz'],
                  '018': ['/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_018_1.fq.gz', '/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_018_2.fq.gz'],
                  '019': ['/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_019_1.fq.gz', '/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_019_2.fq.gz'],
                  '020': ['/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_020_1.fq.gz', '/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_020_2.fq.gz'],
                  '021': ['/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_021_1.fq.gz', '/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_021_2.fq.gz'],
                  '022': ['/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_022_1.fq.gz', '/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_022_2.fq.gz'],
                  '023': ['/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_023_1.fq.gz', '/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_023_2.fq.gz'],
                  '024': ['/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_024_1.fq.gz', '/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_024_2.fq.gz'],
                  '025': ['/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_025_1.fq.gz', '/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_025_2.fq.gz'],
                  '026': ['/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_026_1.fq.gz', '/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_026_2.fq.gz'],
                  '027': ['/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_027_1.fq.gz', '/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_027_2.fq.gz'],
                  '028': ['/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_028_1.fq.gz', '/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_028_2.fq.gz'],
                  '029': ['/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_029_1.fq.gz', '/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_029_2.fq.gz'],
                  '030': ['/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_030_1.fq.gz', '/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_030_2.fq.gz'],
                  '031': ['/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_031_1.fq.gz', '/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_031_2.fq.gz'],
                  '032': ['/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_032_1.fq.gz', '/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_032_2.fq.gz'],
                  '033': ['/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_033_1.fq.gz', '/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_033_2.fq.gz'],
                  '034': ['/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_034_1.fq.gz', '/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_034_2.fq.gz'],
                  '035': ['/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_035_1.fq.gz', '/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_035_2.fq.gz'],
                  '036': ['/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_036_1.fq.gz', '/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_036_2.fq.gz'],
                  '037': ['/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_037_1.fq.gz', '/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_037_2.fq.gz'],
                  '038': ['/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_038_1.fq.gz', '/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_038_2.fq.gz'],
                  '039': ['/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_039_1.fq.gz', '/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_039_2.fq.gz'],
                  '040': ['/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_040_1.fq.gz', '/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_040_2.fq.gz'],
                  '041': ['/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_041_1.fq.gz', '/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_041_2.fq.gz'],
                  '042': ['/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_042_1.fq.gz', '/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_042_2.fq.gz'],
                  '043': ['/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_043_1.fq.gz', '/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_043_2.fq.gz'],
                  '044': ['/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_044_1.fq.gz', '/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_044_2.fq.gz'],
                  '045': ['/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_045_1.fq.gz', '/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_045_2.fq.gz'],
                  '046': ['/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_046_1.fq.gz', '/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_046_2.fq.gz'],
                  '047': ['/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_047_1.fq.gz', '/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_047_2.fq.gz'],
                  '048': ['/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_048_1.fq.gz', '/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_048_2.fq.gz'],
                  '049': ['/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_049_1.fq.gz', '/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_049_2.fq.gz'],
                  '050': ['/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/5samples/raw_data/WC_050_1.fq.gz', '/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/5samples/raw_data/WC_050_2.fq.gz'],
                  '051': ['/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_051_1.fq.gz', '/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_051_2.fq.gz'],
                  '052': ['/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_052_1.fq.gz', '/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_052_2.fq.gz'],
                  '053': ['/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_053_1.fq.gz', '/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_053_2.fq.gz'],
                  '054': ['/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_054_1.fq.gz', '/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_054_2.fq.gz'],
                  '055': ['/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_055_1.fq.gz', '/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_055_2.fq.gz'],
                  '056': ['/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_056_1.fq.gz', '/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_056_2.fq.gz'],
                  '057': ['/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_057_1.fq.gz', '/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_057_2.fq.gz'],
                  '058': ['/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_058_1.fq.gz', '/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_058_2.fq.gz'],
                  '059': ['/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_059_1.fq.gz', '/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_059_2.fq.gz'],
                  '060': ['/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_060_1.fq.gz', '/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_060_2.fq.gz'],
                  '061': ['/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_061_1.fq.gz', '/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_061_2.fq.gz'],
                  '062': ['/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_062_1.fq.gz', '/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_062_2.fq.gz'],
                  '063': ['/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_063_1.fq.gz', '/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_063_2.fq.gz'],
                  '064': ['/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_064_1.fq.gz', '/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_064_2.fq.gz'],
                  '065': ['/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_065_1.fq.gz', '/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_065_2.fq.gz'],
                  '066': ['/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_066_1.fq.gz', '/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_066_2.fq.gz'],
                  '067': ['/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_067_1.fq.gz', '/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_067_2.fq.gz'],
                  '068': ['/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_068_1.fq.gz', '/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_068_2.fq.gz'],
                  '069': ['/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_069_1.fq.gz', '/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_069_2.fq.gz'],
                  '070': ['/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_070_1.fq.gz', '/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_070_2.fq.gz'],
                  '071': ['/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_071_1.fq.gz', '/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_071_2.fq.gz'],
                  '072': ['/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_072_1.fq.gz', '/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_072_2.fq.gz'],
                  '073': ['/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_073_1.fq.gz', '/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_073_2.fq.gz'],
                  '074': ['/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_074_1.fq.gz', '/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_074_2.fq.gz'],
                  '075': ['/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_075_1.fq.gz', '/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_075_2.fq.gz'],
                  '076': ['/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_076_1.fq.gz', '/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_076_2.fq.gz'],
                  '078': ['/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_078_1.fq.gz', '/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_078_2.fq.gz'],
                  '079': ['/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_079_1.fq.gz', '/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_079_2.fq.gz'],
                  '080': ['/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_080_1.fq.gz', '/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_080_2.fq.gz'],
                  '081': ['/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_081_1.fq.gz', '/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_081_2.fq.gz'],
                  '084': ['/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_084_1.fq.gz', '/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_084_2.fq.gz'],
                  '085': ['/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_085_1.fq.gz', '/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_085_2.fq.gz'],
                  '087': ['/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_087_1.fq.gz', '/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_087_2.fq.gz'],
                  '089': ['/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_089_1.fq.gz', '/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_089_2.fq.gz'],
                  '092': ['/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_092_1.fq.gz', '/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_092_2.fq.gz'],
                  '093': ['/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_093_1.fq.gz', '/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_093_2.fq.gz'],
                  '094': ['/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_094_1.fq.gz', '/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_094_2.fq.gz'],
                  '095': ['/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_095_1.fq.gz', '/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_095_2.fq.gz'],
                  '096': ['/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_096_1.fq.gz', '/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_096_2.fq.gz'],
                  '097': ['/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_097_1.fq.gz', '/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_097_2.fq.gz'],
                  '098': ['/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_098_1.fq.gz', '/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_098_2.fq.gz'],
                  '099': ['/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_099_1.fq.gz', '/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_099_2.fq.gz'],
                  '100': ['/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_100_1.fq.gz', '/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_100_2.fq.gz'],
                  '101': ['/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_101_1.fq.gz', '/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_101_2.fq.gz'],
                  '103': ['/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_103_1.fq.gz', '/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_103_2.fq.gz'],
                  '104': ['/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_104_1.fq.gz', '/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_104_2.fq.gz'],
                  '106': ['/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_106_1.fq.gz', '/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_106_2.fq.gz'],
                  '107': ['/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_107_1.fq.gz', '/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_107_2.fq.gz'],
                  '108': ['/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_108_1.fq.gz', '/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_108_2.fq.gz'],
                  '109': ['/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_109_1.fq.gz', '/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_109_2.fq.gz'],
                  '110': ['/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_110_1.fq.gz', '/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_110_2.fq.gz'],
                  '111': ['/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_111_1.fq.gz', '/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_111_2.fq.gz'],
                  '114': ['/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_114_1.fq.gz', '/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_114_2.fq.gz'],
                  '115': ['/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_115_1.fq.gz', '/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_115_2.fq.gz'],
                  '116': ['/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_116_1.fq.gz', '/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_116_2.fq.gz'],
                  '117': ['/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_117_1.fq.gz', '/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_117_2.fq.gz'],
                  '118': ['/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_118_1.fq.gz', '/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_118_2.fq.gz'],
                  '119': ['/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_119_1.fq.gz', '/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_119_2.fq.gz'],
                  '120': ['/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_120_1.fq.gz', '/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_120_2.fq.gz'],
                  '121': ['/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_121_1.fq.gz', '/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_121_2.fq.gz'],
                  '122': ['/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_122_1.fq.gz', '/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_122_2.fq.gz'],
                  '123': ['/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_123_1.fq.gz', '/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_123_2.fq.gz'],
                  '124': ['/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_124_1.fq.gz', '/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_124_2.fq.gz'],
                  '125': ['/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_125_1.fq.gz', '/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_125_2.fq.gz'],
                  '126': ['/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_126_1.fq.gz', '/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_126_2.fq.gz'],
                  '127': ['/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_127_1.fq.gz', '/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_127_2.fq.gz'],
                  '128': ['/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_128_1.fq.gz', '/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_128_2.fq.gz'],
                  '130': ['/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_130_1.fq.gz', '/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_130_2.fq.gz'],
                  '131': ['/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_131_1.fq.gz', '/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_131_2.fq.gz'],
                  '132': ['/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_132_1.fq.gz', '/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_132_2.fq.gz'],
                  '133': ['/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_133_1.fq.gz', '/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_133_2.fq.gz'],
                  '134': ['/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_134_1.fq.gz', '/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_134_2.fq.gz'],
                  '135': ['/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_135_1.fq.gz', '/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_135_2.fq.gz'],
                  '140': ['/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_140_1.fq.gz', '/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_140_2.fq.gz'],
                  '141': ['/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_141_1.fq.gz', '/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_141_2.fq.gz'],
                  '142': ['/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_142_1.fq.gz', '/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_142_2.fq.gz'],
                  '143': ['/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_143_1.fq.gz', '/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_143_2.fq.gz'],
                  '144': ['/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/5samples/raw_data/WC_144_1.fq.gz', '/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/5samples/raw_data/WC_144_2.fq.gz'],
                  '145': ['/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_145_1.fq.gz', '/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_145_2.fq.gz'],
                  '146': ['/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_146_1.fq.gz', '/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_146_2.fq.gz'],
                  '147': ['/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_147_1.fq.gz', '/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_147_2.fq.gz'],
                  '149': ['/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_149_1.fq.gz', '/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_149_2.fq.gz'],
                  '150': ['/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_150_1.fq.gz', '/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_150_2.fq.gz'],
                  '151': ['/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_151_1.fq.gz', '/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_151_2.fq.gz'],
                  '152': ['/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_152_1.fq.gz', '/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_152_2.fq.gz'],
                  '153': ['/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_153_1.fq.gz', '/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_153_2.fq.gz'],
                  '155': ['/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_155_1.fq.gz', '/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_155_2.fq.gz'],
                  '156': ['/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_156_1.fq.gz', '/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_156_2.fq.gz'],
                  '157': ['/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_157_1.fq.gz', '/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_157_2.fq.gz'],
                  '158': ['/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_158_1.fq.gz', '/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_158_2.fq.gz'],
                  '159': ['/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_159_1.fq.gz', '/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_159_2.fq.gz'],
                  '160': ['/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_160_1.fq.gz', '/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_160_2.fq.gz'],
                  '161': ['/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_161_1.fq.gz', '/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_161_2.fq.gz'],
                  '162': ['/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_162_1.fq.gz', '/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_162_2.fq.gz'],
                  '163': ['/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_163_1.fq.gz', '/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RAWDATA/140samples/WC_163_2.fq.gz']}}

all_species_gDNA = ["To", "Tp", "TrR"]
all_species_RNA = ["To", "Tp", "TrR"]
tissues_old = ["floral", "leaf", "stolon", "root"]

tissues = {"To": ["St1", "St2", "St3", "R1", "R2", "R3", "YL1", "YL2", "YL3", "OL1", "OL2", "OL3"],
           "Tp": ["St1", "St2", "St3", "R1", "R2", "R3", "YL1", "YL2", "YL3", "OL1", "OL2", "OL3"],
           "TrR":[       '001','002','003','004','005','006','007','008','009',
                   '010','011','012','013','014','015','016','017','018','019',
                   '020','021','022','023','024','025','026','027','028','029',
                   '030','031','032','033','034','035','036','037','038','039',
                   '040','041','042','043','044','045','046','047','048','049',
                   '050','051','052','053','054','055','056','057','058','059',
                   '060','061','062','063','064','065','066','067','068','069',
                   '070','071','072','073','074','075','076','078','079','080',
                   '081','084','085','087','089','092','093','094','095','096',
                   '097','098','099','100','101','103','104','106','107','108',
                   '109','110','111','114','115','116','117','118','119','120',
                   '121','122','123','124','125','126','127','128','130','131',
                   '132','133','134','135','140','141','142','143','144','145',
                   '146','147','149','150','151','152','153','155','156','157',
                   '158','159','160','161','162','163']}


gwf.target_from_template("ToIndex",
                         bowtie2_index(reference_file, index_ref_file))

def star_mapping(read1, output, output_prefix, genomeDir):
    inputs = [read1, genomeDir+"/genomeParameters.txt"]
    outputs = [output]
    options = {
	"cores": 8,
	"memory": "64g",
	"account": "NChain",
	"walltime": "12:00:00"}
    OFilePrefix = output_prefix

    spec = """
    source /com/extra/STAR/2.5.2b/load.sh
    STAR --runThreadN 8 --genomeDir {dir} --outSAMtype BAM SortedByCoordinate --outFileNamePrefix {of} --readFilesIn {r1} --limitGenomeGenerateRAM=64000000000""".format(r1 = read1, of=output_prefix, dir=genomeDir)

    if read1.split(".")[-1]=="gz":
        spec += " --readFilesCommand zcat"

    return inputs, outputs, options, spec

def bowtie2_mapping(read1, output, genomeDir):
    inputs = [read1, genomeDir+".1.bt2"]
    outputs = [output]
    options = {
	"cores": 8,
	"memory": "24g",
	"account": "NChain",
	"walltime": "12:00:00"}

    spec = """
    source activate HyLiTE

    bowtie2 -N 1 -x {dir} -U {reads} -S {sam} --local -p {cores}
    """.format(reads = read1, sam=output, dir=genomeDir, cores=options["cores"])

    return inputs, outputs, options, spec


for species in all_species_gDNA:
    for i, readfile in enumerate(raw_gDNA[species]):
        gwf.target_from_template(species+"gDNAmap"+str(i+1),
                                 bowtie2_mapping(readfile, "./gDNA/{s}/{s}.{i}.gDNA.sam".format(s=species, i=i+1),
                                                 index_ref_file))

for species in all_species_RNA:
    for tissue in tissues[species]:
        for i, readfile in enumerate(raw_RNA[species][tissue]):
            gwf.target_from_template(species+"_"+tissue+"_"+"RNAmap"+str(i+1),
                                     bowtie2_mapping(readfile, "./RNA/{s}/{s}.{t}.{i}.gDNA.sam".format(s=species, t=tissue, i=i+1),
                                                     index_ref_file))

def samtools_merge(infile1, infile2, outfile):
    inputs = [infile1, infile2]
    outputs = [outfile]
    options = {
        'cores': 1,
        'memory': '2g',
        'account': 'NChain',
        'walltime': '12:00:00'
    }

    spec = """
    source /com/extra/samtools/1.3/load.sh
    samtools merge {} {} {}
    """.format(outfile, infile1, infile2)

    return inputs, outputs, options, spec


for species in all_species_gDNA:
    gwf.target_from_template(species+"gDNAmerge",
                             samtools_merge("./gDNA/{s}/{s}.{i}.gDNA.sam".format(s=species, i=1),
                                            "./gDNA/{s}/{s}.{i}.gDNA.sam".format(s=species, i=2),
                                            "./gDNA/{s}/{s}.gDNA.sam".format(s=species)))

for species in all_species_RNA:
    for tissue in tissues[species]:
        gwf.target_from_template(species+"_"+tissue+"_"+"RNAmerge",
                                 samtools_merge("./RNA/{s}/{s}.{t}.{i}.gDNA.sam".format(s=species, t=tissue, i=1),
                                                "./RNA/{s}/{s}.{t}.{i}.gDNA.sam".format(s=species, t=tissue, i=2),
                                                "./RNA/{s}/{s}.{t}.gDNA.sam".format(s=species, t=tissue)))


def samtools_process(infile, outfile):
    inputs = [infile]
    outputs = [outfile]
    options = {
        'cores': 1,
        'memory': '2g',
        'account': 'NChain',
        'walltime': '12:00:00'
    }

    spec = """
    source /com/extra/samtools/1.3/load.sh

    #samtools view -Sb {infile} -o {infile}.bam
    samtools sort {infile} -o {outfile}
    samtools index {outfile}
    """.format(infile=infile, outfile=outfile)

    return inputs, outputs, options, spec

for species in all_species_gDNA:
    gwf.target_from_template(species+"gDNAsort",
                             samtools_process("./gDNA/{s}/{s}.gDNA.sam".format(s=species),
                                            "./gDNA/{s}.gDNA.s.bam".format(s=species)))

for species in all_species_RNA:
    for tissue in tissues[species]:
        gwf.target_from_template(species+"_"+tissue+"_"+"RNAsort",
                                 samtools_process("./RNA/{s}/{s}.{t}.gDNA.sam".format(s=species, t=tissue),
                                                "./RNA/{s}.{t}.RNA.s.bam".format(s=species, t=tissue)))

def make_protocol_file(in_files, individual, samples, ploidy, filetype, outfile):
    inputs = in_files
    outputs = [outfile]
    options = {
        'cores': 1,
        'memory': '1g',
        'account': 'NChain',
        'walltime': '01:00:00'
    }

    spec = ''
    for inf, ind, sample, p, f in zip(in_files, individual, samples, ploidy, filetype):
        spec += 'echo "'+"\t".join([ind, p, sample, f, inf])+'" >> '+outfile+"\n"

    print(spec)

    return inputs, outputs, options, spec


Name = {"To": "P1", "Tp": "P2", "TrR": "Ch"}
Ploidy = {"To": "1", "Tp": "1", "TrR": "2"}

bamfiles = []
individual = []
ploidy_levels = []
seq_type = []

#### MAKE SAMPLES LIST
samples = []


translate = {"To": "To", "Tp": "Tp",
             "TrR": "TrR"}

samples_map = {"To": [["St1", "R1", "YL1", "OL1"],
                      ["St2", "R2", "YL2", "OL2"],
                      ["St3", "R3", "YL3", "OL3"]],
               "Tp": [["St1", "R1", "YL1", "OL1"],
                      ["St2", "R2", "YL2", "OL2"],
                      ["St3", "R3", "YL3", "OL3"]],
               "TrR": [[     '001','002','003','004','005','006','007','008','009',
                       '010','011','012','013','014','015','016','017','018','019',
                       '020','021','022','023','024','025','026','027','028','029',
                       '030','031','032','033','034','035','036','037','038','039',
                       '040','041','042','043','044','045','046','047','048','049',
                       '050','051','052','053','054','055','056','057','058','059',
                       '060','061','062','063','064','065','066','067','068','069',
                       '070','071','072','073','074','075','076','078','079','080',
                       '081','084','085','087','089','092','093','094','095','096',
                       '097','098','099','100','101','103','104','106','107','108',
                       '109','110','111','114','115','116','117','118','119','120',
                       '121','122','123','124','125','126','127','128','130','131',
                       '132','133','134','135','140','141','142','143','144','145',
                       '146','147','149','150','151','152','153','155','156','157',
                       '158','159','160','161','162','163']]}

#species_included = {"To": False, "Tp": False, "TrR": False}

for species in ["TrR", "To", "Tp"]:
      sample_n = 1
      for current_samples in samples_map[species]:
          for tissue in current_samples:
              bamfiles.append("./RNA/{s}.{t}.RNA.s.bam".format(s=species, t=tissue))
              seq_type.append("RNAseq")
              individual.append(Name[translate[species]])
              ploidy_levels.append(Ploidy[translate[species]])
              samples.append("sample"+str(sample_n))
              sample_n += 1
      bamfiles.append("./gDNA/{s}.gDNA.s.bam".format(s=translate[species]))
      seq_type.append("gDNA")
      individual.append(Name[translate[species]])
      ploidy_levels.append(Ploidy[translate[species]])
      samples.append("sample"+str(sample_n))


#print(bamfiles)
#print(individual)
#print(ploidy_levels)
#print(seq_type)

gwf.target_from_template("HyLiTEprotocolfile",
                         make_protocol_file(bamfiles,
                                            individual,
                                            samples,
                                            ploidy_levels,
                                            seq_type,
                                            "repens_protocol.txt"))

def make_bam_file(in_files,outfile):
    inputs = in_files
    outputs = [outfile]
    options = {
        'cores': 1,
        'memory': '1g',
        'account': 'NChain',
        'walltime': '01:00:00'
    }

    spec = ""
    for ind in in_files:
        spec += 'echo "'+ind+'" >> '+outfile+"\n"

    return inputs, outputs, options, spec

gwf.target_from_template("mpileup_bamfile",
                         make_bam_file(bamfiles, "bamfiles.txt"))

def mpileup(reference, bamfiles, outfile):
    inputs = [reference, bamfiles]
    outputs = [outfile]
    options = {
        'cores': 1,
        'memory': '36g',
        'account': 'NChain',
        'walltime': '72:00:00'
    }

    spec = '''
    source /com/extra/samtools/1.3/load.sh

    samtools mpileup -BQ 0 -d 1000000 -f {ref} -b {bamfiles} > {outfile}
    '''.format(ref=reference, bamfiles=bamfiles, outfile=outfile)

    return inputs, outputs, options, spec


gwf.target_from_template("mpileup",
                         mpileup(reference_file, "bamfiles.txt",
                                 "repens_combo.pileup"))


def fix_mpileup(pileup_file, outfile):
    inputs = [pileup_file]
    outputs = [outfile]
    options = {
        'cores': 1,
        'memory': '2g',
        'account': 'NChain',
        'walltime': '12:00:00'
    }

    spec = '''
    python mpileupfix.py {pileup} > {out}
    '''.format(pileup=pileup_file, out=outfile)

    return inputs, outputs, options, spec

gwf.target_from_template("mpileupfix",
                         fix_mpileup("repens_combo.pileup", "repens_fixed.pileup"))


def run_HyLiTE(reference, protocol, pileup_file):
    inputs = [reference, protocol, pileup_file]
    outputs = ["HyLiTE_output/HyLiTE_output.expression.txt"]
    options = {
        'cores': 1,
        'memory': '124g',
        'account': 'NChain',
        'walltime': '64:00:00'
    }

    spec = '''
    source activate HyLiTE

    HyLiTE -v -f {proto} -p {pileup} -n HyLiTE_output -r {ref}
    '''.format(ref=reference, proto=protocol, pileup=pileup_file)

    return inputs, outputs, options, spec


def HacknHyLiTE(reference, protocol, pileup_file, segments):
    inputs = [reference, protocol, pileup_file]
    outputs = []
    options = {
        'cores': 1,
        'memory': '12g',
        'account': 'NChain',
        'walltime': '48:00:00'
    }

    for i in range(segments):
        outputs.append("./new_subsets/slicing.subset{}.sh".format(i))

    spec = '''
    source activate HyLiTE

    HacknHyLiTE -n {seg} -o new_subsets --name slicing -f {proto} -p {pileup} --options "-r {ref}"
    '''.format(seg=segments, ref=reference, proto=protocol, pileup=pileup_file)

    return inputs, outputs, options, spec


gwf.target_from_template("HyLiTE",
                         run_HyLiTE(reference_file, "repens_protocol.txt", "repens_fixed.pileup"))


n_segments = 80

gwf.target_from_template("HacknHyLiTE",
                         HacknHyLiTE(reference_file, "repens_protocol.txt", "repens_fixed.pileup",
                                     n_segments))

subsets = []
subsets_out = []
for i in range(n_segments):
    subsets.append("./new_subsets/slicing.subset{}.sh".format(i))
    subsets_out.append("./new_subsets/subset{}/slicing.subset{}.expression.txt".format(i, i))


def run_subset(subset, output):
    inputs = [subset]
    outputs = [output]
    options = {
        'cores': 1,
        'memory': '1024g',
        'account': 'NChain',
        'walltime': '192:00:00',
        'queue': 'fat2'
    }

    spec = '''
    source activate HyLiTE

    bash {subset}
    '''.format(subset=subset)

    return inputs, outputs, options, spec


for i, subset in enumerate(zip(subsets, subsets_out)):
    gwf.target_from_template("Subset"+str(i),
                             run_subset(subset[0], subset[1]))

def mergeHyLiTE(subset_out, directory, reference):
    inputs = subset_out
    outputs = []
    options = {
        'cores': 1,
        'memory': '12g',
        'account': 'NChain',
        'walltime': '48:00:00'
    }

    spec = '''
    source activate HyLiTE

    cd {dir}

    HyLiTE-merge -r {ref}
    '''.format(dir=directory, ref=reference)

    return inputs, outputs, options, spec


#gwf.target_from_template("mergeHyLiTE",
#                         mergeHyLiTE(subsets_out, "./slicing_directories/", "."+reference_file))
