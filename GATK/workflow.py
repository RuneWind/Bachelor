from gwf import Workflow
import os

gwf = Workflow()

def add_readgroups(bam, rg_bam):
     inputs = [bam]
     outputs = [rg_bam]
     options = {"cores": 1,
               "memory": "4g",
               "account": "NChain",
               "walltime": "12:00:00"}

     spec = """qx --no-scratch ./picard.sh {bam} {rg_bam} {rg}
     """.format(bam=bam, rg_bam=rg_bam, rg="S9")

     return inputs, outputs, options, spec


def samtools_index(bam_file, indexed_file):
     inputs = [bam_file]
     outputs = [indexed_file]
     options = {"cores": 1,
               "memory": "4g",
               "account": "NChain",
               "walltime": "12:00:00"}

     spec = """
     source activate STAR_clover
     samtools index {input}""".format(input=bam_file)

     return inputs, outputs, options, spec


def unifiedgenotyper(bam_list, output):
     inputs = [bam_list]
     outputs = [output]
     options = {"cores": 12,
               "memory": "32g",
               "account": "NChain",
               "walltime": "72:00:00"}

     spec = """
     ./Unifiedgenotyper.sh {} {}
     """.format(bam_list, output)

     return inputs, outputs, options, spec


###############################################################################################################################################################################
bams = os.listdir(path="/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RG_BAM")

bam_list = []
for files in bams:
     if files.startswith("WC") and files.endswith(".bam"):
          bam_list.append(files)

f = open("/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RG_BAM/bam.list", "w+")
for file in bam_list:
     rg_bam_name = file[:6] + "_RG.bam"
     gwf.target_from_template("rg_"+file[:8], add_readgroups("/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RG_BAM/"+file, "/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RG_BAM/"+rg_bam_name))
     gwf.target_from_template("index_"+file[:8], samtools_index("/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RG_BAM/"+rg_bam_name, "/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RG_BAM/"+rg_bam_name+".bai"))
     f.write("/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RG_BAM/" + file + "\n")
f.close()


gwf.target_from_template("GATK_SNP", unifiedgenotyper("/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/RG_BAM/bam.list", "/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/VARIANT_CALLING/GATK/WC_ALL_SNP.vcf"))


