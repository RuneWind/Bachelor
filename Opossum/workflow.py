from gwf import Workflow
import os

gwf = Workflow()

def generate_MD_tags(bam_file, reference, MD_bam):
	inputs = [bam_file, reference]
	outputs = [MD_bam]
	options = {"cores": 1,
              "memory": "2g",
              "account": "NChain",
              "walltime": "12:00:00"}

	spec = """
   source activate STAR_clover
   samtools calmd {bam_file} {ref_file} > {out}
   """.format(bam_file=bam_file, ref_file=reference, out=MD_bam)

	return inputs, outputs, options, spec

def opossum_preprocessing(bam_file, processed_bam):
	inputs = [bam_file]
	outputs = [processed_bam]
	options = {"cores": 1,
               "memory": "4g",
               "account": "NChain",
               "walltime": "12:00:00"}

	spec = """
	source activate platypus
	opossum --BamFile={inFile} --SoftClipsExist=True --OutFile={outFile}
   """.format(inFile=bam_file, outFile=processed_bam)

	return inputs, outputs, options, spec

################################################################################################################################################
bams = os.listdir(path="/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/MAPPED_READS/S9")

f = open("/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/FIXED_READS/bam.list", "w+")
pwd = "/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/"
bam_list = []
for files in bams:
	if files.endswith(".bam"):
		gwf.target_from_template("MD_tag_"+files[:8], generate_MD_tags(pwd+"MAPPED_READS/S9/"+files, pwd+"REFERENCE/TrR.v5.fasta", pwd+"FIXED_READS/"+files[:8]+"_MD.bam"))
		gwf.target_from_template("opossum_"+files[:8], opossum_preprocessing(pwd+"FIXED_READS/"+files[:8]+"_MD.bam", pwd+"FIXED_READS/"+files[:8]+".bam"))
		f.write(pwd+"FIXED_READS/"+files[:8]+".bam\n")
f.close()
