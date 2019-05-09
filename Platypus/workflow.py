from gwf import Workflow
import os

gwf = Workflow()


def platypus_variant(bam_list, reference, output):
	inputs = [bam_list, reference]
	outputs = [output]
	options = {"cores": 12,
				"memory": "64g",
				"account": "NChain",
				"walltime": "96:00:00"}

	spec = """
	source activate platypus
	platypus callVariants --bamFiles={bam} --refFile={ref} --output={op} --nCPU=12 --filterDuplicates=0 --minMapQual=0 --minGoodQualBases=10 --minBaseQual=20
	""".format(bam=bam_list, ref=reference, op=output)

	return inputs, outputs, options, spec


bams = os.listdir(path="/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/FIXED_READS")

bam_list = []
for files in bams:
     if files.startswith("WC") and files.endswith(".bam"):
          bam_list.append(files)

f = open("/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/FIXED_READS/bam.list", "w+")
for file in bam_list:
     f.write("/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/FIXED_READS/"+file+"\n")
f.close()


gwf.target_from_template("platypus_variant", platypus_variant("/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/FIXED_READS/bam.list", "/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/REFERENCE/TrR.v5.fasta", "./All_Variants.vcf"))