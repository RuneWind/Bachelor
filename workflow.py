from gwf import Workflow
import os

gwf = Workflow()

#################################################################### STAR
def star_index(reference, output):
	inputs = [reference]
	outputs = [output]
	options = {"cores":8, "memory":"16g", "account":"NChain", "walltime": "12:00:00"}

	spec = """
	source /com/extra/STAR/2.5.2b/load.sh
	STAR --runMode genomeGenerate --runThreadN 8 --genomeDir ./REFERENCE --genomeFastaFiles {ref}
	""".format(ref=reference)

	return inputs, outputs, options, spec


def star_mapping(read1, read2, output, identifier):
	inputs = [read1, read2]
	outputs = [output]
	options = {
				"cores": 4,
				"memory": "16g",
				"account": "NChain",
				"walltime": "12:00:00"}
	OFilePrefix = "./MAPPED_READS/S9/" + identifier[:8]

	spec = """
	source /com/extra/STAR/2.5.2b/load.sh
	STAR --runThreadN 8 --genomeDir ./REFERENCE --readFilesCommand zcat --outFilterMultimapNmax 1 --outSAMtype BAM SortedByCoordinate --outSAMunmapped Within KeepPairs --outFileNamePrefix {of} --readFilesIn {r1} {r2}
	""".format(r1 = read1, r2=read2, of=OFilePrefix)

	return inputs, outputs, options, spec


def get_stats(input_file, output_file):
	outputs = [output_file]
	inputs = [input_file]
	options = {
				"cores": 2,
				"memory": "8g",
				"account": "NChain",
				"walltime": "12:00:00"}

	spec = """
	source /com/extra/STAR/2.5.2b/load.sh
	samtools flagstat {file1} > {file2}""".format(file1=input_file, file2=output_file)

	return inputs, outputs, options, spec



def count_matrix(input_files, annotation_file, output_file):
	inputs = []
	outputs = [output_file]
	options = {
				"cores": 8,
				"memory": "4g",
				"account": "NChain",
				"walltime": "12:00:00"}

	spec = """
	source activate STAR_clover
	featureCounts -a {a_file} -o {o_file} -T {cores} -t {feature} -g {attribute}""".format(a_file=annotation_file, o_file=output_file, cores=options["cores"], feature="CDS", attribute="transcript_id")

	for bam_file in input_files:
		spec += " ./MAPPED_READS/{}".format(bam_file)
		inputs.append("./MAPPED_READS/"+bam_file)

	return inputs, outputs, options, spec



############################################################################################################################################################
gwf.target_from_template("star_index",
							star_index("./REFERENCE/TrR.v5.fasta", "REFERENCE/genomeParameters.txt"))



read_files = sorted(os.listdir(path="/faststorage/project/NChain/WHITE_CLOVER/20181211_RNASEQ/new_rnaseq/raw_data"))

raw_read_lst = []
for file in read_files:
	if file.startswith("S9"):
		raw_read_lst.append(file)

for i in range(0,len(raw_read_lst)-1,2):
	split_str = raw_read_lst[i].split("_")
	if len(split_str[1]) == 2:
		split_str[1] = "0" + split_str[1]
	if len(split_str[2]) == 2:
		split_str[2] = split_str[2][0]
	identifier = "_".join(split_str)
	
	frst_rd = "/faststorage/project/NChain/WHITE_CLOVER/20181211_RNASEQ/new_rnaseq/raw_data/" + raw_read_lst[i]
	scnd_rd = "/faststorage/project/NChain/WHITE_CLOVER/20181211_RNASEQ/new_rnaseq/raw_data/" + raw_read_lst[i+1]
	target_name = "star_mapping_" + identifier[:8]
	output_name = "/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/MAPPED_READS/S9/" + identifier[:8] + "Aligned.sortedByCoord.out.bam"
	gwf.target_from_template(target_name,
							star_mapping(frst_rd, scnd_rd, output_name, identifier))



bam_lst = []

for file in sorted(os.listdir(path="/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/MAPPED_READS/S9/"))[:-5]:
	if file.endswith(".bam"):
		bam_lst.append(file)

for file in bam_lst:
	bam_file = "/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/MAPPED_READS/" + file
	stat_file = "/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/MAPPED_READS/" + file + "_stats.txt"
	target_name = "get_stats_" + file[:6]
	gwf.target_from_template(target_name, 
									get_stats(bam_file, stat_file))



gwf.target_from_template("count_features",
							count_matrix(bam_lst, "/faststorage/project/NChain//WHITE_CLOVER/BRAKER_ANNOTATION/pipeline/hrd/annotated/TrR.v5.complete.gtf", "./MAPPED_READS/featureCounts.txt"))




