### trinotate_table_to_tbl.py
### Author: Moran Neuhof
### ---------------------------
### The script converts a table produced by Trinotate (https://github.com/Trinotate/Trinotate.github.io/wiki)
### To the .tbl format requires for transcriptome upload on NCBI website.

### This script performs the conversion based on a long series of trial and error experiments.
### It is tailored to the need of one specific project, uploaded in December 2018.
### As the NCBI submission requirements change from time to time, and the requirements of the 
### user's submissions vary, this script may not work for you or may not do exactly what you want it to do.
### In that case, please use it as inspiration, or copy any useful part of it for your own use.

### I appreciate comments, updates and imprevements to this script any time.

### ---------------------------



import argparse

# converting trinotate table to feature table

# In the end, the structure and feature I would want for transcriptome annotation:
# header line, for each sequence: 
# >Feature (sequence name)
# next line:
# start position	end position	"CDS"
# next line:
# 			product	[annotation]
#			note	"Hypothetical protein"

# the trinotate outfile header:
#gene_id	transcript_id	sprot_Top_BLASTX_hit	RNAMMER	prot_id	prot_coords	sprot_Top_BLASTP_hit	Pfam	SignalP	TmHMM	eggnog	Kegg	gene_ontology_blast	gene_ontology_pfam	transcript	peptide

# write header line: >Feature\t{transcript_id}


def process_blast_hit(blast_hit_line, blast_type):
	"""Receieves a blast hit from trinotate and parses it"""
	# cut and split
	split_blast_hit = blast_hit_line[:blast_hit_line.find(';')].split('^')
	protein_code = split_blast_hit[0]
	st_pos, end_pos = split_blast_hit[2].split(',')[0][2:].split('-')
	record_name = split_blast_hit[5][split_blast_hit[5].find('=')+1:].replace('.', '').replace('_', '-')  # shorter
	# removing brackets
	bracket_pos = record_name.find('{')
	if bracket_pos != -1:
		record_name = record_name[:bracket_pos - 1]
	
	# record_name = split_blast_hit[5][split_blast_hit[5].find('=')+1:]

	annotation = f"{record_name} ({protein_code}) (Trinotate prediction using {blast_type}: {protein_code}, {';'.join(split_blast_hit[2:5])})"
	# annotation = f"{record_name} ({protein_code}) (Trinotate prediction using BLAST: {protein_code}, {';'.join(split_blast_hit[2:5])})"
	annotation_short = f"{record_name} (Trinotate prediction)"

	return annotation, annotation_short, st_pos, end_pos
	# return split_blast_hit, protein_code, st_pos, end_pos, record_name

def replace_positions(st_coords, end_coords):
	"""Replace the order of positions, if needed"""
	if end_coords < st_coords:  # end smaller than start
		return str(end_coords), str(st_coords)  # replace
		# return end_coords, st_coords  # replace
	return str(st_coords), str(end_coords)  # 


def read_seq_len(seq_len_file):
	with open(seq_len_file, 'r') as in_file:
		d = {x[0]: x[1] for x in (line.strip().split('\t') for line in in_file)}
	return d


def read_seqs_to_ignore(seqs_to_ignor_file):
	with open(seqs_to_ignor_file, 'r') as in_file:
		d = {x: True for x in (line.strip() for line in in_file)}
	return d

def get_partial_seq(partial_seq_file):
	d = {}
	with open(partial_seq_file, 'r') as in_file:
		for line in in_file:
			if line.startswith('>'):
				seq_type, no_use_for_this_one, seq_name_pos = line.rstrip().split(' ')[1:]
				seq_type = seq_type[5:]
				seq_name, seq_pos = seq_name_pos.split(':')
				# make sure there's no .1 in the end of the name
				seq_name = seq_name.replace('.1', '')

				st_pos, end_pos = seq_pos.split('(')[0].split('-')  # str
				try:
					d[seq_name][(st_pos, end_pos)] = seq_type  # if a value already exists for this one
				except KeyError:
					d[seq_name] = {(st_pos, end_pos): seq_type}  # set up a pos values for this seq
	return d

def fix_partial_coords(seq_type, st_coords, end_coords):

	if seq_type == '3prime_partial':
		end_coords = f'>{end_coords}'
	elif seq_type == '5prime_partial':
		st_coords = f'<{st_coords}'
	elif seq_type == 'internal':
		st_coords = f'<{st_coords}'
		end_coords = f'>{end_coords}'

	# if complete, do nothing
	return st_coords, end_coords


def check_if_has_mRNA(transcript_id, has_mRNA_dict):
	try:
		return has_mRNA_dict[transcript_id], has_mRNA_dict  # already has a transcript
	except KeyError:
		has_mRNA_dict[transcript_id] = True
		return False, has_mRNA_dict


# parse file
parser = argparse.ArgumentParser()
parser.add_argument('input_file', metavar='Trinotate.txt', type=str, help='The trinotate annotation file')
parser.add_argument('seq_len_file', help='A file containing the sequence lengths', type=str)
parser.add_argument('--seqs_to_ignore_file', help='A file containing a list of sequences to ignore when creating the table', type=str)
parser.add_argument('--transdecoder_file', help='A fasta file containing the transdecoder .pep/.cds output', type=str)

args = parser.parse_args()


if args.seqs_to_ignore_file:
	seqs_to_ignore_dict = read_seqs_to_ignore(args.seqs_to_ignore_file)
else:
	seqs_to_ignore_dict = {}

if args.transdecoder_file:
	partial_seq_dict = get_partial_seq(args.transdecoder_file)
else:
	partial_seq_dict = {}


input_file = args.input_file

out_suffix = input_file.rsplit('.')[-1]
output_file = input_file.replace(out_suffix, 'tbl')

seq_len_dict = read_seq_len(args.seq_len_file)
skipped = 0
has_mRNA_dict = {}  # will store if the transcript has already had an mRNA feature

with open(output_file, 'w') as outfile:
	with open(input_file) as infile:  # tab delimited file
		next(infile)  # skip header
		for split_line in (line.split('\t') for line in infile):
			
			# unpacking and processing line

			transcript_id = split_line[1]  # check if we actually need to ignore this sequence
			try:
				if seqs_to_ignore_dict[transcript_id]:  # if we need to ignore it
					skipped += 1
					continue
			except KeyError:  # no need to ignore, not in the list
				pass

			sprot_Top_BLASTX_hit = split_line[2]  
			sprot_Top_BLASTP_hit = split_line[6]

			# try the BLASTP first
			if sprot_Top_BLASTP_hit != '.':  
				# try blastp
				annotation, annotation_short, hit_st_pos, hit_end_pos = process_blast_hit(sprot_Top_BLASTP_hit, 'BLASTP')
				try:  # use transdecoder CDS annotation
					st_coords, end_coords = split_line[5].split('-')
					end_coords = end_coords.split("[")[0]
					prediction_text = "predicted by Transdecoder"  
					feature_name = "CDS"
				except ValueError:  # CDS unavailble, then don't use CDS
					print(f"{transcript_id}:\tBLASTP hit, no CDS")
					
					tran_has_mRNA, has_mRNA_dict = check_if_has_mRNA(transcript_id, has_mRNA_dict)
					if tran_has_mRNA:
						print(f"{transcript_id}:\tAlready has an mRNA")
						continue

					feature_name = "mRNA"
					# feature_name = "source"  # lets try this one
					st_coords = 1
					end_coords = seq_len_dict[transcript_id]  # all seq len
					prediction_text = "based on blastp hit"  # maybe ignore these lines!!! 

					# use "mRNA" feature

			else:  # no blastp result
				if sprot_Top_BLASTX_hit == '.':  # no protein annotation at all
					continue  # skipping, no blastp or blastx result
				else:
					# check if it has mRNA
					tran_has_mRNA, has_mRNA_dict = check_if_has_mRNA(transcript_id, has_mRNA_dict)
					if tran_has_mRNA:
						print(f"{transcript_id}:\tAlready has an mRNA")
						continue

					# process
					annotation, annotation_short, hit_st_pos, hit_end_pos = process_blast_hit(sprot_Top_BLASTX_hit, 'BLASTX')
					# feature_name = "source"  # lets try this one
					feature_name = "mRNA"
					st_coords = 1
					end_coords = seq_len_dict[transcript_id]  # all seq len
					prediction_text = "based on blastx hit"  # maybe ignore these lines!!! 	
					print(f"{transcript_id}:\tNo BLASTP hit. BLASTX hit")

			
			# make sure these are correct
			st_coords, end_coords = replace_positions(int(st_coords), int(end_coords))
			# now, check if the sequence is partial:
			try:
				seq_type = partial_seq_dict[transcript_id][(st_coords, end_coords)]
			except KeyError:
				seq_type = 'complete'

			st_coords, end_coords = fix_partial_coords(seq_type, st_coords, end_coords)


			# creating feature table record:
			print(f">Feature\t{transcript_id}", file=outfile)  # first row
			print(f"{st_coords}\t{end_coords}\t{feature_name}", file=outfile)  # second row
			
			print(f"\t\t\tproduct\t{annotation_short}", file=outfile)  # short version
			print(f"\t\t\tnote\tHypothetical protein: {annotation}", file=outfile)  # short version

print(f"Skipped {skipped} sequences from list")
print("Done")