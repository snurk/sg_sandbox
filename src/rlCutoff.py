#! /usr/bin/env python

__author__ = "Brandon Pickett"

# ----------- IMPORTS ---------------------------- ||
import sys
import argparse
import gzip
import re

# ----------- CLASSES ---------------------------- ||
class SeqRecord():
	
	def __init__(self, header='', seq='', plus_header='', quals=''):
		self.header = header
		self.seq = seq
		self.plus_header = plus_header
		self.quals = quals
	
# ---------- FUNCTIONS --------------------------- ||
def parseArgs():

	parser = argparse.ArgumentParser(prog="rlCutoff", description="Shorten HiFi sequences longer than a cutoff, optionally compress homopolymers", add_help=True)

	parser.add_argument("-c", "--cutoff", metavar="INT", type=int, action="store", dest="cutoff", help="Maximum read length cutoff. All reads longer than this value will be shorted to that cutoff. The remaining bases are discarded. Cutoff value is assumed to be in uncompressed space, even when -p/--hpc/--homopolymer-compress is used; i.e., the cutoff is applied before the homopolymer compression. [default: do not apply a cutoff]", default=0, required=False)
	parser.add_argument("-f", "--fold", metavar="INT", type=int, action="store", dest="fold", help="When outputting in FASTA format, fold the sequence line every INT nucleotides. Common choices are 60 and 80. [default: do not fold fasta sequence]", default=0, required=False)
	parser.add_argument("-i", "--input", metavar="FILE", type=str, action="store", dest="input", help="Input HiFi reads file in FASTA or FASTQ format, optionally gzipped. Format and compression are determined based on filename.", required=True)
	parser.add_argument("-o", "--output", metavar="FILE", type=str, action="store", dest="output", help="Output HiFi reads file in FASTA or FASTQ format, optionally gzipped. Format and compression are determined based on filename.", required=True)
	parser.add_argument("-p", "--hpc", "--homopolymer-compress", action="store_true", dest="hpc", help="Homopolymer compress (hpc) the output sequences. If the input sequences are already hpc-ed, do not specify this option because it will only waste time. [default: no hpc]", required=False)

	args = parser.parse_args()

	#if ( not re.search(r"\.f(?:ast)?a(?:\.gz)?$", args.input) is None ) and ( not re.search(r"\.f(?:ast)?q(?:\.gz)?$", args.output) is None ):
	if not ( re.search(r"\.f(?:ast)?a(?:\.gz)?$", args.input) is None or re.search(r"\.f(?:ast)?q(?:\.gz)?$", args.output) is None ):
		print("ERROR: Cannot create output FASTQ file from input FASTA file because no quality values exist to copy over.", file=sys.stderr)
		sys.exit(1)
	
	if not re.search(r"\.f(?:ast)?q(?:\.gz)?$", args.output) is None: # if output is FASTQ
		if args.fold > 0: # if folding is requested
			print(f"ERROR: Cannot fold output FASTQ sequences every {args.fold} bases because FASTQ sequences can be on only a single line.", file=sys.stderr)
			sys.exit(1)
		if args.hpc: # if homopolymer compression is requested
			print(f"ERROR: homopolymer compression is unsupported for FASTQ-formatted output files at this time. We haven't implemented a way to handle the quality values in this situation.", file=sys.stderr)
			sys.exit(1)

	return args

def fasta_generator(fd):
	line = fd.readline()
	while line != '':
		header = line.rstrip('\n')[1:]
		seq = []
		line = fd.readline()
		while line != '' and line[0] != '>':
			seq.append(line.rstrip('\n'))
			line = fd.readline()
		yield SeqRecord(header=header, seq=''.join(seq))
		
def fastq_generator(fd):
	line = fd.readline()
	while line != '':
		sr = SeqRecord(header=line.rstrip('\n')[1:])
		sr.seq = ifd.readline().rstrip('\n')
		sr.plus_header = ifd.readline().rstrip('\n')
		sr.quals = ifd.readline().rstrip('\n')
		yield sr
		line = ifd.readline()

def homopolymerCompress(seq):
	seq += '$'
	hpc = []
	i = 0
	c = seq[i]
	i = 1
	while i < len(seq):
		while i < len(seq)-1 and seq[i] == c:
			i += 1
		hpc.append(c)
		c = seq[i]
		i += 1
	return ''.join(hpc)

# ------------- MAIN ----------------------------- ||
if __name__ == "__main__":
	
	# handle the arguments
	args = parseArgs()

	# determine whether to use gunzip on input
	iopencmd = open
	iflag = 'r'
	if not re.search(r"\.gz$", args.input) is None:
		iopencmd = gzip.open
		iflag = "rt"

	# determine whether to use gzip on output
	oopencmd = open
	oflag = 'w'
	if not re.search(r"\.gz$", args.output) is None:
		oopencmd = gzip.open
		oflag = "rt"
	
	# determine whether files are fq and how to parse them
	ifq = not re.search(r"\.f(?:ast)?q(?:\.gz)?$", args.input) is None
	ofq = not re.search(r"\.f(?:ast)?q(?:\.gz)?$", args.output) is None
	record_generator = fastq_generator if ifq else fasta_generator
	header_char = '@' if ofq else '>'

	# actually open the files and do what was requested
	with oopencmd(args.output, oflag) as ofd:
		with iopencmd(args.input, iflag) as ifd:

			# loop through each record in the input file
			for record in record_generator(ifd):
				
				# apply cutoff (if requested)
				if args.cutoff > 0 and len(record.seq) > args.cutoff:
					record.seq = record.seq[:args.cutoff]
					record.quals = record.quals[:args.cutoff]

				# compress homopolymers (if requested)
				if args.hpc:
					record.seq = homopolymerCompress(record.seq)

				# fold fasta output (if requested)
				if not ofq and args.fold > 0:
					record.seq = '\n'.join(record.seq[i:i+args.fold] for i in range(0, len(record.seq), args.fold))

				# write the output
				ofd.write(f"{header_char}{record.header}\n") # write header
				ofd.write(f"{record.seq}\n") # write seq
				if ofq:
					ofd.write(f"{record.plus_header}\n") # write qual header
					ofd.write(f"{record.quals}\n") # write quals

