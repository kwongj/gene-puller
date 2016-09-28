#!/usr/bin/env python
# Script by Jason Kwong
# Retrieve genes and construct alignemnt

# Use modern print function from python 3.x
from __future__ import print_function

# Modules
import argparse
from argparse import RawTextHelpFormatter
import sys
import os
import StringIO
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast import NCBIXML
from Bio.Align.Applications import MuscleCommandline
from Bio import AlignIO
import subprocess
from subprocess import Popen

# Standard functions
# Log a message to stderr
def msg(*args, **kwargs):
	print(*args, file=sys.stderr, **kwargs)

# Log an error to stderr and quit with non-zero error code
def err(*args, **kwargs):
	msg(*args, **kwargs)
	sys.exit(1);

# Other functions
# Check files in FASTA format
def facheck(f):
	if os.path.isfile(f) == False:
		msg('ERROR: Cannot find "{}". Check file exists.'.format(f))
		return 1
	s = open(f, 'r')
	if s.read(1) != '>':
		msg('ERROR: "{}" does not appear to be in FASTA format.'.format(f))
		return 1
	s.close()

# Get sequences from FASTA files
def seqBLAST(f):
	msg('Extracting sequences from {} ...'.format(f))
	result = []
	for seq in SeqIO.parse(args.genes[0], 'fasta'):
		query = str(seq.seq)
		fBLAST = NcbiblastnCommandline(task='blastn', subject=f, evalue=1e-10, perc_identity=90, outfmt=5, dust='no', qcov_hsp_perc=80, culling_limit=1)
		stdout, stderr = fBLAST(stdin=query)
		xmlFILE = StringIO.StringIO()
		xmlFILE.write(stdout)
		xmlFILE.seek(0)
		blastREC = NCBIXML.read(xmlFILE)
		for aln in blastREC.alignments:
			for hsp in aln.hsps:
				result.append(hsp.sbjct)
	return result		

# Get IDs and sequences from SeqRecord
def seqPARSE(file):
	ids = []
	seqs = []
	for seq in SeqIO.parse(file, 'fasta'):
		ids.append(seq.id)
		seqs.append(str(seq.seq))
	return ids, seqs

# Usage
parser = argparse.ArgumentParser(
	formatter_class=RawTextHelpFormatter,
	description='Searches FASTA files for specified query sequences and outputs ClustalW alignment',
	usage='\n  %(prog)s --genes FASTA [OPTIONS] FASTA1 FASTA2 FASTA3 ... FASTAN')
parser.add_argument('fasta', metavar='FASTA', nargs='+', help='FASTA file to search (required)')
parser.add_argument('--genes', metavar='FASTA', required=True, nargs=1, help='File of query genes in FASTA format (required)')
parser.add_argument('--out', metavar='FILE', nargs=1, help='Output file')
parser.add_argument('--version', action='version', version=
	'=====================================\n'
	'%(prog)s v0.1\n'
	'Updated 27-Sep-2016 by Jason Kwong\n'
	'Dependencies: Python 2.x, BioPython, BLAST\n'
	'=====================================')
args = parser.parse_args()

# Set buffer size
buffer = 'N'*300

# Get sequences
geneSEQS = {}
for f in args.fasta:
	f_id = os.path.basename(f)
	f_id = '>' + f_id
	geneLIST = seqBLAST(f)
	geneSEQ = buffer.join(geneLIST)
	geneSEQS[f_id] = geneSEQ

# Add reference sequences
REF = seqPARSE(args.genes[0])
ref_id = '>' + '_'.join(REF[0])
ref_seq = buffer.join(REF[1])
geneSEQS[ref_id] = ref_seq

# Join sequences into string
allSEQS = '\n'.join('{}\n{}'.format(id, seq) for id, seq in sorted(geneSEQS.items()))
allSEQS = allSEQS.rstrip()

# Align sequences with MUSCLE
msg('Aligning sequences with MUSCLE ... ')
muscle = subprocess.Popen(['muscle', '-clwstrict'], stdout=subprocess.PIPE, stderr=subprocess.PIPE, stdin=subprocess.PIPE, universal_newlines=True)
stdout, stderr = muscle.communicate(input=allSEQS)

# Print results to stdout or file
if args.out:
	with open(args.out[0], 'wb') as file:
		file.write(stdout)
else:
	print(stdout)
