#!/usr/bin/env python
# -*- coding: utf-8 -*-

def parse_args():
    '''Parse arguments given to script'''

    import argparse

    parser = argparse.ArgumentParser(description="Given two genomes, get both \
        common and specific regions to both genomes following a kmer-based algorithm")
    parser.add_argument("--max-seeds", metavar="maximum seeds per kmer", 
                        dest="max_seeds", default=20)
    parser.add_argument("-i, --input", dest="cfile", metavar="cluster.similarities file",
    					required = True)
    parser.add_argument("-o, --output", dest="outfile", metavar="output file", required = True)


    args = parser.parse_args()

    return args

# Import modules to obtain some statistics
from numpy import mean
from numpy import median

class Cluster:
	'''This class holds a cluster of sequences with infor about how many
	sequences form the cluster and how similar each sequence is to the ori sequence
	along with methods to obtain max, min, median, and mean similarity values'''
    
    # Constructor needs a list of tuples containing length and similarity of
    # each aligned sequence
	def __init__(self, name, length, simlist, malign=0):
		self.name = name
		# length of cluster reference sequence
		self.length = length
		self.simlist = simlist
		# Build list with lengths disregarding alignments below cutoff
		self.lenabove = [l[0] for l in simlist if l[0] >= malign]
		# Build list with similarities disregarding alignments below cutoff
		self.simabove = [s[1] for s in simlist if s[0] >= malign]

	def set_cutoff(minalign):
		'''Set a minimum alignment length cutoff'''

		# Build list with lengths disregarding alignments below cutoff
		self.lenabove = [l[0] for l in simlist if l[0] >= malign]
		# Build list with similarities disregarding alignments below cutoff
		self.simabove = [s[1] for s in simlist if s[0] >= malign]


	def nseqs(self):
		'''Return number of seqs that form the cluster, above "minalign" cutoff'''

		return len(simabove)

	# Methods for returning stats about lengths and similarities of each cluster
	# alignments
	def minsim(self):
		return min(self.simabove)

	def maxsim(self):
		return max(self.simabove)

	def meansim(self):
		return mean(self.simabove)

	def mediansim(self):
		return median(self.simabove)

	def minlen(self):
		return min(self.lenabove)

	def maxlen(self):
		return max(self.lenabove)

	def meanlen(self):
		return mean(self.lenabove)

	def medianlen(self):
		return median(self.lenabove)


def cluster_parser(cfile, minalign=0):
	'''This function parses a cluster.similarities text yielding Cluster
	instances'''

	with open(cfile) as handle:
		line = handle.readline()
		while line != "": # Read until EOF
			line = line.rstrip().split(":")
			# store cluster id/name
			cname = line[1]
			# First alignment corresponds with reference sequence aligned with
			# itself. Get reference sequence length info
			line = handle.readline()
			line = line.split("/")[1]
			clength = line.split()[0]
			clength = int(clength)
			if clength >= minalign:
				# If reference sequence is above length cutoff read next line
				line = handle.readline()
			# if next line is another alignment
			if line.startswith("#"):
				# list to store lengths and similarities of each alignment
				simlist = []
				# Read each line that starts with # until new cluster is reached
				while line.startswith("#"):
					# clean each line to get relevant info
					line = line.split("/")[1]
					length, sim = line.split()
					# Get rid of parenthesis and percentage
					sim = sim[1:-2]
					# Convert to int and float and store in simlist
					simtup = (int(length), float(sim))
					simlist.append(simtup)
					# Read next line
					line = handle.readline()
				
				else: # If line does not start with <<#>>
					# Yield new Cluster instance
					yield Cluster(cname, clength, simlist)

def write_stats(cparser, outfile):
	'''Get and write cluster stats'''

	with open(outfile, "w") as handle:
		handle.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format("CLUSTER",
			"Min Sim", "Max Sim", "Median Sim", "Min Len", "Max Len", "Median Len"))

		# Total max similarities of all clusters
		max_sims = []
		# Total min similarities of all clusters
		min_sims = []
		# Total similarities of all clusters
		total_sims = []
		# Total min length of all clusters
        min_lengths = []
        # Total max lengths of all clusters
        max_lengths = []
        # Total lengths of all clusters
        total_lengths = []

        cluster_count = 0
		for cluster in cparser:

			C = cluster.name
			minS = cluster.minsim()
			maxS = cluster.maxsim()
			medianS = cluster.mediansim()
			minL = cluster.minlen()
			maxL = cluster.maxlen()
			medianL = cluster.medianlen()
			cluster_count += 1

			max_sims.append(maxS)
			min_sims.append(minS)
			total_sims.extend(cluster.simabove)

			max_lengths.append(maxL)
			min_lengths.append(minL)
			total_lengths.append(cluster.lenabove())

			handle.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(C,minS,maxS,
			medianS, minL, maxL, medianL))

		# Finally write a summary

		handle.write("#"*80)
		handle.write("\n")
		handle.write("Number of clusters:\t{}".format(cluster_count))
		handle.write("Minimum similarity:\t{}".format(min(min_sims)))
		handle.write("Minimum similarity:\t{}".format(min(min_sims)))

		handle.write("Mean similarity:\t{}".format(mean(total.sims)))