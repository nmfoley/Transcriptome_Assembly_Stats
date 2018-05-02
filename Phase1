## Python Homework #3 20180420

###################################################################################
# Some statistics use only the transcriptome itself:

# Total number of basepairs
# Number of transcripts
# Number of "unigenes"
# Mean and median transcript length
# N50 explanation here
###################################################################################

from Bio import SeqIO
import statistics
import math

# Total number of basepairs	
contigsLength = []
total = 0

for seq_record in SeqIO.parse(open("NW-1.Trinity.fasta"), 'fasta'):
	contigsLength.append(len(seq_record.seq))
	total += len(seq_record.seq)	

# Mean, median and mode transcript length	
mean = statistics.mean(contigsLength)
median = statistics.median(contigsLength)

# Number of transcripts
record_dict = SeqIO.to_dict(SeqIO.parse("NW-1.Trinity.fasta", "fasta"))
totalTxn = (len(record_dict))

# Number of "unigenes"
# First split the id into unigene and isoforms, only keep the unigene part
unigenes = []
for key in record_dict:
	unigenes.append(key.split('_i')[0])
# get just the unique unigenes and count them
uni = set(unigenes)
unigenes = len(uni)

# Get N50
all_len = sorted(contigsLength, reverse=True)
index = (len(all_len))/2
ind = math.ceil(index)
n50 = all_len[ind]

# print out all the info
print("The total number of bps: %s" % (total))
print("Mean transcript length: %s" % (mean))
print("Median transcript length: %s" % (median))
print("The total number of transcripts: %s" % (totalTxn))
print("The total number of unigenes: %s" % (unigenes))
print("The N50 value is: %s" % (n50))
