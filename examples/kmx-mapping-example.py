
import pandas as pd
from Bio.Blast.Applications import NcbiblastpCommandline
from Bio.Blast.Applications import NcbimakeblastdbCommandline
from Bio import SeqIO


# download K. marxianus
# https://rest.uniprot.org/uniprotkb/stream?compressed=true&format=fasta&query=accession%3AW0T3G1%20OR%20accession%3AW0T3R8%20OR%20accession%3AW0T4A7%20OR%20accession%3AW0T4C1%20OR%20accession%3AW0T4V8%20OR%20accession%3AW0T5P4%20OR%20accession%3AW0T661%20OR%20accession%3AW0T6B6%20OR%20accession%3AW0T748%20OR%20accession%3AW0T964%20OR%20accession%3AW0T989%20OR%20accession%3AW0T9X4%20OR%20accession%3AW0TA05%20OR%20accession%3AW0TA43%20OR%20accession%3AW0TAA3%20OR%20accession%3AW0TAA4%20OR%20accession%3AW0TAW0%20OR%20accession%3AW0TCD0%20OR%20accession%3AW0TD11%20OR%20accession%3AW0TEF9%20OR%20accession%3AW0TFN1%20OR%20accession%3AW0TGM7%20OR%20accession%3AW0TH64%20OR%20accession%3AW0THY6%20OR%20accession%3AW0TIW1


input_file='./aybrah/yeast_proteome_aybrah.faa'


# create AYbRAH BLASTP DB

cline=NcbimakeblastdbCommandline(cmd='makeblastdb', dbtype='prot', input_file=input_file,out='./aybrah/blastdb/blastdb_prot_aybrah')
cline()

# BLASTP kmx.faa against AYbRAH

cline = NcbiblastpCommandline(query="kmx.faa", db='./aybrah/blastdb/blastdb_prot_aybrah',evalue=1e-5,outfmt=7,out='kmx_1e-25.tsv',num_threads=6)
cline()

# load BLASTP results
blastp=pd.read_csv('kmx_1e-25.tsv',sep='\t',comment='#',header=None)

# df for parsed BLASTP results to map against FOG
kmx_mapped_to_aybrah=pd.DataFrame(columns=['entry_subject','entry_reference','oid_reference','pid','evalue','bitscore','FOG:BLASTP','HOG:BLASTP','FOG:pplacer'])
kmx_no_mappings=[]

records=SeqIO.parse('kmx.faa','fasta')

for record in records:
	print(record.id)
	# find best match
	blastp_matches=blastp[blastp[0]==record.id]
	if len(blastp_matches)==0:
		kmx_no_mappings.append(record.id)
		continue
	row=blastp_matches.iloc[0]
	# protein match
	entry_subject=row[0].split('|')[1]
	entry_reference=row[1].split('|')[1]
	oid_reference=row[1].split('|')[0]
	#
	pid,evalue,bitscore=row[2],row[10],row[11]
	# get AYbRAH details
	try:
		fog_blastp=aybrah.lookup_fog_by_entry.loc[entry_reference].FOG
		hog_blastp=aybrah.get_hog_from_fog(fog_blastp)
	except:
		fog_blastp=''
		hog_blastp=''
	# add to matrix
	kmx_mapped_to_aybrah.loc[len(kmx_mapped_to_aybrah)]=[entry_subject,entry_reference,oid_reference,pid,evalue,bitscore,fog_blastp,hog_blastp,'']


# which proteins are not mapped to K. lactis
kmx_mapped_to_aybrah[kmx_mapped_to_aybrah.oid_reference!='kla']


nonkla=list(filter(None,kmx_mapped_to_aybrah[kmx_mapped_to_aybrah.oid_reference!='kla']['FOG:BLASTP'].tolist()))
fyrment.rxns[fyrment.rxns.FOG.str.contains('|'.join(nonkla))]

# pplacer code, run on silicon rather than locally as it is for Linux


	path_new_sequence='../../pan_proteome/sequences/mafft_add_align/'+str(assembly_uid)+'_new_sequence.fasta'
	SeqIO.write([seq_record], path_new_sequence, "fasta")
	#
	msa_input='~/aybrah_proteomes/aybrah/pplacer/'+hog+'.refpkg/'+hog+'_msa.fasta'
	msa_output='../../pan_proteome/sequences/mafft_add_align/'+seq_record.id.replace('|','_')+'.fasta'
	os.system('mafft --reorder --add '+path_new_sequence+' --auto '+msa_input+' > '+msa_output)
	os.system('pplacer -c ~/aybrah_proteomes/aybrah/pplacer/'+hog+'.refpkg '+msa_output)
	# mafft --add new_sequences --reorder existing_alignment > output












