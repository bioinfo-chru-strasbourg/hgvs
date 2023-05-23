#!/usr/bin/env python
"""
Example usage of HGVS library.

To run the script, first begin by opening a terminal in the root directory
of this software package. Note, the root directory should contain
`setup.py`.

Second, obtain genome sequence in FASTA format, which is required in
example. Genome sequence can be fetched using the following commands:

  cd /tmp
  curl -O http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/chromFa.tar.gz
  tar zxvf chromFa.tar.gz
  cat chr*.fa > hg19.fa
  rm chr*.fa chromFa.tar.gz

This example script can be run using:

  python examples/example1.py

The following output should be displayed:

  chr11 17496508 T C
  NM_000352.3(ABCC8):c.215A>G
  ('NM_000352.3', 'c', '>', CDNACoord(215, -10), CDNACoord(215, -10), 'A', 'G')

"""
from __future__ import print_function

from __future__ import unicode_literals
import pyhgvs as hgvs
import pyhgvs.utils as hgvs_utils
from pyfaidx import Fasta
import duckdb
import pandas as pd


# connexion
conn = duckdb.connect(":memory:", config={"threads":12})

# Databases
genome_file = '/Users/lebechea/HOWARD/databases/genomes/current/hg19.fa'
refgene_file = '/Users/lebechea/HOWARD/databases/refGene/current/refGene.hg19.txt'
refseqlink_file = '/Users/lebechea/HOWARD/databases/refGene/current/ncbiRefSeqLink.hg19.txt'

# Read genome sequence using pyfaidx.
genome = Fasta(genome_file)

# Read RefSeq transcripts into a python dict/model.
with open(refgene_file) as infile:
    transcripts = hgvs_utils.read_transcripts(infile)

# RefGene in database
refgene_structure = {
    "bin": "INTEGER",
    "name": "STRING",
    "chrom": "STRING",
    "strand": "STRING",
    "txStart": "INTEGER",
    "txEnd": "INTEGER",
    "cdsStart": "INTEGER",
    "cdsEnd": "INTEGER",
    "exonCount": "INTEGER",
    "exonStarts": "STRING",
    "exonEnds": "STRING",
    "score": "INTEGER",
    "name2": "STRING",
    "cdsStartStat": "STRING",
    "cdsEndStat": "STRING",
    "exonFrames": "STRING"
}
if refgene_file:
  conn.query(f"CREATE TABLE refgene AS SELECT * FROM read_csv_auto('{refgene_file}',HEADER=False,columns={refgene_structure})")
else:
    sql_structure_list = []
    for col in refgene_structure:
       col_name = col
       col_format = refgene_structure.get(col,"VARCHAR").replace("STRING","VARCHAR")
       sql_structure_list.append(f" {col} {col_format}")
    sql_structure = ",".join(sql_structure_list)
    conn.query(f"CREATE TABLE refgene({sql_structure})")

# RefSeqLink
refseqlink_structure = {
  "id": "STRING",
  "status": "STRING",
  "name": "STRING",
  "product": "STRING",
  "mrnaAcc": "STRING",
  "protAcc": "STRING",
  "locusLinkId": "STRING",
  "omimId": "STRING",
  "hgnc": "STRING",
  "genbank" :"STRING",
  "pseudo": "STRING",
  "gbkey": "STRING",
  "source": "STRING",
  "gene_biotype": "STRING",
  "gene_synonym": "STRING",
  "ncrna_class": "STRING",
  "note": "STRING",
  "description": "STRING",
  "externalId": "STRING",
}
if refseqlink_file:
  conn.query(f"CREATE TABLE refseqlink AS SELECT * FROM read_csv_auto('{refseqlink_file}',HEADER=False,columns={refseqlink_structure})")
else:
    sql_structure_list = []
    for col in refseqlink_structure:
       col_name = col
       col_format = refseqlink_structure.get(col,"VARCHAR").replace("STRING","VARCHAR")
       sql_structure_list.append(f" {col} {col_format}")
    sql_structure = ",".join(sql_structure_list)
    conn.query(f"CREATE TABLE refseqlink({sql_structure})")

# Test databases 
# print(conn.query("SELECT * FROM refgene LIMIT 10"))
# print(conn.query("SELECT * FROM refseqlink LIMIT 10"))


# Provide a callback for fetching a transcript by its name.
def get_transcript(name):
    return transcripts.get(name)

# Get transcript protein by cDNA transcript
def get_transcript_protein(name):
  transcripts_list = list(conn.query(f"SELECT protAcc FROM refseqlink WHERE mrnaAcc='{name}' OR mrnaAcc LIKE '{name}.%'").df()["protAcc"])
  if len(transcripts_list):
     if not pd.isna(transcripts_list[0]):
       return transcripts_list[0]
     else:
        return None
  else:
     return None

def get_transcripts(chr, pos, conn):
  return list(conn.query(f"SELECT name FROM refgene WHERE chrom='{chr}' AND txStart<={pos} AND txEnd>={pos}").df()["name"])


# Parse the HGVS name into genomic coordinates and alleles.
chrom, offset, ref, alt = hgvs.parse_hgvs_name(
    'NM_000352:c.215A>G', genome, get_transcript=get_transcript)
print(chrom, offset, ref, alt)
# Returns variant in VCF style: ('chr11', 17496508, 'T', 'C')
# Notice that since the transcript is on the negative strand, the alleles
# are reverse complemented during conversion.


# Format an HGVS name.
chrom, offset, ref, alt = ('chr11', 17496508, 'T', 'C')
transcripts_list = get_transcripts(chrom, offset, conn)
for transcript_name in transcripts_list:

  transcript = get_transcript(transcript_name)
  transcript_protein = get_transcript_protein(transcript_name)
  exon=transcript.find_exon_number(offset)

  print()
  print(f"Transcript: {transcript_name}")
  print(f"Transcript protein: {transcript_protein}")
  print(f"Exon number: {exon}")

  # cDNA format
  hgvs_name = hgvs.format_hgvs_name(
      chrom, offset, ref, alt, genome=genome, transcript=transcript, transcript_protein=transcript_protein, exon=exon, use_gene=True)
  print(hgvs_name)
  # Protein format
  hgvs_name = hgvs.format_hgvs_name(
      chrom, offset, ref, alt, genome=genome, transcript=transcript, transcript_protein=transcript_protein, exon=exon, use_gene=True, use_protein=True)
  print(hgvs_name)
  # Full format
  hgvs_name = hgvs.format_hgvs_name(
      chrom, offset, ref, alt, genome=genome, transcript=transcript, transcript_protein=transcript_protein, exon=exon, use_gene=True, full_format=True)
  print(hgvs_name)


# hgvs_name = hgvs.HGVSName('NM_000352.3:c.215-10A>G')
# # fields of the HGVS name are available as attributes:
# #
# # hgvs_name.transcript = 'NM_000352.3'
# # hgvs_name.kind = 'c'
# # hgvs_name.mutation_type = '>'
# # hgvs_name.cdna_start = hgvs.CDNACoord(215, -10)
# # hgvs_name.cdna_end = hgvs.CDNACoord(215, -10)
# # hgvs_name.ref_allele = 'A'
# # hgvs_name.alt_allele = 'G'

# print((hgvs_name.transcript,
#        hgvs_name.kind,
#        hgvs_name.mutation_type,
#        hgvs_name.cdna_start,
#        hgvs_name.cdna_end,
#        hgvs_name.ref_allele,
#        hgvs_name.alt_allele))


conn.close()