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
#conn = duckdb.connect(":memory:", config={"threads":12})
conn = duckdb.connect(":memory:")

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
def get_transcript_protein(conn, name):
  transcripts_list = list(conn.query(f"SELECT protAcc FROM refseqlink WHERE mrnaAcc='{name}' OR mrnaAcc LIKE '{name}.%'").df()["protAcc"])
  if len(transcripts_list):
     if not pd.isna(transcripts_list[0]):
       return transcripts_list[0]
     else:
        return None
  else:
     return None

# Get all transcripts of a genomic position
def get_transcripts(conn, chr, pos):
  return list(conn.query(f"SELECT name FROM refgene WHERE chrom='{chr}' AND txStart<={pos} AND txEnd>={pos}").df()["name"])




variants = [
   # SNV
   {
      "chr": "chr7",
      "pos": 55249061,
      "ref": "C",
      "alt": "G",
      "hgvs": "EGFR:NM_001346898:exon20:c.2359C>G:p.Gln787Glu,EGFR:NM_001346900:exon20:c.2200C>G:p.Gln734Glu,EGFR:NM_001346941:exon14:c.1558C>G:p.Gln520Glu,EGFR-AS1:NR_047551:exon2:n.1203G>C:p.Ala401Pro,EGFR:NM_001346899:exon19:c.2224C>G:p.Gln742Glu,EGFR:NM_001346897:exon19:c.2224C>G:p.Gln742Glu,EGFR:NM_005228:exon20:c.2359C>G:p.Gln787Glu"
   },
   {
      "chr": "chr7",
      "pos": 55249062,
      "ref": "A",
      "alt": "T",
      "hgvs": "EGFR:NM_001346898:exon20:c.2360A>T:p.Gln787Leu,EGFR:NM_001346900:exon20:c.2201A>T:p.Gln734Leu,EGFR:NM_001346941:exon14:c.1559A>T:p.Gln520Leu,EGFR-AS1:NR_047551:exon2:n.1202T>A:p.Leu401Gln,EGFR:NM_001346899:exon19:c.2225A>T:p.Gln742Leu,EGFR:NM_001346897:exon19:c.2225A>T:p.Gln742Leu,EGFR:NM_005228:exon20:c.2360A>T:p.Gln787Leu"
   },
   {
      "chr": "chr7",
      "pos": 55249063,
      "ref": "G",
      "alt": "A",
      "hgvs": "EGFR:NM_001346898:exon20:c.2361G>A:p.Gln787Gln,EGFR:NM_001346900:exon20:c.2202G>A:p.Gln734Gln,EGFR:NM_001346941:exon14:c.1560G>A:p.Gln520Gln,EGFR-AS1:NR_047551:exon2:n.1201C>T:p.Ser401Ser,EGFR:NM_001346899:exon19:c.2226G>A:p.Gln742Gln,EGFR:NM_001346897:exon19:c.2226G>A:p.Gln742Gln,EGFR:NM_005228:exon20:c.2361G>A:p.Gln787Gln"
   },
  # INS
   {
      "chr": "chr7",
      "pos": 55249061,
      "ref": "C",
      "alt": "CAAAAAAA",
      "hgvs": "EGFR:NM_001346898:exon20:c.2360_2361ins7:p.Gln787fs,EGFR:NM_001346900:exon20:c.2201_2202ins7:p.Gln734fs,EGFR:NM_001346941:exon14:c.1559_1560ins7:p.Gln520fs,EGFR-AS1:NR_047551:exon2:n.1202_1203ins7:p.Ser401fs,EGFR:NM_001346899:exon19:c.2225_2226ins7:p.Gln742fs,EGFR:NM_001346897:exon19:c.2225_2226ins7:p.Gln742fs,EGFR:NM_005228:exon20:c.2360_2361ins7:p.Gln787fs"
   },
   {
      "chr": "chr7",
      "pos": 55249062,
      "ref": "A",
      "alt": "AAAAAAAA",
      "hgvs": "EGFR:NM_001346898:exon20:c.2360_2361ins7:p.Gln787fs,EGFR:NM_001346900:exon20:c.2201_2202ins7:p.Gln734fs,EGFR:NM_001346941:exon14:c.1559_1560ins7:p.Gln520fs,EGFR-AS1:NR_047551:exon2:n.1202_1203ins7:p.Ser401fs,EGFR:NM_001346899:exon19:c.2225_2226ins7:p.Gln742fs,EGFR:NM_001346897:exon19:c.2225_2226ins7:p.Gln742fs,EGFR:NM_005228:exon20:c.2360_2361ins7:p.Gln787fs"
   },
   {
      "chr": "chr7",
      "pos": 55249063,
      "ref": "G",
      "alt": "GAAAAAAA",
      "hgvs": "EGFR:NM_001346898:exon20:c.2361_2362ins7:p.Leu788fs,EGFR:NM_001346900:exon20:c.2202_2203ins7:p.Leu735fs,EGFR:NM_001346941:exon14:c.1560_1561ins7:p.Leu521fs,EGFR-AS1:NR_047551:exon2:n.1200_1201ins7:p.Glu401fs,EGFR:NM_001346899:exon19:c.2226_2227ins7:p.Leu743fs,EGFR:NM_001346897:exon19:c.2226_2227ins7:p.Leu743fs,EGFR:NM_005228:exon20:c.2361_2362ins7:p.Leu788fs"
   },
   # DEL
    {
      "chr": "chr7",
      "pos": 55249061,
      "ref": "CA",
      "alt": "C",
      "hgvs": "EGFR:NM_001346898:exon20:c.2360delA:p.Gln787fs,EGFR:NM_001346900:exon20:c.2201delA:p.Gln734fs,EGFR:NM_001346941:exon14:c.1559delA:p.Gln520fs,EGFR-AS1:NR_047551:exon2:n.1202delT:p.Leu401fs,EGFR:NM_001346899:exon19:c.2225delA:p.Gln742fs,EGFR:NM_001346897:exon19:c.2225delA:p.Gln742fs,EGFR:NM_005228:exon20:c.2360delA:p.Gln787fs"
   },
  {
      "chr": "chr7",
      "pos": 55249062,
      "ref": "AG",
      "alt": "A",
      "hgvs": "EGFR:NM_001346898:exon20:c.2361delG:p.Gln787fs,EGFR:NM_001346900:exon20:c.2202delG:p.Gln734fs,EGFR:NM_001346941:exon14:c.1560delG:p.Gln520fs,EGFR-AS1:NR_047551:exon2:n.1201delC:p.Ser401fs,EGFR:NM_001346899:exon19:c.2226delG:p.Gln742fs,EGFR:NM_001346897:exon19:c.2226delG:p.Gln742fs,EGFR:NM_005228:exon20:c.2361delG:p.Gln787fs"
   },
    {
      "chr": "chr7",
      "pos": 55249063,
      "ref": "GC",
      "alt": "G",
      "hgvs": "EGFR:NM_001346898:exon20:c.2362delC:p.Leu788fs,EGFR:NM_001346900:exon20:c.2203delC:p.Leu735fs,EGFR:NM_001346941:exon14:c.1561delC:p.Leu521fs,EGFR-AS1:NR_047551:exon2:n.1200delG:p.Ala400fs,EGFR:NM_001346899:exon19:c.2227delC:p.Leu743fs,EGFR:NM_001346897:exon19:c.2227delC:p.Leu743fs,EGFR:NM_005228:exon20:c.2362delC:p.Leu788fs"
   },
   # MNV - DELINS
        {
      "chr": "chr7",
      "pos": 55249061,
      "ref": "CA",
      "alt": "TT",
      "hgvs": "EGFR:NM_001346898:exon20:c.2359_2360delCAinsTT:p.Gln787Leu,EGFR:NM_001346900:exon20:c.2200_2201delCAinsTT:p.Gln734Leu,EGFR:NM_001346941:exon14:c.1558_1559delCAinsTT:p.Gln520Leu,EGFR-AS1:NR_047551:exon2:n.1202_1203delTGinsAA:p.Cys401Asn,EGFR:NM_001346899:exon19:c.2224_2225delCAinsTT:p.Gln742Leu,EGFR:NM_001346897:exon19:c.2224_2225delCAinsTT:p.Gln742Leu,EGFR:NM_005228:exon20:c.2359_2360delCAinsTT:p.Gln787Leu"
   },
    {
      "chr": "chr7",
      "pos": 55249062,
      "ref": "AG",
      "alt": "TT",
      "hgvs": "EGFR:NM_001346898:exon20:c.2360_2361delAGinsTT:p.Gln787Leu,EGFR:NM_001346900:exon20:c.2201_2202delAGinsTT:p.Gln734Leu,EGFR:NM_001346941:exon14:c.1559_1560delAGinsTT:p.Gln520Leu,EGFR-AS1:NR_047551:exon2:n.1201_1202delCTinsAA:p.Ala401Glu,EGFR:NM_001346899:exon19:c.2225_2226delAGinsTT:p.Gln742Leu,EGFR:NM_001346897:exon19:c.2225_2226delAGinsTT:p.Gln742Leu,EGFR:NM_005228:exon20:c.2360_2361delAGinsTT:p.Gln787Leu"
   },
   {
      "chr": "chr7",
      "pos": 55249063,
      "ref": "GC",
      "alt": "TT",
      "hgvs": "EGFR:NM_001346898:exon20:c.2361_2362delGCinsTT:p.GlnLeu787HisPhe,EGFR:NM_001346900:exon20:c.2202_2203delGCinsTT:p.GlnLeu734HisPhe,EGFR:NM_001346941:exon14:c.1560_1561delGCinsTT:p.GlnLeu520HisPhe,EGFR-AS1:NR_047551:exon2:n.1200_1201delGCinsAA:p.GluLeu400GluMet,EGFR:NM_001346899:exon19:c.2226_2227delGCinsTT:p.GlnLeu742HisPhe,EGFR:NM_001346897:exon19:c.2226_2227delGCinsTT:p.GlnLeu742HisPhe,EGFR:NM_005228:exon20:c.2361_2362delGCinsTT:p.GlnLeu787HisPhe"
   },
    # DELINS
       {
      "chr": "chr7",
      "pos": 55249061,
      "ref": "CA",
      "alt": "T",
      "hgvs": "EGFR:NM_001346898:exon20:c.2359_2360delCAinsT:p.Gln787fs,EGFR:NM_001346900:exon20:c.2200_2201delCAinsT:p.Gln734fs,EGFR:NM_001346941:exon14:c.1558_1559delCAinsT:p.Gln520fs,EGFR-AS1:NR_047551:exon2:n.1202_1203delTGinsA:p.Cys401fs,EGFR:NM_001346899:exon19:c.2224_2225delCAinsT:p.Gln742fs,EGFR:NM_001346897:exon19:c.2224_2225delCAinsT:p.Gln742fs,EGFR:NM_005228:exon20:c.2359_2360delCAinsT:p.Gln787fs"
   },
    {
      "chr": "chr7",
      "pos": 55249062,
      "ref": "AG",
      "alt": "T",
      "hgvs": "EGFR:NM_001346898:exon20:c.2360_2361delAGinsT:p.Gln787fs,EGFR:NM_001346900:exon20:c.2201_2202delAGinsT:p.Gln734fs,EGFR:NM_001346941:exon14:c.1559_1560delAGinsT:p.Gln520fs,EGFR-AS1:NR_047551:exon2:n.1201_1202delCTinsA:p.Ala401fs,EGFR:NM_001346899:exon19:c.2225_2226delAGinsT:p.Gln742fs,EGFR:NM_001346897:exon19:c.2225_2226delAGinsT:p.Gln742fs,EGFR:NM_005228:exon20:c.2360_2361delAGinsT:p.Gln787fs"
   },
   {
      "chr": "chr7",
      "pos": 55249063,
      "ref": "GC",
      "alt": "T",
      "hgvs": "EGFR:NM_001346898:exon20:c.2361_2362delGCinsT:p.Gln787fs,EGFR:NM_001346900:exon20:c.2202_2203delGCinsT:p.Gln734fs,EGFR:NM_001346941:exon14:c.1560_1561delGCinsT:p.Gln520fs,EGFR-AS1:NR_047551:exon2:n.1200_1201delGCinsA:p.Glu400fs,EGFR:NM_001346899:exon19:c.2226_2227delGCinsT:p.Gln742fs,EGFR:NM_001346897:exon19:c.2226_2227delGCinsT:p.Gln742fs,EGFR:NM_005228:exon20:c.2361_2362delGCinsT:p.Gln787fs"
   },
]


i = 0
# Perf loop
while i < 1:

  # Format an HGVS name.
  for variant in variants:
    
    chr = variant.get("chr")
    pos = variant.get("pos")
    ref = variant.get("ref")
    alt = variant.get("alt")

    # hgvs list
    hgvs_list = []
    hgvs_full_list = []

    # Transcripts list
    transcripts_list = get_transcripts(conn, chr, pos)
    #transcripts_list = ["NM_001346898"]

    # For each transcipt
    for transcript_name in transcripts_list:

      # Transcript, protien, exon
      transcript = get_transcript(transcript_name)
      transcript_protein = get_transcript_protein(conn, transcript_name)
      exon=transcript.find_exon_number(pos)
      # transcript_protein = "NP_001333827.1"
      # exon=20
      
      # print()
      # print(f"Transcript: {transcript_name}")
      # print(f"Transcript protein: {transcript_protein}")
      # print(f"Exon number: {exon}")

      # cDNA format
      hgvs_name = hgvs.format_hgvs_name(
          chr, pos, ref, alt, genome=genome, transcript=transcript, transcript_protein=transcript_protein, exon=exon, use_gene=True)
      hgvs_list.append(hgvs_name)
      # Protein format
      hgvs_name = hgvs.format_hgvs_name(
          chr, pos, ref, alt, genome=genome, transcript=transcript, transcript_protein=transcript_protein, exon=exon, use_gene=True, use_protein=True)
      hgvs_list.append(hgvs_name)
      # Full format
      hgvs_name = hgvs.format_hgvs_name(
          chr, pos, ref, alt, genome=genome, transcript=transcript, exon=exon, use_gene=True, full_format=True)
      hgvs_full_list.append(hgvs_name)
      # Full format with NP
      # hgvs_name = hgvs.format_hgvs_name(
      #     chr, pos, ref, alt, genome=genome, transcript=transcript, transcript_protein=transcript_protein, exon=exon, use_gene=True, full_format=True)
      # hgvs_full_list.append(hgvs_name)

    if True:
      print()
      print(f"### {chr}:{pos}-{ref}-{alt}")
      print("# HGVS official: "+",".join(hgvs_list))
      print("# HGVS full    : "+",".join(hgvs_full_list))
      #print(variant.get("hgvs"))
      if not ",".join(hgvs_full_list) == variant.get("hgvs"):
        print("failed!!!")
  
  i += 1

conn.close()