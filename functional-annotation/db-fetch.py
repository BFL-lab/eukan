import gzip
import requests
import shutil
import subprocess
import sys
from Bio import SeqIO

uniprot_url = "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.xml.gz"
uniprot_url_r = requests.get(uniprot_url, allow_redirects=True)
open("uniprot_sprot.xml.gz", "wb").write(uniprot_url_r.content)
with open("uniprot_sprot.faa", "w") as outfile:
    for rec in SeqIO.parse(gzip.open("uniprot_sprot.xml.gz", "r"), "uniprot-xml"):
        outfile.write(">%s %s\n%s\n" % (rec.id, rec.description, rec.seq))

pfam_url = "http://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam35.0/Pfam-A.hmm.gz"
pfam_url_r = requests.get(pfam_url, allow_redirects=True)
open("Pfam-A.hmm.gz", "wb").write(pfam_url_r.content)
with gzip.open("Pfam-A.hmm.gz", "rb") as f_in:
    with open("Pfam-A.hmm", "wb") as f_out:
        shutil.copyfileobj(f_in, f_out)

subprocess.call(["hmmpress", "Pfam-A.hmm"], shell=True)
