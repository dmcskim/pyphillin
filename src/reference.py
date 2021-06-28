import pandas as pd
import ftplib
import tarfile

HOSTNAME = 'ftp.ncbi.nlm.nih.gov'
USERNAME = 'anonymous'
PASSWORD = 'dmcskimming@usf.edu'  #make users supply email

def get_assembly_summary():
    #download assembly summary file
    wfile = "genomes/refseq/assembly_summary_refseq.txt"
    #wfile = "refseq/release/RefSeq*catalog.gz"
    ftp = ftplib.FTP(HOSTNAME, USERNAME, PASSWORD)
    ftp.encoding = "utf-8"

    with open(wfile.split('/')[-1], 'wb') as ofile:
        ftp.retrbinary(f"RETR {wfile}", ofile.write)

    ftp.quit()
    return

def get_taxonomy():
    #download taxonomy
    wdir = 'pub/taxonomy/'
    wfile = 'taxdump.tar.gz'

    ftp = ftplib.FTP(HOSTNAME, USERNAME, PASSWORD)
    ftp.encoding = "utf-8"

    with open(wfile, 'wb') as ofile:
        ftp.retrbinary(f"RETR {wdir+wfile}", ofile.write)

    tf = tarfile.open(name=wfile, mode='r:gz')
    tf.extractall('./taxonomy')
    tf.close()
    #extract info
    return

if __name__ == '__main__':
    #scan summary, get reference/representative genomes
    cfile = 'assembly_summary_refseq.txt'
    data = pd.read_csv(cfile, sep='\t', skiprows=1, index_col=0)
