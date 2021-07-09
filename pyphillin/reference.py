import pandas as pd
from ftplib import FTP, error_perm
import tarfile
import zipfile
import subprocess
from os import path
from tqdm import tqdm
from glob import glob

HOSTNAME = 'ftp.ncbi.nlm.nih.gov'
USERNAME = 'anonymous'
PASSWORD = #make users supply email

def get_assembly_summary():
    #download assembly summary file
    wfile = "genomes/refseq/assembly_summary_refseq.txt"
    #wfile = "refseq/release/RefSeq*catalog.gz"
    ftp = FTP(HOSTNAME, USERNAME, PASSWORD)
    ftp.encoding = "utf-8"

    with open(wfile.split('/')[-1], 'wb') as ofile:
        try:
            ftp.retrbinary(f"RETR {wfile}", ofile.write)
        except error_perm as msg:
            print(f"FTP error: {msg}", flush=True)
            pass

    ftp.quit()
    return

def get_taxonomy():
    #download taxonomy
    wdir = 'pub/taxonomy/'
    wfile = 'taxdump.tar.gz'

    ftp = FTP(HOSTNAME, USERNAME, PASSWORD)
    ftp.encoding = "utf-8"

    with open(wfile, 'wb') as ofile:
        try:
            ftp.retrbinary(f"RETR {wdir+wfile}", ofile.write)
        except error_perm as msg:
            print(f"FTP error: {msg}", flush=True)
            pass

    ftp.quit()

    tf = tarfile.open(name=wfile, mode='r:gz')
    tf.extractall('./taxonomy')
    tf.close()
    #extract info
    return

def get_assembly(ftpaddr):
    ttemp = ftpaddr.split('.gov')[1]
    wdir = '/'.join(ttemp.split('/')[:-1])[1:]
    wfile = ttemp.split('/')[-1]
    nfile = '/home/dmcskimming/projects/refseq/'+wfile+'_protein.faa.gz'
    if not path.exists(nfile):

        ftp = FTP(HOSTNAME, USERNAME, PASSWORD)
        ftp.encoding = "utf-8"

        cget = wdir + '/{0}/{0}_protein.faa.gz'.format(wfile)

        with open(nfile, 'wb') as ofile:
            try:
                ftp.retrbinary(f"RETR {cget}", ofile.write)
            except error_perm as msg:
                print(f"FTP error: {msg}", flush=True)
                pass
        ftp.quit()
    return

def download_refseq():
    #scan summary, get reference/representative genomes
    cfile = 'assembly_summary_refseq.txt'
    data = pd.read_csv(cfile, sep='\t', skiprows=1)#, index_col=0)

    keep = []
    for x in tqdm(range(len(data.index))):
        if data.loc[x, 'refseq_category'] != 'na':
            #track wanted
            keep.append(list(data.loc[x,['# assembly_accession',\
                                     'refseq_category',\
                                     'species_taxid',\
                                     'organism_name',\
                                     'ftp_path']]))
            #download data
            get_assembly(data.loc[x, 'ftp_path'])
    return
        
def kofamscan():
    cdir = '/home/dmcskimming/Dropbox/code/pyphillin/tools/kofam_scan/'
    kofam = cdir+'kofam_scan-1.3.0/exec_annotation'

    bdir = '/home/dmcskimming/projects/refseq/'
    ifiles = glob(bdir+'*.faa.gz')

    for tfile in ifiles:
        fname = tfile.split('/')[-1].split('_protein')[0]
        cmd = ['gunzip', tfile]
        res = subprocess.run(cmd)
        if res.returncode != 0:
            print(res)

        cmd = [kofam, '-f', 'mapper', '-o', bdir+fname+'.txt', tfile[:-3]]
        res = subprocess.run(cmd)
        if res.returncode != 0:
            print(res)
    return


if __name__ == '__main__':
    bdir = '/home/dmcskimming/projects/refseq/'
    ifiles = glob(bdir+'*.faa.gz')

    data = pd.read_csv('taxonomy/nodes.dmp', sep="\t|\t", engine='python',\
                      header=None, columns=['tax_id', 'parent_tax_id', 'rank',\
                                           'embl_code', 'division_id', \
                                           'inherited_div_flag', 'gencode_id',\
                                           'inherited_gc_flag', 'MGC_id',\
                                           'inherited_MGC_flag',\
                                            'genbank_hidden_flag',\
                                            'hidden_subtree_root_flag',\
                                            'comments'])


