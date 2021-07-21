import pandas as pd
from ftplib import FTP, error_perm
import tarfile
import zipfile
import subprocess
from os import path
from tqdm import tqdm
from glob import glob
import gzip
from Bio import SeqIO
from collections import defaultdict
from copy import copy

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
    faa = 'protein.faa.gz'
    rrna = 'rna_from_genomic.fna.gz'
    nfile = '/home/dmcskimming/projects/refseq/'+wfile+'_{1}'.format(faa)
    rfile = '/home/dmcskimming/projects/refseq/'+wfile+'_{1}'.format(rrna)

    if not path.exists(nfile) or path.exists(rfile):
        ftp = FTP(HOSTNAME, USERNAME, PASSWORD)
        ftp.encoding = "utf-8"

        cget = wdir + '/{0}/{0}_{1}'.format(wfile, faa)

        with open(nfile, 'wb') as ofile:
            try:
                ftp.retrbinary(f"RETR {cget}", ofile.write)
            except error_perm as msg:
                print(f"FTP error: {msg}", flush=True)
                pass

        cget = wdir + '/{0}/{0}_{1}'.format(wfile, rrna)

        with open(rfile, 'wb') as ofile:
            try:
                ftp.retrbinary(f"RETR {cget}", ofile.write)
            except error_perm as msg:
                print(f"FTP error: {msg}", flush=True)
                pass

        ftp.quit()
    return

def get_rrna_assembly(ftpaddr):
    ttemp = ftpaddr.split('.gov')[1]
    wdir = '/'.join(ttemp.split('/')[:-1])[1:]
    wfile = ttemp.split('/')[-1]
    rrna = 'rna_from_genomic.fna.gz'
    rfile = '/home/dmcskimming/projects/refseq/rrna/'+wfile+'_{0}'.format(rrna)

    if not path.exists(rfile):
        ftp = FTP(HOSTNAME, USERNAME, PASSWORD)
        ftp.encoding = "utf-8"

        cget = wdir + '/{0}/{0}_{1}'.format(wfile, rrna)

        with open(rfile, 'wb') as ofile:
            try:
                ftp.retrbinary(f"RETR {cget}", ofile.write)
            except error_perm as msg:
                print(f"FTP error: {msg}", flush=True)
                pass
        ftp.quit()
    return

def download_refseq(both=True):
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
            if both:
                get_assembly(data.loc[x, 'ftp_path'])
            else:
                #only get rRNA data
                get_rrna_assembly(data.loc[x, 'ftp_path'])
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


def get_taxonomy_nodes():
    bdir = '/home/dmcskimming/projects/refseq/'
    ifiles = glob(bdir+'*.faa.gz')

    data = pd.read_table('taxonomy/nodes.dmp', sep="\t\|\t", engine='python',\
                      header=None, index_col=0,\
                       names=['tax_id', 'parent_tax_id', 'rank',\
                                'embl_code', 'division_id', \
                                'inherited_div_flag', 'gencode_id',\
                                'inherited_gc_flag', 'MGC_id',\
                                'inherited_MGC_flag',\
                                'genbank_hidden_flag',\
                                'hidden_subtree_root_flag',\
                                'comments'])
    return data


if __name__ == '__main__':
    tax_data = get_taxonomy_nodes()

    cfile = 'assembly_summary_refseq.txt'
    data = pd.read_csv(cfile, sep='\t', skiprows=1, index_col=0)

    keep = []
    for x in tqdm(data.index):
        if data.loc[x, 'refseq_category'] != 'na':
            ctax = data.loc[x, 'taxid']
            if ctax in tax_data.index:
                if tax_data.loc[ctax, 'division_id'] == 0:
                    keep.append(x)

    bdir = '/home/dmcskimming/projects/refseq/'
    cdir = bdir+'rrna/'
    ddir = bdir+'profiles/'
    rrna_seqs = []
    rrna_cnts = defaultdict(int)
    kegg_ids = {}
    for x in tqdm(keep):
        #extract 16S sequences
        tfile = glob(cdir+'{0}*'.format(x))
        with gzip.open(tfile[0], 'rt') as handle:
            for record in SeqIO.parse(handle, "fasta"):
                if 'product=16S' in record.description:
                    record.id = x+'|'+record.id
                    rrna_seqs.append(record)
                    rrna_cnts[x] += 1
        #extract K values
        kfile = glob(ddir+'{0}*'.format(x))
        kdata = pd.read_csv(kfile[0], sep='\t', header=None,\
                            names=['native','KEGG'])
        ckegg = defaultdict(int)
        for kx in kdata['KEGG'].dropna().values:
            ckegg[kx] += 1
        kegg_ids[x] = copy(ckegg)
    
    #build tables
    SeqIO.write(rrna_seqs, 'pyphillin_16S_sequences.fasta', 'fasta')
    rrna_res = pd.DataFrame.from_dict(rrna_cnts, orient='index',\
                                      columns=['16S_CN'])
    rrna_res.to_csv('pyphillin_16S_copy_number.tsv', sep='\t', index=True)
    profiles = pd.DataFrame.from_dict(kegg_ids, orient='index')
    profiles.to_csv('pyphillin_KEGG_profiles.tsv', sep='\t', index=True)
