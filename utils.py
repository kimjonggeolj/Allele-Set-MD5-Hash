import pandas as pd
import subprocess
import os
import sys


def shell_do(command, print_cmd=False, log=False, return_log=False, err=False):
    if print_cmd:
        print(f'Executing: {(" ").join(command.split())}', file=sys.stderr)

    res=subprocess.run(command.split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    output = res.stdout.decode('utf-8') + res.stderr.decode('utf-8')
    if log:
        print(output)
    if return_log:
        return output
    if err:
        return res.stderr.decode('utf-8')


def bfiles_to_pfiles(bfile_path=None):
    if not os.path.isfile(f'{bfile_path}.bed'):
        raise FileNotFoundError(f'{bfile_path} does not exist.')

    convert_cmd = f'plink2 --bfile {bfile_path} --make-pgen psam-cols=fid,parents,sex,pheno1,phenos --out {bfile_path}'
    shell_do(convert_cmd)


def create_sampleIDs(geno_path):
    """This generates the sampleIDs from the geno data"""
    if not os.path.isfile(f'{geno_path}.psam'):
        raise FileNotFoundError(f'{geno_path} does not exist.')

    # reformat the .psam file to match the samples in .traw file
    fam_file = pd.read_csv(f'{geno_path}.psam', sep='\s+')
    if '#FID' in fam_file:
        fam_file['sampleID'] = fam_file['#FID'].astype(str) + '_' + fam_file['IID'].astype(str)
    else:
        fam_file['sampleID'] = '0_' + fam_file['#IID'].astype(str)
    return fam_file['sampleID'].tolist()