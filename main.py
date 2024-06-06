import pandas as pd
import os
import argparse
import json

import MD5_plink
import utils


def hash_argparse():
    parser = argparse.ArgumentParser(description="Arguments for MD5 Hashing checksums")

    # i/o arguments
    parser.add_argument('--bfile', type=str, nargs='?', default=None, const=None, help='Genotype: String file path to PLINK1.9 files')
    parser.add_argument('--pfile', type=str, nargs='?', default=None, const=None, help='Genotype: String file path to PLINK2.0 files')
    parser.add_argument('--default_snps', type=bool, nargs='?', default=False, const=True, help='Allele list used for GP2 genotypes')
    parser.add_argument('--snp_list', type=str, nargs='?', default=None, const=None, help='Allele list: String file path')
    parser.add_argument('--out', type=str, nargs='?', default=None, const=None, help='Path prefix for output')

    args = parser.parse_args()
    return args


def get_perfect_callrate_snps(geno_path, out_path):
    # takes in geno_path
    # takes in out_path for outputted snpset
    # returns path to snpset file containing 100% callrate snps by chr:pos

    # extract variants with 100% callrate
    plink_cmd = f'plink2 --pfile {geno_path} --geno 0.0 --make-pgen psam-cols=fid,parents,sex,pheno1,phenos --out {geno_path}_allcall'
    utils.shell_do(plink_cmd)
    # get these 100% callrate snps and compile
    allcalls = pd.read_csv(f'{geno_path}_allcall.pvar', sep='\s+', dtype={'#CHROM':str})
    allcalls['CHR:POS'] = allcalls['#CHROM'].apply(str) + ':' + allcalls['POS'].apply(str)
    allcalls.drop_duplicates(subset='CHR:POS', keep=False, inplace=True)
    allcalls = allcalls['CHR:POS'].str.split(':', expand=True)
    allcalls['POS'] = allcalls.iloc[:,1]
    allcalls['SET'] = 'ALLCALL'

    allcalls.to_csv(f'{out_path}/allcall_snps.set', sep='\t', header=False, index=False)
    return f'{out_path}/allcall_snps.set'


def run_hashing(geno_path, out_path, allele_list='GP2_allcall_snps.set'):
    # takes in geno_path
    # takes in out_path
    # returns tuple:
    # - dictionary of hash:list of samples with same hash (duplicates)
    # - path to file with samples and their hashes

    # should we check if any snps resulted in NaN, and remove all those snps from snpset list?
    # would require to run through the raw_genotype generation twice
    # but necessary to ensure same hash for duplicates with NaNs

    # check validity of allele_list path
    if not os.path.exists(allele_list):
        raise FileNotFoundError(f'{allele_list} does not exist.')

    hasher = MD5_plink.MD5_plink(f'{geno_path}', f'{out_path}', allele_list=allele_list)
    hasher.allele_string_gen()
    hash_out = pd.read_csv(f'{out_path}_MD5_hash.txt', sep='\s+', header=None, names=['IID', 'HASH'])

    # get duplicates from hashes
    grouped_by_hash = hash_out.groupby('HASH').agg(lambda x:list(x)).reset_index()
    grouped_by_hash['NUM_SAMPS'] = grouped_by_hash['IID'].str.len()
    # counts within hashed duplicates
    hashed_dups = grouped_by_hash[grouped_by_hash['NUM_SAMPS'] > 1]
    hashed_dups_out = pd.Series(hashed_dups.IID.values,index=hashed_dups.HASH).to_dict()

    return (hashed_dups_out, f'{out_path}_MD5_hashes.txt')


def handle_main():
    args = hash_argparse()
    args_dict = vars(args)

    # check format of input files
    if (args_dict['bfile'] is None) and (args_dict['pfile'] is None):
        raise KeyError('No bfile of pfile genotypes were provided')
    elif args_dict['bfile'] and (args_dict['pfile'] is None):
        utils.bfiles_to_pfiles(bfile_path=args_dict['bfile'])
        args_dict['geno_path'] = args_dict['bfile']
    else:
        args_dict['geno_path'] = args_dict['pfile']

    geno_path = args_dict['geno_path']
    out_path = args_dict['out']

    # check if allele list should be created
    if args_dict['default_snps']:
        allele_list = 'GP2_allcall_snps.set'
    elif args_dict['snp_list']:
        allele_list = args_dict['snp_list']
    else:
        allele_list = get_perfect_callrate_snps(geno_path, out_path)

    hash_to_ids, hash_file = run_hashing(geno_path, out_path, allele_list)

    json.dump(hash_to_ids, open(f'{out_path}_duplicates.json', 'w'))

    print(f'duplicates contained in {out_path}_duplicates.json')
    print(f'hashes are located in {hash_file}')

    if os.path.getsize(f'{out_path}_missing_alleles.txt') > 0:
        num_samps_with_missing = len(open(f'{out_path}_missing_alleles.txt', 'r').readlines())
        print(f'WARNING: {num_samps_with_missing} sample(s) are missing calls from the allele list. this may result in mismatching hashes for duplicate samples.')
        print(f'missing snps are located in {out_path}_missing_alleles.txt')
    else:
        print('no missing snp calls!')

if __name__ == "__main__":
    handle_main()
