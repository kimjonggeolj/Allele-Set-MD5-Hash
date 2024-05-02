import pandas as pd
import numpy as np
import math
import os
import argparse

import MD5_plink
import utils


def hash_argparse():
    parser = argparse.ArgumentParser(description="Arguments for MD5 Hashing checksums")

    # i/o arguments
    parser.add_argument('--bfile', type=str, nargs='?', default=None, const=None, help='Genotype: String file path to PLINK1.9 files')
    parser.add_argument('--pfile', type=str, nargs='?', default=None, const=None, help='Genotype: String file path to PLINK2.0 files')
    parser.add_argument('--default_snps', type=bool, nargs='?', default=False, const=True, help='Allele list used for GP2 genotypes')
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


def run_hashing(geno_path, out_path, allele_list='GP2_allcall_snps3.set'):
    # takes in geno_path
    # takes in out_path
    # returns tuple:
    # - dictionary of hash:list of samples with same hash (duplicates)
    # - path to file with samples and their hashes

    # should we check if any snps resulted in NaN, and remove all those snps from snpset list?
    # would require to run through the raw_genotype generation twice
    # but necessary to ensure same hash for duplicates with NaNs

    # check num samples
    # if big, split into chunks
    samples = pd.read_csv(f'{geno_path}.psam', sep='\s+', header=None, names=['FID', 'IID', 'PAT', 'MAT', 'SEX', 'PHENO'])
    if samples.shape[0] > 2500:
        # number of chunks
        n = math.ceil(samples.shape[0] / 2500)
        samples_list = np.array_split(samples, n)
        splits = list()
        missings = list()
        for i in range(n):
            # split geno data into chunks
            split = samples_list[i]
            split[['FID', 'IID']].to_csv(f'{geno_path}_split{i+1}.txt', sep='\t', header=False, index=False)
            plink_cmd = f'plink2 --pfile {geno_path} --keep {geno_path}_split{i+1}.txt --make-pgen psam-cols=fid,parents,sex,pheno1,phenos --out {geno_path}_split{i+1}'
            utils.shell_do(plink_cmd)
            # run hash on each chunk
            hasher_split = MD5_plink.MD5_plink(f'{geno_path}_split{i+1}', allele_list=allele_list)
            hashes_split = hasher_split.allele_string_gen()
            # combine all hashes together again
            split = pd.read_csv(f'{geno_path}_split{i+1}_MD5_hash.txt', sep='\s+', header=0, names=['IID', 'HASH'])
            splits.append(split)
            # TODO: handle missing snps files
            # missing = pd.read_csv(f'{geno_path}_split{i+1}_traw_temp_missing_alleles.txt', sep='\s+')
            # missings.append(missing)

            # remove intermediate files (splits)
            for fname in os.listdir(os.path.split(geno_path)[0]):
                if fname.startswith(f'{geno_path}_split'):
                    print(fname)
                    os.remove(f'{geno_path}/{fname}')

        hash_out = pd.concat(splits, axis=0)
        # missing_out = pd.concat(missings, axis=0)

    else:
        hasher = MD5_plink.MD5_plink(f'{geno_path}', allele_list=allele_list)
        hashes = hasher.allele_string_gen()
        hash_out = pd.read_csv(f'{geno_path}_MD5_hash.txt', sep='\s+', header=0, names=['IID', 'HASH'])
        # missing_out = pd.read_csv(f'{geno_path}_traw_temp_missing_alleles.txt', sep='\s+')

    hash_out.to_csv(f'{out_path}_MD5_hashes.txt', sep='\t', header=False, index=False)
    # missing_out.to_csv(f'{out_path}_missing_alleles.txt', sep='\t', header=False, index=False)

    # get duplicates from hashes
    grouped_by_hash = hash_out.groupby('HASH').agg(lambda x:list(x)).reset_index()
    grouped_by_hash['NUM_SAMPS'] = grouped_by_hash['IID'].str.len()
    # counts within hashed duplicates
    hashed_dups = grouped_by_hash[grouped_by_hash['NUM_SAMPS'] > 1]
    hashed_dups_out = pd.Series(hashed_dups.IID.values,index=hashed_dups.HASH).to_dict()

    # TODO: best way to output the duplicates? (currently printing dict to file)
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
    else:
        allele_list = get_perfect_callrate_snps(geno_path, out_path)

    hash_to_ids, hash_file = run_hashing(geno_path, out_path, allele_list)

    with open(f'{out_path}_duplicates.txt', 'w') as f:
        print(hash_to_ids, file=f)
    f.close()

    print(f'duplicates contained in {out_path}_duplicates.txt')
    print(f'hashes are located in {hash_file}')


if __name__ == "__main__":
    handle_main()
