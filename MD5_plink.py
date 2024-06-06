import os
import pandas as pd
import hashlib
import random
import utils

random.seed(42)

class MD5_plink:
    def __init__(self, geno_path, out_path, allele_list='GP2_allcall_snps.set'):
        self.geno_path = geno_path
        self.out_path = out_path
        self.allele_list = allele_list

    # === Following are functions for generating a string list that can be hashed ===
    def traw_gen(self):
        '''
        Generates traw from PLINK using '--export A-transpose'
        Input:
        Output:
        - str: prefix for temporary traw file
        '''

        geno_path = self.geno_path
        allele_list = self.allele_list
        traw_temp = f'{geno_path}_traw_temp'

        # traw is truncated to allele_list
        bashCommand = f'plink2 --pfile {geno_path} --extract bed1 {allele_list} --export A-transpose --out {traw_temp}'
        utils.shell_do(bashCommand)
        return traw_temp

    def allele_string_traw_gen(self, traw_path, samples_to_search):
        '''
        This will convert traw to a single string of alleles. It must be given a gene set file used to generate the traw
        Input:
        - str: prefix to the traw file
        - list of str: list of samples to search for
        Output:
        - list of str: list of alleles_string from traw file
        '''
        allele_list = self.allele_list
        out_path = self.out_path

        # import list of alleles to extract from (allele_list)
        fingerprint_list = pd.read_csv(allele_list, sep='\t', usecols=[0,1], names=['CHR', 'POS'], dtype={'CHR':str})

        # import traw
        traw = pd.read_csv(f'{traw_path}.traw', sep='\t', dtype={'CHR':str})
        # drop any snps that are repetitive
        traw = traw.drop_duplicates(subset=['CHR', 'POS'], keep=False)

        if traw.shape[0] != len(fingerprint_list):
            traw = fingerprint_list.merge(traw.drop(columns='CHR'), how='outer', on='POS').fillna(-1)

        # tranpose df so its long, not wide
        traw = traw.set_index('SNP')
        traw = traw.T

        # prepare list of alleles that are missing
        # will have to remove the warning alleles to ensure identical hashing
        all_missing_alleles = list()
        string_traws_out = list()

        for sample_for_search in samples_to_search:
            # start with blank string, add as it goes
            string_traw = ""
            # warning_alleles = list()
            sampRow = traw.loc[sample_for_search]
            for label, col in traw.items():
                val = sampRow[label]
                if val >= 0 and val < 0.5:
                    allele = col.loc['ALT']+col.loc['ALT']
                elif val >= 0.5 and val < 1.5:
                    allele = col.loc['COUNTED']+col.loc['ALT']
                elif val >= 1.5:
                    allele = col.loc['COUNTED']+col.loc['COUNTED']
                else:
                    allele = 'NN'
                    all_missing_alleles.append((sample_for_search, str(col.loc['CHR'])+':'+str(col.loc['POS'])))
                string_traw = string_traw+allele

            string_traws_out.append(string_traw)

        with open(f'{out_path}_missing_alleles.txt', 'w') as f:
            for snp in all_missing_alleles:
                f.write(f'{snp[0]}\t{snp[1]}\n')
        f.close()

        return string_traws_out


    # Function that generates the hash from input string
    def md5_gen(self, seq):
        """This takes a sequence of alleles as a string and returns md5 hash digest"""
        md5 = hashlib.md5()
        md5.update(seq.encode('utf-8'))
        dig = md5.hexdigest()
        return dig


    # Function that just tapes the two previous functions together
    # presumably what will happen is: plug the output of this function into MD5 function
    def allele_string_gen(self):
        geno_path = self.geno_path
        out_path = self.out_path
        sampleID = utils.create_sampleIDs(geno_path)

        # generate temp traw
        traw_temp = self.traw_gen()

        # take the traw and generate allele string
        allele_hash = []
        allele_strings = self.allele_string_traw_gen(traw_path=traw_temp, samples_to_search=sampleID)
        for i in allele_strings:
            allele_hash.append(self.md5_gen(seq=i))

        # clean up temp
        os.remove(f"{traw_temp}.log")
        os.remove(f"{traw_temp}.traw")

        pd.DataFrame({'sampleID':sampleID, 'hash':allele_hash}, dtype=str).to_csv(f'{out_path}_MD5_hash.txt', sep='\t', header=True, index=False)

        return allele_hash
