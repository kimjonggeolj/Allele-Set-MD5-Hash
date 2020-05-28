import os
import pandas as pd
import subprocess
import hashlib

class MD5_plink:
    def __init__(self, geno_path, sampleID, allele_list='forenseq.hg19.set'):
        self.geno_path = geno_path
        self.sampleID = sampleID
        self.allele_list = allele_list
        self.output_name = "traw_temp"
    
    # === Following are functions for generating a string list that can be hashed ===
    def traw_gen(self):
        geno_path = self.geno_path
        output_name = self.output_name
        allele_list = self. allele_list
        
        """generates traw from PLINK using `--export A-transpose`"""
        # traw is truncated to just allele_list
        bashCommand = "plink --bfile " + geno_path + " --extract range " + allele_list + " --export A-transpose --out " + output_name
        subprocess.run(bashCommand.split(), stdout=subprocess.PIPE)
        
    def allele_string_traw_gen(self, traw, sample_for_search):
        allele_list = self.allele_list
        #sampleID = self.sampleID
        """this will convert traw from PLINK `--export A-transpose` to a single string of alleles. It must be given a gene set file used to generate the traw."""
        # import from traw
        traw = pd.read_csv(traw, sep = '\t')
        # start with blank string, add as it goes
        string_traw = ""
        # import list of alleles to extract from (allele_list)
        forenseq_list = pd.read_csv(allele_list, sep = '\t', names = ['CHR', 'POS', 'remove1', 'remove2'], usecols = ['CHR', 'POS'])
        #   If there are different number of alleles in traw file (presumably truncated to same allele list) and forenseq
        #   ...prepare the warning message because there are alleles missing!
        sampleColumn = traw.columns.get_loc(sample_for_search)
        if len(traw) != len(forenseq_list):
            traw_merged = forenseq_list.merge(traw.drop(columns = 'CHR'), how = "outer", on = "POS")
            warning_alleles = ""
            for index, row in traw_merged.fillna(-1).iterrows():
                if row[sampleColumn] == 0:
                    allele = row[5]+row[5]
                elif row[sampleColumn] == 1:
                    allele = row[4]+row[5]
                elif row[sampleColumn] == -1:
                    allele = "NN"
                    warning_alleles = warning_alleles+" "+str(row[0])+":"+str(row[1])
                else:
                    allele = row[4]+row[4]
                string_traw = string_traw+allele
            print("Warning! Following alleles in the allele list file are missing in the given genotype:", warning_alleles)
        #   If the number of alleles are identical, then no warning message
        else:
            for index, row in traw.iterrows():
                if row[sampleColumn] == 0:
                    allele = row[5]+row[5]
                elif row[sampleColumn] == 1:
                    allele = row[4]+row[5]
                else:
                    allele = row[4]+row[4]
                string_traw = string_traw+allele
            print("All alleles present")
        return string_traw
        
        
        
    def md5_gen(self, seq):
        """This takes a sequence of alleles as a string and returns md5 hash digest"""
        md5 = hashlib.md5()
        md5.update(seq.encode('utf-8'))
        dig = md5.hexdigest()
        print(dig)
        return dig
    
    # Function that just tapes the two previous functions together
    # presumably what will happen is: plug the output of this function into MD5 function
    def allele_string_gen(self):
        geno_path = self.geno_path
        sampleID = self.sampleID
        allele_list = self.allele_list
        output_name = self.output_name
        # generate temp traw
        self.traw_gen()
        # take the traw and generate allele string
        # IF THERE ARE MULTIPLE IDs
        if isinstance(sampleID, list):
            allele_hash = []
            for i in sampleID:
                allele_string = self.allele_string_traw_gen(traw="traw_temp.traw", sample_for_search = i)
                # Now create hash
                allele_hash.append(self.md5_gen(allele_string))
        # if there is only one
        else:
            allele_string = self.allele_string_traw_gen(traw="traw_temp.traw", sample_for_search = sampleID)
            # Now create hash
            allele_hash = self.md5_gen(allele_string)
        # clean up temp
        os.remove("traw_temp.log")
        os.remove("traw_temp.traw")
        return allele_hash
        
#test_sample = "001_10_001_10"
with open('test_sample.txt', 'r') as file:
    test_sample = file.read().replace('\n', '')
hasher = MD5_plink('/data/CARD/PD/imputed_data/PLINK/PLINK_HARD.NEUROX_DBGAP', test_sample)
hash_example = hasher.allele_string_gen()
hash_example

# if you want to generate a df of multiple hashes from multiple samples...say the first 20 samples
fam = pd.read_csv("/data/CARD/PD/imputed_data/PLINK/PLINK_HARD.NEUROX_DBGAP.fam", sep="\s", header=None, names=["FID","IID","r1","r2","r3","r4"], usecols=["FID","IID"], nrows=20)
test_sample_list = fam['FID'].values+'_'+fam['IID'].values
hasher_multiple = MD5_plink('/data/CARD/PD/imputed_data/PLINK/PLINK_HARD.NEUROX_DBGAP', test_sample_list.tolist())
df = pd.DataFrame({'sampleID':test_sample_list, 'hash':hasher_multiple.allele_string_gen()})
df.head()