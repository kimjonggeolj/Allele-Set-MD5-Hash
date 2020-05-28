import os
import pandas as pd
import subprocess

class MD5_plink():
    # === Following are functions for generating a string list that can be hashed ===
    def traw_gen(input_name, output_name = "traw_temp", allele_list = 'forenseq.hg19.set'):
        """generates traw from PLINK using `--export A-transpose`"""
        # traw is truncated to just allele_list
        bashCommand = "plink --bfile " + input_name + " --extract range " + allele_list + " --export A-transpose --out " + output_name
        subprocess.run(bashCommand.split(), stdout=subprocess.PIPE)
        
    def allele_string_traw_gen(traw, sampleID, allele_list = 'forenseq.hg19.set'):
        """this will convert traw from PLINK `--export A-transpose` to a single string of alleles. It must be given a gene set file used to generate the traw."""
        # import from traw
        traw = pd.read_csv(traw, sep = '\t')
        # start with blank string, add as it goes
        string_traw = ""
        # import list of alleles to extract from (allele_list)
        forenseq_list = pd.read_csv(allele_list, sep = '\t', names = ['CHR', 'POS', 'remove1', 'remove2'], usecols = ['CHR', 'POS'])
        #   If there are different number of alleles in traw file (presumably truncated to same allele list) and forenseq
        #   ...prepare the warning message because there are alleles missing!
        sampleColumn = traw.columns.get_loc(sampleID)
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
            return string_traw
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
    # Function that just tapes the two previous functions together
    # presumably what will happen is: plug the output of this function into MD5 function
    def allele_string_gen(input_name, allele_list = 'forenseq.hg19.set'):
        # generate temp traw
        MD5_plink.traw_gen(input_name = input_name)
        # take the traw and generate allele string
        df=pd.DataFrame(columns = ["sampleID", "allele_list"])
        sampleID=pd.read_csv(input_name+".fam", sep="\s", header=None, names=["FID","IID","r1","r2","r3","r4"], usecols=["FID","IID"])
        sampleID=sampleID['FID']+"_"+sampleID['IID']
        # instead of "allele_list" it would be "hash" column instead in the final script
        # NOTE: Maybe make another column with the missing SNPs instead of a printed warning???
        for i in sampleID:
            new_row={'sampleID':i, "allele_list":MD5_plink.allele_string_traw_gen(traw = "traw_temp.traw", sampleID = i)}
            df=df.append(new_row, ignore_index=True)
        # clean up temp
        os.remove("traw_temp.log")
        os.remove("traw_temp.traw")
        return df