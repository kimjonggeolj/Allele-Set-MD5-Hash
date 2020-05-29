# ForenSeq MD5 Hasher

ForenSeq MD5 Hasher is a tool that generates a MD5 hash of selected genetic variants in a PLINK binary file. It defaults to six ForenSeq SNPs, but a specific variant set can be provided to generate hashes of different variants. It requires participant IDs of interest and the PLINK binary files.

Quick-Start
---
Import the class MD5_plink. Provide the class relevant attributes `geno_path` and `sampleID`, where `geno_path` is the directory of the PLINK binary file and `sampleID` is the list or a string of participant ID(s) of interest. Use MD5_plink.allele_string_gen() to generate the a list of hash from a given list of participant IDs.

```
test_sample = "FID_IID"
hasher = MD5_plink(geno_path='PLINK_geno', sampleID=test_sample)
hash_example = hasher.allele_string_gen()
hash_example
```
