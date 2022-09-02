myrepo
======
Suite of tools I developed to facilitate my genomics data analysis work or for practice.

1. mpileupModalDepthCalc.cc
   Calculates frnome-wide modal depth for aligned sequence data in pileup format

2. extract_parents_vcf.cc

   Extracts variant data for specific genome regions in specific individuals from large variant data

3. random_dnaSeq.cc 

   Generates 3 random sequences in fasta format from a single large fasta file input

4. my_Assembler_concatemer1.cc

   Ongoing development of tool which applies deBruijn graph to the de-novo assemembly of short
   read sequence data. Attempts resolution of repeats and assembly into haplotypes. 

5. jsonloader5.cc

   Forms part of a genome variation database software still under development.
   Jsonloader5 implements a bespoke semaphore-controlled de-serialisation of json database file
   into a custom shared memory created and allocated memory space equivalent to json size.

5. menu_check.cc
  
   My very first c++ code. Takes user's menu choice from stdin and checks against available menu.
   Persistent menu request until selected menu is in chef's menu list 
