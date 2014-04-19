XphasedMultiBridging
====================

This is a tutorial of how to run X-phased Multibridging. 

python batchProcessing.py 

If you issue the command above, you will be guided through a sample run using synthetic genome with synthetic reads. The reads are corrupted by indel noise. The pipeline will start from data generation to finishing the assembly and comparing with the reference.

Remark 1: There are various unit test available in debugging.py that allow you to test each modules of the code[but manual hacking is needed to use the debugging.py]. 

Remark 2: By uncommentiing appropriate lines towards the end in batchProcessing.py, one can change the synthetic genome with real genome by providing a reference geneome. For this simulation code, it was tested that it will take several hours to a day to run on a genome of size ~1.5M. 
