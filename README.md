XphasedMultiBridging
====================

This is a tutorial of how to run X-phased Multibridging. 

python batchProcessing.py 

If you issue the command above, you will be guided through a sample run using synthetic genome with synthetic reads. The reads are corrupted by indel noise. The pipeline will start from data generation to finishing the assembly and comparing with the reference.

To cross-check, you may want to compare the UnitTest_motherGen.fasta(the reference) and rec.fasta(the reconstructed genome) .

Remark 1: There are various unit test available in debugging.py that allow you to test each modules of the code[but manual hacking is needed to use the debugging.py]. For the run-down and system architecture of the assembler, one can refer to "Figure 13: Pipeline of the prototype assembler" on P. 26 in the paper of "Near-optimal Assembly for Shotgun Sequencing with Noisy Reads" available at http://arxiv.org/pdf/1402.6971v1.pdf

Remark 2: By uncommentiing appropriate lines towards the end in batchProcessing.py, one can change the synthetic genome with real genome by providing a reference geneome. For this simulation code, it was tested that it will take several hours to a day to run on a genome of size ~1.5M. 
