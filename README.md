# Psip1 data and scripts

Files in this directory include HMM models for Psip1, for eukaryotic PhoX, eukaryotic PhoAaty, 
Synechococcus and cyanobacterial PhoAaty, as well as the fasta files used to create this models.

Additional files include [PMEs_Psip1_AphoAty_sequences.faa](Psip1Paper_SuppFiles%2FPMEs_Psip1_AphoAty_sequences.faa) which contains all the 
sequences of known phosphatases used to create the dendrogram. [RefSe_plusAndrew_-20.txt](Psip1Paper_SuppFiles%2FRefSe_plusAndrew_-20.txt) and
subsequent derivatives include the sequences used to create the phylogenetic tree in Figure 5 which contains the putative
domain of Psip1.

----------------------------------------------------------------
Scripts folder contains some of the scripts used to analyze and produce the figures in the paper. Raw data is not 
provided, but it can be retrieved if requested.

- [Psip1_activity_PaperfigFig1.py](Psip1Paper_SuppFiles%2FScripts%2FPsip1_activity_Paperfig1.py) - Output the activity of Psip1 in the presence of different metals, pH, substrates and 
concentrations of PNPP represented in Figure 1.
- [BlastP_clustergram.py](Psip1Paper_SuppFiles%2FScripts%2FBlastP_clustergram.py) - Output the dendrogram of all phosphatases and the grouping based on agglomerative hyrarchical 
clustering in Figure 2.
- [MapCruise.py](Psip1Paper_SuppFiles%2FScripts%2FMapCruise.py)- Output the map of the cruises for AMT22 and AMT23 in Figure 4.
- [MapGraph.py](Psip1Paper_SuppFiles%2FScripts%2FMapGraph.py)- Output the map of the stations and the pie charts with the expression of the different phosphatases based 
on relative abundance in Figure 4.
- [Reference_Genome_Indexing.py](Psip1Paper_SuppFiles%2FScripts%2FReference_Genome_Indexing.py) - Creates the reference genome file needed to mapp the reads using BBMAP.

----------------------------------------------------------------
[ManifestFiles](Psip1Paper_SuppFiles%2FManifestFiles) contains the manifest files used to upload the reads to the repository

