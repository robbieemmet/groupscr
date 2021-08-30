# groupscr
Supporting Information Code for Emmet et al. 2021, "A spatial capture-recapture model for group-living species", Ecology.

This is code for running the simulation portion of Robert L. Emmet, Ben C. Augustine, Briana Abrahms, Lindsey Rich, and Beth Gardner. 2021. "A spatial capture-recapture model for group-living species." Ecology.

Any scripts listed in MetadataS1.docx is also included in the online Supplementary Information for the paper.

Other scripts contain code we used to conduct and analyze simulations.

"runallscens_knownid.R" simulates data and fits models for the cluster SCR model fit to data generated from the same model.

"runallscens_groupscr.R" simulates data from the cluster CSR model and fits a standard group-level SCR model (found in "SCR0_group.txt").

"runallscens_bcrw_cluster.R" simulates data assuming hidden Markov model group dynamic movement and fits the cluster SCR model.

The files containing "Process Sims" in the name can be used to process the simulation results from the above scripts.

Questions about the code can be directed to Robert L. Emmet at robert.l.emmet@gmail.com.
