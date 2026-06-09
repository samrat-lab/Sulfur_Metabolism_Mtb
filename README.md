## Sulfur limitation rewires metabolic state and defines a therapeutic vulnerability in Mycobacterium tuberculosis

# Matlab codes description
**Data_pre_processing.m:** Code for pre-processing the RNA seq data

**Building_metabolic_models_Eflux.m:** Code for building Eflux models for control case (+ sulfur) and test case (- sulfur)

**Mean_match_flux_evaluation.m:** Code for evaluating flux solutions using GP sampler

**Final_flux_values.m:** Filtering the optimal flux profile for both control and test case

**Pathway_enrichment_of_pert_rxns.m:** Finding the subsystems associated with altered reactions, and within each subsystem organizing the number of upregulated and downregulated reactions
