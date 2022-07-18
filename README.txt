==================================================================================

  BARONESS: BAsis foR the OrthogonaNal complEment of the Stoichiometric Subspace

==================================================================================

Matlab was used to develop the function used here.The function basisOrthComp returns details about the network and its stoichiometric subspace, its orthogonal complement of the stoichiometric subspace, and the basis for the orthogonal complement of the stoichiometric subspace with labels for the nonpivot species. Furthermore, the output variables 'model', 'basis', and 'nonpivot' allow the user to view the following, respectively:

   - Complete network with all the species listed in the 'species' field of the structure 'model'
   - Table showing the basis for the orthogonal complement of the stoichiometric subspace
   - Species representing nonpivot columns



====
Note
====

The idea for the first part of the algorithm comes from [2].



=================================
How to fill out 'model' structure
=================================

'model' is the input for the function basisOrthComp. It is a structure, representing the CRN, with the following fields:

   - id: name of the model
   - species: a list of all species in the network; this is left blank since incorporated into the function is a step which compiles all species used in the model
   - reaction: a list of all reactions in the network, each with the following subfields:
        - id: a string representing the reaction
        - reactant: has the following further subfields:
             - species: a list of strings representing the species in the reactant complex
             - stoichiometry: a list of numbers representing the stoichiometric coefficient of each species in the reactant complex (listed in the same order of the species)
        - product: has the following further subfields:
             - species: a list of strings representing the species in the product complex
             - stoichiometry: a list of numbers representing the stoichiometric coefficient of each species in the product complex (listed in the same order of the species)
        - reversible: has the value true or false indicating if the reaction is reversible or not, respectively
        - kinetic: has the following further subfields:
             - reactant1: a list of numbers representing the kinetic order of each species in the reactant complex in the left to right direction (listed in the same order of the species)
             - reactant2: a list of numbers representing the kinetic order of each species in the reactant complex in the right to left direction (listed in the same order of the species) (empty if the reaction is not reversible)

To fill out the 'model' structure, write a string for 'model.id': this is just to put a name to the network. To add the reactions to the network, use the function addReaction where the output is 'model'. addReaction is developed to make the input of reactions of the CRN easier than the input in [7]:

   addReaction
      - OUTPUT: Returns a structure called 'model' with added field 'reaction' with subfields 'id', 'reactant', 'product', 'reversible', and 'kinetic'. The output variable 'model' allows the user to view the network with the added reaction.
      - INPUTS
           - model: a structure, representing the CRN
           - id: visual representation of the reaction, e.g., reactant -> product (string)
           - reactant_species: species of the reactant complex (cell)
           - reactant_stoichiometry: stoichiometry of the species of the reactant complex (cell)
           - reactant_kinetic: kinetic orders of the species of the reactant complex (array)
           - product_species: species of the product complex (cell)
           - product_stoichiometry: stoichiometry of the species of the product complex (cell)
           - product_kinetic: "kinetic orders" of the species of the product complex, if the reaction is reversible (array); if the reaction in NOT reversible, leave blank
           - reversible: logical; whether the reaction is reversible or not (true or false)
      * Make sure the function addReaction is in the same folder/path being used as the current working directory.



========
Examples
========

2 examples are included in this folder:

   - Example 1: Deficiency 0 subnetwork of deficiency-oriented coarsening of the insulin signaling system [1]

   - Example 2: Deficiency 7 subnetwork of deficiency-oriented coarsening of the insulin signaling system [1]



===================
Contact Information
===================

For questions, comments, and suggestions, feel free to contact me at pvnlubenia@yahoo.co.uk.


- Patrick Lubenia (18 July 2022)



==========
References
==========

   [1] Lubenia P, Mendoza E, Lao A (2022) Reaction network analysis of metabolic insulin signaling (submitted)

   [2] Soranzo N, Altafini C (2009) ERNEST: a toolbox for chemical reaction network theory. Bioinform 25(21):2853Ð2854. https://doi.org/10.1093/bioinformatics/btp513