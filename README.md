This program implements a model for determining phase behavior for a 
binary mixture of two charged colloid particles phase separates.
The target application is a mixture of globular RNA with negative 
charge and proteins with positive charge.

*** Compilation:

    g++ -O3 -o phasesep phasesep.c

*** Requirements: 
  
    radial distribution functions 
    provided are RDFs for RNA-RNA, RNA-protein, protein-protein interactions from CG simulations

*** Usage:

    phasesep [options]
      -verbose                            // provide extra output
      -highres                            // output with more digits
      -tag <name>                         // add tag to output
      -temp <value>                       // temperature in K
      -n <rnavalue> <posvalue>            // number of molecules (can be fractional)
      -nrna <value> -npos <value>        
      -c <rnavalue> <posvalue>            // concentrations in mM   
      -crna <value> -cpos <value> 
      -len <value>                        // system size: len*len*len
      -kappa <value>                      // salt screening; typical: 0.5-2.0 

      -fcharge <value>                    // effective charge screening parameter, model 1, typical: 10-30
      -abcharge <value> <value>           // effective charge screening parameter, model 2, typical: 2 0.6 
      -maxvol <value>                     // maximum allowable macromolecular volume fraction in condensate, 
                                          // typical: 0.3-0.5
     
      -q <rnavalue> <posvalue>            // charges
      -qrna <value> -qpos <value>
      -r <rnavalue> <posvalue>            // radii
      -rrna <value> -rpos <value> 
      -epsilon <value>                    // model parameter; default: 4.0
      -a0 <value>                         // model parameter; default: 3.0

      -rnarna <name>                      // name for alternative RNA-RNA RDF profile
      -rnapos <name>                      // name for alternative RNA-POS RDF profile
      -posrna <name>                      // name for alternative POS-RNA RDF profile
      -pospos <name>                      // name for alternative POS-POS RDF profile
      -rdfcut <value>                     // cutoff for RDF beyond which it is set to 1
      -maxrad <value>                     // radial integration limit 
    
      -thresh <value>                     // numerical error tolerance; default: 0.0005
      -vfac <value>                       // multiplier when scanning for volume; default: 0.9995
      -nrscale <value>                    // step size when scanning for concentration; default: 0.02

*** Example:
 
    1. RNA with trypsin at sufficiently high concentration:
       ./phasesep -qrna -46 -rrna 1.47 -qpos 6 -rpos 1.81 -crna 0.4 -cpos 0.3 -kappa 1.3
          phasesep       :[rna_mM]: 0.243405 16.25 :[pos_mM]: 0.059332 24.66 :[clusterrad_nm]: 13.27 :[volfrac]: 0.13 0.37


    2. RNA with trypsin at a concentration that is too low:     
       ./phasesep -qrna -46 -rrna 1.47 -qpos 6 -rpos 1.81 -crna 0.4 -cpos 0.05 -kappa 1.3
          disperse       :[rna_mM]: 0.400000  0.40 :[pos_mM]: 0.050000  0.05 :[clusterrad_nm]: 62.04 :[volfrac]: 0.00 0.00

    3. RNA with trypsin with alternative effective charge model and using a reduced maximum macromolecular volume fraction:
       ./phasesep -qrna -46 -rrna 1.47 -qpos 6 -rpos 1.81 -crna 0.4 -cpos 0.3 -kappa 1.17 -abcharge 2.0 0.6 -maxvol 0.3
          phasesep       :[rna_mM]: 0.306871  8.62 :[pos_mM]: 0.133601 14.98 :[clusterrad_nm]: 13.88 :[volfrac]: 0.07 0.22
       

*** Citation:

    Bercem Dutagaci, Grzegorz Nawrocki, Joyce Goodluck, Ali Akbar Ashkarran, Charles G. Hoogstraten, Lisa Lapidus, Michael Feig: 
    Charge-Driven Condensation of RNA and Proteins Suggests Broad Role of Phase Separation in Cytoplasmic Environments
    bioRxiv(2020)



