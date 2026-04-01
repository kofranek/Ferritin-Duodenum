within ;
package FeMetabolism
   //#########################################################################
  // based on Ensculescu model modified by Pugachov FeMetabolism, version 1.63
  //##########################################################################

  model FeMetabolismModel
    "Implementation of Enculescu et al. PLOS Comput. Biol. 2017 Fe metabolism model: final version"

    //+Constants

    //################################################
    // Main constants in ODES
    //################################################

    // Hepcidin expression
    //################################################
    constant Real hep_norm(unit="ug") = 0.664938 "Hepcidin amount (s16 ) - NORM";
    constant Real Bmp6_norm(unit="ug") = 15.8291 "Bmp6 amount (s18) - NORM";
    constant Real LPS_norm(fixed=true,unit="ug") = 0.000 "LPS (Lipopolysaccharide) amount (s23) - NORM";
    constant Real Il6mRNA_norm(unit="ug") = 0.000 "IL-6 mRNA amount (s21) - NORM";
    constant Real Il6_norm(unit="ug") = 0.000 "IL-L amount (s22) - NORM";

    // Ferroportin regulation
    //################################################
    constant Real Fpn_liv_mRNA_norm(unit="ug") = 0.922424 "Fpn mRNA amount in liver (s7) - NORM";
    constant Real Fpn_spl_mRNA_norm(unit="ug") = 0.922424 "Fpn mRNA amount in spleen (s9) - NORM";
    constant Real Fpn_duo_mRNA_norm(unit="ug") = 0.922424 "Fpn mRNA amount in duodenum (s8) - NORM";
    constant Real Fpn_res_mRNA_norm(unit="ug") = 0.922424 "Fpn mRNA amount in other organs (s11) - NORM";

    constant Real Fpn_liv_norm(unit="ug") = 1.00036 "Fpn amount in liver (s12) - NORM";
    constant Real Fpn_spl_norm(unit="ug") = 1.00049 "Fpn amount in spleen (s14) - NORM";
    constant Real Fpn_duo_norm(unit="ug") = 1.00017 "Fpn amount in duodenum (s13) - NORM";
    constant Real Fpn_res_norm(unit="ug") = 1.00020 "Fpn amount in other organs (s15) - NORM";

    // Dynamics of the iron pools
    //################################################
    constant Real Fe_liv_norm(unit="ug") = 76.9562 "Fe amount in liver (s2) - NORM";
    constant Real Fe_spl_norm(unit="ug") = 17.7295 "Fe amount in spleen (s5) - NORM";
    constant Real Fe_bm_norm(unit="ug") = 63.1596 "Fe amount in bones (s3) - NORM";
    constant Real Fe_RBC_norm(unit="ug") = 1016.720 "Fe amount in red blood cells (s6) - NORM";
    constant Real Fe_duo_norm(unit="ug") = 2.97118 "Fe amount in duodenum (s4) - NORM";
    constant Real Fe_res_norm(unit="ug") = 466.811 "Fe amount in other organs (s10) - NORM";
    constant Real Fe_ser_norm(unit="ug") = 1.51304 "Fe amount in serum (s1) - NORM";

    constant Real Fe_total_norm(unit="ug") = Fe_liv_norm + Fe_spl_norm + Fe_bm_norm + Fe_RBC_norm + Fe_duo_norm + Fe_res_norm + Fe_ser_norm "Total Fe amount - NORM";

    constant Real Fe_ser_input_norm(unit="ug.h-1") = 20.897783 "Fe amount in serum: input - NORM";
    constant Real Fe_ser_output_norm(unit="ug.h-1") = 20.897783 "Fe amount in serum: output - NORM";

    // ###############################################
    // Auxiliary constants
    // ###############################################

    // Duodenum

    constant Real Fe_duo_2_norm(unit="ug") = 0.7427 "Fe++ in duodenum - NORM";
    constant Real Fe_duo_3_norm(unit="ug") = 2.2283 "Fe+++ in duodenum, basically in ferritin - NORM";

    constant Real Fe_duo_to_ferritin_norm = 2.2283;
    constant Real Fe_duo_from_ferritin_norm = 2.2283;

    constant Real Fe_duo_intake_norm = 5.4441032 "Fe intake from food - NORM";

    constant Real Fe_duo_in_food_norm = 1.8323029 "Fe input to duodenum from intestines, prev:in_1 - NORM";
    constant Real Fe_duo_in_ser_norm = 1.066469 "Fe input to duodenum from serum, prev:in_2 - NORM";
    constant Real Fe_duo_out_ser_norm = 2.625664 "Fe output from duodenum to serum, prev:out_1 - NORM";
    constant Real Fe_duo_out_loss_norm = 0.27310795 "Fe output from duodenum to intestines (Fe loss), prev:out_2 - NORM";
    constant Real Fe_duo_unused_norm = Fe_duo_intake_norm * (1 - min(1 / Fe_duo_norm, 1)) "Unused Fe going through intestines - NORM";

    constant Real Fpn_duo_in_1_norm = 0.027948529 "Fpn in duodenum creation: standard from mRNA - NORM";
    constant Real Fpn_duo_in_2_norm = 0.013292233 "Fpn in duodenum creation by activation of Fe_duo - NORM";
    constant Real Fpn_duo_out_1_norm = 0.021124585 "Degradation (internalization) of Fpn in duodenum: standard degradation rate - NORM";
    constant Real Fpn_duo_out_2_norm = 0.020116176 "Degradation (internalization) of Fpn in duodenum by hepcidine - NORM";

    constant Real Fpn_duo_in_norm = Fpn_duo_in_1_norm + Fpn_duo_in_2_norm "Fpn in duodenum synthesis - NORM";
    constant Real Fpn_duo_out_norm = Fpn_duo_out_1_norm + Fpn_duo_out_2_norm "Fpn in duodenum degradation - NORM";

    constant Real Fpn_duo_mRNA_in_norm = 1 "Fpn mRNA in duodenum synthesis - NORM";
    constant Real Fpn_duo_mRNA_out_norm = 1 "Fpn mRNA in duodenum degradation - NORM";
    constant Real Fpn_duo_mRNA_inhib_norm = 0 "Fpn mRNA in duodenum production inhibition by IL-6 (or indirectly by LPS) - NORM";

    // Liver

    constant Real Fe_liv_2_norm(unit="ug") = 19.2487 "Fe++ in liver - NORM";
    constant Real Fe_liv_3_norm(unit="ug") = 57.7461 "Fe+++ in liver, basically in ferritin - NORM";

    constant Real Fe_liv_to_ferritin_norm = 57.7461;
    constant Real Fe_liv_from_ferritin_norm = 57.7461;

    constant Real Fe_liv_in_ser_norm = 5.992713 "Fe input to liver from serum - NORM";
    constant Real Fe_liv_out_ser_norm = 5.992713 "Fe output from liver to serum - NORM";

    constant Real Fpn_liv_in_1_norm = 0.11963841 "Fpn in liver creation: standard from mRNA - NORM";
    constant Real Fpn_liv_in_2_norm = 0.030545777 "Fpn in liver creation by activation of Fe_liv - NORM";
    constant Real Fpn_liv_out_1_norm = 0.05538275 "Degradation (internalization) of Fpn in liver: standard degradation rate - NORM";
    constant Real Fpn_liv_out_2_norm = 0.09480143 "Degradation (internalization) of Fpn in liver by hepcidine - NORM";

    constant Real Fpn_liv_in_norm = Fpn_liv_in_1_norm + Fpn_liv_in_2_norm "Fpn in liver synthesis - NORM";
    constant Real Fpn_liv_out_norm = Fpn_liv_out_1_norm + Fpn_liv_out_2_norm "Fpn in liver degradation - NORM";

    constant Real Fpn_liv_mRNA_in_norm = 1 "Fpn mRNA in liver synthesis - NORM";
    constant Real Fpn_liv_mRNA_out_norm = 1 "Fpn mRNA in liver degradation - NORM";
    constant Real Fpn_liv_mRNA_inhib_norm = 0 "Fpn mRNA in liver production inhibition by IL-6 (or indirectly by LPS) - NORM";

    // Serum

    constant Real Fe_ser_in_liv_norm = 5.992713 "Fe flow from liver to serum - NORM";
    constant Real Fe_ser_in_spl_norm = 4.2752566 "Fe flow from spleen to serum - NORM";
    constant Real Fe_ser_in_duo_norm = 2.625664 "Fe flow from duodenum to serum - NORM";
    constant Real Fe_ser_in_res_norm = 8.00415 "Fe flow from other organs to serum - NORM";

    constant Real Fe_ser_out_liv_norm = 5.992713 "Fe flow from serum to liver - NORM";
    constant Real Fe_ser_out_bm_norm = 4.2752566 "Fe flow from serum to bone marrow - NORM";
    constant Real Fe_ser_out_duo_norm = 1.066469 "Fe flow from serum to duodenum - NORM";
    constant Real Fe_ser_out_res_norm = 9.563345 "Fe flow from serum to other organs - NORM";

    // Hepcidine related stuff

    constant Real hep_in_norm = 0.04654567 "Production of hepcidin - NORM";
    constant Real hep_out_norm = 0.04654567 "Degradation of hepcidin - NORM";

    constant Real Il6mRNA_in_norm = 0 "Production of IL-6 mRNA - NORM";
    constant Real Il6mRNA_out_norm = 0 "Degradation of IL-6 mRNA - NORM";

    constant Real Il6_in_norm = 0 "Production of IL-6 - NORM";
    constant Real Il6_out_norm = 0 "Degradation of IL-6 - NORM";

    constant Real Bmp6_in_norm = 37.929714 "Bmp6 production rate - NORM";
    constant Real Bmp6_out_norm = 37.929714 "Bmp6 degradation rate - NORM";

    // Spleen

    constant Real Fe_spl_2_norm(unit="ug") = 4.43525 "Fe++ in spleen - NORM";
    constant Real Fe_spl_3_norm(unit="ug") = 13.30575 "Fe++ in spleen - NORM";

    constant Real Fe_spl_to_ferritin_norm = 13.30575;
    constant Real Fe_spl_from_ferritin_norm = 13.30575;

    constant Real Fe_spl_in_RBC_norm = 3.663764 "Fe flow to spleen from RBC - NORM";
    constant Real Fe_spl_in_bm_norm = 0.61149263 "Fe flow to spleen from bone marrow - NORM";
    constant Real Fe_spl_out_ser_norm = 4.2752566 "Fe flow from spleen to serum - NORM";

    constant Real Fpn_spl_in_1_norm = 0.020959321 "Fpn in spleen creation: standard from mRNA - NORM";
    constant Real Fpn_spl_in_2_norm = 0.005212414 "Fpn in spleen creation by activation of Fe_spl - NORM";
    constant Real Fpn_spl_out_1_norm = 0.0030254605 "Degradation (internalization) of Fpn in spleen: standard degradation rate - NORM";
    constant Real Fpn_spl_out_2_norm = 0.023146275 "Degradation (internalization) of Fpn in spleen by hepcidine - NORM";

    constant Real Fpn_spl_in_norm = Fpn_spl_in_1_norm + Fpn_spl_in_2_norm "Fpn in spleen synthesis - NORM";
    constant Real Fpn_spl_out_norm = Fpn_spl_out_1_norm + Fpn_spl_out_2_norm "Fpn in spleen degradation - NORM";

    constant Real Fpn_spl_mRNA_in_norm = 1 "Fpn mRNA in spleen synthesis - NORM";
    constant Real Fpn_spl_mRNA_out_norm = 1 "Fpn mRNA in spleen degradation - NORM";
    constant Real Fpn_spl_mRNA_inhib_norm = 0 "Fpn mRNA in spleen production inhibition by IL-6 (or indirectly by LPS) - NORM";

    // Other organs

    constant Real Fe_res_2_norm(unit="ug") = 116.74625 "Fe++ in other organs - NORM";
    constant Real Fe_res_3_norm(unit="ug") = 350.23875 "Fe++ in other organs - NORM";

    constant Real Fe_res_to_ferritin_norm = 350.23875;
    constant Real Fe_res_from_ferritin_norm = 350.23875;

    constant Real Fe_res_in_ser_norm = 9.563345 "Fe flow to other organs from serum - NORM";
    constant Real Fe_res_out_ser_norm = 8.00415 "Fe flow from other organs to serum - NORM";
    constant Real Fe_res_out_loss_norm = 1.559195 "Fe loss from other organs (skin peeling, etc.) - NORM";

    constant Real Fpn_res_in_1_norm = 0.0046672816 "Fpn in other organs creation: standard from mRNA - NORM";
    constant Real Fpn_res_in_2_norm = 0.2478532 "Fpn in other organs creation by activation of Fe_res - NORM";
    constant Real Fpn_res_out_1_norm = 0.028923023 "Degradation (internalization) of Fpn in other organs: standard degradation rate - NORM";
    constant Real Fpn_res_out_2_norm = 0.22359747 "Degradation (internalization) of Fpn in other organs by hepcidine - NORM";

    constant Real Fpn_res_in_norm = Fpn_res_in_1_norm + Fpn_res_in_2_norm "Fpn in spleen synthesis - NORM";
    constant Real Fpn_res_out_norm = Fpn_res_out_1_norm + Fpn_res_out_2_norm "Fpn in spleen degradation - NORM";

    constant Real Fpn_res_mRNA_in_norm = 1 "Fpn mRNA in other organs synthesis - NORM";
    constant Real Fpn_res_mRNA_out_norm = 1 "Fpn mRNA in other organs degradation - NORM";
    constant Real Fpn_res_mRNA_inhib_norm = 0 "Fpn mRNA in other organs production inhibition by IL-6 (or indirectly by LPS) - NORM";

    // Bone marrow

    constant Real Fe_bm_in_ser_norm = 4.2752566 "Fe flow to bone marrow from serum - NORM";
    constant Real Fe_bm_out_RBC_norm = 3.663764 "Fe flow from bone marrow to RBC - NORM";
    constant Real Fe_bm_out_spl_norm = 0.61149263 "Fe flow from bone marrow to spleen - NORM";

    // Red blood cells (RBC)

    constant Real Fe_RBC_in_bm_norm = 3.663764 "Fe flow to RBC from bone marrow - NORM";
    constant Real Fe_RBC_out_spl_norm = 3.663764 "Fe flow from RBC to spleen - NORM";

    //-Constants

    //+Relative indicators

    //################################################
    // Main relative indicators in ODES
    //################################################

    // Hepcidin expression
    //################################################
    Real hep_rel "Hepcidin amount (s16 ) - REL";
    Real Bmp6_rel "Bmp6 amount (s18) - REL";

    // Ferroportin regulation
    //################################################
    Real Fpn_liv_mRNA_rel "Fpn mRNA amount in liver (s7) - REL";
    Real Fpn_spl_mRNA_rel "Fpn mRNA amount in spleen (s9) - REL";
    Real Fpn_duo_mRNA_rel "Fpn mRNA amount in duodenum (s8) - REL";
    Real Fpn_res_mRNA_rel "Fpn mRNA amount in other organs (s11) - REL";

    Real Fpn_liv_rel "Fpn amount in liver (s12) - REL";
    Real Fpn_spl_rel "Fpn amount in spleen (s14) - REL";
    Real Fpn_duo_rel "Fpn amount in duodenum (s13) - REL";
    Real Fpn_res_rel "Fpn amount in other organs (s15) - REL";

    // Dynamics of the iron pools
    //################################################
    Real Fe_liv_rel "Fe amount in liver (s2) - REL";
    Real Fe_spl_rel "Fe amount in spleen (s5) - REL";
    Real Fe_bm_rel "Fe amount in bones (s3) - REL";
    Real Fe_RBC_rel "Fe amount in red blood cells (s6) - REL";
    Real Fe_duo_rel "Fe amount in duodenum (s4) - REL";
    Real Fe_res_rel "Fe amount in other organs (s10) - REL";
    Real Fe_ser_rel "Fe amount in serum (s1) - REL";

    Real Fe_total_rel "Total Fe amount - REL";

    Real Fe_ser_input_rel "Fe amount in serum: input - REL";
    Real Fe_ser_output_rel "Fe amount in serum: output - REL";

    // ###############################################
    // Auxiliary relative indicators
    // ###############################################

    // Duodenum

    Real Fe_duo_2_rel "Fe++ in duodenum - REL";
    Real Fe_duo_3_rel "Fe+++ in duodenum, basically in ferritin - REL";

    Real Fe_duo_to_ferritin_rel;
    Real Fe_duo_from_ferritin_rel;

    Real Fe_duo_intake_rel "Fe intake from food - REL";

    Real Fe_duo_in_food_rel "Fe input to duodenum from intestines, prev:in_1 - REL";
    Real Fe_duo_in_ser_rel "Fe input to duodenum from serum, prev:in_2 - REL";
    Real Fe_duo_out_ser_rel "Fe output from duodenum to serum, prev:out_1 - REL";
    Real Fe_duo_out_loss_rel "Fe output from duodenum to intestines (Fe loss), prev:out_2 - REL";
    Real Fe_duo_unused_rel "Unused Fe going through intestines - REL";

    Real Fpn_duo_in_1_rel "Fpn in duodenum creation: standard from mRNA - REL";
    Real Fpn_duo_in_2_rel "Fpn in duodenum creation by activation of Fe_duo - REL";
    Real Fpn_duo_out_1_rel "Degradation (internalization) of Fpn in duodenum: standard degradation rate - REL";
    Real Fpn_duo_out_2_rel "Degradation (internalization) of Fpn in duodenum by hepcidine - REL";

    Real Fpn_duo_in_rel "Fpn in duodenum synthesis - REL";
    Real Fpn_duo_out_rel "Fpn in duodenum degradation - REL";

    Real Fpn_duo_mRNA_in_rel "Fpn mRNA in duodenum synthesis - REL";
    Real Fpn_duo_mRNA_out_rel "Fpn mRNA in duodenum degradation - REL";

    // Liver

    Real Fe_liv_2_rel "Fe++ in liver - REL";
    Real Fe_liv_3_rel "Fe+++ in liver, basically in ferritin - REL";

    Real Fe_liv_to_ferritin_rel;
    Real Fe_liv_from_ferritin_rel;

    Real Fe_liv_in_ser_rel "Fe input to liver from serum - REL";
    Real Fe_liv_out_ser_rel "Fe output from liver to serum - REL";

    Real Fpn_liv_in_1_rel "Fpn in liver creation: standard from mRNA - REL";
    Real Fpn_liv_in_2_rel "Fpn in liver creation by activation of Fe_liv - REL";
    Real Fpn_liv_out_1_rel "Degradation (internalization) of Fpn in liver: standard degradation rate - REL";
    Real Fpn_liv_out_2_rel "Degradation (internalization) of Fpn in liver by hepcidine - REL";

    Real Fpn_liv_in_rel "Fpn in liver synthesis - REL";
    Real Fpn_liv_out_rel "Fpn in liver degradation - REL";

    Real Fpn_liv_mRNA_in_rel "Fpn mRNA in liver synthesis - REL";
    Real Fpn_liv_mRNA_out_rel "Fpn mRNA in liver degradation - REL";

    // Serum

    Real Fe_ser_in_liv_rel "Fe flow from liver to serum - REL";
    Real Fe_ser_in_spl_rel "Fe flow from spleen to serum - REL";
    Real Fe_ser_in_duo_rel "Fe flow from duodenum to serum - REL";
    Real Fe_ser_in_res_rel "Fe flow from other organs to serum - REL";

    Real Fe_ser_out_liv_rel "Fe flow from serum to liver - REL";
    Real Fe_ser_out_bm_rel "Fe flow from serum to bone marrow - REL";
    Real Fe_ser_out_duo_rel "Fe flow from serum to duodenum - REL";
    Real Fe_ser_out_res_rel "Fe flow from serum to other organs - REL";

    // Hepcidine related stuff

    Real hep_in_rel "Production of hepcidin - REL";
    Real hep_out_rel "Degradation of hepcidin - REL";

    Real Bmp6_in_rel "Bmp6 production rate - REL";
    Real Bmp6_out_rel "Bmp6 degradation rate - REL";

    // Spleen

    Real Fe_spl_2_rel "Fe++ in spleen - REL";
    Real Fe_spl_3_rel "Fe++ in spleen - REL";

    Real Fe_spl_to_ferritin_rel;
    Real Fe_spl_from_ferritin_rel;

    Real Fe_spl_in_RBC_rel "Fe flow to spleen from RBC - REL";
    Real Fe_spl_in_bm_rel "Fe flow to spleen from bone marrow - REL";
    Real Fe_spl_out_ser_rel "Fe flow from spleen to serum - REL";

    Real Fpn_spl_in_1_rel "Fpn in spleen creation: standard from mRNA - REL";
    Real Fpn_spl_in_2_rel "Fpn in spleen creation by activation of Fe_spl - REL";
    Real Fpn_spl_out_1_rel "Degradation (internalization) of Fpn in spleen: standard degradation rate - REL";
    Real Fpn_spl_out_2_rel "Degradation (internalization) of Fpn in spleen by hepcidine - REL";

    Real Fpn_spl_in_rel "Fpn in spleen synthesis - REL";
    Real Fpn_spl_out_rel "Fpn in spleen degradation - REL";

    Real Fpn_spl_mRNA_in_rel "Fpn mRNA in spleen synthesis - REL";
    Real Fpn_spl_mRNA_out_rel "Fpn mRNA in spleen degradation - REL";

    // Other organs

    Real Fe_res_2_rel "Fe++ in other organs - REL";
    Real Fe_res_3_rel "Fe++ in other organs - REL";

    Real Fe_res_to_ferritin_rel;
    Real Fe_res_from_ferritin_rel;

    Real Fe_res_in_ser_rel "Fe flow to other organs from serum - REL";
    Real Fe_res_out_ser_rel "Fe flow from other organs to serum - REL";
    Real Fe_res_out_loss_rel "Fe loss from other organs (skin peeling, etc.) - REL";

    Real Fpn_res_in_1_rel "Fpn in other organs creation: standard from mRNA - REL";
    Real Fpn_res_in_2_rel "Fpn in other organs creation by activation of Fe_res - REL";
    Real Fpn_res_out_1_rel "Degradation (internalization) of Fpn in other organs: standard degradation rate - REL";
    Real Fpn_res_out_2_rel "Degradation (internalization) of Fpn in other organs by hepcidine - REL";

    Real Fpn_res_in_rel "Fpn in spleen synthesis - REL";
    Real Fpn_res_out_rel "Fpn in spleen degradation - REL";

    Real Fpn_res_mRNA_in_rel "Fpn mRNA in other organs synthesis - REL";
    Real Fpn_res_mRNA_out_rel "Fpn mRNA in other organs degradation - REL";

    // Bone marrow

    Real Fe_bm_in_ser_rel "Fe flow to bone marrow from serum - REL";
    Real Fe_bm_out_RBC_rel "Fe flow from bone marrow to RBC - REL";
    Real Fe_bm_out_spl_rel "Fe flow from bone marrow to spleen - REL";

    // Red blood cells (RBC)

    Real Fe_RBC_in_bm_rel "Fe flow to RBC from bone marrow - REL";
    Real Fe_RBC_out_spl_rel "Fe flow from RBC to spleen - REL";

    //-Relative indicators

    //+1.06

    parameter Real hep_regulation = 1 "Hepcidin regulation: 1 - auto, 0 - manual";
    parameter Real hep_manual(unit="ug") = 0.664938;

    Real hep_auto(start=hep_norm,unit="ug");
    Real hep_regulation_switcher;

    //-1.06

    //################################################
    // Main varibles in ODES
    //################################################

    // Hepcidin expression
    //################################################
    Real hep(start=0.664938,unit="ug") "Hepcidin amount (s16 )";
    Real Bmp6(start=15.8291,unit="ug") "Bmp6 amount (s18)";
    Real LPS(start=0.000,fixed=true,unit="ug") "LPS (Lipopolysaccharide) amount (s23)"; //1.00 hodnota z clanku
    Real Il6mRNA(start=0.000,unit="ug") "IL-6 mRNA amount (s21)"; //NESEDI JEDNOTKY !!!
    Real Il6(start=0.000,unit="ug") "IL-L amount (s22)";

    // Ferroportin regulation
    //################################################
    Real Fpn_liv_mRNA(start=0.922424,unit="ug") "Fpn mRNA amount in liver (s7)";
    Real Fpn_spl_mRNA(start=0.922424,unit="ug") "Fpn mRNA amount in spleen (s9)";
    Real Fpn_duo_mRNA(start=0.922424,unit="ug") "Fpn mRNA amount in duodenum (s8)";
    Real Fpn_res_mRNA(start=0.922424,unit="ug") "Fpn mRNA amount in other organs (s11)";

    Real Fpn_liv(start=1.00036,unit="ug") "Fpn amount in liver (s12)";
    Real Fpn_spl(start=1.00049,unit="ug") "Fpn amount in spleen (s14)";
    Real Fpn_duo(start=1.00017,unit="ug") "Fpn amount in duodenum (s13)";
    Real Fpn_res(start=1.00020,unit="ug") "Fpn amount in other organs (s15)";

    // Dynamics of the iron pools
    //################################################
    Real Fe_liv(start=76.9562,unit="ug") "Fe amount in liver (s2)";
    Real Fe_spl(start=17.7295,unit="ug") "Fe amount in spleen (s5)";
    Real Fe_bm(start=63.1596,unit="ug") "Fe amount in bones (s3)";
    Real Fe_RBC(start=1016.720,unit="ug") "Fe amount in red blood cells (s6)";
    Real Fe_duo(start=2.97118,unit="ug") "Fe amount in duodenum (s4)";
    Real Fe_res(start=466.811,unit="ug") "Fe amount in other organs (s10)";
    Real Fe_ser(start=1.51304,unit="ug") "Fe amount in serum (s1)";

    Real Fe_total(unit="ug") = Fe_liv + Fe_spl + Fe_bm + Fe_RBC + Fe_duo + Fe_res + Fe_ser "Total Fe amount";

    Real Fe_ser_input(unit="ug.h-1") "Fe amount in serum: input";
    Real Fe_ser_output(unit="ug.h-1") "Fe amount in serum: output";

    //################################################
    // Model parameters
    // Format: rate constants (k): k^{X}_{Y} -->> k_{X}_{Y}; (K): K^{X}_{Y} -->> K_{X}_{Y}; MIN/MAX value -->> MIN/MAX in description
    //################################################

    // Hepcidin expression
    //################################################
    parameter Real k_hep_deg(unit="h-1") = 0.07 "Hepcidin degradation rate (k20), 0.067/0.070";
    parameter Real K_Bmp6(unit="ug") = 19.65 "Michaelis-Menten constant Bmp6 synthesis (k37), 16.5/55.7";
    parameter Real v_Bmp6_max(unit="h-1") = 1.6015*K_Bmp6 "Bmp6 maximal synthesis rate (k21*k37 = k21*K_Bmp6) k21 = 1.6015, K_Bmp6 = 19.65, 14.2/126.5";
    parameter Real Tf(unit="ug") = 1000.0 "Paremeter determining the maximal amount of iron that can be bound to transferrin (k32)";
    parameter Real k_Bmp6_deg(unit="h-1") = 2.3962 "Bmp6 degradation rate (k22), 1.0/9.5";
    parameter Real k_LPS_deg(unit="h-1") = 5.8560 "LPS degradation rate, 5.9/5.9";
    parameter Real K_Il6mRNA(unit="ug") = 2.6e-6 "Michaelis-Menten constant Il6mRNA synthesis, 2.6-e6/2.6e-6";
    parameter Real k_Il6mRNA_deg(unit="h-1") = 0.2814 "Il6mRNA degradation rate, 0.28/0.28";
    parameter Real k_Il6_syn(unit="h-1") = 4.1067*157.4 "Il6 synthesis rate, 4.1067*k23, k23 = 157.4, 136/872"; // NESEDI JEDNOTKY!!!
    parameter Real k_Il6_deg(unit="h-1") = 4.4465 "Il6 degradation rate, 4.45/4.45";

    // Ferroportin regulation
    //################################################
    parameter Real K_2(unit="a.u.") = 0.0012836 "Constant FpnmRNA production (k41), 3.0e-4/2.0e-3";
    parameter Real k_FpnmRNA_deg(unit="h-1") = 1.0841 "FpnmRNA degradation rate (k25), 1.02/1.20";
    parameter Real K_liv_1(unit="a.u.") = 30.66 "Constant Fpn_liv_mRNA production (k24), 28.0/32.5";
    parameter Real K_spl_1(unit="a.u.") = 1.0867*K_liv_1 "Constant Fpn_spl_mRNA production (k36*K_liv_1, k36 = 1.086700), 33.1/34.5";
    parameter Real K_duo_1(unit="a.u.") = 0.020001*K_liv_1 "Constant Fpn_duo_mRNA production (k35*K_liv_1, k35 = 0.020001), 0.56/0.79";
    parameter Real K_res_1(unit="a.u.") = 0.36629*K_liv_1 "Constant Fpn_res_mRNA production (k42*K_liv_1, k42 = 0.366290), 7.80/43.7";
    parameter Real k_Fpnliv_syn(unit="h-1") = 0.1297 "Fpnliv synthesis rate (k28), 0.07/0.14";
    parameter Real k_Fpnspl_syn(unit="h-1") = 0.022722 "Fpnspl synthesis rate (k30), 0.015/0.027";
    parameter Real k_Fpnduo_syn(unit="h-1") = 0.030299 "Fpnduo synthesis rate (k29), 0.01/0.25";
    parameter Real k_Fpnres_syn(unit="h-1") = 0.0050598 "Fpnres synthesis rate (k40), 0.004/0.108";
    parameter Real k_liv_1(unit="ug-1") = 0.0033177 "Constant Fpn_liv production (k17), 0.002/0.006";
    parameter Real k_spl_1(unit="ug-1") = 0.014027 "Constant Fpn_spl production (k27), 0.005/0.028";
    parameter Real k_duo_1(unit="ug-1") = 0.16007 "Constant Fpn_duo production (k46), 0.06/0.49";
    parameter Real k_res_1(unit="ug-1") = 0.11376 "Constant Fpn_res production (k47), 0.004/0.152";
    parameter Real k_Fpnliv_deg(unit="h-1") = 0.055363 "Fpnliv degradation rate (k12), 0.01/0.06";
    parameter Real k_Fpnspl_deg(unit="h-1") = k_Fpnliv_deg*0.054621 "Fpnspl degradation rate (k12*k18, k18 = 0.054621), 0.0007/0.0038";
    parameter Real k_Fpnduo_deg(unit="h-1") = k_Fpnliv_deg*0.3815 "Fpnduo degradation rate (k12*k15, k15 = 0.3815), 0.0056/0.147";
    parameter Real k_Fpnres_deg(unit="h-1") = k_Fpnliv_deg*0.52232 "Fpnres degradation rate (k12*k38, k38 = 0.52232), 0.025/0.129";
    parameter Real k_liv_2(unit="a.u.") = 2.5743 "Constant Fpnliv degradation (k13), 2.11/12.93";
    parameter Real k_spl_2(unit="a.u.") = k_liv_2*4.4694 "Constant Fpnspl degradation (k13*k19, k19 = 4.4694), 9.2/73.7";
    parameter Real k_duo_2(unit="a.u.") = k_liv_2*0.55631 "Constant Fpnduo degradation (k13*k16, k16 = 0.55631), 0.78/4.16";
    parameter Real k_res_2(unit="a.u.") = k_liv_2*4.5163 "Constant Fpnres degradation (k13*k39, k39 = 4.5163), 3.7/38.1";

    // Dynamics of the iron pools
    //################################################
    parameter Real v_liv_1(unit="ug.h-1") = 3.9607 "Low liver iron uptake (k4), 2.78/9.81";
    parameter Real v_liv_2(unit="ug.h-1") = 14.3810 "High liver iron uptake (k45)"; //!!! Inconsistecy xml model vs. paper !!!, xml model value adopted;
    parameter Real th(unit="ug") = 2.6870 "Threshold serum iron value (k26), 2.08/3.00";
    parameter Real u_liv(unit="h-1") = 0.077844 "Liver iron export rate (k1), 0.05/0.19";
    parameter Real Fe_liv_max(unit="uq") = 119.55 "Threshold value liver iron export (k43), 100.0/159.0";
    parameter Real v_spl_1(unit="h-1") = 0.0036035 "Spleen iron uptake rate from RBC (k10), 0.002 = value for trace exp., 0.004/0.005";
    parameter Real v_spl_2(unit="h-1") = 0.0096817 "Spleen iron uptake rate from bones (k7), 0.008/0.019";
    parameter Real u_spl(unit="h-1") = 0.24102 "Spleen export rate (k3), 0.21/0.36";
    parameter Real Fe_spl_max(unit="ug") = 88.216 "Threshold value spleen iron export (k44), 57.0/95.0";
    parameter Real v_bm(unit="h-1") = 2.8256 "Bone marrow uptake rate (k5), 2.66/4.16";
    parameter Real v_RBC(unit="h-1") = 0.058008 "RBC uptake rate (k8), 0.055/0.075"; // !!! Problem u rovnice (14) - doplnit prvni term na prave strane
    parameter Real k14(unit="1") = 0.19419 "k14 constant";
    parameter Real k31(unit="1") = 33.333 "k31 constant";
    parameter Real K_duo(unit="ug.h-1") = 177.34 "Saturation parameter duodenal uptake, corresponds to k11, units??"; // !!! Problem u rovnice (14) - doplnit prvni term na prave strane
    parameter Real v_duo(unit="h-1") = 0.70485 "Duodenal uptake rate from blood (k6), 0.6/1.3";
    parameter Real u_duo(unit="h-1") = 0.88356 "Duodenal export rate (k2), 0.66/1.38";
    parameter Real v_res(unit="h-1") = 6.3206 "Other organs uptake rate (k34), 5.4/10.0";
    parameter Real u_res(unit="h-1") = 0.017143 "Other organs export rate (k33), 0.014/0.030";
    parameter Real v_duo_lost(unit="h-1") = 0.091919 "Iron lost rate duodenum (k48), 0.001/0.320";
    parameter Real u_res_lost(unit="h-1") = 0.0033401 "Iron lost rate rest (k9), 0.002/0.004";
    parameter Real Fe_res_max(unit="ug") = 510.68 "Limit value, iron lost rest (k49), 510/947";

    // Food income parameters
    //################################################
    parameter Real v_duo_max(unit="uq.h-1") = 9.86 "Maximal duodenal uptake from food (SI. Eq.14)";
    parameter Real Fe_food(unit="ug.h-1") = 218.6322 "Food iron content (SI. Eq.14), default: 218.6322";


    // Pathophysiology parameters
    //################################################

     type exprese = enumeration(
        Nulova
           "Nulova",
        Polovicni
           "Polovicni",
        Fyziologicka
           "Fyziologicka",
        Zvysena
           "Zvysena");

    // Iron deficiency (anemia)
     parameter Real bleeding(unit="ug.h-1") = 0 "Bleeding, to be set";
     parameter Real malabsorption = 1 "iron malabsorption coefficient: <0;1>; 0 = no absorption, 1 = physiologic";
     parameter Real loss_factor = 1 "iron loss factor, 1 = physiologic, >1 enhanced loss";

    // Iron overload (hemochromatosis)
     parameter Real transfusion(unit="ug.h-1") = 0 "transfusion rate, to be set";
     parameter Boolean unregulated_absorption = false;

    // Hepcidine-related
     parameter Real hep_knockout = 1;
    // parameter exprese hep_exp = exprese.Fyziologicka; //defaultni hodnota, mozne nastaveni: Fyziologicka, Nulova, Polovicni, Zvysena
    // parameter Real hep_knockout=
    //   if hep_exp == exprese.Nulova then 0
    //   elseif hep_exp == exprese.Polovicni then 0.5
    //   elseif hep_exp == exprese.Fyziologicka then 1
    //   else 2;

    // Fpn-related
     parameter Real Fpn_duo_knockout = 1;
    // parameter exprese Fpn_duo_exp = exprese.Fyziologicka; //defaultni hodnota, mozne nastaveni: Fyziologicka, Nulova, Polovicni, Zvysena
    // parameter Real Fpn_duo_knockout=
    //   if Fpn_duo_exp == exprese.Nulova then 0
    //   elseif Fpn_duo_exp == exprese.Polovicni then 0.5
    //   elseif Fpn_duo_exp == exprese.Fyziologicka then 1
    //   else 2;

     parameter Real Fpn_liv_knockout = 1;
    // parameter exprese Fpn_liv_exp = exprese.Fyziologicka; //defaultni hodnota, mozne nastaveni: Fyziologicka, Nulova, Polovicni, Zvysena
    // parameter Real Fpn_liv_knockout=
    //   if Fpn_liv_exp == exprese.Nulova then 0
    //   elseif Fpn_liv_exp == exprese.Polovicni then 0.5
    //   elseif Fpn_liv_exp == exprese.Fyziologicka then 1
    //   else 2;

     parameter Real Fpn_spl_knockout = 1;
    // parameter exprese Fpn_spl_exp = exprese.Fyziologicka; //defaultni hodnota, mozne nastaveni: Fyziologicka, Nulova, Polovicni, Zvysena
    // parameter Real Fpn_spl_knockout=
    //   if Fpn_spl_exp == exprese.Nulova then 0
    //   elseif Fpn_spl_exp == exprese.Polovicni then 0.5
    //   elseif Fpn_spl_exp == exprese.Fyziologicka then 1
    //   else 2;

     parameter Real Fpn_res_knockout = 1;
    // parameter exprese Fpn_res_exp = exprese.Fyziologicka; //defaultni hodnota, mozne nastaveni: Fyziologicka, Nulova, Polovicni, Zvysena
    // parameter Real Fpn_res_knockout=
    //   if Fpn_res_exp == exprese.Nulova then 0
    //   elseif Fpn_res_exp == exprese.Polovicni then 0.5
    //   elseif Fpn_res_exp == exprese.Fyziologicka then 1
    //   else 2;


    // Auxiliary variables
    // ###############################################

    // General

    parameter Real to_ferritin_rate = 3 "rate of flow Fe++ -> Fe+++";
    parameter Real from_ferritin_rate = to_ferritin_rate/3 "rate of flow Fe+++ -> Fe++";

    // Duodenum

    Real Fe_duo_2(start=0.7427,unit="ug") "Fe++ in duodenum";
    Real Fe_duo_3(start=2.2283,unit="ug") "Fe+++ in duodenum, basically in ferritin";

    Real Fe_duo_to_ferritin = to_ferritin_rate * Fe_duo_2;
    Real Fe_duo_from_ferritin = from_ferritin_rate * Fe_duo_3;

    Real Fe_duo_intake = v_duo_max*Fe_food/(Fe_food+K_duo) "Fe intake from food";

    Real Fe_duo_in_food =   (if unregulated_absorption then Fe_duo_intake else Fe_duo_intake*min(1/Fe_duo, 1)*malabsorption) "Fe input to duodenum from intestines, prev:in_1";
    Real Fe_duo_in_ser =   v_duo*Fe_ser "Fe input to duodenum from serum, prev:in_2";
    Real Fe_duo_out_ser =  u_duo*Fe_duo*Fpn_duo "Fe output from duodenum to serum, prev:out_1";
    Real Fe_duo_out_loss =  v_duo_lost*Fe_duo "Fe output from duodenum to intestines (Fe loss), prev:out_2";
    Real Fe_duo_unused = Fe_duo_intake*(1-min(1/Fe_duo, 1)) "Unused Fe going through intestines";

    Real Fpn_duo_in_1 = k_Fpnduo_syn*Fpn_duo_mRNA "Fpn in duodenum creation: standard from mRNA";
    Real Fpn_duo_in_2 = k_Fpnduo_syn*k_duo_1*Fe_duo*Fpn_duo_mRNA "Fpn in duodenum creation by activation of Fe_duo";
    Real Fpn_duo_out_1 = k_Fpnduo_deg*Fpn_duo "Degradation (internalization) of Fpn in duodenum: standard degradation rate";
    Real Fpn_duo_out_2 = k_Fpnduo_deg*k_duo_2*hep*Fpn_duo "Degradation (internalization) of Fpn in duodenum by hepcidine";

    Real Fpn_duo_in =    Fpn_duo_in_1 + Fpn_duo_in_2 "Fpn in duodenum synthesis";
    Real Fpn_duo_out =   Fpn_duo_out_1 + Fpn_duo_out_2 "Fpn in duodenum degradation";

    Real Fpn_duo_mRNA_in =  Fpn_duo_knockout/(1 + Fpn_duo_mRNA_inhib) "Fpn mRNA in duodenum synthesis";
    Real Fpn_duo_mRNA_out = k_FpnmRNA_deg*Fpn_duo_mRNA "Fpn mRNA in duodenum degradation";
    Real Fpn_duo_mRNA_inhib = (K_duo_1*Il6)/(K_2 + Il6) "Fpn mRNA in duodenum production inhibition by IL-6 (or indirectly by LPS)";

    // Liver

    Real Fe_liv_2(start=19.2487,unit="ug") "Fe++ in liver";
    Real Fe_liv_3(start=57.7461,unit="ug") "Fe+++ in liver, basically in ferritin";

    Real Fe_liv_to_ferritin = to_ferritin_rate * Fe_liv_2;
    Real Fe_liv_from_ferritin = from_ferritin_rate * Fe_liv_3;

    Real Fe_liv_in_ser = v_liv_1*min(Fe_ser, th) + v_liv_2*max(Fe_ser - th, 0) "Fe input to liver from serum";
    Real Fe_liv_out_ser = u_liv*min(Fe_liv, Fe_liv_max)*Fpn_liv "Fe output from liver to serum";
    Real Fpn_liv_in_1 = k_Fpnliv_syn*Fpn_liv_mRNA "Fpn in liver creation: standard from mRNA";
    Real Fpn_liv_in_2 = k_Fpnliv_syn*k_liv_1*Fe_liv*Fpn_liv_mRNA "Fpn in liver creation by activation of Fe_liv";
    Real Fpn_liv_out_1 = k_Fpnliv_deg*Fpn_liv "Degradation (internalization) of Fpn in liver: standard degradation rate";
    Real Fpn_liv_out_2 = k_Fpnliv_deg*k_liv_2*hep*Fpn_liv "Degradation (internalization) of Fpn in liver by hepcidine";

    Real Fpn_liv_in = Fpn_liv_in_1 + Fpn_liv_in_2 "Fpn in liver synthesis";
    Real Fpn_liv_out = Fpn_liv_out_1 + Fpn_liv_out_2 "Fpn in liver degradation";

    Real Fpn_liv_mRNA_in = Fpn_liv_knockout/(1 + Fpn_liv_mRNA_inhib) "Fpn mRNA in liver synthesis";
    Real Fpn_liv_mRNA_out = k_FpnmRNA_deg*Fpn_liv_mRNA "Fpn mRNA in liver degradation";
    Real Fpn_liv_mRNA_inhib = (K_liv_1*Il6)/(K_2 + Il6) "Fpn mRNA in liver production inhibition by IL-6 (or indirectly by LPS)";

    // Serum

    Real Fe_ser_in_liv = u_liv*MIN(Fe_liv, Fe_liv_max)*Fpn_liv "Fe flow from liver to serum";
    Real Fe_ser_in_spl = u_spl*MIN(Fe_spl,Fe_spl_max)*Fpn_spl "Fe flow from spleen to serum";
    Real Fe_ser_in_duo = u_duo*Fe_duo*Fpn_duo "Fe flow from duodenum to serum";
    Real Fe_ser_in_res = u_res*Fe_res*Fpn_res "Fe flow from other organs to serum";

    Real Fe_ser_out_liv = v_liv_1*MIN(Fe_ser, th) + v_liv_2*MAX(Fe_ser - th, 0) "Fe flow from serum to liver";
    Real Fe_ser_out_bm =  v_bm*Fe_ser "Fe flow from serum to bone marrow";
    Real Fe_ser_out_duo = v_duo*Fe_ser "Fe flow from serum to duodenum";
    Real Fe_ser_out_res = v_res*Fe_ser "Fe flow from serum to other organs";

    // Hepcidine related stuff

    Real hep_in =  Promoter(Il6, Bmp6)*hep_knockout "Production of hepcidin";
    Real hep_out = k_hep_deg*hep "Degradation of hepcidin";

    Real Il6mRNA_in =  LPS/(LPS + K_Il6mRNA) "Production of IL-6 mRNA";
    Real Il6mRNA_out = k_Il6mRNA_deg*Il6mRNA "Degradation of IL-6 mRNA";

    Real Il6_in =  k_Il6_syn*(Il6mRNA^4) "Production of IL-6";
    Real Il6_out = k_Il6_deg*Il6 "Degradation of IL-6";

    Real Bmp6_in = v_Bmp6_max*(Fe_liv/(K_Bmp6 + Fe_liv))*MIN(Fe_ser, Tf) "Bmp6 production rate";
    Real Bmp6_out = k_Bmp6_deg*Bmp6 "Bmp6 degradation rate";

    // Spleen

    Real Fe_spl_2(start=4.43525,unit="ug") "Fe++ in spleen";
    Real Fe_spl_3(start=13.30575,unit="ug") "Fe++ in spleen";

    Real Fe_spl_to_ferritin = to_ferritin_rate * Fe_spl_2;
    Real Fe_spl_from_ferritin = from_ferritin_rate * Fe_spl_3;

    Real Fe_spl_in_RBC = v_spl_1*Fe_RBC "Fe flow to spleen from RBC";
    Real Fe_spl_in_bm = v_spl_2*Fe_bm "Fe flow to spleen from bone marrow";
    Real Fe_spl_out_ser = u_spl*MIN(Fe_spl, Fe_spl_max)*Fpn_spl "Fe flow from spleen to serum";

    Real Fpn_spl_in_1 = k_Fpnspl_syn*Fpn_spl_mRNA "Fpn in spleen creation: standard from mRNA";
    Real Fpn_spl_in_2 = k_Fpnspl_syn*k_spl_1*Fe_spl*Fpn_spl_mRNA "Fpn in spleen creation by activation of Fe_spl";
    Real Fpn_spl_out_1 = k_Fpnspl_deg*Fpn_spl "Degradation (internalization) of Fpn in spleen: standard degradation rate";
    Real Fpn_spl_out_2 = k_Fpnspl_deg*k_spl_2*hep*Fpn_spl "Degradation (internalization) of Fpn in spleen by hepcidine";

    Real Fpn_spl_in = Fpn_spl_in_1 + Fpn_spl_in_2 "Fpn in spleen synthesis";
    Real Fpn_spl_out = Fpn_spl_out_1 + Fpn_spl_out_2 "Fpn in spleen degradation";

    Real Fpn_spl_mRNA_in = Fpn_spl_knockout/(1 + Fpn_spl_mRNA_inhib) "Fpn mRNA in spleen synthesis";
    Real Fpn_spl_mRNA_out = k_FpnmRNA_deg*Fpn_spl_mRNA "Fpn mRNA in spleen degradation";
    Real Fpn_spl_mRNA_inhib = (K_spl_1*Il6)/(K_2 + Il6) "Fpn mRNA in spleen production inhibition by IL-6 (or indirectly by LPS)";

    // Other organs

    Real Fe_res_2(start=116.74625,unit="ug") "Fe++ in other organs";
    Real Fe_res_3(start=350.23875,unit="ug") "Fe++ in other organs";

    Real Fe_res_to_ferritin = to_ferritin_rate * Fe_res_2;
    Real Fe_res_from_ferritin = from_ferritin_rate * Fe_res_3;

    Real Fe_res_in_ser = v_res*Fe_ser "Fe flow to other organs from serum";
    Real Fe_res_out_ser = u_res*Fe_res*Fpn_res "Fe flow from other organs to serum";
    Real Fe_res_out_loss = u_res_lost*MIN(Fe_res, Fe_res_max)*loss_factor "Fe loss from other organs (skin peeling, etc.)";

    Real Fpn_res_in_1 = k_Fpnres_syn*Fpn_res_mRNA "Fpn in other organs creation: standard from mRNA";
    Real Fpn_res_in_2 = k_Fpnres_syn*k_res_1*Fe_res*Fpn_res_mRNA "Fpn in other organs creation by activation of Fe_res";
    Real Fpn_res_out_1 = k_Fpnres_deg*Fpn_res "Degradation (internalization) of Fpn in other organs: standard degradation rate";
    Real Fpn_res_out_2 = k_Fpnres_deg*k_res_2*hep*Fpn_res "Degradation (internalization) of Fpn in other organs by hepcidine";

    Real Fpn_res_in = Fpn_res_in_1 + Fpn_res_in_2 "Fpn in spleen synthesis";
    Real Fpn_res_out = Fpn_res_out_1 + Fpn_res_out_2 "Fpn in spleen degradation";

    Real Fpn_res_mRNA_in = Fpn_res_knockout/(1 + Fpn_res_mRNA_inhib) "Fpn mRNA in other organs synthesis";
    Real Fpn_res_mRNA_out = k_FpnmRNA_deg*Fpn_res_mRNA "Fpn mRNA in other organs degradation";
    Real Fpn_res_mRNA_inhib = (K_res_1*Il6)/(K_2 + Il6) "Fpn mRNA in other organs production inhibition by IL-6 (or indirectly by LPS)";

    // Bone marrow

    Real Fe_bm_in_ser = v_bm*Fe_ser "Fe flow to bone marrow from serum";
    Real Fe_bm_out_RBC = v_RBC*Fe_bm "Fe flow from bone marrow to RBC";
    Real Fe_bm_out_spl = v_spl_2*Fe_bm "Fe flow from bone marrow to spleen";

    // Red blood cells (RBC)

    Real Fe_RBC_in_bm = v_RBC*Fe_bm "Fe flow to RBC from bone marrow";
    Real Fe_RBC_out_spl = v_spl_1*Fe_RBC "Fe flow from RBC to spleen";

    //################################################
    // INITIAL EQUATIONS
    //################################################

  initial equation

    // Hepcidine-related stuff

    //+1.06

    //der(hep) = 0; // eq. 1

    der(hep_auto) = 0; // eq. 1

    //-1.06

    der(Bmp6) = 0; // eq. 2
    //der(LPS) = 0; // eq. 3
    der(Il6mRNA) = 0; // eq 4.
    der(Il6) = 0; // eq 5.

    // Duodenum

    der(Fe_duo_2) = 0; // eq 14.
    der(Fe_duo_3) = 0;
    der(Fpn_duo) = 0; // eq 7.3
    der(Fpn_duo_mRNA) = 0; // eq 6.3

    // Serum

    der(Fe_ser) = 0; // eq 16.

    // Liver

    der(Fe_liv_2) = 0;
    der(Fe_liv_3) = 0;
    der(Fpn_liv) = 0; // eq 7.1
    der(Fpn_liv_mRNA) = 0; // eq 6.1

    // Spleen

    der(Fe_spl_2) = 0;
    der(Fe_spl_3) = 0;
    der(Fpn_spl) = 0; // eq 7.2
    der(Fpn_spl_mRNA) = 0; // eq 6.2

    // Other organs

    der(Fe_res_2) = 0;
    der(Fe_res_3) = 0;
    der(Fpn_res) = 0; // eq 7.4
    der(Fpn_res_mRNA) = 0; // eq 6.4

    // Bone marrow

    der(Fe_bm) = 0; //eq 12.

    // Red blood cells

    der(Fe_RBC) = 0; //eq 13.

    //################################################
    // EQUATIONS
    //################################################

  equation

    // Hepcidine-related stuff

    //+1.06

    //der(hep) = hep_in - hep_out; // eq. 1

    hep_regulation_switcher = smooth(1, noEvent(if hep_regulation > 0.5 then 1.0 else 0));
    der(hep_auto) = hep_regulation_switcher * (hep_in - hep_out);
    hep = hep_regulation_switcher * hep_auto + (1 - hep_regulation_switcher) * hep_manual;

    //-1.06

    der(Bmp6) = Bmp6_in - Bmp6_out; // eq. 2
    der(LPS) = -k_LPS_deg*LPS; // eq. 3
    der(Il6mRNA) = Il6mRNA_in - Il6mRNA_out; // eq 4.
    der(Il6) = Il6_in - Il6_out; // eq 5.

    // Duodenum

    Fe_duo = Fe_duo_2 + Fe_duo_3;
    der(Fe_duo_2) = Fe_duo_in_food + Fe_duo_in_ser - Fe_duo_out_ser - Fe_duo_out_loss + Fe_duo_from_ferritin - Fe_duo_to_ferritin; // eq 14.
    der(Fe_duo_3) = Fe_duo_to_ferritin - Fe_duo_from_ferritin;
    der(Fpn_duo) = Fpn_duo_in - Fpn_duo_out; // eq 7.3
    der(Fpn_duo_mRNA) = Fpn_duo_mRNA_in - Fpn_duo_mRNA_out; // eq 6.3

    // Serum

    Fe_ser_input = Fe_ser_in_liv + Fe_ser_in_spl + Fe_ser_in_duo + Fe_ser_in_res + transfusion;
    Fe_ser_output = Fe_ser_out_liv + Fe_ser_out_bm + Fe_ser_out_duo + Fe_ser_out_res + bleeding;

    der(Fe_ser) = Fe_ser_input - Fe_ser_output; // eq 16.

    // Liver

    Fe_liv = Fe_liv_2 + Fe_liv_3;
    der(Fe_liv_2) = Fe_liv_in_ser - Fe_liv_out_ser + Fe_liv_from_ferritin - Fe_liv_to_ferritin;
    der(Fe_liv_3) = Fe_liv_to_ferritin - Fe_liv_from_ferritin;
    der(Fpn_liv) = Fpn_liv_in - Fpn_liv_out; // eq 7.1
    der(Fpn_liv_mRNA) = Fpn_liv_mRNA_in - Fpn_liv_mRNA_out; // eq 6.1

    // Spleen

    Fe_spl = Fe_spl_2 + Fe_spl_3;
    der(Fe_spl_2) = Fe_spl_in_RBC + Fe_spl_in_bm - Fe_spl_out_ser + Fe_spl_from_ferritin - Fe_spl_to_ferritin;
    der(Fe_spl_3) = Fe_spl_to_ferritin - Fe_spl_from_ferritin;
    der(Fpn_spl) = Fpn_spl_in - Fpn_spl_out; // eq 7.2
    der(Fpn_spl_mRNA) = Fpn_spl_mRNA_in - Fpn_spl_mRNA_out; // eq 6.2

    // Other organs

    Fe_res = Fe_res_2 + Fe_res_3;
    der(Fe_res_2) = Fe_res_in_ser - Fe_res_out_ser - Fe_res_out_loss + Fe_res_from_ferritin - Fe_res_to_ferritin;
    der(Fe_res_3) = Fe_res_to_ferritin - Fe_res_from_ferritin;
    der(Fpn_res) = Fpn_res_in - Fpn_res_out; // eq 7.4
    der(Fpn_res_mRNA) = Fpn_res_mRNA_in - Fpn_res_mRNA_out; // eq 6.4

    // Bone marrow

    der(Fe_bm) = Fe_bm_in_ser - Fe_bm_out_RBC - Fe_bm_out_spl; //eq 12.

    // Red blood cells

    der(Fe_RBC) = Fe_RBC_in_bm - Fe_RBC_out_spl - bleeding + transfusion; //eq 13.

    //+Relative indicators

    // Hepcidin expression
    //################################################
    hep_rel = hep / hep_norm;
    Bmp6_rel = Bmp6 / Bmp6_norm;

    // Ferroportin regulation
    //################################################
    Fpn_liv_mRNA_rel = Fpn_liv_mRNA / Fpn_liv_mRNA_norm;
    Fpn_spl_mRNA_rel = Fpn_spl_mRNA / Fpn_spl_mRNA_norm;
    Fpn_duo_mRNA_rel = Fpn_duo_mRNA / Fpn_duo_mRNA_norm;
    Fpn_res_mRNA_rel = Fpn_res_mRNA / Fpn_res_mRNA_norm;
    Fpn_liv_rel = Fpn_liv / Fpn_liv_norm;
    Fpn_spl_rel = Fpn_spl / Fpn_spl_norm;
    Fpn_duo_rel = Fpn_duo / Fpn_duo_norm;
    Fpn_res_rel = Fpn_res / Fpn_res_norm;

    // Dynamics of the iron pools
    //################################################
    Fe_liv_rel = Fe_liv / Fe_liv_norm;
    Fe_spl_rel = Fe_spl / Fe_spl_norm;
    Fe_bm_rel = Fe_bm / Fe_bm_norm;
    Fe_RBC_rel = Fe_RBC / Fe_RBC_norm;
    Fe_duo_rel = Fe_duo / Fe_duo_norm;
    Fe_res_rel = Fe_res / Fe_res_norm;
    Fe_ser_rel = Fe_ser / Fe_ser_norm;
    Fe_total_rel = Fe_total / Fe_total_norm;
    Fe_ser_input_rel = Fe_ser_input / Fe_ser_input_norm;
    Fe_ser_output_rel = Fe_ser_output / Fe_ser_output_norm;

    // ###############################################
    // Auxiliary relative indicators
    // ###############################################

    // Duodenum

    Fe_duo_2_rel = Fe_duo_2 / Fe_duo_2_norm;
    Fe_duo_3_rel = Fe_duo_3 / Fe_duo_3_norm;
    Fe_duo_to_ferritin_rel = Fe_duo_to_ferritin / Fe_duo_to_ferritin_norm;
    Fe_duo_from_ferritin_rel = Fe_duo_from_ferritin / Fe_duo_from_ferritin_norm;
    Fe_duo_intake_rel = Fe_duo_intake / Fe_duo_intake_norm;
    Fe_duo_in_food_rel = Fe_duo_in_food / Fe_duo_in_food_norm;
    Fe_duo_in_ser_rel = Fe_duo_in_ser / Fe_duo_in_ser_norm;
    Fe_duo_out_ser_rel = Fe_duo_out_ser / Fe_duo_out_ser_norm;
    Fe_duo_out_loss_rel = Fe_duo_out_loss / Fe_duo_out_loss_norm;
    Fe_duo_unused_rel = Fe_duo_unused / Fe_duo_unused_norm;
    Fpn_duo_in_1_rel = Fpn_duo_in_1 / Fpn_duo_in_1_norm;
    Fpn_duo_in_2_rel = Fpn_duo_in_2 / Fpn_duo_in_2_norm;
    Fpn_duo_out_1_rel = Fpn_duo_out_1 / Fpn_duo_out_1_norm;
    Fpn_duo_out_2_rel = Fpn_duo_out_2 / Fpn_duo_out_2_norm;
    Fpn_duo_in_rel = Fpn_duo_in / Fpn_duo_in_norm;
    Fpn_duo_out_rel = Fpn_duo_out / Fpn_duo_out_norm;
    Fpn_duo_mRNA_in_rel = Fpn_duo_mRNA_in / Fpn_duo_mRNA_in_norm;
    Fpn_duo_mRNA_out_rel = Fpn_duo_mRNA_out / Fpn_duo_mRNA_out_norm;

    // Liver

    Fe_liv_2_rel = Fe_liv_2 / Fe_liv_2_norm;
    Fe_liv_3_rel = Fe_liv_3 / Fe_liv_3_norm;
    Fe_liv_to_ferritin_rel = Fe_liv_to_ferritin / Fe_liv_to_ferritin_norm;
    Fe_liv_from_ferritin_rel = Fe_liv_from_ferritin / Fe_liv_from_ferritin_norm;
    Fe_liv_in_ser_rel = Fe_liv_in_ser / Fe_liv_in_ser_norm;
    Fe_liv_out_ser_rel = Fe_liv_out_ser / Fe_liv_out_ser_norm;
    Fpn_liv_in_1_rel = Fpn_liv_in_1 / Fpn_liv_in_1_norm;
    Fpn_liv_in_2_rel = Fpn_liv_in_2 / Fpn_liv_in_2_norm;
    Fpn_liv_out_1_rel = Fpn_liv_out_1 / Fpn_liv_out_1_norm;
    Fpn_liv_out_2_rel = Fpn_liv_out_2 / Fpn_liv_out_2_norm;
    Fpn_liv_in_rel = Fpn_liv_in / Fpn_liv_in_norm;
    Fpn_liv_out_rel = Fpn_liv_out / Fpn_liv_out_norm;
    Fpn_liv_mRNA_in_rel = Fpn_liv_mRNA_in / Fpn_liv_mRNA_in_norm;
    Fpn_liv_mRNA_out_rel = Fpn_liv_mRNA_out / Fpn_liv_mRNA_out_norm;

    // Serum

    Fe_ser_in_liv_rel = Fe_ser_in_liv / Fe_ser_in_liv_norm;
    Fe_ser_in_spl_rel = Fe_ser_in_spl / Fe_ser_in_spl_norm;
    Fe_ser_in_duo_rel = Fe_ser_in_duo / Fe_ser_in_duo_norm;
    Fe_ser_in_res_rel = Fe_ser_in_res / Fe_ser_in_res_norm;
    Fe_ser_out_liv_rel = Fe_ser_out_liv / Fe_ser_out_liv_norm;
    Fe_ser_out_bm_rel = Fe_ser_out_bm / Fe_ser_out_bm_norm;
    Fe_ser_out_duo_rel = Fe_ser_out_duo / Fe_ser_out_duo_norm;
    Fe_ser_out_res_rel = Fe_ser_out_res / Fe_ser_out_res_norm;

    // Hepcidine related stuff

    hep_in_rel = hep_in / hep_in_norm;
    hep_out_rel = hep_out / hep_out_norm;
    Bmp6_in_rel = Bmp6_in / Bmp6_in_norm;
    Bmp6_out_rel = Bmp6_out / Bmp6_out_norm;

    // Spleen

    Fe_spl_2_rel = Fe_spl_2 / Fe_spl_2_norm;
    Fe_spl_3_rel = Fe_spl_3 / Fe_spl_3_norm;
    Fe_spl_to_ferritin_rel = Fe_spl_to_ferritin / Fe_spl_to_ferritin_norm;
    Fe_spl_from_ferritin_rel = Fe_spl_from_ferritin / Fe_spl_from_ferritin_norm;
    Fe_spl_in_RBC_rel = Fe_spl_in_RBC / Fe_spl_in_RBC_norm;
    Fe_spl_in_bm_rel = Fe_spl_in_bm / Fe_spl_in_bm_norm;
    Fe_spl_out_ser_rel = Fe_spl_out_ser / Fe_spl_out_ser_norm;
    Fpn_spl_in_1_rel = Fpn_spl_in_1 / Fpn_spl_in_1_norm;
    Fpn_spl_in_2_rel = Fpn_spl_in_2 / Fpn_spl_in_2_norm;
    Fpn_spl_out_1_rel = Fpn_spl_out_1 / Fpn_spl_out_1_norm;
    Fpn_spl_out_2_rel = Fpn_spl_out_2 / Fpn_spl_out_2_norm;
    Fpn_spl_in_rel = Fpn_spl_in / Fpn_spl_in_norm;
    Fpn_spl_out_rel = Fpn_spl_out / Fpn_spl_out_norm;
    Fpn_spl_mRNA_in_rel = Fpn_spl_mRNA_in / Fpn_spl_mRNA_in_norm;
    Fpn_spl_mRNA_out_rel = Fpn_spl_mRNA_out / Fpn_spl_mRNA_out_norm;

    // Other organs

    Fe_res_2_rel = Fe_res_2 / Fe_res_2_norm;
    Fe_res_3_rel = Fe_res_3 / Fe_res_3_norm;
    Fe_res_to_ferritin_rel = Fe_res_to_ferritin / Fe_res_to_ferritin_norm;
    Fe_res_from_ferritin_rel = Fe_res_from_ferritin / Fe_res_from_ferritin_norm;
    Fe_res_in_ser_rel = Fe_res_in_ser / Fe_res_in_ser_norm;
    Fe_res_out_ser_rel = Fe_res_out_ser / Fe_res_out_ser_norm;
    Fe_res_out_loss_rel = Fe_res_out_loss / Fe_res_out_loss_norm;
    Fpn_res_in_1_rel = Fpn_res_in_1 / Fpn_res_in_1_norm;
    Fpn_res_in_2_rel = Fpn_res_in_2 / Fpn_res_in_2_norm;
    Fpn_res_out_1_rel = Fpn_res_out_1 / Fpn_res_out_1_norm;
    Fpn_res_out_2_rel = Fpn_res_out_2 / Fpn_res_out_2_norm;
    Fpn_res_in_rel = Fpn_res_in / Fpn_res_in_norm;
    Fpn_res_out_rel = Fpn_res_out / Fpn_res_out_norm;
    Fpn_res_mRNA_in_rel = Fpn_res_mRNA_in / Fpn_res_mRNA_in_norm;
    Fpn_res_mRNA_out_rel = Fpn_res_mRNA_out / Fpn_res_mRNA_out_norm;

    // Bone marrow

    Fe_bm_in_ser_rel = Fe_bm_in_ser / Fe_bm_in_ser_norm;
    Fe_bm_out_RBC_rel = Fe_bm_out_RBC / Fe_bm_out_RBC_norm;
    Fe_bm_out_spl_rel = Fe_bm_out_spl / Fe_bm_out_spl_norm;

    //Red blood cells (RBC)

    Fe_RBC_in_bm_rel = Fe_RBC_in_bm / Fe_RBC_in_bm_norm;
    Fe_RBC_out_spl_rel = Fe_RBC_out_spl / Fe_RBC_out_spl_norm;

    //-Relative indicators

    annotation (
      Icon(coordinateSystem(preserveAspectRatio=false)),
      Diagram(coordinateSystem(preserveAspectRatio=false)),
      experiment(
        StopTime=1000000,
        __Dymola_NumberOfIntervals=1000000,
        Tolerance=1e-12,
        __Dymola_Algorithm="Dassl"));
  end FeMetabolismModel;

  function facB1 "facB1 function (id f9)"
    input Real facB1_i1;
    input Real facB1_i2;
    output Real facB1_o1 "result";
  algorithm
    facB1_o1 := ((pSMAD(pSMAD_i1=facB1_i1, pSMAD_i2=facB1_i2))^1.7807)/0.4391;
  end facB1;

  function facB2 "facB2 function (id f10)"
    input Real facB2_i1;
    input Real facB2_i2;
    output Real facB2_o1 "result";
  algorithm
    facB2_o1 := ((pSMAD(pSMAD_i1=facB2_i1, pSMAD_i2=facB2_i2))^1.7807)/16.8738;
  end facB2;

  function facSig "facSig function (id f5)"
    input Real facSig_i1;
    input Real facSig_i2;
    output Real facSig_o1 "result";
  algorithm
    facSig_o1 := (0.1285*2.852*hill(
        hill_i1=facSig_i1,
        hill_i2=1.0242,
        hill_i3=7.7388) - (1 + 0.4135*(0.0583 + 1.949*hill(
        hill_i1=facSig_i2,
        hill_i2=1.4481,
        hill_i3=140.244))))/(1 + 0.4135*0.0583);
  end facSig;

  function facSig1 "facSig1 function (id f6)"
    input Real facSig1_i1;
    output Real facSig1_o1 "result";
  algorithm
    facSig1_o1 := 0.1285*2.852*hill(
        hill_i1=facSig1_i1,
        hill_i2=1.0242,
        hill_i3=7.7388)/(1 + 0.4135*0.0583);
  end facSig1;

  function facST "facST function (id f11)"
    input Real facST_i1;
    input Real facST_i2;
    output Real facST_o1 "result";
  algorithm
    facST_o1 := pSTAT(pSTAT_i1=facST_i1, pSTAT_i2=facST_i2)/206.3988;
  end facST;

  function Freg "Freg function (id f14)"
    input Real Freg_i1;
    input Real Freg_i2;
    output Real Freg_o1 "result";
  algorithm
    Freg_o1 := Freg1(Freg1_i1=Freg_i1, Freg1_i2=Freg_i2)/Freg2(Freg2_i1=Freg_i1,
      Freg2_i2=Freg_i2);
  end Freg;

  function Freg1 "Freg1 function (id f12)"
    input Real Freg1_i1;
    input Real Freg1_i2;
    output Real Freg1_o1 "result";
  algorithm
    Freg1_o1 := 1 + FeMetabolism.facB1(facB1_i1=Freg1_i1, facB1_i2=Freg1_i2)*
      537.649 + FeMetabolism.facB2(facB2_i1=Freg1_i1, facB2_i2=Freg1_i2)*4972.6
       + FeMetabolism.facST(facST_i1=Freg1_i1, facST_i2=Freg1_i2)*584.75 +
      FeMetabolism.facB1(facB1_i1=Freg1_i1, facB1_i2=Freg1_i2)*
      FeMetabolism.facB2(facB2_i1=Freg1_i1, facB2_i2=Freg1_i2)*537.649*4972.6
       + FeMetabolism.facB1(facB1_i1=Freg1_i1, facB1_i2=Freg1_i2)*
      FeMetabolism.facST(facST_i1=Freg1_i1, facST_i2=Freg1_i2)*537.649*584.75*
      5.3869 + FeMetabolism.facB2(facB2_i1=Freg1_i1, facB2_i2=Freg1_i2)*
      FeMetabolism.facST(facST_i1=Freg1_i1, facST_i2=Freg1_i2)*4972.6*584.75 +
      FeMetabolism.facB1(facB1_i1=Freg1_i1, facB1_i2=Freg1_i2)*
      FeMetabolism.facB2(facB2_i1=Freg1_i1, facB2_i2=Freg1_i2)*
      FeMetabolism.facST(facST_i1=Freg1_i1, facST_i2=Freg1_i2)*537.649*4972.6*
      584.75*5.3869;
  end Freg1;

  function Freg2 "Freg2 function (id f13)"
    input Real Freg2_i1;
    input Real Freg2_i2;
    output Real Freg2_o1 "result";
  algorithm
    Freg2_o1 := 1 + FeMetabolism.facB1(facB1_i1=Freg2_i1, facB1_i2=Freg2_i2) +
      FeMetabolism.facB2(facB2_i1=Freg2_i1, facB2_i2=Freg2_i2) +
      FeMetabolism.facST(facST_i1=Freg2_i1, facST_i2=Freg2_i2) +
      FeMetabolism.facB1(facB1_i1=Freg2_i1, facB1_i2=Freg2_i2)*
      FeMetabolism.facB2(facB2_i1=Freg2_i1, facB2_i2=Freg2_i2) +
      FeMetabolism.facB1(facB1_i1=Freg2_i1, facB1_i2=Freg2_i2)*
      FeMetabolism.facST(facST_i1=Freg2_i1, facST_i2=Freg2_i2)*5.3869 +
      FeMetabolism.facB2(facB2_i1=Freg2_i1, facB2_i2=Freg2_i2)*
      FeMetabolism.facST(facST_i1=Freg2_i1, facST_i2=Freg2_i2) +
      FeMetabolism.facB1(facB1_i1=Freg2_i1, facB1_i2=Freg2_i2)*
      FeMetabolism.facB2(facB2_i1=Freg2_i1, facB2_i2=Freg2_i2)*
      FeMetabolism.facST(facST_i1=Freg2_i1, facST_i2=Freg2_i2)*5.3869;
  end Freg2;

  function hill "Hill function (id f4)"
    input Real hill_i1;
    input Real hill_i2;
    input Real hill_i3;
    output Real hill_o1 "result";
    constant Real eps=1e-9;
    Real i1_lim;
    Real i3_lim;
  algorithm
    if hill_i1 < eps then
      i1_lim := eps;
    else
      i1_lim := hill_i1;
    end if;

    if hill_i3 < eps then
      i3_lim := eps;
    else
      i3_lim := hill_i3;
    end if;

    hill_o1 := (i1_lim^hill_i2)/(i1_lim^hill_i2 + i3_lim^hill_i2);
  end hill;

  function Promoter "Promoter occupancy function (id f15)"
    input Real Promoter_i1;
    input Real Promoter_i2;
    output Real Promoter_o1 "result";
  algorithm
    Promoter_o1 := FeMetabolism.Freg(Freg_i1=Promoter_i1, Freg_i2=Promoter_i2)/(
      FeMetabolism.Freg(Freg_i1=Promoter_i1, Freg_i2=Promoter_i2) + 6804.7);
  end Promoter;

  function pSMAD "pSMAD function (id f8)"
    input Real pSMAD_i1;
    input Real pSMAD_i2;
    output Real pSMAD_o1 "result";
  algorithm
    pSMAD_o1 := 0.0583 + 1.949*FeMetabolism.hill(
        hill_i1=pSMAD_i2,
        hill_i2=1.4481,
        hill_i3=140.244)/(1 + 0.1285*pSTAT(pSTAT_i1=pSMAD_i1, pSTAT_i2=pSMAD_i2));
  end pSMAD;

  function pSTAT "pSTAT function (id f7)"
    input Real pSTAT_i1;
    input Real pSTAT_i2;
    output Real pSTAT_o1 "result";
  algorithm
    pSTAT_o1 := (0.5/0.1285)*(FeMetabolism.facSig(facSig_i1=pSTAT_i1, facSig_i2=
      pSTAT_i2) + sqrt((FeMetabolism.facSig(facSig_i1=pSTAT_i1, facSig_i2=
      pSTAT_i2))^2 + 4*FeMetabolism.facSig1(facSig1_i1=pSTAT_i1)));
  end pSTAT;

  function MIN "MINimum function (id f1)"
    input Real MIN_i1;
    input Real MIN_i2;
    output Real MIN_o1 "result";
  algorithm
    MIN_o1 := (MIN_i1 + MIN_i2 - abs(MIN_i1 - MIN_i2))/2;
  end MIN;

  function MAX "MAXimum function (id f2)"
    input Real MAX_i1;
    input Real MAX_i2;
    output Real MAX_o1 "result";
  algorithm
    MAX_o1 := (MAX_i1 + MAX_i2 + abs(MAX_i1 - MAX_i2))/2;
  end MAX;

  model Serum
    Bodylight.Types.RealIO.MassFlowRateInput Fe_ser_in_liv_SI annotation (
        Placement(transformation(extent={{-188,-10},{-148,30}}),
          iconTransformation(extent={{-114,76},{-92,98}})));
    Bodylight.Types.RealIO.MassFlowRateInput Fe_ser_in_spl_SI annotation (
        Placement(transformation(extent={{-188,-10},{-148,30}}),
          iconTransformation(extent={{-114,54},{-92,76}})));
    Bodylight.Types.RealIO.MassFlowRateInput Fe_ser_in_duo_SI annotation (
        Placement(transformation(extent={{-188,-10},{-148,30}}),
          iconTransformation(extent={{-114,32},{-92,54}})));
    Bodylight.Types.RealIO.MassFlowRateInput Fe_ser_in_res_SI annotation (
        Placement(transformation(extent={{-188,-10},{-148,30}}),
          iconTransformation(extent={{-114,10},{-92,32}})));
    Bodylight.Types.RealIO.MassFlowRateInput transfusion_SI annotation (Placement(
          transformation(extent={{-188,-10},{-148,30}}), iconTransformation(
          extent={{-11,-11},{11,11}},
          rotation=90,
          origin={-45,-111})));
    Bodylight.Types.RealIO.MassFlowRateInput Fe_ser_out_liv_SI annotation (
        Placement(transformation(extent={{-188,-10},{-148,30}}),
          iconTransformation(extent={{-114,-12},{-92,10}})));
    Bodylight.Types.RealIO.MassFlowRateInput Fe_ser_out_duo_SI annotation (
        Placement(transformation(extent={{-188,-10},{-148,30}}),
          iconTransformation(extent={{-114,-34},{-92,-12}})));
    Bodylight.Types.RealIO.MassFlowRateInput Fe_ser_out_res_SI annotation (
        Placement(transformation(extent={{-188,-10},{-148,30}}),
          iconTransformation(extent={{-114,-56},{-92,-34}})));
    Bodylight.Types.RealIO.MassFlowRateInput bleeding_SI annotation (Placement(
          transformation(extent={{-188,-10},{-148,30}}), iconTransformation(
          extent={{-11,-11},{11,11}},
          rotation=90,
          origin={19,-111})));
    Bodylight.Types.RealIO.MassFlowRateInput Fe_ser_out_bm_SI annotation (
        Placement(transformation(extent={{-188,-10},{-148,30}}),
        iconTransformation(extent={{-114,-78},{-92,-56}})));
    //Fe_ser(start=1.51304,unit="ug") "Fe amount in serum (s1)";
    //microgram to kilogram conversion: 1.51304 ug=1.51304e-12 kg;_
    Bodylight.Types.RealIO.MassOutput Fe_ser_SI( start=1.51304e-12) annotation (
        Placement(transformation(extent={{-178,36},{-158,56}}),
          iconTransformation(
          extent={{-10,-10},{10,10}},
          rotation=0,
          origin={110,84})));
    Bodylight.Types.MassFlowRate Fe_ser_input;
    Bodylight.Types.MassFlowRate Fe_ser_output;
  equation
    Fe_ser_input = Fe_ser_in_liv_SI + Fe_ser_in_spl_SI + Fe_ser_in_duo_SI + Fe_ser_in_res_SI + transfusion_SI;
    Fe_ser_output = Fe_ser_out_liv_SI + Fe_ser_out_bm_SI + Fe_ser_out_duo_SI + Fe_ser_out_res_SI + bleeding_SI;
    der(Fe_ser_SI) = Fe_ser_input - Fe_ser_output;
    annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
            Rectangle(
            extent={{100,-100},{-100,100}},
            lineColor={28,108,200},
            fillColor={255,255,0},
            fillPattern=FillPattern.Solid), Text(
            extent={{-50,42},{46,-10}},
            textColor={28,108,200},
            textString="Serum")}), Diagram(coordinateSystem(preserveAspectRatio
            =false)));

                                                   // eq 16.

  end Serum;

  model Other_Organ
    Bodylight.Types.RealIO.MassFlowRateOutput Fe_ser_in_res_SI annotation (
        Placement(transformation(extent={{-170,-14},{-150,6}}),
          iconTransformation(extent={{100,80},{120,100}})));
    Bodylight.Types.RealIO.MassFlowRateOutput Fe_ser_out_res_SI annotation (
        Placement(transformation(extent={{-170,-14},{-150,6}}),
          iconTransformation(extent={{100,60},{120,80}})));
    annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
            Rectangle(
            extent={{100,-100},{-100,100}},
            lineColor={28,108,200},
            fillColor={255,255,0},
            fillPattern=FillPattern.Solid), Text(
            extent={{-56,48},{46,-10}},
            textColor={28,108,200},
            textString="Other_Organ")}), Diagram(coordinateSystem(
            preserveAspectRatio=false)));
  end Other_Organ;

  model Duodenum
    Bodylight.Types.RealIO.MassFlowRateOutput Fe_ser_in_duo_SI annotation (
        Placement(transformation(extent={{-178,58},{-158,78}}),
          iconTransformation(extent={{100,74},{120,94}})));
    Bodylight.Types.RealIO.MassFlowRateOutput Fe_ser_out_duo_SI annotation (
        Placement(transformation(extent={{-178,58},{-158,78}}),
          iconTransformation(extent={{100,54},{120,74}})));
    annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
            Rectangle(
            extent={{100,-100},{-100,100}},
            lineColor={28,108,200},
            fillColor={255,255,0},
            fillPattern=FillPattern.Solid), Text(
            extent={{-56,48},{46,-10}},
            textColor={28,108,200},
            textString="Duodenum")}), Diagram(coordinateSystem(
            preserveAspectRatio=false)));
  end Duodenum;

  model Bone_Marrow
    Bodylight.Types.RealIO.MassFlowRateOutput Fe_ser_out_bm_SI annotation (
        Placement(transformation(extent={{-176,-44},{-156,-24}}),
          iconTransformation(extent={{100,76},{120,96}})));
    Bodylight.Types.RealIO.MassFlowRateOutput Fe_bm_out_RBC_SI annotation (
        Placement(transformation(extent={{-176,-44},{-156,-24}}),
          iconTransformation(extent={{100,56},{120,76}})));
    Bodylight.Types.RealIO.MassFlowRateOutput Fe_bm_out_spl_SI annotation (
        Placement(transformation(extent={{-176,-44},{-156,-24}}),
          iconTransformation(extent={{100,36},{120,56}})));
    annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
            Rectangle(
            extent={{100,-100},{-100,100}},
            lineColor={28,108,200},
            fillColor={255,255,0},
            fillPattern=FillPattern.Solid), Text(
            extent={{-56,48},{46,-10}},
            textColor={28,108,200},
            textString="Bone_Marrow")}), Diagram(coordinateSystem(
            preserveAspectRatio=false)));
  end Bone_Marrow;

  model Red_Blood_Cells
    Bodylight.Types.RealIO.MassFlowRateInput Fe_RBC_in_bm_SI annotation (
        Placement(transformation(extent={{-156,-36},{-116,4}}),
          iconTransformation(extent={{-130,48},{-96,82}})));
    Bodylight.Types.RealIO.MassFlowRateInput bleeding_SI annotation (Placement(
          transformation(extent={{-156,8},{-116,48}}),   iconTransformation(
            extent={{-132,-16},{-98,18}})));
    Bodylight.Types.RealIO.MassFlowRateInput transfusion_SI annotation (Placement(
          transformation(extent={{-160,44},{-120,84}}),  iconTransformation(
            extent={{-132,-74},{-98,-40}})));
    Bodylight.Types.RealIO.MassOutput Fe_RBC_SI( start=1016.720) annotation (Placement(
          transformation(extent={{104,22},{124,42}}),     iconTransformation(
          extent={{100,72},{120,92}})));
    //Real Fe_RBC(start=1016.720,unit="ug") "Fe amount in red blood cells (s6)";
    Bodylight.Types.RealIO.MassFlowRateOutput Fe_RBC_out_spl_SI annotation (
        Placement(transformation(extent={{104,-10},{124,10}}),
          iconTransformation(extent={{100,50},{120,70}})));
    parameter Real v_spl_1(unit="h-1") = 0.0036035
       "Spleen iron uptake rate from RBC (k10), 
      0.002 = value for trace exp., 0.004/0.005";
  equation
    der(Fe_RBC_SI) = Fe_RBC_in_bm_SI - Fe_RBC_out_spl_SI - bleeding_SI + transfusion_SI;
    annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
            Rectangle(
            extent={{100,-100},{-100,100}},
            lineColor={28,108,200},
            fillColor={255,255,0},
            fillPattern=FillPattern.Solid), Text(
            extent={{-56,48},{46,-10}},
            textColor={28,108,200},
            textString="Red_Blood_Cells")}), Diagram(coordinateSystem(
            preserveAspectRatio=false)));

    // Real Fe_RBC_in_bm = v_RBC*Fe_bm "Fe flow to RBC from bone marrow";
    // Real Fe_bm_out_RBC = v_RBC*Fe_bm "Fe flow from bone marrow to RBC";
    // Fe_bm_out_RBC = Real Fe_RBC_in_bm calculated in bone marrow



   // Fe_RBC_out_spl_SI = v_spl_1/3600*Fe_RBC_SI "Fe flow from RBC to spleen";
                                                                                         //eq 13.

  end Red_Blood_Cells;

  model Spleen
    Bodylight.Types.RealIO.MassFlowRateOutput Fe_ser_in_spl_SI annotation (
        Placement(transformation(extent={{-180,26},{-160,46}}),
          iconTransformation(extent={{100,74},{120,94}})));
    annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
            Rectangle(
            extent={{100,-100},{-100,100}},
            lineColor={28,108,200},
            fillColor={255,255,0},
            fillPattern=FillPattern.Solid), Text(
            extent={{-56,48},{46,-10}},
            textColor={28,108,200},
            textString="Spleen")}), Diagram(coordinateSystem(
            preserveAspectRatio=false)));
  end Spleen;

  model Liver
    Bodylight.Types.RealIO.MassFlowRateOutput Fe_ser_in_liv_SI annotation (
        Placement(transformation(extent={{-180,26},{-160,46}}),
          iconTransformation(extent={{100,78},{120,98}})));
    Bodylight.Types.RealIO.MassFlowRateOutput Fe_ser_out_liv annotation (
        Placement(transformation(extent={{-172,8},{-152,28}}),
          iconTransformation(extent={{100,58},{120,78}})));
    annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
            Rectangle(
            extent={{100,-100},{-100,100}},
            lineColor={28,108,200},
            fillColor={255,255,0},
            fillPattern=FillPattern.Solid), Text(
            extent={{-56,48},{46,-10}},
            textColor={28,108,200},
            textString="Liver")}),  Diagram(coordinateSystem(
            preserveAspectRatio=false)));
  end Liver;

  model Iron_Metabolism
    Liver liver
      annotation (Placement(transformation(extent={{-76,60},{-56,80}})));
    Duodenum duodenum
      annotation (Placement(transformation(extent={{-76,-8},{-56,12}})));
    Bone_Marrow bone_Marrow
      annotation (Placement(transformation(extent={{20,-36},{40,-16}})));
    Red_Blood_Cells red_Blood_Cells
      annotation (Placement(transformation(extent={{62,-10},{82,10}})));
    Spleen spleen
      annotation (Placement(transformation(extent={{38,18},{58,38}})));
    Serum serum annotation (Placement(transformation(extent={{-12,22},{8,42}})));
    Other_Organ other_Organ
      annotation (Placement(transformation(extent={{-52,-40},{-32,-20}})));
    Bodylight.Types.Constants.MassFlowRateConst Transfusion
      annotation (Placement(transformation(extent={{-18,8},{-10,16}})));
    Bodylight.Types.Constants.MassFlowRateConst Bleeding
      annotation (Placement(transformation(extent={{16,10},{8,18}})));
  equation
    connect(liver.Fe_ser_in_liv_SI, serum.Fe_ser_in_liv_SI) annotation (Line(
          points={{-55,78.8},{-38,78.8},{-38,40.7},{-12.3,40.7}}, color={0,0,
            127}));
    connect(serum.Fe_ser_out_liv_SI, liver.Fe_ser_out_liv) annotation (Line(
          points={{-12.3,31.9},{-42,31.9},{-42,76.8},{-55,76.8}}, color={0,0,
            127}));
    connect(spleen.Fe_ser_in_spl_SI, serum.Fe_ser_in_spl_SI) annotation (Line(
          points={{59,36.4},{62,36.4},{62,60},{-34,60},{-34,38.5},{-12.3,38.5}},
          color={0,0,127}));
    connect(serum.Fe_ser_in_duo_SI, duodenum.Fe_ser_in_duo_SI) annotation (Line(
          points={{-12.3,36.3},{-34,36.3},{-34,10.4},{-55,10.4}}, color={0,0,
            127}));
    connect(serum.Fe_ser_out_duo_SI, duodenum.Fe_ser_out_duo_SI) annotation (
        Line(points={{-12.3,29.7},{-32,29.7},{-32,8.4},{-55,8.4}}, color={0,0,
            127}));
    connect(serum.Fe_ser_in_res_SI, other_Organ.Fe_ser_in_res_SI) annotation (
        Line(points={{-12.3,34.1},{-28,34.1},{-28,-21},{-31,-21}}, color={0,0,
            127}));
    connect(other_Organ.Fe_ser_out_res_SI, serum.Fe_ser_out_res_SI) annotation
      (Line(points={{-31,-23},{-24,-23},{-24,27.5},{-12.3,27.5}}, color={0,0,
            127}));
    connect(bone_Marrow.Fe_ser_out_bm_SI, serum.Fe_ser_out_bm_SI) annotation (
        Line(points={{41,-17.4},{46,-17.4},{46,-10},{-20,-10},{-20,25.3},{-12.3,
            25.3}}, color={0,0,127}));
    connect(Transfusion.y, serum.transfusion_SI) annotation (Line(points={{-9,
            12},{-6.5,12},{-6.5,20.9}}, color={0,0,127}));
    connect(Bleeding.y, serum.bleeding_SI) annotation (Line(points={{7,14},{
            -0.1,14},{-0.1,20.9}}, color={0,0,127}));
    connect(red_Blood_Cells.bleeding_SI, Bleeding.y) annotation (Line(points={{
            60.5,-1.1},{0,-1.1},{0,14},{7,14}}, color={0,0,127}));
    connect(red_Blood_Cells.transfusion_SI, Transfusion.y) annotation (Line(
          points={{60.5,-5.7},{-6,-5.7},{-6,12},{-9,12}}, color={0,0,127}));
    connect(bone_Marrow.Fe_bm_out_RBC_SI, red_Blood_Cells.Fe_RBC_in_bm_SI)
      annotation (Line(points={{41,-19.4},{52,-19.4},{52,3.9},{60.5,3.9}},
          color={0,0,127}));
    annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
            Rectangle(
            extent={{100,-100},{-100,100}},
            lineColor={28,108,200},
            fillColor={255,255,0},
            fillPattern=FillPattern.Solid), Text(
            extent={{-56,48},{46,-10}},
            textColor={28,108,200},
            textString="Iron Metabolism")}), Diagram(coordinateSystem(
            preserveAspectRatio=false)));
  end Iron_Metabolism;
  annotation (uses(Bodylight(version="1.0")));
end FeMetabolism;
