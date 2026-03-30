within ;
package EnterocyteMucosalBlock "Enterocyte mucosal block"

  package models
    model CellularFerritinIronStorageModel "Cellular ferritin iron storage"

      //Species
      Bodylight.Types.Concentration LIP(
        displayUnit = "mol/L",
        start = 1e-05 * 1e3) "labile iron pool";
      Bodylight.Types.Concentration FT_cage(
        displayUnit = "mol/L",
        start = 5e-09 * 1e3) "FT-cage";
      Bodylight.Types.Concentration core(
        displayUnit = "mol/L",
        start = 7.5e-06 * 1e3) "core";
      Bodylight.Types.Concentration DFP(
        displayUnit = "mol/L",
        start = 0) "diferric peroxo complex";

      Real atoms_per_cage_transient "Transient number of Fe atoms that are stored inside the core of a ferritin cage";

      parameter Integer H = 4 "H subunits";
      parameter Integer L = 24 - H "L subunits";

      parameter Integer rN = 50;
      parameter Integer rO = 2;

      parameter BodylightExtension.Types.ReactionRateFirstOrder k_FTlysis = 1.203e-05;

      parameter BodylightExtension.Types.MolarReactionRate FT_Expression = 6.015e-14 * 1e3;

      //FT degradation
      BodylightExtension.Types.MolarReactionRate FT_Degradation;

      //FT degradation core release
      BodylightExtension.Types.MolarReactionRate CoreRelease;

      //Oxidation (2 LIP -> DFP)
      BodylightExtension.Types.MolarReactionRate Oxidation;
      parameter BodylightExtension.Types.ReactionRateFirstOrder k_cat_oxidation = 591 "catalytic turnover number";
      parameter Bodylight.Types.Concentration K_m_oxidation(
        displayUnit = "mol/L") = 0.35 "Michaelis constant";
      parameter Real n_oxidation = 1.3 "Hill coefficient)";

      //Reduction (DFP -> 2 LIP)
      BodylightExtension.Types.MolarReactionRate Reduction;
      parameter BodylightExtension.Types.ReactionRateFirstOrder k_deg = 0.2605 "rate constant";

      //Nucleation (2 DFP -> 4 core)
      BodylightExtension.Types.MolarReactionRate Nucleation;
      parameter BodylightExtension.Types.ReactionRateThirdOrder k_cat_nucleation = 5e7 * 1e-6 "catalytic turnover number";
      parameter Bodylight.Types.Concentration K_i_nucleation(
        displayUnit = "mol/L") = 0.461598e-3 * 1e3 "inhibition constant";
      parameter Integer n_nucleation = 4 "Hill coefficient";

      //Mineralization (DFP -> 2 core)
      BodylightExtension.Types.MolarReactionRate Mineralization;
      parameter BodylightExtension.Types.ReactionRateFirstOrder k_cat_mineralization = 0.101564 "catalytic turnover number";
      parameter Bodylight.Types.Concentration K_m_mineralization(
        displayUnit = "mol/L") = 5e-06 * 1e3 "Michaelis constant";
      parameter Bodylight.Types.Concentration K_i_mineralization(
        displayUnit = "mol/L") = 4.6458e-3 * 1e3 "inhibition constant";
      parameter Integer n_mineralization = 4 "Hill coefficient";
      parameter Integer m_mineralization = 8 "Hill coefficient";

    equation

      atoms_per_cage_transient = core / FT_cage;

      FT_Degradation = k_FTlysis * FT_cage;

      CoreRelease = k_FTlysis * core;

      Oxidation = (k_cat_oxidation * (H + rO) / (24 + rO) * FT_cage * LIP ^ n_oxidation)
        / (K_m_oxidation ^ n_oxidation + LIP ^ n_oxidation);

      Reduction = k_deg * DFP;

      Nucleation = k_cat_nucleation * DFP ^ 2 * FT_cage * (L + rN) / (24 + rN)
        * K_i_nucleation ^ n_nucleation / (K_i_nucleation ^ n_nucleation + core ^ n_nucleation);

      Mineralization = (k_cat_mineralization * DFP * core) / (K_m_mineralization + DFP)
        * K_i_mineralization ^ n_mineralization / (K_i_mineralization ^ n_mineralization + core ^ n_mineralization)
        * (4300 ^ m_mineralization - atoms_per_cage_transient ^ m_mineralization) / 4300 ^ m_mineralization;

      der(LIP) = -2 * Oxidation + 2 * Reduction + CoreRelease;

      der(FT_cage) = -FT_Degradation + FT_Expression;

      der(core) = 2 * Mineralization + 4 * Nucleation - CoreRelease;

      der(DFP) = Oxidation - Mineralization - Reduction - 2 * Nucleation;

      annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
        coordinateSystem(preserveAspectRatio=false)),
        experiment(
          StopTime=500,
          Tolerance=1e-07,
          __Dymola_Algorithm="Dassl"));

    end CellularFerritinIronStorageModel;

    model EnterocyteMucosalBlockModel "Enterocyte mucosal block"

      //Compartments

      constant Bodylight.Types.Volume cell(
        displayUnit = "L") = 1.4e-12 * 1e-3;
      constant Bodylight.Types.Volume upper(
        displayUnit = "L") = 6.67e-10 * 1e-3;
      constant Bodylight.Types.Volume lower(
        displayUnit = "L") = 8.57e-10 * 1e-3;
      constant Bodylight.Types.Area apical_mem(
        displayUnit = "cm2") = 1.5e-5 * 1e-4;
      constant Bodylight.Types.Area BLM(
        displayUnit = "cm2") = 1.5e-7 * 1e-4;

      //Species

      //cell

      Bodylight.Types.Concentration FT_cage(
        displayUnit = "mol/L",
        start = 2.375189822e-9 * 1e3) "FT-cage";
      Bodylight.Types.Concentration core(
        displayUnit = "mol/L",
        start = 3.682217017e-6 * 1e3) "core";
      Bodylight.Types.Concentration DFP(
        displayUnit = "mol/L",
        start = 1.344769304e-10 * 1e3) "diferric peroxo complex";
      Bodylight.Types.Concentration LIP(
        displayUnit = "mol/L",
        start = 1.223884748e-7 * 1e3) "labile iron pool";
      Bodylight.Types.Concentration IRPs_active(
        displayUnit = "mol/L",
        start = 6.889335935e-11 * 1e3) "iron regulatory proteins (active)";
      Bodylight.Types.Concentration IRPs_inactive(
        displayUnit = "mol/L",
        start = 7.264345126e-12 * 1e3) "iron regulatory proteins (inactive)";

      Real atoms_per_cage_transient "Transient number of Fe atoms that are stored inside the core of a ferritin cage";

      //lower

      Bodylight.Types.Concentration Fe_blood(
        displayUnit = "mol/L",
        start = 4.958456433e-9 * 1e3);
      Bodylight.Types.Concentration body_fe(
        displayUnit = "mol/L",
        start = 0);

      //upper

      Bodylight.Types.Concentration Fe_lumen(
        displayUnit = "mol/L",
        start = 1.25e-8 * 1e3);

      //apical_mem

      BodylightExtension.Types.SurfaceConcentration DMT1(
        start = 3.123079828e-12 * 1e4);
      BodylightExtension.Types.SurfaceConcentration DMT1_vesicular(
        start = 1.876958455e-12 * 1e4);

      //BLM

      BodylightExtension.Types.SurfaceConcentration FPN_active(
        start = 9.979749648e-14 * 1e4);
      BodylightExtension.Types.SurfaceConcentration FPN_internalized(
        start = 2.02503521e-16 * 1e4);

      //Reactions

      //FT expression ( -> FT-cage;  IRPs_active)
      BodylightExtension.Types.MolarReactionRate FT_Expression;
      parameter BodylightExtension.Types.MolarReactionRate k_cat_FT_Expression = 7.68e-14 * 1e3;
      parameter Integer n_FT_Expression = 1;
      parameter Bodylight.Types.Concentration K_FT_Expression(
        displayUnit = "mol/L") = 1.4e-11 * 1e3;

      //FT degradation (FT-cage -> )
      BodylightExtension.Types.MolarReactionRate FT_Degradation;
      parameter BodylightExtension.Types.ReactionRateFirstOrder k_FT_Degradation = 5.461499585e-6;

      //FT degradation core release (core -> LIP; FT-cage core)
      BodylightExtension.Types.MolarReactionRate FT_Degradation_Core_Release;

      //FT Fe oxidation (2 * LIP -> DFP; FT-cage)
      BodylightExtension.Types.MolarReactionRate FT_Fe_Oxidation;
      parameter BodylightExtension.Types.ReactionRateFirstOrder k_cat_FT_Fe_Oxidation = 591;
      parameter Bodylight.Types.Concentration K_m_FT_Fe_Oxidation(
        displayUnit = "mol/L") = 0.35e-3 * 1e3;
      parameter Real n_FT_Fe_Oxidation = 1.3;

      //FT Fe Reduction (DFP -> 2 * LIP)
      BodylightExtension.Types.MolarReactionRate FT_Fe_Reduction;
      parameter BodylightExtension.Types.ReactionRateFirstOrder k_FT_Fe_Reduction = 0.2605;

      //FT Nucleation (2 * DFP -> 4 * core; FT-cage core)
      BodylightExtension.Types.MolarReactionRate FT_Nucleation;
      parameter BodylightExtension.Types.ReactionRateThirdOrder k_cat_FT_Nucleation = 5e7 * 1e-6;
      parameter Bodylight.Types.Concentration K_i_FT_Nucleation(
        displayUnit = "mol/L") = 0.461598e-3 * 1e3;
      parameter Integer n_FT_Nucleation = 4;

      //FT core formation (Mineralization) (DFP -> 2 * core; core)
      BodylightExtension.Types.MolarReactionRate FT_Core_Formation;
      parameter BodylightExtension.Types.ReactionRateFirstOrder k_cat_FT_Core_Formation = 0.101564;
      parameter Bodylight.Types.Concentration K_m_FT_Core_Formation(
        displayUnit = "mol/L") = 5e-06 * 1e3;
      parameter Bodylight.Types.Concentration K_i_FT_Core_Formation(
        displayUnit = "mol/L") = 4.6458e-3 * 1e3;
      parameter Integer n_FT_Core_Formation = 4;
      parameter Integer m_FT_Core_Formation = 8;

      //IRPs degradation (IRPs_active -> IRPs_inactive;  LIP)
      BodylightExtension.Types.MolarReactionRate IRPs_Degradation;
      parameter BodylightExtension.Types.ReactionRateSecondOrder
        k_cat_IRPs_Degradation = 3.99474 * 1e-3;

      //IRPs activation (IRPs_inactive -> IRPs_active)
      BodylightExtension.Types.MolarReactionRate IRPs_Activation;
      parameter BodylightExtension.Types.ReactionRateFirstOrder k_cat_IRPs_Activation = 4.63671e-6;

      //Body Sequestration (Fe_blood -> body_fe)
      BodylightExtension.Types.MolarReactionRate Body_Sequestration;
      parameter BodylightExtension.Types.ReactionRateFirstOrder k_cat_Body_Sequestration = 0.329e-3;

      //Fe Basal Uptake (Fe_blood -> LIP)
      Bodylight.Types.MolarFlowRate Fe_Basal_Uptake;
      parameter Bodylight.Types.VolumeFlowRate k_cat_Fe_Basal_Uptake(
        displayUnit="l/s")= 2.22055e-16 * 1e-3;

      //paracellular Fe transport (Fe_lumen = Fe_blood)
      Bodylight.Types.MolarFlowRate Paracellular_Fe_Transport;
      parameter Bodylight.Types.VolumeFlowRate k_for_Paracellular_Fe_Transport(
        displayUnit="l/s")= 3.87951e-22 * 1e-3;
      parameter Bodylight.Types.VolumeFlowRate k_rev_Paracellular_Fe_Transport(
        displayUnit="l/s")= 3.1746e-15 * 1e-3;

      //DMT1 endocytosis free (DMT1 -> DMT1_vesicular)
      BodylightExtension.Types.MolarFluxPerArea DMT1_Endocytosis_Free;
      parameter BodylightExtension.Types.ReactionRateFirstOrder k_cat_DMT1_Endocytosis_Free = 29.4233;

      //DMT1 endocytosis LIP modified (DMT1 -> DMT1_vesicular;  LIP)
      Bodylight.Types.MolarFlowRate DMT1_Endocytosis_Modified;
      parameter BodylightExtension.Types.DiffusionCoefficient k_cat_DMT1_Endocytosis_Modified(
        displayUnit="cm2/s") = 14.516 * 1e-4;
      parameter Bodylight.Types.Concentration K_m_DMT1_Endocytosis_Modified(
        displayUnit = "mol/L") = 2.80591 * 1e3;
      parameter Real n_DMT1_Endocytosis_Modified = 1.03128;

      //DMT1 fusion (DMT1_vesicular -> DMT1;  DMT1)
      BodylightExtension.Types.MolarFluxPerArea DMT1_Fusion;
      parameter BodylightExtension.Types.ReactionRateFirstOrder k_cat_DMT1_Fusion = 48.9989;

      //DMT1 iron transport (Fe_lumen -> LIP;  DMT1)
      BodylightExtension.Types.MolarFluxPerArea DMT1_Iron_Transport;
      parameter BodylightExtension.Types.ReactionRateFirstOrder k_cat_DMT1_Iron_Transport = 6844.7;
      parameter Bodylight.Types.Concentration K_m_DMT1_Iron_Transport(
        displayUnit = "mol/L") = 2.835 * 1e3;
      parameter Integer n_DMT1_Iron_Transport = 1;

      //FPN-inactivation (FPN_active -> FPN_internalized;  Fe_blood)
      Bodylight.Types.MolarFlowRate FPN_Inactivation;
      parameter BodylightExtension.Types.DiffusionCoefficient k_cat_FPN_Inactivation = 1.44264e-6 * 1e-4;
      parameter Bodylight.Types.Concentration K_m_FPN_Inactivation(
        displayUnit = "mol/L") = 1.22073e-5 * 1e3;
      parameter Real n_FPN_Inactivation = 2.71609;

      //FPN-activation (FPN_internalized -> FPN_active)
      Bodylight.Types.MolarFlowRate FPN_Activation;
      parameter BodylightExtension.Types.DiffusionCoefficient k_cat_FPN_Activation(
        displayUnit="cm2/s") = 4.37363e-13 * 1e-4;

      //FPN iron transport (LIP -> Fe_blood;  FPN_active)
      BodylightExtension.Types.MolarFluxPerArea FPN_Iron_Transport;
      parameter BodylightExtension.Types.ReactionRateFirstOrder k_cat_FPN_Iron_Transport = 1.88317;
      parameter Bodylight.Types.Concentration K_m_FPN_Iron_Transport(
        displayUnit = "mol/L") = 2.31608e-6 * 1e3;
      parameter Integer n_FPN_Iron_Transport = 1;

      //Global Quantities

      parameter Integer H = 24 "H subunits";
      parameter Integer L = 24 - H "L subunits";

      parameter Integer rN = 50;
      parameter Integer rO = 2;

    equation

      atoms_per_cage_transient = core / FT_cage;

      FT_Expression = k_cat_FT_Expression * (1 - IRPs_active ^ n_FT_Expression / (K_FT_Expression ^ n_FT_Expression + IRPs_active ^ n_FT_Expression));

      FT_Degradation = k_FT_Degradation * FT_cage;

      FT_Degradation_Core_Release = k_FT_Degradation * core;

      FT_Fe_Oxidation = (k_cat_FT_Fe_Oxidation * (H + rO) / (24 + rO) * FT_cage * LIP ^ n_FT_Fe_Oxidation)
        / (K_m_FT_Fe_Oxidation ^ n_FT_Fe_Oxidation + LIP ^ n_FT_Fe_Oxidation);

      FT_Fe_Reduction = k_FT_Fe_Reduction * DFP;

      FT_Nucleation = k_cat_FT_Nucleation * DFP ^ 2 * FT_cage * (L + rN) / (24 + rN)
        * K_i_FT_Nucleation ^ n_FT_Nucleation / (K_i_FT_Nucleation ^ n_FT_Nucleation + core ^ n_FT_Nucleation);

      FT_Core_Formation = (k_cat_FT_Core_Formation * DFP * core) / (K_m_FT_Core_Formation + DFP)
        * K_i_FT_Core_Formation ^ n_FT_Core_Formation / (K_i_FT_Core_Formation ^ n_FT_Core_Formation + core ^ n_FT_Core_Formation)
        * (4300 ^ m_FT_Core_Formation - atoms_per_cage_transient ^ m_FT_Core_Formation) / 4300 ^ m_FT_Core_Formation;

      DMT1_Endocytosis_Free = k_cat_DMT1_Endocytosis_Free * DMT1;

      DMT1_Endocytosis_Modified = k_cat_DMT1_Endocytosis_Modified * DMT1 * LIP ^ n_DMT1_Endocytosis_Modified /
        (K_m_DMT1_Endocytosis_Modified ^ n_DMT1_Endocytosis_Modified + LIP ^ n_DMT1_Endocytosis_Modified);

      DMT1_Fusion = k_cat_DMT1_Fusion * DMT1_vesicular;

      DMT1_Iron_Transport = k_cat_DMT1_Iron_Transport * DMT1 * Fe_lumen ^ n_DMT1_Iron_Transport /
        (K_m_DMT1_Iron_Transport ^ n_DMT1_Iron_Transport + Fe_lumen ^ n_DMT1_Iron_Transport);

      FPN_Inactivation = k_cat_FPN_Inactivation * FPN_active * Fe_blood ^ n_FPN_Inactivation / (K_m_FPN_Inactivation ^ n_FPN_Inactivation + Fe_blood ^ n_FPN_Inactivation);

      FPN_Activation = k_cat_FPN_Activation * FPN_internalized;

      FPN_Iron_Transport = k_cat_FPN_Iron_Transport * FPN_active * LIP ^ n_FPN_Iron_Transport / (K_m_FPN_Iron_Transport ^ n_FPN_Iron_Transport + LIP ^ n_FPN_Iron_Transport);

      IRPs_Degradation = k_cat_IRPs_Degradation * IRPs_active * LIP;

      IRPs_Activation = k_cat_IRPs_Activation * IRPs_inactive;

      Body_Sequestration = k_cat_Body_Sequestration * Fe_blood;

      Fe_Basal_Uptake = k_cat_Fe_Basal_Uptake * Fe_blood;

      Paracellular_Fe_Transport = k_for_Paracellular_Fe_Transport * Fe_lumen - k_rev_Paracellular_Fe_Transport * Fe_blood;

      //cell

      der(FT_cage) * cell = -cell * FT_Degradation
        + cell * FT_Expression;

      der(core) * cell = 2 * cell * FT_Core_Formation
        + 4 * cell * FT_Nucleation
        - cell * FT_Degradation_Core_Release;

      der(DFP) * cell = -cell * FT_Core_Formation
        - cell * FT_Fe_Reduction
        + cell * FT_Fe_Oxidation
        - 2 * cell * FT_Nucleation;

      der(LIP) * cell = 2 * cell * FT_Fe_Reduction
        - 2 * cell * FT_Fe_Oxidation
        + apical_mem * DMT1_Iron_Transport
        + Fe_Basal_Uptake
        + cell * FT_Degradation_Core_Release
        - BLM * FPN_Iron_Transport;

      der(IRPs_active) * cell = cell * IRPs_Activation
        - cell * IRPs_Degradation;

      der(IRPs_inactive) * cell = -cell * IRPs_Activation
        + cell * IRPs_Degradation;

      //lower

      der(Fe_blood) * lower = -Fe_Basal_Uptake
        + Paracellular_Fe_Transport
        - lower * Body_Sequestration
        + BLM * FPN_Iron_Transport;

      der(body_fe) * lower = lower * Body_Sequestration;

      //upper

      der(Fe_lumen) * upper = -apical_mem * DMT1_Iron_Transport
        - Paracellular_Fe_Transport;

      //apical_mem

      der(DMT1) * apical_mem = -apical_mem * DMT1_Endocytosis_Free
        + apical_mem * DMT1_Fusion
        - DMT1_Endocytosis_Modified;

      der(DMT1_vesicular) * apical_mem = apical_mem * DMT1_Endocytosis_Free
        - apical_mem * DMT1_Fusion
        + DMT1_Endocytosis_Modified;

      //BLM

      der(FPN_active) * BLM = -FPN_Inactivation
        + FPN_Activation;

      der(FPN_internalized) * BLM = FPN_Inactivation
        - FPN_Activation;

      annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
            coordinateSystem(preserveAspectRatio=false)));
    end EnterocyteMucosalBlockModel;

    model EnterocyteMucosalBlockShortModel "Enterocyte mucosal block (short)"

      //Species

      Bodylight.Types.Concentration FT_cage(
        displayUnit = "mol/L",
        start = 2.375189822e-9 * 1e3) "FT-cage";
      Bodylight.Types.Concentration core(
        displayUnit = "mol/L",
        start = 3.682217017e-6 * 1e3) "core";
      Bodylight.Types.Concentration DFP(
        displayUnit = "mol/L",
        start = 1.344769304e-10 * 1e3) "diferric peroxo complex";
      Bodylight.Types.Concentration LIP(
        displayUnit = "mol/L",
        start = 1.223884748e-7 * 1e3) "labile iron pool";
      Bodylight.Types.Concentration IRPs_active(
        displayUnit = "mol/L",
        start = 6.889335935e-11 * 1e3) "iron regulatory proteins (active)";
      Bodylight.Types.Concentration IRPs_inactive(
        displayUnit = "mol/L",
        start = 7.264345126e-12 * 1e3) "iron regulatory proteins (inactive)";

      Real atoms_per_cage_transient "Transient number of Fe atoms that are stored inside the core of a ferritin cage";

      //Reactions

      //FT expression ( -> FT-cage;  IRPs_active)
      BodylightExtension.Types.MolarReactionRate FT_Expression;
      parameter BodylightExtension.Types.MolarReactionRate k_cat_FT_Expression = 7.68e-14 * 1e3;
      parameter Integer n_FT_Expression = 1;
      parameter Bodylight.Types.Concentration K_FT_Expression(
        displayUnit = "mol/L") = 1.4e-11 * 1e3;

      //FT degradation (FT-cage -> )
      BodylightExtension.Types.MolarReactionRate FT_Degradation;
      parameter BodylightExtension.Types.ReactionRateFirstOrder k_FT_Degradation = 5.461499585e-6;

      //FT degradation core release (core -> LIP; FT-cage core)
      BodylightExtension.Types.MolarReactionRate FT_Degradation_Core_Release;

      //FT Fe oxidation (2 * LIP -> DFP; FT-cage)
      BodylightExtension.Types.MolarReactionRate FT_Fe_Oxidation;
      parameter BodylightExtension.Types.ReactionRateFirstOrder k_cat_FT_Fe_Oxidation = 591;
      parameter Bodylight.Types.Concentration K_m_FT_Fe_Oxidation(
        displayUnit = "mol/L") = 0.35e-3 * 1e3;
      parameter Real n_FT_Fe_Oxidation = 1.3;

      //FT Fe Reduction (DFP -> 2 * LIP)
      BodylightExtension.Types.MolarReactionRate FT_Fe_Reduction;
      parameter BodylightExtension.Types.ReactionRateFirstOrder k_FT_Fe_Reduction = 0.2605;

      //FT Nucleation (2 * DFP -> 4 * core; FT-cage core)
      BodylightExtension.Types.MolarReactionRate FT_Nucleation;
      parameter BodylightExtension.Types.ReactionRateThirdOrder k_cat_FT_Nucleation = 5e7 * 1e-6;
      parameter Bodylight.Types.Concentration K_i_FT_Nucleation(
        displayUnit = "mol/L") = 0.461598e-3 * 1e3;
      parameter Integer n_FT_Nucleation = 4;

      //FT core formation (Mineralization) (DFP -> 2 * core; core)
      BodylightExtension.Types.MolarReactionRate FT_Core_Formation;
      parameter BodylightExtension.Types.ReactionRateFirstOrder k_cat_FT_Core_Formation = 0.101564;
      parameter Bodylight.Types.Concentration K_m_FT_Core_Formation(
        displayUnit = "mol/L") = 5e-06 * 1e3;
      parameter Bodylight.Types.Concentration K_i_FT_Core_Formation(
        displayUnit = "mol/L") = 4.6458e-3 * 1e3;
      parameter Integer n_FT_Core_Formation = 4;
      parameter Integer m_FT_Core_Formation = 8;

      //IRPs degradation (IRPs_active -> IRPs_inactive;  LIP)
      BodylightExtension.Types.MolarReactionRate IRPs_Degradation;
      parameter BodylightExtension.Types.ReactionRateSecondOrder
        k_cat_IRPs_Degradation = 3.99474 * 1e-3;

      //IRPs activation (IRPs_inactive -> IRPs_active)
      BodylightExtension.Types.MolarReactionRate IRPs_Activation;
      parameter BodylightExtension.Types.ReactionRateFirstOrder k_cat_IRPs_Activation = 4.63671e-6;

      //Global Quantities

      parameter Integer H = 24 "H subunits";
      parameter Integer L = 24 - H "L subunits";

      parameter Integer rN = 50;
      parameter Integer rO = 2;

    equation

      atoms_per_cage_transient = core / FT_cage;

      FT_Core_Formation = (k_cat_FT_Core_Formation * DFP * core) / (K_m_FT_Core_Formation + DFP)
        * K_i_FT_Core_Formation ^ n_FT_Core_Formation / (K_i_FT_Core_Formation ^ n_FT_Core_Formation + core ^ n_FT_Core_Formation)
        * (4300 ^ m_FT_Core_Formation - atoms_per_cage_transient ^ m_FT_Core_Formation) / 4300 ^ m_FT_Core_Formation;

      FT_Degradation = k_FT_Degradation * FT_cage;

      FT_Degradation_Core_Release = k_FT_Degradation * core;

      FT_Expression = k_cat_FT_Expression * (1 - IRPs_active ^ n_FT_Expression / (K_FT_Expression ^ n_FT_Expression + IRPs_active ^ n_FT_Expression));

      FT_Fe_Oxidation = (k_cat_FT_Fe_Oxidation * (H + rO) / (24 + rO) * FT_cage * LIP ^ n_FT_Fe_Oxidation)
        / (K_m_FT_Fe_Oxidation ^ n_FT_Fe_Oxidation + LIP ^ n_FT_Fe_Oxidation);

      FT_Fe_Reduction = k_FT_Fe_Reduction * DFP;

      FT_Nucleation = k_cat_FT_Nucleation * DFP ^ 2 * FT_cage * (L + rN) / (24 + rN)
        * K_i_FT_Nucleation ^ n_FT_Nucleation / (K_i_FT_Nucleation ^ n_FT_Nucleation + core ^ n_FT_Nucleation);

      IRPs_Activation = k_cat_IRPs_Activation * IRPs_inactive;

      IRPs_Degradation = k_cat_IRPs_Degradation * IRPs_active * LIP;

      der(FT_cage) = -FT_Degradation
        + FT_Expression;

      der(core) = 2 * FT_Core_Formation
        + 4 * FT_Nucleation
        - FT_Degradation_Core_Release;

      der(DFP) = -FT_Core_Formation
        - FT_Fe_Reduction
        + FT_Fe_Oxidation
        - 2 * FT_Nucleation;

      der(LIP) = 2 * FT_Fe_Reduction
        - 2 * FT_Fe_Oxidation
        + FT_Degradation_Core_Release;

      der(IRPs_active) = IRPs_Activation
        - IRPs_Degradation;

      der(IRPs_inactive) = -IRPs_Activation
        + IRPs_Degradation;

      annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
            coordinateSystem(preserveAspectRatio=false)));
    end EnterocyteMucosalBlockShortModel;

    model EnterocyteMucosalBlockShortModelCellular
      "Enterocyte mucosal block (short) + Cellular (=start values from CellularFerritinIronStorageModel)"

      //Species

      Bodylight.Types.Concentration FT_cage(
        displayUnit = "mol/L",
        start = 5e-09 * 1e3) "FT-cage";
      Bodylight.Types.Concentration core(
        displayUnit = "mol/L",
        start = 7.5e-06 * 1e3) "core";
      Bodylight.Types.Concentration DFP(
        displayUnit = "mol/L",
        start = 0) "diferric peroxo complex";
      Bodylight.Types.Concentration LIP(
        displayUnit = "mol/L",
        start = 1e-05 * 1e3) "labile iron pool";
      Bodylight.Types.Concentration IRPs_active(
        displayUnit = "mol/L",
        start = 6.889335935e-11 * 1e3) "iron regulatory proteins (active)";
      Bodylight.Types.Concentration IRPs_inactive(
        displayUnit = "mol/L",
        start = 7.264345126e-12 * 1e3) "iron regulatory proteins (inactive)";

      Real atoms_per_cage_transient "Transient number of Fe atoms that are stored inside the core of a ferritin cage";

      //Reactions

      //FT expression ( -> FT-cage;  IRPs_active)
      BodylightExtension.Types.MolarReactionRate FT_Expression;
      parameter BodylightExtension.Types.MolarReactionRate k_cat_FT_Expression = 7.68e-14 * 1e3;
      parameter Integer n_FT_Expression = 1;
      parameter Bodylight.Types.Concentration K_FT_Expression(
        displayUnit = "mol/L") = 1.4e-11 * 1e3;

      //FT degradation (FT-cage -> )
      BodylightExtension.Types.MolarReactionRate FT_Degradation;
      parameter BodylightExtension.Types.ReactionRateFirstOrder k_FT_Degradation = 5.461499585e-6;

      //FT degradation core release (core -> LIP; FT-cage core)
      BodylightExtension.Types.MolarReactionRate FT_Degradation_Core_Release;

      //FT Fe oxidation (2 * LIP -> DFP; FT-cage)
      BodylightExtension.Types.MolarReactionRate FT_Fe_Oxidation;
      parameter BodylightExtension.Types.ReactionRateFirstOrder k_cat_FT_Fe_Oxidation = 591;
      parameter Bodylight.Types.Concentration K_m_FT_Fe_Oxidation(
        displayUnit = "mol/L") = 0.35e-3 * 1e3;
      parameter Real n_FT_Fe_Oxidation = 1.3;

      //FT Fe Reduction (DFP -> 2 * LIP)
      BodylightExtension.Types.MolarReactionRate FT_Fe_Reduction;
      parameter BodylightExtension.Types.ReactionRateFirstOrder k_FT_Fe_Reduction = 0.2605;

      //FT Nucleation (2 * DFP -> 4 * core; FT-cage core)
      BodylightExtension.Types.MolarReactionRate FT_Nucleation;
      parameter BodylightExtension.Types.ReactionRateThirdOrder k_cat_FT_Nucleation = 5e7 * 1e-6;
      parameter Bodylight.Types.Concentration K_i_FT_Nucleation(
        displayUnit = "mol/L") = 0.461598e-3 * 1e3;
      parameter Integer n_FT_Nucleation = 4;

      //FT core formation (Mineralization) (DFP -> 2 * core; core)
      BodylightExtension.Types.MolarReactionRate FT_Core_Formation;
      parameter BodylightExtension.Types.ReactionRateFirstOrder k_cat_FT_Core_Formation = 0.101564;
      parameter Bodylight.Types.Concentration K_m_FT_Core_Formation(
        displayUnit = "mol/L") = 5e-06 * 1e3;
      parameter Bodylight.Types.Concentration K_i_FT_Core_Formation(
        displayUnit = "mol/L") = 4.6458e-3 * 1e3;
      parameter Integer n_FT_Core_Formation = 4;
      parameter Integer m_FT_Core_Formation = 8;

      //IRPs degradation (IRPs_active -> IRPs_inactive;  LIP)
      BodylightExtension.Types.MolarReactionRate IRPs_Degradation;
      parameter BodylightExtension.Types.ReactionRateSecondOrder k_cat_IRPs_Degradation = 3.99474 * 1e-3;

      //IRPs activation (IRPs_inactive -> IRPs_active)
      BodylightExtension.Types.MolarReactionRate IRPs_Activation;
      parameter BodylightExtension.Types.ReactionRateFirstOrder k_cat_IRPs_Activation = 4.63671e-6;

      //Global Quantities

      parameter Integer H = 4 "H subunits";
      parameter Integer L = 24 - H "L subunits";

      parameter Integer rN = 50;
      parameter Integer rO = 2;

    equation

      atoms_per_cage_transient = core / FT_cage;

      FT_Core_Formation = (k_cat_FT_Core_Formation * DFP * core) / (K_m_FT_Core_Formation + DFP)
        * K_i_FT_Core_Formation ^ n_FT_Core_Formation / (K_i_FT_Core_Formation ^ n_FT_Core_Formation + core ^ n_FT_Core_Formation)
        * (4300 ^ m_FT_Core_Formation - atoms_per_cage_transient ^ m_FT_Core_Formation) / 4300 ^ m_FT_Core_Formation;

      FT_Degradation = k_FT_Degradation * FT_cage;

      FT_Degradation_Core_Release = k_FT_Degradation * core;

      FT_Expression = k_cat_FT_Expression * (1 - IRPs_active ^ n_FT_Expression / (K_FT_Expression ^ n_FT_Expression + IRPs_active ^ n_FT_Expression));

      FT_Fe_Oxidation = (k_cat_FT_Fe_Oxidation * (H + rO) / (24 + rO) * FT_cage * LIP ^ n_FT_Fe_Oxidation)
        / (K_m_FT_Fe_Oxidation ^ n_FT_Fe_Oxidation + LIP ^ n_FT_Fe_Oxidation);

      FT_Fe_Reduction = k_FT_Fe_Reduction * DFP;

      FT_Nucleation = k_cat_FT_Nucleation * DFP ^ 2 * FT_cage * (L + rN) / (24 + rN)
        * K_i_FT_Nucleation ^ n_FT_Nucleation / (K_i_FT_Nucleation ^ n_FT_Nucleation + core ^ n_FT_Nucleation);

      IRPs_Activation = k_cat_IRPs_Activation * IRPs_inactive;

      IRPs_Degradation = k_cat_IRPs_Degradation * IRPs_active * LIP;

      der(FT_cage) = -FT_Degradation + FT_Expression;

      der(core) = 2 * FT_Core_Formation
        + 4 * FT_Nucleation
        - FT_Degradation_Core_Release;

      der(DFP) = -FT_Core_Formation
        - FT_Fe_Reduction
        + FT_Fe_Oxidation
        - 2 * FT_Nucleation;

      der(LIP) = 2 * FT_Fe_Reduction
        - 2 * FT_Fe_Oxidation
        + FT_Degradation_Core_Release;

      der(IRPs_active) = IRPs_Activation
        - IRPs_Degradation;

      der(IRPs_inactive) = -IRPs_Activation
        + IRPs_Degradation;

      annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
            coordinateSystem(preserveAspectRatio=false)));
    end EnterocyteMucosalBlockShortModelCellular;

    model FerritinIronStorage

      Bodylight.Types.RealIO.ConcentrationInput Fe_total_set annotation (Placement(
            transformation(extent={{-258,-32},{-218,8}}), iconTransformation(extent={{-120,56},
                {-92,84}})));

      Bodylight.Types.RealIO.ConcentrationOutput Fe_in_FT annotation (Placement(
            transformation(extent={{-248,26},{-228,46}}), iconTransformation(extent
              ={{100,64},{120,84}})));
      Bodylight.Types.RealIO.ConcentrationOutput LIP(start = 6.15243e-7 * 1e3) annotation (Placement(
            transformation(extent={{-238,50},{-218,70}}), iconTransformation(extent
          ={{100,32},{120,52}})));

      Bodylight.Types.RealIO.FractionOutput Fract_Fe_in_Ft annotation (Placement(
            transformation(extent={{-248,-26},{-228,-6}}), iconTransformation(
              extent={{100,-28},{120,-8}})));
      Bodylight.Types.RealIO.FractionOutput Fract_LIP annotation (Placement(
            transformation(extent={{-248,-26},{-228,-6}}), iconTransformation(
          extent={{102,-64},{122,-44}})));
      Bodylight.Types.Concentration Fe_total;
      Bodylight.Types.Concentration Fe_total_need = Fe_total_set - Fe_total;

      Bodylight.Types.Concentration core(
        start = 7.5e-06 * 1e3) "core";
      Bodylight.Types.Concentration DFP(
        start = 0) "diferric peroxo complex";
      //Bodylight.Types.Concentration FT_cage;

      Real atoms_per_cage_transient "Transient number of Fe atoms that are stored inside the core of a ferritin cage";

      parameter Integer H = 4 "H subunits";
      parameter Integer L = 24 - H "L subunits";

      parameter Integer rN = 50;
      parameter Integer rO = 2;

      parameter Bodylight.Types.Frequency k_FTlysis = 1.203e-05;
      parameter Bodylight.Types.Frequency k_Fe_total_set_achieve_time = 1e-2;

     // BodylightExtension.Types.MolarReactionRate FT_Expression( start = 16.015e-14 * 1000);

    /*    
    parameter Real FT_Expression(
     quantity = "ReactionRate",
     unit = "mol/(m3.s)",
     displayUnit = "mol/(l.s)")
    
    = 16.015e-14 * 1000;
 */

     /*  
  //FT degradation
  Real FT_Degradation(
    quantity = "ReactionRate",
    unit = "mol/(m3.s)",
    displayUnit = "mol/(l.s)");
 */

      //FT degradation core release
      BodylightExtension.Types.MolarReactionRate CoreRelease;

      //Oxidation (2 LIP -> DFP)
      BodylightExtension.Types.MolarReactionRate Oxidation;

      parameter BodylightExtension.Types.ReactionRateFirstOrder k_cat_oxidation = 591 "catalytic turnover number";
      parameter Bodylight.Types.Concentration K_m_oxidation(
        displayUnit = "mol/L") = 0.35 "Michaelis constant";
      parameter Real n_oxidation = 1.3 "Hill coefficient)";

      //Reduction (DFP -> 2 LIP)
      BodylightExtension.Types.MolarReactionRate Reduction;

      parameter BodylightExtension.Types.ReactionRateFirstOrder k_deg = 0.2605 "rate constant";

      //Nucleation (2 DFP -> 4 core)
      BodylightExtension.Types.MolarReactionRate Nucleation;

      parameter BodylightExtension.Types.ReactionRateThirdOrder k_cat_nucleation = 5e07 * 1e-6 "catalytic turnover number";
      parameter Bodylight.Types.Concentration K_i_nucleation(
        displayUnit = "mol/L") = 0.461598 "inhibition constant";
      parameter Integer n_nucleation = 4 "Hill coefficient";

      //Mineralization (DFP -> 2 core)
      BodylightExtension.Types.MolarReactionRate Mineralization;

      parameter BodylightExtension.Types.ReactionRateFirstOrder k_cat_mineralization = 0.101564 "catalytic turnover number";
      parameter Bodylight.Types.Concentration K_m_mineralization(
        displayUnit = "mol/L") = 5e-03 "Michaelis constant";
      parameter Bodylight.Types.Concentration K_i_mineralization(
        displayUnit = "mol/L") = 4.6458 "inhibition constant";
      parameter Integer n_mineralization = 4 "Hill coefficient";
      parameter Integer m_mineralization = 8 "Hill coefficient";

    // parameter Bodylight.Types.Concentration FT_cage_norm(
    //    displayUnit = "mol/L") = 1.33125e-05 "FT cage (norm)";


      //How to recaltulato Fe_total
      //Fe_total_need = Fe_total_set - Fe_total;
      //core_init = core+Fe_total_need*Fract_Fe_in_Ft;
      //LIP_init=LIP+Fe_total_need*Fract_LIP;

      Bodylight.Types.RealIO.ConcentrationInput FT_cage annotation (Placement(
            transformation(extent={{-258,-32},{-218,8}}), iconTransformation(extent
              ={{-120,-14},{-92,14}})));
    initial equation
     // FT_cage = FT_cage_norm;

      //Fe_total_need = Fe_total_set - Fe_total;
      //LIP = Fe_total_set - (core + DFP);

    equation
      //FT_Expression = 16.015e-14 * 1000;
      // FT_Expression=Ft_expressionIn;

      Fe_in_FT = core + DFP;
      Fe_total = Fe_in_FT + LIP;
      //Fe_total_need = Fe_total_set - Fe_total;

      Fract_Fe_in_Ft = Fe_in_FT / Fe_total;
      Fract_LIP = 1 - Fract_Fe_in_Ft;

      atoms_per_cage_transient = core / FT_cage;

      //FT_Degradation = k_FTlysis * FT_cage;

      CoreRelease = k_FTlysis * core;

      Oxidation = (k_cat_oxidation * (H + rO) / (24 + rO) * FT_cage * LIP ^ n_oxidation)
        / (K_m_oxidation ^ n_oxidation + LIP ^ n_oxidation);

      Reduction = k_deg * DFP;

      Nucleation = k_cat_nucleation * DFP ^ 2 * FT_cage * (L + rN) / (24 + rN)
        * K_i_nucleation ^ n_nucleation / (K_i_nucleation ^ n_nucleation + core ^ n_nucleation);

      Mineralization = (k_cat_mineralization * DFP * core) / (K_m_mineralization + DFP)
        * K_i_mineralization ^ n_mineralization / (K_i_mineralization ^ n_mineralization + core ^ n_mineralization)
        * (4300 ^ m_mineralization - atoms_per_cage_transient ^ m_mineralization) / 4300 ^ m_mineralization;

      //der(FT_cage) = -FT_Degradation + FT_Expression;
      //-FT_Degradation + FT_Expression=0;

      der(LIP) = -2 * Oxidation + 2 * Reduction + CoreRelease
        + Fe_total_need * k_Fe_total_set_achieve_time * Fract_LIP;

      der(core) = 2 * Mineralization + 4 * Nucleation - CoreRelease
        + Fe_total_need * k_Fe_total_set_achieve_time * Fract_Fe_in_Ft;

      der(DFP) = Oxidation - Mineralization - Reduction - 2 * Nucleation;

      annotation (Diagram(coordinateSystem(extent={{-100,-100},{100,100}})), Icon(
            coordinateSystem(extent={{-100,-100},{100,100}}), graphics={
            Rectangle(
              extent={{100,-100},{-102,100}},
              lineColor={28,108,200},
              fillColor={255,255,0},
              fillPattern=FillPattern.Solid),
            Text(
              extent={{-96,-112},{100,-134}},
              textColor={28,108,200},
              textString="%name"),
            Text(
              extent={{-24,42},{-88,96}},
              textColor={28,108,200},
              textString="Fe_total_set"),
            Text(
              extent={{-24,-24},{-88,30}},
              textColor={28,108,200},
              textString="FT_cage"),
            Text(
              extent={{22,54},{92,30}},
              textColor={28,108,200},
              textString="LIP")}));
    end FerritinIronStorage;

    model Test_FT_storage
        extends Modelica.Icons.Example;
      Bodylight.Types.Constants.ConcentrationConst Fe_total(k(displayUnit=
              "mmol/l") = 0.0038)
        annotation (Placement(transformation(extent={{-94,70},{-86,78}})));
      Bodylight.Types.Constants.ConcentrationConst FT_cage(k(displayUnit=
              "mmol/l") = 1.33125e-05)
        annotation (Placement(transformation(extent={{-106,20},{-98,28}})));
      FerritinIronStorage ferritinIronStorage
        annotation (Placement(transformation(extent={{-24,-4},{32,52}})));
    equation
      connect(Fe_total.y, ferritinIronStorage.Fe_total_set) annotation (Line(
            points={{-85,74},{-36,74},{-36,43.6},{-25.68,43.6}}, color={0,0,127}));
      connect(ferritinIronStorage.FT_cage, FT_cage.y)
        annotation (Line(points={{-25.68,24},{-97,24}}, color={0,0,127}));
      annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
            coordinateSystem(preserveAspectRatio=false)));
    end Test_FT_storage;

    model FerritinLysis
      BodylightExtension.Types.RealIO.MolarReactionRateOutput FT_lysis annotation (
          Placement(transformation(extent={{106,-130},{126,-110}}),
            iconTransformation(extent={{100,-10},{120,10}})));
      Bodylight.Types.RealIO.ConcentrationInput FT_cage annotation (Placement(
            transformation(extent={{-334,-58},{-294,-18}}), iconTransformation(
          extent={{-126,-12},{-100,14}})));
      parameter Bodylight.Types.Frequency k_FTlysis = 1.203e-05;
    equation
     FT_lysis = k_FTlysis * FT_cage;
      annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
              Rectangle(
              extent={{-100,100},{100,-100}},
              lineColor={28,108,200},
              fillColor={255,255,0},
              fillPattern=FillPattern.Solid), Text(
              extent={{-110,-106},{110,-124}},
              textColor={28,108,200},
              textString="%name")}), Diagram(coordinateSystem(preserveAspectRatio=false)));

    end FerritinLysis;

    model FT_cage_regulation
      "Control loop or production and degradation FT_cage through LIP"
      Bodylight.Types.RealIO.ConcentrationInput LIP annotation (Placement(
            transformation(extent={{-250,62},{-210,102}}), iconTransformation(
              extent={{-124,42},{-84,82}})));
      Bodylight.Types.RealIO.ConcentrationOutput FT_cage(start = 2.375189822e-9 * 1e3) annotation (Placement(
            transformation(extent={{-224,58},{-204,78}}), iconTransformation(extent
              ={{94,60},{114,80}})));
      Bodylight.Types.RealIO.ConcentrationOutput IRPs_active(start = 6.889335935e-11 * 1e3) annotation (Placement(
            transformation(extent={{-240,52},{-220,72}}), iconTransformation(extent
              ={{94,-24},{114,-4}})));
      Bodylight.Types.RealIO.ConcentrationOutput IRPs_inactive(start = 7.264345126e-12 * 1e3) annotation (
          Placement(transformation(extent={{-234,12},{-214,32}}),
            iconTransformation(extent={{94,-62},{114,-42}})));
      BodylightExtension.Types.MolarReactionRate FT_Expression;
      BodylightExtension.Types.MolarReactionRate k_cat_FT_Expression = 7.68e-14 * 1e3;

      //FT expression ( -> FT-cage;  IRPs_active)
      parameter Integer n_FT_Expression = 1;
      parameter Bodylight.Types.Concentration K_FT_Expression(
        displayUnit = "mol/L") = 1.4e-11 * 1e3;

      //FT degradation (FT-cage -> )
      BodylightExtension.Types.MolarReactionRate FT_Degradation;
      parameter BodylightExtension.Types.ReactionRateFirstOrder k_FT_Degradation = 5.461499585e-6;

      //IRPs degradation (IRPs_active -> IRPs_inactive;  LIP)
      BodylightExtension.Types.MolarReactionRate IRPs_Degradation;
      parameter BodylightExtension.Types.ReactionRateSecondOrder k_cat_IRPs_Degradation = 3.99474 * 1e-3;

      //IRPs activation (IRPs_inactive -> IRPs_active)
      BodylightExtension.Types.MolarReactionRate IRPs_Activation;
      parameter BodylightExtension.Types.ReactionRateFirstOrder k_cat_IRPs_Activation = 4.63671e-6;

    equation
    //  FT_Degradation = k_FT_Degradation * FT_cage;

    //  FT_Expression = k_cat_FT_Expression * (1 -
    //    IRPs_active ^ n_FT_Expression
    //    / (K_FT_Expression ^ n_FT_Expression + IRPs_active ^ n_FT_Expression));

    //  IRPs_Activation = k_cat_IRPs_Activation * IRPs_inactive;
    //  IRPs_Degradation = k_cat_IRPs_Degradation * IRPs_active * LIP;

    //  der(FT_cage) = FT_Degradation + FT_Expression;
    //  der(IRPs_active) = IRPs_Activation
    //    - IRPs_Degradation;
    //  der(IRPs_inactive) = -IRPs_Activation
    //    + IRPs_Degradation;

      FT_Degradation = k_FT_Degradation * FT_cage;

      FT_Expression = k_cat_FT_Expression * (1 -
        IRPs_active ^ n_FT_Expression
        / (K_FT_Expression ^ n_FT_Expression + IRPs_active ^ n_FT_Expression));

      IRPs_Activation = k_cat_IRPs_Activation * IRPs_inactive;

      IRPs_Degradation = k_cat_IRPs_Degradation * IRPs_active * LIP;

      der(FT_cage) = -FT_Degradation
        + FT_Expression;

      der(IRPs_active) = IRPs_Activation
        - IRPs_Degradation;

      der(IRPs_inactive) = -IRPs_Activation
        + IRPs_Degradation;

      annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
              Rectangle(
              extent={{100,-100},{-100,100}},
              lineColor={28,108,200},
              fillColor={255,255,0},
              fillPattern=FillPattern.Solid),
            Text(
              extent={{-78,-110},{68,-130}},
              textColor={28,108,200},
              textString="%name"),
            Text(
              extent={{-78,72},{-24,48}},
              textColor={28,108,200},
              textString="LIP"),
            Text(
              extent={{-4,80},{84,54}},
              textColor={28,108,200},
              textString="FT_cage")}),          Diagram(coordinateSystem(
              preserveAspectRatio=false)));

    end FT_cage_regulation;

    model Test_FT_cageRegulation
        extends Modelica.Icons.Example;
      Bodylight.Types.Constants.ConcentrationConst Fe_total(k(displayUnit=
              "mmol/l") = 0.00380474)
        annotation (Placement(transformation(extent={{-94,70},{-86,78}})));
      FerritinIronStorage ferritinIronStorage
        annotation (Placement(transformation(extent={{-24,-4},{32,52}})));
      FT_cage_regulation fT_cage_regulation
        annotation (Placement(transformation(extent={{-40,-82},{-4,-46}})));
      Bodylight.Types.Constants.ConcentrationConst LIP_in(k(displayUnit="mol/l")
           = 0.0001223884748)
        annotation (Placement(transformation(extent={{-98,-58},{-90,-50}})));
      Bodylight.Types.Constants.ConcentrationConst FT_cage_in(k(displayUnit=
              "mol/l") = 2.375189822e-06)
        annotation (Placement(transformation(extent={{-92,20},{-84,28}})));
    equation
      connect(Fe_total.y, ferritinIronStorage.Fe_total_set) annotation (Line(
            points={{-85,74},{-36,74},{-36,43.6},{-25.68,43.6}}, color={0,0,127}));
      connect(fT_cage_regulation.FT_cage, ferritinIronStorage.FT_cage)
        annotation (Line(points={{-3.28,-51.4},{34,-51.4},{34,-98},{-70,-98},{
              -70,24},{-25.68,24}}, color={0,0,127}));
      connect(ferritinIronStorage.LIP, fT_cage_regulation.LIP) annotation (Line(
            points={{34.8,35.76},{62,35.76},{62,-30},{-58,-30},{-58,-52.84},{
              -40.72,-52.84}}, color={0,0,127}));
      annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
            coordinateSystem(preserveAspectRatio=false)));
    end Test_FT_cageRegulation;

    model FerritinCageBlockShortModel "Enterocyte mucosal block (short)"

      //cell

      Bodylight.Types.Concentration FT_cage(
        displayUnit = "mol/L",
        start = 2.375189822e-9 * 1e3) "FT-cage";
      Bodylight.Types.Concentration core(
        displayUnit = "mol/L",
        start = 3.682217017e-6 * 1e3) "core";
      Bodylight.Types.Concentration DFP(
        displayUnit = "mol/L",
        start = 1.344769304e-10 * 1e3) "diferric peroxo complex";
      Bodylight.Types.Concentration LIP(
        displayUnit = "mol/L",
        start = 1.223884748e-7 * 1e3) "labile iron pool";
      Bodylight.Types.Concentration IRPs_active(
        displayUnit = "mol/L",
        start = 6.889335935e-11 * 1e3) "iron regulatory proteins (active)";
      Bodylight.Types.Concentration IRPs_inactive(
        displayUnit = "mol/L",
        start = 7.264345126e-12 * 1e3) "iron regulatory proteins (inactive)";

      Real atoms_per_cage_transient "Transient number of Fe atoms that are stored inside the core of a ferritin cage";

      //Reactions

      //FT expression ( -> FT-cage;  IRPs_active)
      BodylightExtension.Types.MolarReactionRate FT_Expression;
      BodylightExtension.Types.MolarReactionRate k_cat_FT_Expression=
          7.68e-14*1e3;
      parameter Integer n_FT_Expression = 1;
      parameter Bodylight.Types.Concentration K_FT_Expression(
        displayUnit = "mol/L") = 1.4e-11 * 1e3;

      //FT degradation (FT-cage -> )
      BodylightExtension.Types.MolarReactionRate FT_Degradation;
      parameter Bodylight.Types.Frequency k_FT_Degradation = 5.461499585e-6;

      //FT degradation core release (core -> LIP; FT-cage core)
      BodylightExtension.Types.MolarReactionRate
        FT_Degradation_Core_Release;

      //FT Fe oxidation (2 * LIP -> DFP; FT-cage)
      BodylightExtension.Types.MolarReactionRate FT_Fe_Oxidation;
      parameter Bodylight.Types.Frequency k_cat_FT_Fe_Oxidation = 591;
      parameter Bodylight.Types.Concentration K_m_FT_Fe_Oxidation(
        displayUnit = "mol/L") = 0.35e-3 * 1e3;
      parameter Real n_FT_Fe_Oxidation = 1.3;

      //FT Fe Reduction (DFP -> 2 * LIP)
      BodylightExtension.Types.MolarReactionRate FT_Fe_Reduction;
      parameter Bodylight.Types.Frequency k_FT_Fe_Reduction = 0.2605;

      //FT Nucleation (2 * DFP -> 4 * core; FT-cage core)
      BodylightExtension.Types.MolarReactionRate FT_Nucleation;
      parameter BodylightExtension.Types.ReactionRateThirdOrder k_cat_FT_Nucleation = 5e7 * 1e-6;
      parameter Bodylight.Types.Concentration K_i_FT_Nucleation(
        displayUnit = "mol/L") = 0.461598e-3 * 1e3;
      parameter Integer n_FT_Nucleation = 4;

      //FT core formation (Mineralization) (DFP -> 2 * core; core)
      BodylightExtension.Types.MolarReactionRate FT_Core_Formation;
      parameter Bodylight.Types.Frequency k_cat_FT_Core_Formation = 0.101564;
      parameter Bodylight.Types.Concentration K_m_FT_Core_Formation(
        displayUnit = "mol/L") = 5e-06 * 1e3;
      parameter Bodylight.Types.Concentration K_i_FT_Core_Formation(
        displayUnit = "mol/L") = 4.6458e-3 * 1e3;
      parameter Integer n_FT_Core_Formation = 4;
      parameter Integer m_FT_Core_Formation = 8;

      //IRPs degradation (IRPs_active -> IRPs_inactive;  LIP)
      BodylightExtension.Types.MolarReactionRate IRPs_Degradation;
      parameter BodylightExtension.Types.ReactionRateSecondOrder
        k_cat_IRPs_Degradation=3.99474*1e-3;

      //IRPs activation (IRPs_inactive -> IRPs_active)
      BodylightExtension.Types.MolarReactionRate IRPs_Activation;
      parameter Bodylight.Types.Frequency k_cat_IRPs_Activation = 4.63671e-6;

      //added parameter
      parameter Bodylight.Types.Frequency k_Fe_total_set_achieve_time=1e-2;

      //Global Quantities

      parameter Integer H = 24 "H subunits";
      parameter Integer L = 24 - H "L subunits";

      parameter Integer rN = 50;
      parameter Integer rO = 2;

      // Added variables
      Bodylight.Types.RealIO.ConcentrationInput Fe_total_set annotation (Placement(
            transformation(extent={{-258,-32},{-218,8}}), iconTransformation(extent={{-122,
                -12},{-94,16}})));
      Bodylight.Types.Concentration Fe_total;
      Bodylight.Types.Concentration Fe_total_need=Fe_total_set - Fe_total;
      Bodylight.Types.Concentration Fe_in_FT;
      Bodylight.Types.Fraction Fract_Fe_in_Ft;
      Bodylight.Types.Fraction Fract_LIP;

    equation

      Fe_in_FT = core + DFP;
      Fe_total = Fe_in_FT + LIP;
      Fract_Fe_in_Ft = Fe_in_FT / Fe_total;
      Fract_LIP = 1 - Fract_Fe_in_Ft;

      atoms_per_cage_transient = core / FT_cage;

      FT_Core_Formation = (k_cat_FT_Core_Formation * DFP * core) / (K_m_FT_Core_Formation + DFP)
        * K_i_FT_Core_Formation ^ n_FT_Core_Formation / (K_i_FT_Core_Formation ^ n_FT_Core_Formation + core ^ n_FT_Core_Formation)
        * (4300 ^ m_FT_Core_Formation - atoms_per_cage_transient ^ m_FT_Core_Formation) / 4300 ^ m_FT_Core_Formation;

      FT_Degradation = k_FT_Degradation * FT_cage;

      FT_Degradation_Core_Release = k_FT_Degradation * core;

      FT_Expression = k_cat_FT_Expression * (1 - IRPs_active ^ n_FT_Expression / (K_FT_Expression ^ n_FT_Expression + IRPs_active ^ n_FT_Expression));

      FT_Fe_Oxidation = (k_cat_FT_Fe_Oxidation * (H + rO) / (24 + rO) * FT_cage * LIP ^ n_FT_Fe_Oxidation)
        / (K_m_FT_Fe_Oxidation ^ n_FT_Fe_Oxidation + LIP ^ n_FT_Fe_Oxidation);

      FT_Fe_Reduction = k_FT_Fe_Reduction * DFP;

      FT_Nucleation = k_cat_FT_Nucleation * DFP ^ 2 * FT_cage * (L + rN) / (24 + rN)
        * K_i_FT_Nucleation ^ n_FT_Nucleation / (K_i_FT_Nucleation ^ n_FT_Nucleation + core ^ n_FT_Nucleation);

      IRPs_Activation = k_cat_IRPs_Activation * IRPs_inactive;

      IRPs_Degradation = k_cat_IRPs_Degradation * IRPs_active * LIP;

      der(FT_cage) = - FT_Degradation +  FT_Expression;

      der(core) = 2 *  FT_Core_Formation
        + 4  * FT_Nucleation
        - FT_Degradation_Core_Release
        + Fe_total_need * k_Fe_total_set_achieve_time * Fract_Fe_in_Ft;


      //der(core) = 2 * Mineralization + 4 * Nucleation - CoreRelease
      //  + Fe_total_need*k_Fe_total_set_achieve_time*Fract_Fe_in_Ft;

      der(DFP)  = - FT_Core_Formation
        -  FT_Fe_Reduction
        +  FT_Fe_Oxidation
        - 2 * FT_Nucleation;

      der(LIP) = 2 * FT_Fe_Reduction
        - 2 * FT_Fe_Oxidation
        + FT_Degradation_Core_Release
        + Fe_total_need * k_Fe_total_set_achieve_time * Fract_LIP;

     // der(LIP) = -2 * Oxidation + 2 * Reduction + CoreRelease
     //   +Fe_total_need*k_Fe_total_set_achieve_time*Fract_LIP;

      der(IRPs_active) = IRPs_Activation
        - IRPs_Degradation;

      der(IRPs_inactive) = -IRPs_Activation
        + IRPs_Degradation;

      annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
              Rectangle(
              extent={{-100,100},{100,-100}},
              lineColor={28,108,200},
              fillColor={255,255,0},
              fillPattern=FillPattern.Solid)}),                      Diagram(
            coordinateSystem(preserveAspectRatio=false)));
    end FerritinCageBlockShortModel;

    model Test_FerritinCageBlockShortModel
        extends Modelica.Icons.Example;
      Bodylight.Types.Constants.ConcentrationConst Fe_total(k(displayUnit=
              "mmol/l") = 0.00380474)
        annotation (Placement(transformation(extent={{-86,38},{-78,46}})));
      FerritinCageBlockShortModel ferritinCageBlockShortModel
        annotation (Placement(transformation(extent={{-44,24},{-2,64}})));
    equation
      connect(Fe_total.y, ferritinCageBlockShortModel.Fe_total_set) annotation
        (Line(points={{-77,42},{-54,42},{-54,44.4},{-45.68,44.4}}, color={0,0,
              127}));
      annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
            coordinateSystem(preserveAspectRatio=false)),
        experiment(
          StopTime=4000000,
          __Dymola_NumberOfIntervals=5000,
          __Dymola_Algorithm="Dassl"));
    end Test_FerritinCageBlockShortModel;

    model FerritinCageBlockShortModel_withOutputs
      "Enterocyte mucosal block (short)"

      //cell

      Bodylight.Types.Concentration FT_cage(
        displayUnit = "mol/L",
        start = 2.375189822e-9 * 1e3) "FT-cage";
      Bodylight.Types.Concentration core(
        displayUnit = "mol/L",
        start = 3.682217017e-6 * 1e3) "core";
      Bodylight.Types.Concentration DFP(
        displayUnit = "mol/L",
        start = 1.344769304e-10 * 1e3) "diferric peroxo complex";
      Bodylight.Types.Concentration IRPs_active(
        displayUnit = "mol/L",
        start = 6.889335935e-11 * 1e3) "iron regulatory proteins (active)";
      Bodylight.Types.Concentration IRPs_inactive(
        displayUnit = "mol/L",
        start = 7.264345126e-12 * 1e3) "iron regulatory proteins (inactive)";

      Real atoms_per_cage_transient "Transient number of Fe atoms that are stored inside the core of a ferritin cage";

      //Reactions

      //FT expression ( -> FT-cage;  IRPs_active)
      BodylightExtension.Types.MolarReactionRate FT_Expression;
      BodylightExtension.Types.MolarReactionRate k_cat_FT_Expression=
          7.68e-14*1e3;
      parameter Integer n_FT_Expression = 1;
      parameter Bodylight.Types.Concentration K_FT_Expression(
        displayUnit = "mol/L") = 1.4e-11 * 1e3;

      //FT degradation (FT-cage -> )
      BodylightExtension.Types.MolarReactionRate FT_Degradation;
      parameter Bodylight.Types.Frequency k_FT_Degradation = 5.461499585e-6;

      //FT degradation core release (core -> LIP; FT-cage core)
      BodylightExtension.Types.MolarReactionRate
        FT_Degradation_Core_Release;

      //FT Fe oxidation (2 * LIP -> DFP; FT-cage)
      BodylightExtension.Types.MolarReactionRate FT_Fe_Oxidation;
      parameter Bodylight.Types.Frequency k_cat_FT_Fe_Oxidation = 591;
      parameter Bodylight.Types.Concentration K_m_FT_Fe_Oxidation(
        displayUnit = "mol/L") = 0.35e-3 * 1e3;
      parameter Real n_FT_Fe_Oxidation = 1.3;

      //FT Fe Reduction (DFP -> 2 * LIP)
      BodylightExtension.Types.MolarReactionRate FT_Fe_Reduction;
      parameter Bodylight.Types.Frequency k_FT_Fe_Reduction = 0.2605;

      //FT Nucleation (2 * DFP -> 4 * core; FT-cage core)
      BodylightExtension.Types.MolarReactionRate FT_Nucleation;
      parameter BodylightExtension.Types.ReactionRateThirdOrder k_cat_FT_Nucleation = 5e7 * 1e-6;
      parameter Bodylight.Types.Concentration K_i_FT_Nucleation(
        displayUnit = "mol/L") = 0.461598e-3 * 1e3;
      parameter Integer n_FT_Nucleation = 4;

      //FT core formation (Mineralization) (DFP -> 2 * core; core)
      BodylightExtension.Types.MolarReactionRate FT_Core_Formation;
      parameter Bodylight.Types.Frequency k_cat_FT_Core_Formation = 0.101564;
      parameter Bodylight.Types.Concentration K_m_FT_Core_Formation(
        displayUnit = "mol/L") = 5e-06 * 1e3;
      parameter Bodylight.Types.Concentration K_i_FT_Core_Formation(
        displayUnit = "mol/L") = 4.6458e-3 * 1e3;
      parameter Integer n_FT_Core_Formation = 4;
      parameter Integer m_FT_Core_Formation = 8;

      //IRPs degradation (IRPs_active -> IRPs_inactive;  LIP)
      BodylightExtension.Types.MolarReactionRate IRPs_Degradation;
      parameter BodylightExtension.Types.ReactionRateSecondOrder
        k_cat_IRPs_Degradation=3.99474*1e-3;

      //IRPs activation (IRPs_inactive -> IRPs_active)
      BodylightExtension.Types.MolarReactionRate IRPs_Activation;
      parameter Bodylight.Types.Frequency k_cat_IRPs_Activation = 4.63671e-6;

      //added parameter
      parameter Bodylight.Types.Frequency k_Fe_total_set_achieve_time=1e-2;

      //Global Quantities

      parameter Integer H = 24 "H subunits";
      parameter Integer L = 24 - H "L subunits";

      parameter Integer rN = 50;
      parameter Integer rO = 2;

      Bodylight.Types.RealIO.ConcentrationInput Fe_total_set;

      // Added variables

      Bodylight.Types.Concentration Fe_total;
      Bodylight.Types.Concentration Fe_total_need=Fe_total_set - Fe_total;
      Bodylight.Types.RealIO.ConcentrationOutput LIP(start = 6.15243e-7 * 1e3) annotation (Placement(
            transformation(extent={{112,16},{132,36}}),   iconTransformation(extent
          ={{100,32},{120,52}})));
      Bodylight.Types.RealIO.ConcentrationOutput Fe_in_FT annotation (Placement(
            transformation(extent={{124,44},{144,64}}),   iconTransformation(extent={{100,58},
                {120,78}})));
      Bodylight.Types.RealIO.FractionOutput Fract_Fe_in_Ft annotation (Placement(
            transformation(extent={{114,-46},{134,-26}}), iconTransformation(extent
              ={{100,-34},{120,-14}})));
      Bodylight.Types.RealIO.FractionOutput Fract_LIP annotation (Placement(
            transformation(extent={{132,-68},{152,-48}}), iconTransformation(extent
              ={{100,-70},{120,-50}})));
      Bodylight.Types.RealIO.ConcentrationInput Fe_total_norm annotation (Placement(
            transformation(extent={{-258,-32},{-218,8}}), iconTransformation(extent={{-118,10},
                {-90,38}})));
      Bodylight.Types.RealIO.FractionInput Fe_total_fract annotation (Placement(
            transformation(extent={{-316,-106},{-276,-66}}), iconTransformation(
              extent={{-116,-46},{-92,-22}})));
    equation

      Fe_total_set= Fe_total_norm*Fe_total_fract;
      Fe_in_FT = core + DFP;
      Fe_total = Fe_in_FT + LIP;
      Fract_Fe_in_Ft = Fe_in_FT / Fe_total;
      Fract_LIP = 1 - Fract_Fe_in_Ft;

      atoms_per_cage_transient = core / FT_cage;

      FT_Core_Formation = (k_cat_FT_Core_Formation * DFP * core) / (K_m_FT_Core_Formation + DFP)
        * K_i_FT_Core_Formation ^ n_FT_Core_Formation / (K_i_FT_Core_Formation ^ n_FT_Core_Formation + core ^ n_FT_Core_Formation)
        * (4300 ^ m_FT_Core_Formation - atoms_per_cage_transient ^ m_FT_Core_Formation) / 4300 ^ m_FT_Core_Formation;

      FT_Degradation = k_FT_Degradation * FT_cage;

      FT_Degradation_Core_Release = k_FT_Degradation * core;

      FT_Expression = k_cat_FT_Expression * (1 - IRPs_active ^ n_FT_Expression / (K_FT_Expression ^ n_FT_Expression + IRPs_active ^ n_FT_Expression));

      FT_Fe_Oxidation = (k_cat_FT_Fe_Oxidation * (H + rO) / (24 + rO) * FT_cage * LIP ^ n_FT_Fe_Oxidation)
        / (K_m_FT_Fe_Oxidation ^ n_FT_Fe_Oxidation + LIP ^ n_FT_Fe_Oxidation);

      FT_Fe_Reduction = k_FT_Fe_Reduction * DFP;

      FT_Nucleation = k_cat_FT_Nucleation * DFP ^ 2 * FT_cage * (L + rN) / (24 + rN)
        * K_i_FT_Nucleation ^ n_FT_Nucleation / (K_i_FT_Nucleation ^ n_FT_Nucleation + core ^ n_FT_Nucleation);

      IRPs_Activation = k_cat_IRPs_Activation * IRPs_inactive;

      IRPs_Degradation = k_cat_IRPs_Degradation * IRPs_active * LIP;

      der(FT_cage) = - FT_Degradation +  FT_Expression;

      der(core) = 2 *  FT_Core_Formation
        + 4  * FT_Nucleation
        - FT_Degradation_Core_Release
        + Fe_total_need * k_Fe_total_set_achieve_time * Fract_Fe_in_Ft;

      //der(core) = 2 * Mineralization + 4 * Nucleation - CoreRelease
      //  + Fe_total_need*k_Fe_total_set_achieve_time*Fract_Fe_in_Ft;

      der(DFP)  = - FT_Core_Formation
        -  FT_Fe_Reduction
        +  FT_Fe_Oxidation
        - 2 * FT_Nucleation;

      der(LIP) = 2 * FT_Fe_Reduction
        - 2 * FT_Fe_Oxidation
        + FT_Degradation_Core_Release
        + Fe_total_need * k_Fe_total_set_achieve_time * Fract_LIP;

     // der(LIP) = -2 * Oxidation + 2 * Reduction + CoreRelease
     //   +Fe_total_need*k_Fe_total_set_achieve_time*Fract_LIP;

      der(IRPs_active) = IRPs_Activation
        - IRPs_Degradation;

      der(IRPs_inactive) = -IRPs_Activation
        + IRPs_Degradation;

      annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
              Rectangle(
              extent={{-100,100},{100,-100}},
              lineColor={28,108,200},
              fillColor={255,255,0},
              fillPattern=FillPattern.Solid),
            Text(
              extent={{54,70},{94,60}},
              textColor={28,108,200},
              textString="Fe_in_FT"),
            Text(
              extent={{42,44},{82,34}},
              textColor={28,108,200},
              textString="LIP"),
            Text(
              extent={{58,-18},{88,-24}},
              textColor={28,108,200},
              textString="Fract_Fe_inFt"),
            Text(
              extent={{58,-50},{80,-62}},
              textColor={28,108,200},
              textString="Fract_LIP1"),
            Text(
              extent={{-80,44},{-46,8}},
              textColor={28,108,200},
              textString="Fe_total"),
            Text(
              extent={{-82,-22},{-22,-48}},
              textColor={28,108,200},
              textString="Fe_total_fract")}),                        Diagram(
            coordinateSystem(preserveAspectRatio=false)));
    end FerritinCageBlockShortModel_withOutputs;

    model Test_FerritinCageBlockShortModel_withOutputs
        extends Modelica.Icons.Example;
      Bodylight.Types.Constants.ConcentrationConst Fe_total(k(displayUnit=
              "mmol/l") = 0.0380474)
        annotation (Placement(transformation(extent={{-92,24},{-84,32}})));
      FerritinCageBlockShortModel_withOutputs
        ferritinCageBlockShortModel_withOutputs
        annotation (Placement(transformation(extent={{-44,-32},{56,66}})));
    equation
      connect(Fe_total.y, ferritinCageBlockShortModel_withOutputs.Fe_total_set)
        annotation (Line(points={{-83,28},{-64.5,28},{-64.5,28.76},{-46,28.76}},
            color={0,0,127}));
      annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
            coordinateSystem(preserveAspectRatio=false)),
        experiment(
          StopTime=4000000,
          __Dymola_NumberOfIntervals=5000,
          __Dymola_Algorithm="Dassl"));
    end Test_FerritinCageBlockShortModel_withOutputs;

    model Test_FerritinCageBlockShortModel_withOutputs_k
        extends Modelica.Icons.Example;
      Bodylight.Types.Constants.ConcentrationConst Fe_total(k(displayUnit=
              "mmol/l") = 0.00380474)
        annotation (Placement(transformation(extent={{-92,24},{-84,32}})));
      FerritinCageBlockShortModel_withOutputs
        ferritinCageBlockShortModel_withOutputs
        annotation (Placement(transformation(extent={{-44,-32},{56,66}})));
      Bodylight.Types.Constants.FractionConst fraction(k=1)
        annotation (Placement(transformation(extent={{-92,-6},{-84,2}})));
    equation
      connect(Fe_total.y, ferritinCageBlockShortModel_withOutputs.Fe_total_set)
        annotation (Line(points={{-83,28},{-64.5,28},{-64.5,28.76},{-46,28.76}},
            color={0,0,127}));
      connect(fraction.y, ferritinCageBlockShortModel_withOutputs.Fe_total_fract)
        annotation (Line(points={{-83,-2},{-58,-2},{-58,0.34},{-46,0.34}},
            color={0,0,127}));
      annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
            coordinateSystem(preserveAspectRatio=false)),
        experiment(
          StopTime=4000000,
          __Dymola_NumberOfIntervals=5000,
          __Dymola_Algorithm="Dassl"));
    end Test_FerritinCageBlockShortModel_withOutputs_k;
  end models;
  annotation (uses(Modelica(version="4.0.0"), Bodylight(version="1.0")));
end EnterocyteMucosalBlock;
