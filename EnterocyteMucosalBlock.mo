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
        start = 1.2239676e-7 * 1e3) "labile iron pool";
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
      Bodylight.Types.RealIO.ConcentrationOutput LIP(start = 1.2239676e-7 * 1e3) annotation (Placement(
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
              extent={{38,70},{94,60}},
              textColor={28,108,200},
              textString="Fe_in_FT",
              horizontalAlignment=TextAlignment.Right),
            Text(
              extent={{32,44},{94,34}},
              textColor={28,108,200},
              textString="LIP",
              horizontalAlignment=TextAlignment.Right),
            Text(
              extent={{4,-22},{90,-32}},
              textColor={28,108,200},
              textString="Fract_Fe_inFt",
              horizontalAlignment=TextAlignment.Right),
            Text(
              extent={{26,-56},{92,-68}},
              textColor={28,108,200},
              horizontalAlignment=TextAlignment.Right,
              textString="Fract_LIP"),
            Text(
              extent={{-86,30},{22,16}},
              textColor={28,108,200},
              horizontalAlignment=TextAlignment.Left,
              textString="Fe_total_norm"),
            Text(
              extent={{-90,-28},{14,-40}},
              textColor={28,108,200},
              horizontalAlignment=TextAlignment.Left,
              textString="Fe_total_fract")}),                        Diagram(
            coordinateSystem(preserveAspectRatio=false)));
    end FerritinCageBlockShortModel_withOutputs;

    model Test_FerritinCageBlockShortModel_withOutputs
        extends Modelica.Icons.Example;
      Bodylight.Types.Constants.ConcentrationConst Fe_total_norm(k(displayUnit=
              "mmol/l") = 0.00380474)
        annotation (Placement(transformation(extent={{-92,24},{-84,32}})));
      FerritinCageBlockShortModel_withOutputs
        ferritinCageBlockShortModel_withOutputs
        annotation (Placement(transformation(extent={{-44,-32},{56,66}})));
      Bodylight.Types.Constants.FractionConst fraction(k=1)
        annotation (Placement(transformation(extent={{-92,-6},{-84,2}})));
    equation
      connect(Fe_total_norm.y, ferritinCageBlockShortModel_withOutputs.Fe_total_set)
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
    end Test_FerritinCageBlockShortModel_withOutputs;
  end models;

  package FeMetabolism

   package algorithms

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
        Freg1_o1 := 1 + facB1(facB1_i1=Freg1_i1, facB1_i2=Freg1_i2)*
          537.649 + facB2(facB2_i1=Freg1_i1, facB2_i2=Freg1_i2)*4972.6
           + facST(facST_i1=Freg1_i1, facST_i2=Freg1_i2)*584.75 +
          facB1(facB1_i1=Freg1_i1, facB1_i2=Freg1_i2)*
          facB2(facB2_i1=Freg1_i1, facB2_i2=Freg1_i2)*537.649*4972.6
           + facB1(facB1_i1=Freg1_i1, facB1_i2=Freg1_i2)*
          facST(facST_i1=Freg1_i1, facST_i2=Freg1_i2)*537.649*584.75*
          5.3869 + facB2(facB2_i1=Freg1_i1, facB2_i2=Freg1_i2)*
          facST(facST_i1=Freg1_i1, facST_i2=Freg1_i2)*4972.6*584.75 +
          facB1(facB1_i1=Freg1_i1, facB1_i2=Freg1_i2)*
          facB2(facB2_i1=Freg1_i1, facB2_i2=Freg1_i2)*
          facST(facST_i1=Freg1_i1, facST_i2=Freg1_i2)*537.649*4972.6*
          584.75*5.3869;
      end Freg1;

      function Freg2 "Freg2 function (id f13)"
        input Real Freg2_i1;
        input Real Freg2_i2;
        output Real Freg2_o1 "result";
      algorithm
        Freg2_o1 := 1 + facB1(facB1_i1=Freg2_i1, facB1_i2=Freg2_i2) +
          facB2(facB2_i1=Freg2_i1, facB2_i2=Freg2_i2) +
          facST(facST_i1=Freg2_i1, facST_i2=Freg2_i2) +
          facB1(facB1_i1=Freg2_i1, facB1_i2=Freg2_i2)*
          facB2(facB2_i1=Freg2_i1, facB2_i2=Freg2_i2) +
          facB1(facB1_i1=Freg2_i1, facB1_i2=Freg2_i2)*
          facST(facST_i1=Freg2_i1, facST_i2=Freg2_i2)*5.3869 +
          facB2(facB2_i1=Freg2_i1, facB2_i2=Freg2_i2)*
          facST(facST_i1=Freg2_i1, facST_i2=Freg2_i2) +
          facB1(facB1_i1=Freg2_i1, facB1_i2=Freg2_i2)*
          facB2(facB2_i1=Freg2_i1, facB2_i2=Freg2_i2)*
          facST(facST_i1=Freg2_i1, facST_i2=Freg2_i2)*5.3869;
      end Freg2;

      function hill "Hill function (id f4)"
        input Real hill_i1;
        input Real hill_i2;
        input Real hill_i3;
        output Real hill_o1 "result";

      protected
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
        Promoter_o1 := Freg(Freg_i1=Promoter_i1, Freg_i2=Promoter_i2)/(
          Freg(Freg_i1=Promoter_i1, Freg_i2=Promoter_i2) + 6804.7);
      end Promoter;

      function pSMAD "pSMAD function (id f8)"
        input Real pSMAD_i1;
        input Real pSMAD_i2;
        output Real pSMAD_o1 "result";
      algorithm
        pSMAD_o1 := 0.0583 + 1.949*hill(
            hill_i1=pSMAD_i2,
            hill_i2=1.4481,
            hill_i3=140.244)/(1 + 0.1285*pSTAT(pSTAT_i1=pSMAD_i1, pSTAT_i2=pSMAD_i2));
      end pSMAD;

      function pSTAT "pSTAT function (id f7)"
        input Real pSTAT_i1;
        input Real pSTAT_i2;
        output Real pSTAT_o1 "result";
      algorithm
        pSTAT_o1 := (0.5/0.1285)*(facSig(facSig_i1=pSTAT_i1, facSig_i2=
          pSTAT_i2) + sqrt((facSig(facSig_i1=pSTAT_i1, facSig_i2=
          pSTAT_i2))^2 + 4*facSig1(facSig1_i1=pSTAT_i1)));
      end pSTAT;
   end algorithms;

    model Serum
      Bodylight.Types.RealIO.MassFlowRateInput Liver_in "Fe input from liver" annotation (
          Placement(transformation(extent={{-120,72},{-92,100}}),
            iconTransformation(extent={{-120,72},{-92,100}})));
      Bodylight.Types.RealIO.MassFlowRateInput Spleen_in "Fe input from spleen" annotation
        (Placement(transformation(extent={{-132,50},{-92,90}}), iconTransformation(
              extent={{-120,42},{-92,70}})));
      Bodylight.Types.RealIO.MassFlowRateInput Duodenum_in "Fe input from duodenum" annotation (Placement(transformation(extent={{-132,50},{-92,90}}),
            iconTransformation(extent={{-120,12},{-92,40}})));
      Bodylight.Types.RealIO.MassFlowRateInput OtherOrgans_in "Fe input from other organs" annotation (Placement(transformation(extent={{-132,50},{-92,90}}),
            iconTransformation(extent={{-120,-18},{-92,10}})));

      Bodylight.Types.RealIO.MassFlowRateInput Transfusion "transfusion rate"
        annotation (Placement(transformation(extent={{-132,50},{-92,90}}),
            iconTransformation(extent={{-118,-68},{-90,-40}})));
      Bodylight.Types.RealIO.MassFlowRateInput Bleeding "bleeding rate" annotation
        (Placement(transformation(extent={{-132,50},{-92,90}}), iconTransformation(
          extent={{-118,-100},{-90,-72}})));

      Bodylight.Types.RealIO.MassFlowRateOutput Liver_out "Fe output to liver" annotation (
          Placement(transformation(extent={{96,76},{116,96}}), iconTransformation(
              extent={{96,76},{116,96}})));
      Bodylight.Types.RealIO.MassFlowRateOutput BoneMarrow_out "Fe output to bone marrow" annotation (Placement(transformation(extent={{96,70},{116,90}}),
            iconTransformation(extent={{96,50},{116,70}})));
      Bodylight.Types.RealIO.MassFlowRateOutput Duodenum_out "Fe output to bone duodenum" annotation (Placement(transformation(extent={{96,70},{116,90}}),
            iconTransformation(extent={{96,24},{116,44}})));
      Bodylight.Types.RealIO.MassFlowRateOutput OtherOrgans_out "Fe output to other organs" annotation (Placement(transformation(extent={{96,70},{116,90}}),
          iconTransformation(extent={{96,-2},{116,18}})));

      Bodylight.Types.RealIO.MassOutput Fe(
        start = 1.51304 * 1e-9,
        displayUnit = "ug") "Fe amount in serum" annotation (Placement(transformation(extent={{92,74},{112,
          94}}), iconTransformation(extent={{92,-86},{112,-66}})));

      Bodylight.Types.MassFlowRate Fe_input(
        displayUnit = "ug.h-1") "Fe ammount in serum: input";
      Bodylight.Types.MassFlowRate Fe_output(
        displayUnit = "ug.h-1") "Fe ammount in serum: output";

      parameter Bodylight.Types.Mass th(
        displayUnit = "ug") = 2.6870 * 1e-9 "Threshold serum iron value (k26), 2.08/3.00";
      parameter Bodylight.Types.Frequency v_liv_1(
        displayUnit = "h-1") = 3.9607 / (60 * 60) "Low liver iron uptake (k4), 2.78/9.81";
      parameter Bodylight.Types.Frequency v_liv_2(
        displayUnit="h-1") = 14.3810 / (60 * 60) "High liver iron uptake (k45)"; //!!! Inconsistecy xml model vs. paper !!!, xml model value adopted;
      parameter Bodylight.Types.Frequency v_bm(
        displayUnit="h-1") = 2.8256 / (60 * 60) "Bone marrow uptake rate (k5), 2.66/4.16";
      parameter Bodylight.Types.Frequency v_duo(
        displayUnit="h-1") = 0.70485 / (60 * 60) "Duodenal uptake rate from blood (k6), 0.6/1.3";
      parameter Bodylight.Types.Frequency v_res(
        displayUnit="h-1") = 6.3206 / (60 * 60) "Other organs uptake rate (k34), 5.4/10.0";

    initial equation

      der(Fe) = 0; // eq 16.

    equation

      Liver_out = v_liv_1 * min(Fe, th) + v_liv_2 * max(Fe - th, 0);

      BoneMarrow_out = v_bm * Fe;

      Duodenum_out = v_duo * Fe;

      OtherOrgans_out = v_res * Fe;

      Fe_input = Liver_in + Spleen_in + Duodenum_in + OtherOrgans_in + Transfusion;
      Fe_output = Liver_out + BoneMarrow_out + Duodenum_out + OtherOrgans_out + Bleeding;

      der(Fe) = Fe_input - Fe_output;

      annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
            Rectangle(
              extent={{-100,100},{100,-100}},
              lineColor={28,108,200},
              fillColor={255,255,0},
              fillPattern=FillPattern.Solid),
            Text(
              extent={{-48,8},{72,-58}},
              textColor={28,108,200},
              textString="%name"),
            Text(
              extent={{-100,94},{-48,78}},
              textColor={28,108,200},
              textString="Liver"),
            Text(
              extent={{-96,62},{-40,46}},
              textColor={28,108,200},
              textString="Spleen"),
            Text(
              extent={{-88,36},{-26,12}},
              textColor={28,108,200},
              textString="Duodenum"),
            Text(
              extent={{-104,8},{-36,-20}},
              textColor={28,108,200},
              textString="Other
organs"),   Text(
              extent={{46,90},{96,74}},
              textColor={28,108,200},
              textString="Liver"),
            Text(
              extent={{34,74},{106,46}},
              textColor={28,108,200},
              textString="Bone
marrow"),   Text(
              extent={{26,42},{98,28}},
              textColor={28,108,200},
              textString="Duodenum"),
            Text(
              extent={{32,22},{104,-8}},
              textColor={28,108,200},
              textString="Other
organs"),   Text(
              extent={{-88,-44},{-28,-64}},
              textColor={28,108,200},
              textString="transfusion"),
            Text(
              extent={{-92,-78},{-38,-92}},
              textColor={28,108,200},
              textString="bleeding"),
            Text(
              extent={{54,-68},{102,-82}},
              textColor={28,108,200},
              textString="Fe")}),
                Diagram(coordinateSystem(preserveAspectRatio=false)));
    end Serum;

    model Test_Serum
                     extends Modelica.Icons.Example;
      Serum serum
        annotation (Placement(transformation(extent={{-34,-18},{28,40}})));
      Bodylight.Types.Constants.MassFlowRateConst Fe_liver(k(displayUnit=
              "ug.h-1") = 1.6646425e-12)
        annotation (Placement(transformation(extent={{-96,54},{-82,66}})));
      Bodylight.Types.Constants.MassFlowRateConst Fe_spleen(k(displayUnit=
              "ug.h-1") = 1.1875712777778e-12)
        annotation (Placement(transformation(extent={{-96,36},{-82,48}})));
      Bodylight.Types.Constants.MassFlowRateConst Fe_duodenum(k(displayUnit=
              "ug.h-1") = 7.2935111111111e-13)
        annotation (Placement(transformation(extent={{-96,18},{-82,30}})));
      Bodylight.Types.Constants.MassFlowRateConst Fe_otherorgans(k(displayUnit=
              "ug.h-1") = 2.223375e-12)
        annotation (Placement(transformation(extent={{-96,2},{-82,14}})));
      Bodylight.Types.Constants.MassFlowRateConst transfusion(k(displayUnit=
              "ug.h-1") = 0)
        annotation (Placement(transformation(extent={{-96,-16},{-82,-4}})));
      Bodylight.Types.Constants.MassFlowRateConst bleeding(k(displayUnit=
              "ug.h-1") = 0)
        annotation (Placement(transformation(extent={{-96,-34},{-82,-22}})));
    equation
      connect(Fe_spleen.y, serum.Spleen_in) annotation (Line(points={{-80.25,42},
              {-46,42},{-46,27.24},{-35.86,27.24}}, color={0,0,127}));
      connect(Fe_duodenum.y, serum.Duodenum_in) annotation (Line(points={{
              -80.25,24},{-46,24},{-46,18.54},{-35.86,18.54}}, color={0,0,127}));
      connect(Fe_otherorgans.y, serum.OtherOrgans_in) annotation (Line(points={
              {-80.25,8},{-44,8},{-44,9.84},{-35.86,9.84}}, color={0,0,127}));
      connect(Fe_liver.y, serum.Liver_in) annotation (Line(points={{-80.25,60},
              {-44,60},{-44,35.94},{-35.86,35.94}}, color={0,0,127}));
      connect(transfusion.y, serum.Transfusion) annotation (Line(points={{
              -80.25,-10},{-46,-10},{-46,-4.66},{-35.24,-4.66}}, color={0,0,127}));
      connect(bleeding.y, serum.Bleeding) annotation (Line(points={{-80.25,-28},
              {-44,-28},{-44,-13.94},{-35.24,-13.94}}, color={0,0,127}));
      annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
            coordinateSystem(preserveAspectRatio=false)));
    end Test_Serum;

    model BoneMarrow
      Bodylight.Types.RealIO.MassFlowRateInput Serum_in "Fe input from serum" annotation (
          Placement(transformation(extent={{-120,-14},{-92,14}}),
            iconTransformation(extent={{-120,-14},{-92,14}})));
      Bodylight.Types.RealIO.MassFlowRateOutput RBC_out "Fe output to RBC"
        annotation (Placement(transformation(extent={{96,76},{116,96}}),
            iconTransformation(extent={{96,10},{116,30}})));
      Bodylight.Types.RealIO.MassFlowRateOutput Spleen_out
        "Fe output to spleen"
        annotation (Placement(transformation(extent={{96,76},{116,96}}),
            iconTransformation(extent={{96,-30},{116,-10}})));

      Bodylight.Types.Mass Fe(
        start = 63.1596 * 1e-9,
        displayUnit = "ug") "Fe ammount in bones (s3)";
      Bodylight.Types.MassFlowRate Fe_input(
        displayUnit = "ug.h-1") "Fe ammount in bones: input";
      Bodylight.Types.MassFlowRate Fe_output(
        displayUnit = "ug.h-1") "Fe ammount in bones: output";

      parameter Bodylight.Types.Frequency v_RBC(
        displayUnit="h-1") = 0.058008 / (60 * 60) "RBC uptake rate (k8), 0.055/0.075"; // !!! Problem u rovnice (14) - doplnit prvni term na prave strane
      parameter Bodylight.Types.Frequency v_spl_2(
        displayUnit="h-1") = 0.0096817 / (60 * 60) "Spleen iron uptake rate from bones (k7), 0.008/0.019";

    initial equation

      der(Fe) = 0; // eq 12.

    equation

      RBC_out = v_RBC * Fe;
      Spleen_out = v_spl_2 * Fe;

      Fe_input = Serum_in;
      Fe_output = RBC_out + Spleen_out;

      der(Fe) = Fe_input - Fe_output;

      annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
            Rectangle(
              extent={{-100,100},{100,-100}},
              lineColor={28,108,200},
              fillColor={255,255,0},
              fillPattern=FillPattern.Solid),
            Text(
              extent={{-62,-24},{58,-90}},
              textColor={28,108,200},
              textString="%name"),
            Text(
              extent={{-92,8},{-40,-8}},
              textColor={28,108,200},
              textString="Serum"),
            Text(
              extent={{46,28},{98,12}},
              textColor={28,108,200},
              textString="RBC"),
            Text(
              extent={{42,-12},{94,-28}},
              textColor={28,108,200},
              textString="Spleen")}), Diagram(coordinateSystem(preserveAspectRatio=false)));
    end BoneMarrow;

    model Test_BoneMarrow
                          extends Modelica.Icons.Example
      annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
            coordinateSystem(preserveAspectRatio=false)));
      BoneMarrow boneMarrow
        annotation (Placement(transformation(extent={{-28,-28},{30,32}})));
      Bodylight.Types.Constants.MassFlowRateConst Fe_serum(k(displayUnit=
              "ug.h-1") = 1.1875712777778e-12)
        annotation (Placement(transformation(extent={{-92,-4},{-78,8}})));
    equation
      connect(Fe_serum.y, boneMarrow.Serum_in)
        annotation (Line(points={{-76.25,2},{-29.74,2}}, color={0,0,127}));
    end Test_BoneMarrow;

    model RBC
      Bodylight.Types.RealIO.MassFlowRateInput BoneMarrow_in
        "Fe input from bone marrow"                                                      annotation (Placement(transformation(extent={{-120,
                72},{-92,100}}), iconTransformation(extent={{-120,58},{-92,86}})));
      Bodylight.Types.RealIO.MassFlowRateInput Transfusion "transfusion rate" annotation (Placement(transformation(extent={{-132,50},{-92,90}}),
            iconTransformation(extent={{-120,-52},{-92,-24}})));
      Bodylight.Types.RealIO.MassFlowRateInput Bleeding "bleeding rate" annotation
        (Placement(transformation(extent={{-132,50},{-92,90}}), iconTransformation(
          extent={{-120,-90},{-92,-62}})));
      Bodylight.Types.RealIO.MassFlowRateOutput Spleen_out "Fe output to spleen" annotation (Placement(transformation(extent={{96,76},{116,96}}),
            iconTransformation(extent={{94,-10},{114,10}})));

      Bodylight.Types.Mass Fe(
        start = 63.1596 * 1e-9,
        displayUnit = "ug") "Fe ammount in bones (s3)";
      Bodylight.Types.MassFlowRate Fe_input(
        displayUnit = "ug.h-1") "Fe ammount in bones: input";
      Bodylight.Types.MassFlowRate Fe_output(
        displayUnit = "ug.h-1") "Fe ammount in bones: output";

      parameter Bodylight.Types.Frequency v_spl_1(
        displayUnit = "h-1") = 0.0036035 / (60 * 60) "Spleen iron uptake rate from RBC (k10), 0.002 = value for trace exp., 0.004/0.005";

    initial equation

      der(Fe) = 0; // eq 13.

    equation

      Spleen_out = v_spl_1 * Fe;

      Fe_input = BoneMarrow_in + Transfusion;
      Fe_output = Spleen_out + Bleeding;

      der(Fe) = Fe_input - Fe_output;

      annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
            Rectangle(
              extent={{-100,100},{100,-100}},
              lineColor={28,108,200},
              fillColor={255,255,0},
              fillPattern=FillPattern.Solid),
            Text(
              extent={{-60,70},{60,4}},
              textColor={28,108,200},
              textString="%name"),
            Text(
              extent={{-88,80},{-18,62}},
              textColor={28,108,200},
              textString="Bone marrow"),
            Text(
              extent={{-86,-28},{-26,-48}},
              textColor={28,108,200},
              textString="transfusion"),
            Text(
              extent={{-88,-68},{-34,-82}},
              textColor={28,108,200},
              textString="bleeding"),
            Text(
              extent={{40,8},{92,-8}},
              textColor={28,108,200},
              textString="Spleen")}), Diagram(coordinateSystem(preserveAspectRatio=false)));
    end RBC;

    model Test_RBC
                   extends Modelica.Icons.Example
      annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
            coordinateSystem(preserveAspectRatio=false)));
      RBC rBC annotation (Placement(transformation(extent={{-34,-20},{26,40}})));
      Bodylight.Types.Constants.MassFlowRateConst transfusion(k(displayUnit=
              "ug.h-1") = 0)
        annotation (Placement(transformation(extent={{-98,-6},{-84,6}})));
      Bodylight.Types.Constants.MassFlowRateConst bleeding(k(displayUnit=
              "ug.h-1") = 0)
        annotation (Placement(transformation(extent={{-98,-22},{-84,-10}})));
      Bodylight.Types.Constants.MassFlowRateConst Fe_bm(k(displayUnit="ug.h-1")
           = 1.0177121666667e-12)
        annotation (Placement(transformation(extent={{-98,26},{-84,38}})));
    equation
      connect(transfusion.y, rBC.Transfusion) annotation (Line(points={{-82.25,
              0},{-30,0},{-30,-1.4},{-35.8,-1.4}}, color={0,0,127}));
      connect(bleeding.y, rBC.Bleeding) annotation (Line(points={{-82.25,-16},{
              -40,-16},{-40,-12.8},{-35.8,-12.8}}, color={0,0,127}));
      connect(Fe_bm.y, rBC.BoneMarrow_in) annotation (Line(points={{-82.25,32},
              {-32.025,32},{-32.025,31.6},{-35.8,31.6}}, color={0,0,127}));
    end Test_RBC;

    model Spleen
      Bodylight.Types.RealIO.MassFlowRateInput BoneMarrow_in
        "Fe input from bone marrow"                                                      annotation (Placement(transformation(extent={{-120,
                72},{-92,100}}), iconTransformation(extent={{-120,72},{-92,100}})));
      Bodylight.Types.RealIO.MassFlowRateInput RBC_in "Fe input from RBC" annotation (Placement(transformation(extent={{-120,72},{-92,100}}),
            iconTransformation(extent={{-120,42},{-92,70}})));
      Bodylight.Types.RealIO.MassFlowRateOutput Serum_out "Fe output to serum" annotation (Placement(transformation(extent={{98,78},{118,98}}),
          iconTransformation(extent={{98,78},{118,98}})));

      Bodylight.Types.RealIO.MassInput hep "Hepcidin" annotation (Placement(
            transformation(extent={{-120,72},{-92,100}}), iconTransformation(extent
              ={{-120,-70},{-92,-42}})));
      Bodylight.Types.RealIO.MassInput il6 "IL6" annotation (Placement(
            transformation(extent={{-120,72},{-92,100}}), iconTransformation(extent
              ={{-120,-100},{-92,-72}})));

      Bodylight.Types.Mass Fe(
        start = 17.7295 * 1e-9,
        displayUnit = "ug") "Fe amount in spleen (s5)";
      Bodylight.Types.Mass FPN(
        start = 1.00049 * 1e-9,
        displayUnit = "ug") "FPN amount in spleen (s14)";
      Bodylight.Types.Mass FPN_mRNA(
        start = 0.922424 * 1e-9,
        displayUnit = "ug") "FPN mRNA amount in spleen (s9)";

      parameter Bodylight.Types.Mass Fe_max(
        displayUnit = "ug") = 88.216 * 1e-9 "Threshold value spleen iron export (k44), 57.0/95.0";
      parameter BodylightExtension.Types.MassSpecificRate u = 0.24102 * 1e9 / (60 * 60) "Spleen export rate (k3), 0.21/0.36";

      Bodylight.Types.MassFlowRate FPN_synthesis "FPN synthesis";
      Bodylight.Types.MassFlowRate FPN_degradation "FPN degradation";
      parameter Bodylight.Types.Frequency k_synth(
        displayUnit = "h-1") = 0.022722 / (60 * 60) "Fpnspl synthesis rate (k30), 0.015/0.027";
      parameter Bodylight.Types.Frequency k_deg(
        displayUnit="h-1") = 0.055363 * 0.054621 / (60 * 60) "Fpnspl degradation rate (k12*k18, k18 = 0.054621), 0.0007/0.0038";
      parameter BodylightExtension.Types.ReverseMass k_Fe = 0.014027 * 1e9 "Constant Fpn_spl production (k27), 0.005/0.028";
      parameter BodylightExtension.Types.ReverseMass k_hep = 2.5743 * 4.4694 * 1e9 "Constant Fpnspl degradation (k13*k19, k19 = 4.4694), 9.2/73.7";

      Bodylight.Types.MassFlowRate FPN_mRNA_synthesis "FPN synthesis";
      Bodylight.Types.MassFlowRate FPN_mRNA_degradation "FPN degradation";
      parameter Bodylight.Types.MassFlowRate K_mRNA = 1 / (60 * 60) * 1e-9;
      parameter Real K1_mRNA = 30.66 * 1.0867 "Constant Fpn_spl_mRNA production (k36*K_liv_1, k36 = 1.086700), 33.1/34.5";
      parameter Bodylight.Types.Mass K2_mRNA(
        displayUnit = "ug") = 0.0012836 * 1e-9 "Constant FpnmRNA production (k41), 3.0e-4/2.0e-3";
      parameter Bodylight.Types.Frequency k_mRNA_deg(
        displayUnit = "h-1") = 1.0841 / (60 * 60) "FpnmRNA degradation rate (k25), 1.02/1.20";

    initial equation

      der(Fe) = 0;
      der(FPN) = 0;
      der(FPN_mRNA) = 0;

    equation

      Serum_out = u * min(Fe, Fe_max) * FPN;

      der(Fe) = RBC_in + BoneMarrow_in - Serum_out;

      FPN_synthesis = k_synth * (1 + k_Fe * Fe) * FPN_mRNA;
      FPN_degradation = k_deg * (1 + k_hep * hep) * FPN;

      der(FPN) = FPN_synthesis - FPN_degradation;

      FPN_mRNA_synthesis = K_mRNA / (1 + (K1_mRNA * il6 / (K2_mRNA + il6)));
      FPN_mRNA_degradation = k_mRNA_deg * FPN_mRNA;

      der(FPN_mRNA) = FPN_mRNA_synthesis - FPN_mRNA_degradation;

      annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
              Rectangle(
              extent={{-100,100},{100,-100}},
              lineColor={28,108,200},
              fillColor={255,255,0},
              fillPattern=FillPattern.Solid),
            Text(
              extent={{-88,94},{-18,76}},
              textColor={28,108,200},
              textString="Bone marrow"),
            Text(
              extent={{-106,62},{-40,48}},
              textColor={28,108,200},
              textString="RBC"),
            Text(
              extent={{50,94},{98,80}},
              textColor={28,108,200},
              textString="Serum"),
            Text(
              extent={{-58,38},{62,-28}},
              textColor={28,108,200},
              textString="%name"),
            Text(
              extent={{-96,-50},{-30,-64}},
              textColor={28,108,200},
              textString="Hepcidin"),
            Text(
              extent={{-110,-80},{-44,-94}},
              textColor={28,108,200},
              textString="IL6")}),              Diagram(coordinateSystem(
              preserveAspectRatio=false)));
    end Spleen;

    model Test_Spleen
                      extends Modelica.Icons.Example
      annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
            coordinateSystem(preserveAspectRatio=false)));
      Spleen spleen
        annotation (Placement(transformation(extent={{-20,-20},{20,20}})));
      Bodylight.Types.Constants.MassFlowRateConst Fe_bm(k(displayUnit="ug.h-1")
           = 1.6985906388889e-13)
        annotation (Placement(transformation(extent={{-98,10},{-84,22}})));
      Bodylight.Types.Constants.MassFlowRateConst Fe_RBC(k(displayUnit="ug.h-1")
           = 1.0177122222222e-12)
        annotation (Placement(transformation(extent={{-98,-6},{-84,6}})));
      Bodylight.Types.Constants.MassConst hep(k(displayUnit="ug") =
          6.6493815e-10)
        annotation (Placement(transformation(extent={{-96,-20},{-88,-12}})));
      Bodylight.Types.Constants.MassConst il6(k(displayUnit="ug") = 0)
        annotation (Placement(transformation(extent={{-96,-32},{-88,-24}})));
    equation
      connect(Fe_bm.y, spleen.BoneMarrow_in) annotation (Line(points={{-82.25,
              16},{-24,16},{-24,17.2},{-21.2,17.2}}, color={0,0,127}));
      connect(Fe_RBC.y, spleen.RBC_in) annotation (Line(points={{-82.25,0},{-30,
              0},{-30,11.2},{-21.2,11.2}}, color={0,0,127}));
      connect(hep.y, spleen.hep) annotation (Line(points={{-87,-16},{-32,-16},{
              -32,-11.2},{-21.2,-11.2}}, color={0,0,127}));
      connect(il6.y, spleen.il6) annotation (Line(points={{-87,-28},{-30,-28},{
              -30,-17.2},{-21.2,-17.2}}, color={0,0,127}));
    end Test_Spleen;

    model Liver
      Bodylight.Types.RealIO.MassFlowRateInput Serum_in "Fe input from serum"
        annotation (Placement(transformation(extent={{-120,72},{-92,100}}),
            iconTransformation(extent={{-120,72},{-92,100}})));
      Bodylight.Types.RealIO.MassFlowRateOutput Serum_out "Fe output to serum" annotation (Placement(transformation(extent={{98,78},{118,98}}),
          iconTransformation(extent={{98,78},{118,98}})));

      Bodylight.Types.RealIO.MassInput hep "Hepcidin" annotation (Placement(
            transformation(extent={{-120,72},{-92,100}}), iconTransformation(extent
              ={{-120,-70},{-92,-42}})));
      Bodylight.Types.RealIO.MassInput il6 "IL6" annotation (Placement(
            transformation(extent={{-120,72},{-92,100}}), iconTransformation(extent
              ={{-120,-100},{-92,-72}})));

      Bodylight.Types.RealIO.MassOutput Fe(
        start = 76.9562 * 1e-9,
        displayUnit="ug") "Fe amount in liver" annotation (Placement(transformation(extent={{92,74},{112,
          94}}), iconTransformation(extent={{92,-80},{112,-60}})));

      Bodylight.Types.Mass FPN(
        start = 1.00036 * 1e-9,
        displayUnit = "ug") "FPN amount in liver (s12)";
      Bodylight.Types.Mass FPN_mRNA(
        start = 0.922424 * 1e-9,
        displayUnit = "ug") "FPN mRNA amount in liver (s7)";

      parameter Bodylight.Types.Mass Fe_max(
        displayUnit = "ug") = 119.55 * 1e-9 "Threshold value liver iron export (k43), 100.0/159.0";
      parameter BodylightExtension.Types.MassSpecificRate u = 0.077844 * 1e9 / (60 * 60) "Liver iron export rate (k1), 0.05/0.19";

      Bodylight.Types.MassFlowRate FPN_synthesis "FPN synthesis";
      Bodylight.Types.MassFlowRate FPN_degradation "FPN degradation";
      parameter Bodylight.Types.Frequency k_synth(
        displayUnit = "h-1") = 0.1297 / (60 * 60) "Fpnliv synthesis rate (k28), 0.07/0.14";
      parameter Bodylight.Types.Frequency k_deg(
        displayUnit="h-1") = 0.055363 / (60 * 60) "Fpnliv degradation rate (k12), 0.01/0.06";
      parameter BodylightExtension.Types.ReverseMass k_Fe = 0.0033177 * 1e9 "Constant Fpn_liv production (k17), 0.002/0.006";
      parameter BodylightExtension.Types.ReverseMass k_hep = 2.5743 * 1e9 "Constant Fpnliv degradation (k13), 2.11/12.93";

      Bodylight.Types.MassFlowRate FPN_mRNA_synthesis "FPN synthesis";
      Bodylight.Types.MassFlowRate FPN_mRNA_degradation "FPN degradation";
      parameter Bodylight.Types.MassFlowRate K_mRNA = 1 / (60 * 60) * 1e-9;
      parameter Real K1_mRNA = 30.66 "Constant Fpn_liv_mRNA production (k24), 28.0/32.5";
      parameter Bodylight.Types.Mass K2_mRNA(
        displayUnit = "ug") = 0.0012836 * 1e-9 "Constant FpnmRNA production (k41), 3.0e-4/2.0e-3";
      parameter Bodylight.Types.Frequency k_mRNA_deg(
        displayUnit = "h-1") = 1.0841 / (60 * 60) "FpnmRNA degradation rate (k25), 1.02/1.20";

    initial equation

      der(Fe) = 0;
      der(FPN) = 0;
      der(FPN_mRNA) = 0;

    equation

      Serum_out = u * min(Fe, Fe_max) * FPN;

      der(Fe) = Serum_in - Serum_out;

      FPN_synthesis = k_synth * (1 + k_Fe * Fe) * FPN_mRNA;
      FPN_degradation = k_deg * (1 + k_hep * hep) * FPN;

      der(FPN) = FPN_synthesis - FPN_degradation;

      FPN_mRNA_synthesis = K_mRNA / (1 + (K1_mRNA * il6 / (K2_mRNA + il6)));
      FPN_mRNA_degradation = k_mRNA_deg * FPN_mRNA;

      der(FPN_mRNA) = FPN_mRNA_synthesis - FPN_mRNA_degradation;

      annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
              Rectangle(
              extent={{-100,100},{100,-100}},
              lineColor={28,108,200},
              fillColor={255,255,0},
              fillPattern=FillPattern.Solid),
            Text(
              extent={{-94,92},{-40,78}},
              textColor={28,108,200},
              textString="Serum"),
            Text(
              extent={{50,94},{98,80}},
              textColor={28,108,200},
              textString="Serum"),
            Text(
              extent={{-58,38},{62,-28}},
              textColor={28,108,200},
              textString="%name"),
            Text(
              extent={{-96,-50},{-30,-64}},
              textColor={28,108,200},
              textString="Hepcidin"),
            Text(
              extent={{-110,-80},{-44,-94}},
              textColor={28,108,200},
              textString="IL6"),
            Text(
              extent={{54,-62},{102,-76}},
              textColor={28,108,200},
              textString="Fe")}),               Diagram(coordinateSystem(
              preserveAspectRatio=false)));
    end Liver;

    model Test_Liver
                     extends Modelica.Icons.Example
      annotation (

                  Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
            coordinateSystem(preserveAspectRatio=false)));
      Liver liver
        annotation (Placement(transformation(extent={{-20,-20},{20,20}})));
      Bodylight.Types.Constants.MassFlowRateConst Fe_serum(k(displayUnit=
              "ug.h-1") = 1.6646425e-12)
        annotation (Placement(transformation(extent={{-90,10},{-76,22}})));
      Bodylight.Types.Constants.MassConst hep(k(displayUnit="ug") =
          6.6493815e-10)
        annotation (Placement(transformation(extent={{-92,-14},{-84,-6}})));
      Bodylight.Types.Constants.MassConst il6(k(displayUnit="ug") =
          6.6320815e-41)
        annotation (Placement(transformation(extent={{-92,-28},{-84,-20}})));
    equation
      connect(Fe_serum.y, liver.Serum_in) annotation (Line(points={{-74.25,16},
              {-30,16},{-30,17.2},{-21.2,17.2}}, color={0,0,127}));
      connect(hep.y, liver.hep) annotation (Line(points={{-83,-10},{-30,-10},{
              -30,-11.2},{-21.2,-11.2}}, color={0,0,127}));
      connect(il6.y, liver.il6) annotation (Line(points={{-83,-24},{-32,-24},{
              -32,-17.2},{-21.2,-17.2}}, color={0,0,127}));
    end Test_Liver;

    model Duodenum
      Bodylight.Types.RealIO.MassFlowRateInput Serum_in "Fe input from serum"
        annotation (Placement(transformation(extent={{-120,72},{-92,100}}),
            iconTransformation(extent={{-120,72},{-92,100}})));
      Bodylight.Types.RealIO.MassFlowRateInput Food_in "Fe input from food"
        annotation (Placement(transformation(extent={{-120,72},{-92,100}}),
            iconTransformation(extent={{-120,42},{-92,70}})));
      Bodylight.Types.RealIO.MassFlowRateOutput Serum_out "Fe output to serum" annotation (Placement(transformation(extent={{98,78},{118,98}}),
          iconTransformation(extent={{98,78},{118,98}})));

      Bodylight.Types.RealIO.MassInput hep "Hepcidin" annotation (Placement(
            transformation(extent={{-120,72},{-92,100}}), iconTransformation(extent
              ={{-120,-70},{-92,-42}})));
      Bodylight.Types.RealIO.MassInput il6 "IL6" annotation (Placement(
            transformation(extent={{-120,72},{-92,100}}), iconTransformation(extent
              ={{-120,-100},{-92,-72}})));

      Bodylight.Types.Mass Fe(
        start = 2.97118 * 1e-9,
        displayUnit = "ug") "Fe amount in duodenum (s4)";
      Bodylight.Types.Mass FPN(
        start = 1.00017 * 1e-9,
        displayUnit = "ug") "Fpn amount in duodenum (s13)";
      Bodylight.Types.Mass FPN_mRNA(
        start = 0.922424 * 1e-9,
        displayUnit = "ug") "Fpn mRNA amount in duodenum (s8)";

      parameter BodylightExtension.Types.MassSpecificRate u = 0.88356 * 1e9 / (60 * 60) "Duodenal export rate (k2), 0.66/1.38";

      Bodylight.Types.MassFlowRate Fe_intake "Fe intake from food";
      parameter Bodylight.Types.MassFlowRate v_max = 9.86 / (60 * 60) * 1e-9 "Maximal duodenal uptake from food (SI. Eq.14)";
      parameter Bodylight.Types.MassFlowRate K_sat = 177.34 / (60 * 60) * 1e-9 "Saturation parameter duodenal uptake, corresponds to k11";

      Bodylight.Types.MassFlowRate Fe_absorption "Fe input to duodenum from intestines, prev:in_1";
      parameter Boolean unregulated_absorption = false;
      parameter Real malabsorption = 1 "iron malabsorption coefficient: <0;1>; 0 = no absorption, 1 = physiologic";
      parameter Bodylight.Types.Mass K_absorption(
        displayUnit = "ug") = 1 * 1e-9;

      Bodylight.Types.MassFlowRate Fe_loss "Fe output from duodenum to intestines (Fe loss), prev:out_2";
      parameter Bodylight.Types.Frequency v_lost = 0.091919 / (60 * 60) "Iron lost rate duodenum (k48), 0.001/0.320";

      Bodylight.Types.MassFlowRate FPN_synthesis "FPN synthesis";
      Bodylight.Types.MassFlowRate FPN_degradation "FPN degradation";
      parameter Bodylight.Types.Frequency k_synth(
        displayUnit = "h-1") = 0.030299 / (60 * 60) "Fpnduo synthesis rate (k29), 0.01/0.25";
      parameter Bodylight.Types.Frequency k_deg(
        displayUnit="h-1") = 0.055363 * 0.3815 / (60 * 60) "Fpnduo degradation rate (k12*k15, k15 = 0.3815), 0.0056/0.147";
      parameter BodylightExtension.Types.ReverseMass k_Fe = 0.16007 * 1e9 "Constant Fpn_duo production (k46), 0.06/0.49";
      parameter BodylightExtension.Types.ReverseMass k_hep = 2.5743 * 0.55631 * 1e9 "Constant Fpnduo degradation (k13*k16, k16 = 0.55631), 0.78/4.16";

      Bodylight.Types.MassFlowRate FPN_mRNA_synthesis "FPN synthesis";
      Bodylight.Types.MassFlowRate FPN_mRNA_degradation "FPN degradation";
      parameter Bodylight.Types.MassFlowRate K_mRNA =  1 / (60 * 60) * 1e-9;
      parameter Real K1_mRNA = 30.66;
      parameter Bodylight.Types.Mass K2_mRNA(
        displayUnit = "ug") = 0.0012836 * 1e-9 "Constant FpnmRNA production (k41), 3.0e-4/2.0e-3";
      parameter Bodylight.Types.Frequency k_mRNA_deg(
        displayUnit = "h-1") = 1.0841 / (60 * 60) "FpnmRNA degradation rate (k25), 1.02/1.20";

    initial equation

      der(Fe) = 0;
      der(FPN) = 0;
      der(FPN_mRNA) = 0;

    equation

      Serum_out = u * Fe * FPN;
      Fe_loss = v_lost * Fe;

      Fe_intake = v_max * Food_in / (Food_in + K_sat);
      Fe_absorption =
        if unregulated_absorption
          then Fe_intake
        else
          Fe_intake * min(K_absorption / Fe, 1) * malabsorption;

      der(Fe) = Serum_in - Serum_out + Fe_absorption - Fe_loss;

      FPN_synthesis = k_synth * (1 + k_Fe * Fe) * FPN_mRNA;
      FPN_degradation = k_deg * (1 + k_hep * hep) * FPN;

      der(FPN) = FPN_synthesis - FPN_degradation;

      FPN_mRNA_synthesis = K_mRNA / (1 + (K1_mRNA * il6 / (K2_mRNA + il6)));
      FPN_mRNA_degradation = k_mRNA_deg * FPN_mRNA;

      der(FPN_mRNA) = FPN_mRNA_synthesis - FPN_mRNA_degradation;

      annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
              Rectangle(
              extent={{-100,100},{100,-100}},
              lineColor={28,108,200},
              fillColor={255,255,0},
              fillPattern=FillPattern.Solid),
            Text(
              extent={{-94,92},{-40,78}},
              textColor={28,108,200},
              textString="Serum"),
            Text(
              extent={{50,94},{98,80}},
              textColor={28,108,200},
              textString="Serum"),
            Text(
              extent={{-58,38},{62,-28}},
              textColor={28,108,200},
              textString="%name"),
            Text(
              extent={{-96,-50},{-30,-64}},
              textColor={28,108,200},
              textString="Hepcidin"),
            Text(
              extent={{-110,-80},{-44,-94}},
              textColor={28,108,200},
              textString="IL6"),
            Text(
              extent={{-94,62},{-40,48}},
              textColor={28,108,200},
              textString="Food")}),             Diagram(coordinateSystem(
              preserveAspectRatio=false)));
    end Duodenum;

    model Test_Duodenum
                        extends Modelica.Icons.Example
      annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
            coordinateSystem(preserveAspectRatio=false)));
      Duodenum duodenum
        annotation (Placement(transformation(extent={{-20,-20},{20,20}})));
      Bodylight.Types.Constants.MassFlowRateConst Fe_serum(k(displayUnit="ug.h-1")
           = 2.9624138888889e-13)
        annotation (Placement(transformation(extent={{-98,20},{-84,32}})));
      Bodylight.Types.Constants.MassConst hep(k(displayUnit="ug") = 6.6493815e-10)
        annotation (Placement(transformation(extent={{-98,-14},{-90,-6}})));
      Bodylight.Types.Constants.MassConst il6(k(displayUnit="ug") = 0)
        annotation (Placement(transformation(extent={{-98,-26},{-90,-18}})));
      Bodylight.Types.Constants.MassFlowRateConst Fe_food(k(displayUnit="ug.h-1")
           = 6.0731166666667e-11)
        annotation (Placement(transformation(extent={{-98,2},{-84,14}})));
    equation
      connect(Fe_serum.y, duodenum.Serum_in) annotation (Line(points={{-82.25,26},
              {-51.725,26},{-51.725,17.2},{-21.2,17.2}}, color={0,0,127}));
      connect(hep.y, duodenum.hep) annotation (Line(points={{-89,-10},{-30,-10},{-30,
              -11.2},{-21.2,-11.2}}, color={0,0,127}));
      connect(il6.y, duodenum.il6) annotation (Line(points={{-89,-22},{-32,-22},{-32,
              -17.2},{-21.2,-17.2}}, color={0,0,127}));
      connect(Fe_food.y, duodenum.Food_in) annotation (Line(points={{-82.25,8},{-32,
              8},{-32,11.2},{-21.2,11.2}}, color={0,0,127}));
    end Test_Duodenum;

    model OtherOrgans
      Bodylight.Types.RealIO.MassFlowRateInput Serum_in "Fe input from serum"
        annotation (Placement(transformation(extent={{-120,72},{-92,100}}),
            iconTransformation(extent={{-120,72},{-92,100}})));
      Bodylight.Types.RealIO.MassFlowRateOutput Serum_out "Fe output to serum" annotation (Placement(transformation(extent={{98,78},{118,98}}),
          iconTransformation(extent={{98,78},{118,98}})));

      Bodylight.Types.RealIO.MassInput hep "Hepcidin" annotation (Placement(
            transformation(extent={{-120,72},{-92,100}}), iconTransformation(extent
              ={{-120,-70},{-92,-42}})));
      Bodylight.Types.RealIO.MassInput il6 "IL6" annotation (Placement(
            transformation(extent={{-120,72},{-92,100}}), iconTransformation(extent
              ={{-120,-100},{-92,-72}})));

      Bodylight.Types.Mass Fe(
        start = 466.811 * 1e-9,
        displayUnit = "ug") "Fe amount in other organs (s10)";
      Bodylight.Types.Mass FPN(
        start = 1.00020 * 1e-9,
        displayUnit = "ug") "Fpn amount in other organs (s15)";
      Bodylight.Types.Mass FPN_mRNA(
        start = 0.922424 * 1e-9,
        displayUnit = "ug") "Fpn mRNA amount in other organs (s11)";

      parameter Bodylight.Types.Mass Fe_max(
        displayUnit = "ug") = 510.68 * 1e-9 "Limit value, iron lost rest (k49), 510/947";
      parameter BodylightExtension.Types.MassSpecificRate u = 0.017143 * 1e9 / (60 * 60) "Other organs export rate (k33), 0.014/0.030";

      Bodylight.Types.MassFlowRate Fe_loss;
      parameter Real loss_factor = 1 "iron loss factor, 1 = physiologic, >1 enhanced loss";
      parameter Bodylight.Types.Frequency u_lost = 0.0033401 / (60 * 60) "Iron lost rate rest (k9), 0.002/0.004";

      Bodylight.Types.MassFlowRate FPN_synthesis "FPN synthesis";
      Bodylight.Types.MassFlowRate FPN_degradation "FPN degradation";
      parameter Bodylight.Types.Frequency k_synth(
        displayUnit = "h-1") = 0.0050598 / (60 * 60) "Fpnres synthesis rate (k40), 0.004/0.108";
      parameter Bodylight.Types.Frequency k_deg(
        displayUnit="h-1") = 0.055363 * 0.52232 / (60 * 60) "Fpnres degradation rate (k12*k38, k38 = 0.52232), 0.025/0.129";
      parameter BodylightExtension.Types.ReverseMass k_Fe = 0.11376 * 1e9 "Constant Fpn_res production (k47), 0.004/0.152";
      parameter BodylightExtension.Types.ReverseMass k_hep = 2.5743 * 4.5163 * 1e9 "Constant Fpnres degradation (k13*k39, k39 = 4.5163), 3.7/38.1";

      Bodylight.Types.MassFlowRate FPN_mRNA_synthesis "FPN synthesis";
      Bodylight.Types.MassFlowRate FPN_mRNA_degradation "FPN degradation";
      parameter Bodylight.Types.MassFlowRate K_mRNA = 1 / (60 * 60) * 1e-9;
      parameter Real K1_mRNA = 30.66 * 0.36629 "Constant Fpn_res_mRNA production (k42*K_liv_1, k42 = 0.366290), 7.80/43.7";
      parameter Bodylight.Types.Mass K2_mRNA(
        displayUnit = "ug") = 0.0012836 * 1e-9 "Constant FpnmRNA production (k41), 3.0e-4/2.0e-3";
      parameter Bodylight.Types.Frequency k_mRNA_deg(
        displayUnit = "h-1") = 1.0841 / (60 * 60) "FpnmRNA degradation rate (k25), 1.02/1.20";

    initial equation

      der(Fe) = 0;
      der(FPN) = 0;
      der(FPN_mRNA) = 0;

    equation

      Serum_out = u * Fe * FPN;
      Fe_loss = u_lost * min(Fe, Fe_max) * loss_factor;

      der(Fe) = Serum_in - Serum_out - Fe_loss;

      FPN_synthesis = k_synth * (1 + k_Fe * Fe) * FPN_mRNA;
      FPN_degradation = k_deg * (1 + k_hep * hep) * FPN;

      der(FPN) = FPN_synthesis - FPN_degradation;

      FPN_mRNA_synthesis = K_mRNA / (1 + (K1_mRNA * il6 / (K2_mRNA + il6)));
      FPN_mRNA_degradation = k_mRNA_deg * FPN_mRNA;

      der(FPN_mRNA) = FPN_mRNA_synthesis - FPN_mRNA_degradation;

      annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
              Rectangle(
              extent={{-100,100},{100,-100}},
              lineColor={28,108,200},
              fillColor={255,255,0},
              fillPattern=FillPattern.Solid),
            Text(
              extent={{-94,92},{-40,78}},
              textColor={28,108,200},
              textString="Serum"),
            Text(
              extent={{50,94},{98,80}},
              textColor={28,108,200},
              textString="Serum"),
            Text(
              extent={{-58,38},{62,-28}},
              textColor={28,108,200},
              textString="%name"),
            Text(
              extent={{-96,-50},{-30,-64}},
              textColor={28,108,200},
              textString="Hepcidin"),
            Text(
              extent={{-110,-80},{-44,-94}},
              textColor={28,108,200},
              textString="IL6")}),              Diagram(coordinateSystem(
              preserveAspectRatio=false)));
    end OtherOrgans;

    model Test_OtherOrgans
      extends Modelica.Icons.Example
      annotation (

                  Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
            coordinateSystem(preserveAspectRatio=false)));
      OtherOrgans otherOrgans
        annotation (Placement(transformation(extent={{-20,-20},{40,40}})));
      Bodylight.Types.Constants.MassFlowRateConst Fe_serum(k(displayUnit=
              "ug.h-1") = 2.6564847222222e-12)
        annotation (Placement(transformation(extent={{-96,30},{-82,42}})));
      Bodylight.Types.Constants.MassConst hep(k(displayUnit="ug") =
          6.6493815e-10)
        annotation (Placement(transformation(extent={{-92,-10},{-84,-2}})));
      Bodylight.Types.Constants.MassConst il6(k(displayUnit="ug") = 0)
        annotation (Placement(transformation(extent={{-92,-22},{-84,-14}})));
    equation
      connect(Fe_serum.y, otherOrgans.Serum_in) annotation (Line(points={{
              -80.25,36},{-51.025,36},{-51.025,35.8},{-21.8,35.8}}, color={0,0,
              127}));
      connect(hep.y, otherOrgans.hep) annotation (Line(points={{-83,-6},{-21.8,
              -6},{-21.8,-6.8}}, color={0,0,127}));
      connect(il6.y, otherOrgans.il6) annotation (Line(points={{-83,-18},{-24,
              -18},{-24,-15.8},{-21.8,-15.8}}, color={0,0,127}));
    end Test_OtherOrgans;

    model Hepcidin
      Bodylight.Types.RealIO.MassInput il6 annotation (Placement(transformation(
              extent={{-120,60},{-80,100}}), iconTransformation(extent={{-120,60},{-80,
                100}})));
      Bodylight.Types.RealIO.MassInput Bmp6 annotation (Placement(transformation(
              extent={{-120,60},{-80,100}}), iconTransformation(extent={{-120,18},{-80,
                58}})));
      Bodylight.Types.RealIO.MassOutput hep(
        start = 0.664938 * 1e-9,
        displayUnit = "ug") "Hepcidin amount (s16)" annotation (Placement(transformation(
              extent={{96,70},{116,90}}), iconTransformation(extent={{96,70},{116,90}})));

      parameter Bodylight.Types.Frequency k_deg = 0.07 / (60 * 60) "Hepcidin degradation rate (k20), 0.067/0.070";
      parameter Real knockout = 1;

      parameter BodylightExtension.Types.ReverseMass K_mass = 1e9;
      parameter Bodylight.Types.Frequency K_freq = 1 / (60 * 60);

      Bodylight.Types.MassFlowRate hep_in "Production of hepcidin";
      Bodylight.Types.MassFlowRate hep_out "Degradation of hepcidin";

    initial equation

      der(hep) = 0; // eq. 1

    equation

      hep_in = algorithms.Promoter(il6 * K_mass, Bmp6 * K_mass) / K_mass * K_freq * knockout;

      hep_out = k_deg * hep;

      der(hep) = hep_in - hep_out; // eq. 1

      annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
            Rectangle(
              extent={{-100,100},{100,-100}},
              lineColor={28,108,200},
              fillColor={255,255,0},
              fillPattern=FillPattern.Solid),
            Text(
              extent={{-54,34},{66,-32}},
              textColor={28,108,200},
              textString="%name"),
            Text(
              extent={{-96,92},{-20,68}},
              textColor={28,108,200},
              textString="IL-6"),
            Text(
              extent={{-80,48},{-20,24}},
              textColor={28,108,200},
              textString="Bmp6"),
            Text(
              extent={{34,88},{92,70}},
              textColor={28,108,200},
              textString="Hepcidin")}), Diagram(coordinateSystem(
              preserveAspectRatio=false)));
    end Hepcidin;

    model Test_Hepcidin
      extends Modelica.Icons.Example;
      Hepcidin hepcidin
        annotation (Placement(transformation(extent={{-40,-40},{40,40}})));
      Bodylight.Types.Constants.MassConst il6(k(displayUnit="ug") = 0)
        annotation (Placement(transformation(extent={{-94,34},{-86,42}})));
      Bodylight.Types.Constants.MassConst Bmp6(k(displayUnit="ug") =
          1.5829111e-08)
        annotation (Placement(transformation(extent={{-94,16},{-86,24}})));
    equation
      connect(il6.y, hepcidin.il6) annotation (Line(points={{-85,38},{-54,38},{
              -54,32},{-40,32}}, color={0,0,127}));
      connect(Bmp6.y, hepcidin.Bmp6) annotation (Line(points={{-85,20},{-54,20},
              {-54,15.2},{-40,15.2}}, color={0,0,127}));
    end Test_Hepcidin;

    model IL6 "IL-6"
      Bodylight.Types.RealIO.MassInput LPS annotation (Placement(transformation(
              extent={{-120,60},{-80,100}}), iconTransformation(extent={{-120,60},{-80,
                100}})));
      Bodylight.Types.RealIO.MassOutput il6(
        start = 0*1e-9,
        displayUnit = "ug") "IL-6 amount (s22)"
                            annotation (Placement(transformation(
              extent={{94,74},{114,94}}), iconTransformation(extent={{94,74},{114,94}})));

      Bodylight.Types.Mass il6mRNA(
        start = 0 * 1e-9,
        displayUnit = "ug") "IL-6 mRNA ammount (s21)";

      Bodylight.Types.MassFlowRate il6mRNA_in "Production of IL-6 mRNA";
      Bodylight.Types.MassFlowRate il6mRNA_out "Degradation of IL-6 mRNA";
      parameter Bodylight.Types.Mass K_mRNA = 2.6e-6 * 1e-9 "Michaelis-Menten constant Il6mRNA synthesis, 2.6-e6/2.6e-6";
      parameter Bodylight.Types.Frequency k_deg_mRNA = 0.2814 / (60 * 60) "Il6mRNA degradation rate, 0.28/0.28";

      Bodylight.Types.MassFlowRate il6_in "Production of IL-6";
      Bodylight.Types.MassFlowRate il6_out "Degradation of IL-6";
      parameter Bodylight.Types.Frequency k_syn = 4.1067 * 157.4 / (60 * 60) "Il6 synthesis rate, 4.1067*k23, k23 = 157.4, 136/872";
      parameter Bodylight.Types.Frequency k_deg = 4.4465 / (60 * 60) "Il6 degradation rate, 4.45/4.45";
      parameter BodylightExtension.Types.ReverseMass K_mass = 1e9;

    initial equation

      der(il6mRNA) = 0; // eq 4.
      der(il6) = 0; // eq 5.

    equation

      il6mRNA_in = LPS / (LPS + K_mRNA);
      il6mRNA_out = k_deg_mRNA * il6mRNA;

      der(il6mRNA) = il6mRNA_in - il6mRNA_out;

      il6_in = k_syn * K_mass ^ 3 * il6mRNA ^ 4;
      il6_out = k_deg * il6;

      der(il6) = il6_in - il6_out;

      annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
              Rectangle(
              extent={{-100,100},{100,-100}},
              lineColor={28,108,200},
              fillColor={255,255,0},
              fillPattern=FillPattern.Solid),
            Text(
              extent={{-52,38},{68,-28}},
              textColor={28,108,200},
              textString="%name"),
            Text(
              extent={{28,94},{104,70}},
              textColor={28,108,200},
              textString="IL-6"),
            Text(
              extent={{-90,90},{-14,66}},
              textColor={28,108,200},
              textString="LPS")}),   Diagram(coordinateSystem(preserveAspectRatio=false)));
    end IL6;

    model Test_il6
      extends Modelica.Icons.Example;
      IL6 iL6 annotation (Placement(transformation(extent={{-40,-40},{40,40}})));
      Bodylight.Types.Constants.MassConst LPS(k(displayUnit="ug") = 0)
        annotation (Placement(transformation(extent={{-98,32},{-90,40}})));
    equation
      connect(LPS.y, iL6.LPS) annotation (Line(points={{-89,36},{-54,36},{-54,
              32},{-40,32}}, color={0,0,127}));
    end Test_il6;

    model LPS
      Bodylight.Types.RealIO.MassOutput LPS(
        start = 0 * 1e-9,
        displayUnit = "ug") "LPS (Lipopolysaccharide) amount (s23)"
                            annotation (Placement(transformation(
          extent={{92,74},{112,94}}), iconTransformation(extent={{92,74},{112,94}})));

      parameter Bodylight.Types.Frequency k_deg = 5.8560 / (60 * 60) "LPS degradation rate, 5.9/5.9";

    initial equation

      der(LPS) = 0;

    equation

        der(LPS) = -k_deg * LPS;

      annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
            Rectangle(
              extent={{-100,100},{100,-100}},
              lineColor={28,108,200},
              fillColor={255,255,0},
              fillPattern=FillPattern.Solid),
            Text(
              extent={{-54,38},{66,-28}},
              textColor={28,108,200},
              textString="%name"),
            Text(
              extent={{24,94},{100,70}},
              textColor={28,108,200},
              textString="LPS")}), Diagram(coordinateSystem(preserveAspectRatio=false)));
    end LPS;

    model Bmp6
      Bodylight.Types.RealIO.MassOutput Bmp6(
        start = 15.8291 * 1e-9,
        displayUnit="ug") "Bmp6 amount (s18)" annotation (Placement(transformation(extent={{92,76},{112,
                96}}), iconTransformation(extent={{92,76},{112,96}})));
      Bodylight.Types.RealIO.MassInput Fe_liv "Fe amount in liver" annotation (
          Placement(transformation(extent={{-120,60},{-80,100}}),
            iconTransformation(extent={{-120,60},{-80,100}})));
      Bodylight.Types.RealIO.MassInput Fe_ser "Fe amount in serum" annotation (
          Placement(transformation(extent={{-120,60},{-80,100}}),
          iconTransformation(extent={{-120,16},{-80,56}})));

      Bodylight.Types.MassFlowRate Bmp6_in "Bmp6 production rate";
      Bodylight.Types.MassFlowRate Bmp6_out "Bmp6 production rate";
      Bodylight.Types.Frequency v_max "Bmp6 maximal synthesis rate";
      parameter Bodylight.Types.Mass Tf = 1000.0 * 1e-9 "Paremeter determining the maximal amount of iron that can be bound to transferrin (k32)";
      parameter Bodylight.Types.Mass K = 19.65 * 1e-9 "Michaelis-Menten constant Bmp6 synthesis (k37), 16.5/55.7";
      parameter BodylightExtension.Types.MassSpecificRate k_max = 1.6015 * 1e9 / (60 * 60);
      parameter Bodylight.Types.Frequency k_deg = 2.3962 / (60 * 60) "Bmp6 degradation rate (k22), 1.0/9.5";

    initial equation

      der(Bmp6) = 0;

    equation

      v_max = k_max * K;
      Bmp6_in = v_max * Fe_liv / (K + Fe_liv) * min(Fe_ser, Tf);
      Bmp6_out = k_deg * Bmp6;

      der(Bmp6) = Bmp6_in - Bmp6_out;

      annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
            Rectangle(
              extent={{-100,100},{100,-100}},
              lineColor={28,108,200},
              fillColor={255,255,0},
              fillPattern=FillPattern.Solid),
            Text(
              extent={{-74,90},{-4,68}},
              textColor={28,108,200},
              textString="Fe Liver"),
            Text(
              extent={{-76,54},{4,18}},
              textColor={28,108,200},
              textString="Fe Serum"),
            Text(
              extent={{28,96},{88,72}},
              textColor={28,108,200},
              textString="Bmp6"),
            Text(
              extent={{-56,20},{64,-46}},
              textColor={28,108,200},
              textString="%name")}), Diagram(coordinateSystem(preserveAspectRatio=false)));
    end Bmp6;

    model Test_Bmp6
      extends Modelica.Icons.Example;
      Bmp6 bmp6
        annotation (Placement(transformation(extent={{-40,-40},{40,40}})));
      Bodylight.Types.Constants.MassConst Fe_liv(k(displayUnit="ug") =
          7.695616e-08)
        annotation (Placement(transformation(extent={{-98,32},{-86,42}})));
      Bodylight.Types.Constants.MassConst Fe_ser(k(displayUnit="ug") =
          1.5130439e-09)
        annotation (Placement(transformation(extent={{-98,4},{-86,14}})));
    equation
      connect(Fe_liv.y, bmp6.Fe_liv) annotation (Line(points={{-84.5,37},{-54,
              37},{-54,32},{-40,32}}, color={0,0,127}));
      connect(Fe_ser.y, bmp6.Fe_ser) annotation (Line(points={{-84.5,9},{-54,9},
              {-54,14.4},{-40,14.4}}, color={0,0,127}));
    end Test_Bmp6;

    model Test_Full
      extends Modelica.Icons.Example;
      Serum serum
        annotation (Placement(transformation(extent={{-10,-14},{10,6}})));
      Liver liver
        annotation (Placement(transformation(extent={{-52,2},{-32,22}})));
      Spleen spleen
        annotation (Placement(transformation(extent={{22,16},{42,36}})));
      Duodenum duodenum
        annotation (Placement(transformation(extent={{-52,-28},{-32,-8}})));
      OtherOrgans otherOrgans
        annotation (Placement(transformation(extent={{-10,-42},{10,-22}})));
      BoneMarrow boneMarrow
        annotation (Placement(transformation(extent={{26,-40},{46,-20}})));
      RBC rBC annotation (Placement(transformation(extent={{60,-16},{80,4}})));
      Hepcidin hepcidin
        annotation (Placement(transformation(extent={{-18,36},{2,56}})));
      IL6 iL6 annotation (Placement(transformation(extent={{-50,36},{-30,56}})));
      Bmp6 bmp6
        annotation (Placement(transformation(extent={{-50,60},{-30,80}})));
      LPS lPS annotation (Placement(transformation(extent={{-86,36},{-66,56}})));
      Bodylight.Types.Constants.MassFlowRateConst transfusion(k(displayUnit=
              "ug.h-1") = 0)
        annotation (Placement(transformation(extent={{-40,-60},{-26,-48}})));
      Bodylight.Types.Constants.MassFlowRateConst bleeding(k(displayUnit=
              "ug.h-1") = 0)
        annotation (Placement(transformation(extent={{-40,-76},{-26,-64}})));
      Bodylight.Types.Constants.MassFlowRateConst Fe_food(k(displayUnit=
              "ug.h-1") = 6.0731166666667e-11)
        annotation (Placement(transformation(extent={{-98,-22},{-84,-10}})));
    equation
      connect(spleen.Serum_out, serum.Spleen_in) annotation (Line(points={{42.8,
              34.8},{50,34.8},{50,10},{-20,10},{-20,1.6},{-10.6,1.6}}, color={0,
              0,127}));
      connect(duodenum.Serum_out, serum.Duodenum_in) annotation (Line(points={{
              -31.2,-9.2},{-26,-9.2},{-26,-1.4},{-10.6,-1.4}}, color={0,0,127}));
      connect(otherOrgans.Serum_out, serum.OtherOrgans_in) annotation (Line(
            points={{10.8,-23.2},{16,-23.2},{16,-16},{-20,-16},{-20,-4.4},{
              -10.6,-4.4}}, color={0,0,127}));
      connect(serum.BoneMarrow_out, boneMarrow.Serum_in) annotation (Line(
            points={{10.6,2},{20,2},{20,-30},{25.4,-30}}, color={0,0,127}));
      connect(boneMarrow.RBC_out, rBC.BoneMarrow_in) annotation (Line(points={{
              46.6,-28},{48,-28},{48,1.2},{59.4,1.2}}, color={0,0,127}));
      connect(rBC.Spleen_out, spleen.RBC_in) annotation (Line(points={{80.4,-6},
              {86,-6},{86,6},{14,6},{14,31.6},{21.4,31.6}}, color={0,0,127}));
      connect(boneMarrow.Spleen_out, spleen.BoneMarrow_in) annotation (Line(
            points={{46.6,-32},{50,-32},{50,-4},{22,-4},{22,12},{8,12},{8,34.6},
              {21.4,34.6}}, color={0,0,127}));
      connect(serum.Liver_out, liver.Serum_in) annotation (Line(points={{10.6,
              4.6},{12,4.6},{12,28},{-60,28},{-60,20.6},{-52.6,20.6}}, color={0,
              0,127}));
      connect(liver.Serum_out, serum.Liver_in) annotation (Line(points={{-31.2,
              20.8},{-26,20.8},{-26,4.6},{-10.6,4.6}}, color={0,0,127}));
      connect(serum.Duodenum_out, duodenum.Serum_in) annotation (Line(points={{
              10.6,-0.6},{14,-0.6},{14,-18},{-28,-18},{-28,-4},{-58,-4},{-58,
              -9.4},{-52.6,-9.4}}, color={0,0,127}));
      connect(serum.OtherOrgans_out, otherOrgans.Serum_in) annotation (Line(
            points={{10.6,-3.2},{18,-3.2},{18,-20},{-16,-20},{-16,-23.4},{-10.6,
              -23.4}}, color={0,0,127}));
      connect(hepcidin.hep, spleen.hep) annotation (Line(points={{2.6,54},{4,54},
              {4,20.4},{21.4,20.4}}, color={0,0,127}));
      connect(hepcidin.hep, liver.hep) annotation (Line(points={{2.6,54},{4,54},
              {4,30},{-60,30},{-60,6.4},{-52.6,6.4}}, color={0,0,127}));
      connect(liver.hep, duodenum.hep) annotation (Line(points={{-52.6,6.4},{
              -60,6.4},{-60,-23.6},{-52.6,-23.6}}, color={0,0,127}));
      connect(hepcidin.hep, duodenum.hep) annotation (Line(points={{2.6,54},{4,
              54},{4,30},{-60,30},{-60,-23.6},{-52.6,-23.6}}, color={0,0,127}));
      connect(duodenum.hep, otherOrgans.hep) annotation (Line(points={{-52.6,
              -23.6},{-60,-23.6},{-60,-37.6},{-10.6,-37.6}}, color={0,0,127}));
      connect(iL6.il6, hepcidin.il6) annotation (Line(points={{-29.6,54.4},{-18,
              54.4},{-18,54}}, color={0,0,127}));
      connect(iL6.il6, spleen.il6) annotation (Line(points={{-29.6,54.4},{-24,
              54.4},{-24,24},{-10,24},{-10,17.4},{21.4,17.4}}, color={0,0,127}));
      connect(iL6.il6, liver.il6) annotation (Line(points={{-29.6,54.4},{-24,
              54.4},{-24,24},{-64,24},{-64,3.4},{-52.6,3.4}}, color={0,0,127}));
      connect(iL6.il6, duodenum.il6) annotation (Line(points={{-29.6,54.4},{-24,
              54.4},{-24,24},{-64,24},{-64,-26.6},{-52.6,-26.6}}, color={0,0,
              127}));
      connect(iL6.il6, otherOrgans.il6) annotation (Line(points={{-29.6,54.4},{
              -24,54.4},{-24,24},{-64,24},{-64,-40.6},{-10.6,-40.6}}, color={0,
              0,127}));
      connect(bmp6.Bmp6, hepcidin.Bmp6) annotation (Line(points={{-29.8,78.6},{
              -22,78.6},{-22,49.8},{-18,49.8}}, color={0,0,127}));
      connect(lPS.LPS, iL6.LPS) annotation (Line(points={{-65.8,54.4},{-50,54.4},
              {-50,54}}, color={0,0,127}));
      connect(liver.Fe, bmp6.Fe_liv) annotation (Line(points={{-31.8,5},{-28,5},
              {-28,-2},{-88,-2},{-88,78},{-50,78}}, color={0,0,127}));
      connect(serum.Fe, bmp6.Fe_ser) annotation (Line(points={{10.2,-11.6},{46,
              -11.6},{46,82},{-60,82},{-60,73.6},{-50,73.6}}, color={0,0,127}));
      connect(transfusion.y, serum.Transfusion) annotation (Line(points={{
              -24.25,-54},{-18,-54},{-18,-9.4},{-10.4,-9.4}}, color={0,0,127}));
      connect(transfusion.y, rBC.Transfusion) annotation (Line(points={{-24.25,
              -54},{52,-54},{52,-9.8},{59.4,-9.8}}, color={0,0,127}));
      connect(bleeding.y, serum.Bleeding) annotation (Line(points={{-24.25,-70},
              {-14,-70},{-14,-12.6},{-10.4,-12.6}}, color={0,0,127}));
      connect(bleeding.y, rBC.Bleeding) annotation (Line(points={{-24.25,-70},{
              54,-70},{54,-13.6},{59.4,-13.6}}, color={0,0,127}));
      connect(Fe_food.y, duodenum.Food_in) annotation (Line(points={{-82.25,-16},
              {-66.425,-16},{-66.425,-12.4},{-52.6,-12.4}}, color={0,0,127}));
    end Test_Full;
  end FeMetabolism;
  annotation (uses(Modelica(version="4.0.0"), Bodylight(version="1.0")));
end EnterocyteMucosalBlock;
