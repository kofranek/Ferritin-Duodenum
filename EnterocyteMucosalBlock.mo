within ;
package EnterocyteMucosalBlock "Enterocyte mucosal block"

  package types
    type SurfaceConcentration = Real(
      quantity="AmountOfSubstancePerArea",
      unit="mol/m2",
      displayUnit="mol/cm2"
      );
    type AmountFluxPerArea = Real(
      quantity    = "AmountFluxPerArea",
      unit        = "mol/(m2.s)",
      displayUnit = "mol/(cm2.s)"
      );
    type AmountRatePerVolume = Real(
      quantity    = "AmountRatePerVolume",
      unit        = "mol/(m3.s)",
      displayUnit = "mol/(l.s)"
      );
    type SecondOrderRateConstant = Real(
      quantity    = "SecondOrderRateConstant",
      unit        = "m3/(mol.s)",
      displayUnit = "l/(mol.s)"
      );
  end types;

  package models
    model CellularFerritinIronStorageModel "Cellular ferritin iron storage"

      constant Modelica.Units.SI.Volume V_cell(
        displayUnit = "L") = 1.4e-12 / 1000 "cell volume";

      //Species
      Modelica.Units.SI.MolarConcentration LIP(
        displayUnit = "mol/L",
        start = 1e-05 * 1000) "labile iron pool";
      Modelica.Units.SI.MolarConcentration FT_cage(
        displayUnit = "mol/L",
        start = 5e-09 * 1000) "FT-cage";
      Modelica.Units.SI.MolarConcentration core(
        displayUnit = "mol/L",
        start = 7.5e-06 * 1000) "core";
      Modelica.Units.SI.MolarConcentration DFP(
        displayUnit = "mol/L",
        start = 0) "diferric peroxo complex";

      Real atoms_per_cage_transient "Transient number of Fe atoms that are stored inside the core of a ferritin cage";

      parameter Integer H = 4 "H subunits";
      parameter Integer L = 24 - H "L subunits";
      //parameter Integer atoms_per_cage_transient = 1500 "Transient number of Fe atoms that are stored inside the core of a ferritin cage";
      parameter Integer rN = 50;
      parameter Integer rO = 2;

      parameter Modelica.Units.SI.Frequency k_FTlysis = 1.203e-05;

      parameter Real FT_Expression(
        quantity = "ReactionRate",
        unit = "mol/(m3.s)",
        displayUnit = "mol/(l.s)")
         = 6.015e-14 * 1000;

      //FT degradation
      Real FT_Degradation(
        quantity = "ReactionRate",
        unit = "mol/(m3.s)",
        displayUnit = "mol/(l.s)");

      //FT degradation core release
      Real CoreRelease(
        quantity = "ReactionRate",
        unit = "mol/(m3.s)",
          displayUnit = "mol/(l.s)");

      //Oxidation (2 LIP -> DFP)
      Real Oxidation(
        quantity = "ReactionRate",
        unit = "mol/(m3.s)",
        displayUnit = "mol/(l.s)");

      parameter Modelica.Units.SI.Frequency k_cat_oxidation = 591 "catalytic turnover number";
      parameter Modelica.Units.SI.MolarConcentration K_m_oxidation(
        displayUnit = "mol/L") = 0.35 "Michaelis constant";
      parameter Real n_oxidation = 1.3 "Hill coefficient)";

      //Reduction (DFP -> 2 LIP)
      Real Reduction(
        quantity = "ReactionRate",
        unit = "mol/(m3.s)",
        displayUnit = "mol/(l.s)");

      parameter Modelica.Units.SI.Frequency k_deg = 0.2605 "rate constant";

      //Nucleation (2 DFP -> 4 core)
      Real Nucleation(
        quantity = "ReactionRate",
        unit = "mol/(m3.s)",
        displayUnit = "mol/(l.s)");

      parameter Modelica.Units.SI.Frequency k_cat_nucleation = 5e07 "catalytic turnover number";
      parameter Modelica.Units.SI.MolarConcentration K_i_nucleation(
        displayUnit = "mol/L") = 0.461598 "inhibition constant";
      parameter Integer n_nucleation = 4 "Hill coefficient";

      //Mineralization (DFP -> 2 core)
      Real Mineralization(
        quantity = "ReactionRate",
        unit = "mol/(m3.s)",
        displayUnit = "mol/(l.s)");

      parameter Modelica.Units.SI.Frequency k_cat_mineralization = 0.101564 "catalytic turnover number";
      parameter Modelica.Units.SI.MolarConcentration K_m_mineralization(
        displayUnit = "mol/L") = 5e-03 "Michaelis constant";
      parameter Modelica.Units.SI.MolarConcentration K_i_mineralization(
        displayUnit = "mol/L") = 4.6458 "inhibition constant";
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
        coordinateSystem(preserveAspectRatio=false)));

    end CellularFerritinIronStorageModel;

    model EnterocyteMucosalBlockModel "Enterocyte mucosal block"

      //Compartments

      constant Modelica.Units.SI.Volume cell(
        displayUnit = "L") = 1.4e-12 * 1e-3;
      constant Modelica.Units.SI.Volume upper(
        displayUnit = "L") = 6.67e-10 * 1e-3;
      constant Modelica.Units.SI.Volume lower(
        displayUnit = "L") = 8.57e-10 * 1e-3;
      constant Modelica.Units.SI.Area apical_mem(
        displayUnit = "cm2") = 1.5e-5 * 1e-4;
      constant Modelica.Units.SI.Area BLM(
        displayUnit = "cm2") = 1.5e-7 * 1e-4;

      //Species

      //cell

      Modelica.Units.SI.MolarConcentration FT_cage(
        displayUnit = "mol/L",
        start = 2.375189822e-9 * 1e3) "FT-cage";
      Modelica.Units.SI.MolarConcentration core(
        displayUnit = "mol/L",
        start = 3.682217017e-6 * 1e3) "core";
      Modelica.Units.SI.MolarConcentration DFP(
        displayUnit = "mol/L",
        start = 1.344769304e-10 * 1e3) "diferric peroxo complex";
      Modelica.Units.SI.MolarConcentration LIP(
        displayUnit = "mol/L",
        start = 1.223884748e-7 * 1e3) "labile iron pool";
      Modelica.Units.SI.MolarConcentration IRPs_active(
        displayUnit = "mol/L",
        start = 6.889335935e-11 * 1e3) "iron regulatory proteins (active)";
      Modelica.Units.SI.MolarConcentration IRPs_inactive(
        displayUnit = "mol/L",
        start = 7.264345126e-12 * 1e3) "iron regulatory proteins (inactive)";

      //lower

      Modelica.Units.SI.MolarConcentration Fe_blood(
        displayUnit = "mol/L",
        start = 4.958456433e-9 * 1e3);
      Modelica.Units.SI.MolarConcentration body_fe(
        displayUnit = "mol/L",
        start = 0);

      //upper

      Modelica.Units.SI.MolarConcentration Fe_lumen(
        displayUnit = "mol/L",
        start = 1.25e-8 * 1e3);

      //apical_mem

      EnterocyteMucosalBlock.types.SurfaceConcentration DMT1(start=
            3.123079828e-12*1e4);
      EnterocyteMucosalBlock.types.SurfaceConcentration DMT1_vesicular(start=
            1.876958455e-12*1e4);

      //BLM

      EnterocyteMucosalBlock.types.SurfaceConcentration FPN_active(start=
            9.979749648e-14*1e4);
      EnterocyteMucosalBlock.types.SurfaceConcentration FPN_internalized(start=
            2.02503521e-16*1e4);

      //Reactions

      //FT expression ( -> FT-cage;  IRPs_active)
      EnterocyteMucosalBlock.types.AmountRatePerVolume FT_Expression;
      EnterocyteMucosalBlock.types.AmountRatePerVolume k_cat_FT_Expression=
          7.68e-14*1e3;
      parameter Integer n_FT_Expression = 1;
      parameter Modelica.Units.SI.MolarConcentration K_FT_Expression(
        displayUnit = "mol/L") = 1.4e-11 * 1e3;

      //FT degradation (FT-cage -> )
      EnterocyteMucosalBlock.types.AmountRatePerVolume FT_Degradation;
      parameter Modelica.Units.SI.Frequency k_FT_Degradation = 5.461499585e-6;

      //FT degradation core release (core -> LIP; FT-cage core)
      EnterocyteMucosalBlock.types.AmountRatePerVolume
        FT_Degradation_Core_Release;

      //FT Fe oxidation (2 * LIP -> DFP; FT-cage)
      EnterocyteMucosalBlock.types.AmountRatePerVolume FT_Fe_Oxidation;
      parameter Modelica.Units.SI.Frequency k_cat_FT_Fe_Oxidation = 591;
      parameter Modelica.Units.SI.MolarConcentration K_m_FT_Fe_Oxidation(
        displayUnit = "mol/L") = 0.35e-3 * 1e3;
      parameter Real n_FT_Fe_Oxidation = 1.3;

      //FT Fe Reduction (DFP -> 2 * LIP)
      EnterocyteMucosalBlock.types.AmountRatePerVolume FT_Fe_Reduction;
      parameter Modelica.Units.SI.Frequency k_FT_Fe_Reduction = 0.2605;

      //FT Nucleation (2 * DFP -> 4 * core; FT-cage core)
      EnterocyteMucosalBlock.types.AmountRatePerVolume FT_Nucleation;
      parameter Modelica.Units.SI.Frequency k_cat_FT_Nucleation = 5e7;
      parameter Modelica.Units.SI.MolarConcentration K_i_FT_Nucleation(
        displayUnit = "mol/L") = 0.461598e-3 * 1e3;
      parameter Integer n_FT_Nucleation = 4;

      //FT core formation (Mineralization) (DFP -> 2 * core; core)
      EnterocyteMucosalBlock.types.AmountRatePerVolume FT_Core_Formation;
      parameter Modelica.Units.SI.Frequency k_cat_FT_Core_Formation = 0.101564;
      parameter Modelica.Units.SI.MolarConcentration K_m_FT_Core_Formation(
        displayUnit = "mol/L") = 5e-06 * 1e3;
      parameter Modelica.Units.SI.MolarConcentration K_i_FT_Core_Formation(
        displayUnit = "mol/L") = 4.6458e-3 * 1e3;
      parameter Integer n_FT_Core_Formation = 4;
      parameter Integer m_FT_Core_Formation = 8;

      //IRPs degradation (IRPs_active -> IRPs_inactive;  LIP)
      EnterocyteMucosalBlock.types.AmountRatePerVolume IRPs_Degradation;
      parameter EnterocyteMucosalBlock.types.SecondOrderRateConstant
        k_cat_IRPs_Degradation=3.99474*1e-3;

      //IRPs activation (IRPs_inactive -> IRPs_active)
      EnterocyteMucosalBlock.types.AmountRatePerVolume IRPs_Activation;
      parameter Modelica.Units.SI.Frequency k_cat_IRPs_Activation = 4.63671e-6;

      //Body Sequestration (Fe_blood -> body_fe)
      EnterocyteMucosalBlock.types.AmountRatePerVolume Body_Sequestration;
      parameter Modelica.Units.SI.Frequency k_cat_Body_Sequestration = 0.329e-3;

      //Fe Basal Uptake (Fe_blood -> LIP)
      Modelica.Units.SI.MolarFlowRate Fe_Basal_Uptake;
      parameter Modelica.Units.SI.VolumeFlowRate k_cat_Fe_Basal_Uptake(
        displayUnit="l/s")= 2.22055e-16 * 1e-3;

      //paracellular Fe transport (Fe_lumen = Fe_blood)
      Modelica.Units.SI.MolarFlowRate Paracellular_Fe_Transport;
      parameter Modelica.Units.SI.VolumeFlowRate k_for_Paracellular_Fe_Transport(
        displayUnit="l/s")= 3.87951e-22 * 1e-3;
      parameter Modelica.Units.SI.VolumeFlowRate k_rev_Paracellular_Fe_Transport(
        displayUnit="l/s")= 3.1746e-15 * 1e-3;

      //DMT1 endocytosis free (DMT1 -> DMT1_vesicular)
      EnterocyteMucosalBlock.types.AmountFluxPerArea DMT1_Endocytosis_Free;
      parameter Modelica.Units.SI.Frequency k_cat_DMT1_Endocytosis_Free = 29.4233;

      //DMT1 endocytosis LIP modified (DMT1 -> DMT1_vesicular;  LIP)
      Modelica.Units.SI.MolarFlowRate DMT1_Endocytosis_Modified;
      parameter Modelica.Units.SI.DiffusionCoefficient k_cat_DMT1_Endocytosis_Modified(
        displayUnit="cm2/s") = 14.516 * 1e-4;
      parameter Modelica.Units.SI.MolarConcentration K_m_DMT1_Endocytosis_Modified(
        displayUnit = "mol/L") = 2.80591 * 1e3;
      parameter Real n_DMT1_Endocytosis_Modified = 1.03128;

      //DMT1 fusion (DMT1_vesicular -> DMT1;  DMT1)
      EnterocyteMucosalBlock.types.AmountFluxPerArea DMT1_Fusion;
      parameter Modelica.Units.SI.Frequency k_cat_DMT1_Fusion = 48.9989;

      //DMT1 iron transport (Fe_lumen -> LIP;  DMT1)
      EnterocyteMucosalBlock.types.AmountFluxPerArea DMT1_Iron_Transport;
      parameter Modelica.Units.SI.Frequency k_cat_DMT1_Iron_Transport = 6844.7;
      parameter Modelica.Units.SI.MolarConcentration K_m_DMT1_Iron_Transport(
        displayUnit = "mol/L") = 2.835 * 1e3;
      parameter Integer n_DMT1_Iron_Transport = 1;

      //FPN-inactivation (FPN_active -> FPN_internalized;  Fe_blood)
      Modelica.Units.SI.MolarFlowRate FPN_Inactivation;
      parameter Modelica.Units.SI.Frequency k_cat_FPN_Inactivation = 1.44264e-6;
      parameter Modelica.Units.SI.MolarConcentration K_m_FPN_Inactivation(
        displayUnit = "mol/L") = 1.22073e-5 * 1e3;
      parameter Real n_FPN_Inactivation = 2.71609;

      //FPN-activation (FPN_internalized -> FPN_active)
      Modelica.Units.SI.MolarFlowRate FPN_Activation;
      parameter Modelica.Units.SI.DiffusionCoefficient k_cat_FPN_Activation(
        displayUnit="cm2/s") = 4.37363e-13 * 1e-4;

      //FPN iron transport (LIP -> Fe_blood;  FPN_active)
      EnterocyteMucosalBlock.types.AmountFluxPerArea FPN_Iron_Transport;
      parameter Modelica.Units.SI.Frequency k_cat_FPN_Iron_Transport = 1.88317;
      parameter Modelica.Units.SI.MolarConcentration K_m_FPN_Iron_Transport(
        displayUnit = "mol/L") = 2.31608e-6 * 1e3;
      parameter Integer n_FPN_Iron_Transport = 1;

      //Global Quantities

      Real atoms_per_cage_transient "Transient number of Fe atoms that are stored inside the core of a ferritin cage";

      parameter Integer H = 24 "H subunits";
      parameter Integer L = 24 - H "L subunits";
      //parameter Real atoms_per_cage_transient = 1550.28325867654 "Transient number of Fe atoms that are stored inside the core of a ferritin cage";
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

      //e-25, must be e-29
      FPN_Inactivation = k_cat_FPN_Inactivation * FPN_active * Fe_blood ^ n_FPN_Inactivation / (K_m_FPN_Inactivation ^ n_FPN_Inactivation + Fe_blood ^ n_FPN_Inactivation);

      FPN_Activation = k_cat_FPN_Activation * FPN_internalized;

      FPN_Iron_Transport = k_cat_FPN_Iron_Transport * FPN_active * LIP ^ n_FPN_Iron_Transport / (K_m_FPN_Iron_Transport ^ n_FPN_Iron_Transport + LIP ^ n_FPN_Iron_Transport);

      IRPs_Degradation = k_cat_IRPs_Degradation * IRPs_active * LIP;

      IRPs_Activation = k_cat_IRPs_Activation * IRPs_inactive;

      Body_Sequestration = k_cat_Body_Sequestration * Fe_blood;

      Fe_Basal_Uptake = k_cat_Fe_Basal_Uptake * Fe_blood;

      Paracellular_Fe_Transport = k_for_Paracellular_Fe_Transport * Fe_lumen - k_rev_Paracellular_Fe_Transport * Fe_blood;

      //cell

      der(FT_cage) * cell = -cell * FT_Degradation + cell * FT_Expression;

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

      //Compartments

      constant Modelica.Units.SI.Volume cell(
        displayUnit = "L") = 1.4e-12 * 1e-3;

      //Species

      //cell

      Modelica.Units.SI.MolarConcentration FT_cage(
        displayUnit = "mol/L",
        start = 2.375189822e-9 * 1e3) "FT-cage";
      Modelica.Units.SI.MolarConcentration core(
        displayUnit = "mol/L",
        start = 3.682217017e-6 * 1e3) "core";
      Modelica.Units.SI.MolarConcentration DFP(
        displayUnit = "mol/L",
        start = 1.344769304e-10 * 1e3) "diferric peroxo complex";
      Modelica.Units.SI.MolarConcentration LIP(
        displayUnit = "mol/L",
        start = 1.223884748e-7 * 1e3) "labile iron pool";
      Modelica.Units.SI.MolarConcentration IRPs_active(
        displayUnit = "mol/L",
        start = 6.889335935e-11 * 1e3) "iron regulatory proteins (active)";
      Modelica.Units.SI.MolarConcentration IRPs_inactive(
        displayUnit = "mol/L",
        start = 7.264345126e-12 * 1e3) "iron regulatory proteins (inactive)";

      //Reactions

      //FT expression ( -> FT-cage;  IRPs_active)
      EnterocyteMucosalBlock.types.AmountRatePerVolume FT_Expression;
      EnterocyteMucosalBlock.types.AmountRatePerVolume k_cat_FT_Expression=
          7.68e-14*1e3;
      parameter Integer n_FT_Expression = 1;
      parameter Modelica.Units.SI.MolarConcentration K_FT_Expression(
        displayUnit = "mol/L") = 1.4e-11 * 1e3;

      //FT degradation (FT-cage -> )
      EnterocyteMucosalBlock.types.AmountRatePerVolume FT_Degradation;
      parameter Modelica.Units.SI.Frequency k_FT_Degradation = 5.461499585e-6;

      //FT degradation core release (core -> LIP; FT-cage core)
      EnterocyteMucosalBlock.types.AmountRatePerVolume
        FT_Degradation_Core_Release;

      //FT Fe oxidation (2 * LIP -> DFP; FT-cage)
      EnterocyteMucosalBlock.types.AmountRatePerVolume FT_Fe_Oxidation;
      parameter Modelica.Units.SI.Frequency k_cat_FT_Fe_Oxidation = 591;
      parameter Modelica.Units.SI.MolarConcentration K_m_FT_Fe_Oxidation(
        displayUnit = "mol/L") = 0.35e-3 * 1e3;
      parameter Real n_FT_Fe_Oxidation = 1.3;

      //FT Fe Reduction (DFP -> 2 * LIP)
      EnterocyteMucosalBlock.types.AmountRatePerVolume FT_Fe_Reduction;
      parameter Modelica.Units.SI.Frequency k_FT_Fe_Reduction = 0.2605;

      //FT Nucleation (2 * DFP -> 4 * core; FT-cage core)
      EnterocyteMucosalBlock.types.AmountRatePerVolume FT_Nucleation;
      parameter Modelica.Units.SI.Frequency k_cat_FT_Nucleation = 5e7;
      parameter Modelica.Units.SI.MolarConcentration K_i_FT_Nucleation(
        displayUnit = "mol/L") = 0.461598e-3 * 1e3;
      parameter Integer n_FT_Nucleation = 4;

      //FT core formation (Mineralization) (DFP -> 2 * core; core)
      EnterocyteMucosalBlock.types.AmountRatePerVolume FT_Core_Formation;
      parameter Modelica.Units.SI.Frequency k_cat_FT_Core_Formation = 0.101564;
      parameter Modelica.Units.SI.MolarConcentration K_m_FT_Core_Formation(
        displayUnit = "mol/L") = 5e-06 * 1e3;
      parameter Modelica.Units.SI.MolarConcentration K_i_FT_Core_Formation(
        displayUnit = "mol/L") = 4.6458e-3 * 1e3;
      parameter Integer n_FT_Core_Formation = 4;
      parameter Integer m_FT_Core_Formation = 8;

      //IRPs degradation (IRPs_active -> IRPs_inactive;  LIP)
      EnterocyteMucosalBlock.types.AmountRatePerVolume IRPs_Degradation;
      parameter EnterocyteMucosalBlock.types.SecondOrderRateConstant
        k_cat_IRPs_Degradation=3.99474*1e-3;

      //IRPs activation (IRPs_inactive -> IRPs_active)
      EnterocyteMucosalBlock.types.AmountRatePerVolume IRPs_Activation;
      parameter Modelica.Units.SI.Frequency k_cat_IRPs_Activation = 4.63671e-6;

      //Global Quantities

      Real atoms_per_cage_transient "Transient number of Fe atoms that are stored inside the core of a ferritin cage";

      parameter Integer H = 24 "H subunits";
      parameter Integer L = 24 - H "L subunits";

      parameter Integer rN = 50;
      parameter Integer rO = 2;

    equation

      atoms_per_cage_transient = core / FT_cage;

      //10
      FT_Core_Formation = (k_cat_FT_Core_Formation * DFP * core) / (K_m_FT_Core_Formation + DFP)
        * K_i_FT_Core_Formation ^ n_FT_Core_Formation / (K_i_FT_Core_Formation ^ n_FT_Core_Formation + core ^ n_FT_Core_Formation)
        * (4300 ^ m_FT_Core_Formation - atoms_per_cage_transient ^ m_FT_Core_Formation) / 4300 ^ m_FT_Core_Formation;

      //11
      FT_Degradation = k_FT_Degradation * FT_cage;

      //12
      FT_Degradation_Core_Release = k_FT_Degradation * core;

      //13
      FT_Expression = k_cat_FT_Expression * (1 - IRPs_active ^ n_FT_Expression / (K_FT_Expression ^ n_FT_Expression + IRPs_active ^ n_FT_Expression));

      //14
      FT_Fe_Oxidation = (k_cat_FT_Fe_Oxidation * (H + rO) / (24 + rO) * FT_cage * LIP ^ n_FT_Fe_Oxidation)
        / (K_m_FT_Fe_Oxidation ^ n_FT_Fe_Oxidation + LIP ^ n_FT_Fe_Oxidation);

      //15
      FT_Fe_Reduction = k_FT_Fe_Reduction * DFP;

      //16
      FT_Nucleation = k_cat_FT_Nucleation * DFP ^ 2 * FT_cage * (L + rN) / (24 + rN)
        * K_i_FT_Nucleation ^ n_FT_Nucleation / (K_i_FT_Nucleation ^ n_FT_Nucleation + core ^ n_FT_Nucleation);

      //17
      IRPs_Activation = k_cat_IRPs_Activation * IRPs_inactive;

      //18
      IRPs_Degradation = k_cat_IRPs_Degradation * IRPs_active * LIP;

      //cell

      der(FT_cage) * cell = -cell * FT_Degradation + cell * FT_Expression;

      der(core) * cell = 2 * cell * FT_Core_Formation
        + 4 * cell * FT_Nucleation
        - cell * FT_Degradation_Core_Release;

      der(DFP) * cell = -cell * FT_Core_Formation
        - cell * FT_Fe_Reduction
        + cell * FT_Fe_Oxidation
        - 2 * cell * FT_Nucleation;

      der(LIP) * cell = 2 * cell * FT_Fe_Reduction
        - 2 * cell * FT_Fe_Oxidation
        + cell * FT_Degradation_Core_Release;

      der(IRPs_active) * cell = cell * IRPs_Activation
        - cell * IRPs_Degradation;

      der(IRPs_inactive) * cell = -cell * IRPs_Activation
        + cell * IRPs_Degradation;

      annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
            coordinateSystem(preserveAspectRatio=false)));
    end EnterocyteMucosalBlockShortModel;

    model EnterocyteMucosalBlockShortDeviationModel
      "Enterocyte mucosal block (short) + Fe cell total Deviation"

      //Compartments

      constant Modelica.Units.SI.Volume cell(
        displayUnit = "L") = 1.4e-12 * 1e-3;

      //Species

      //cell

      Modelica.Units.SI.MolarConcentration FT_cage(
        displayUnit = "mol/L",
        start = 2.375189822e-9 * 1e3) "FT-cage";
      Modelica.Units.SI.MolarConcentration core(
        displayUnit = "mol/L",
        start = 3.682217017e-6 * 1e3) "core";
      Modelica.Units.SI.MolarConcentration DFP(
        displayUnit = "mol/L",
        start = 1.344769304e-10 * 1e3) "diferric peroxo complex";
      Modelica.Units.SI.MolarConcentration LIP(
        displayUnit = "mol/L",
        start = 1.223884748e-7 * 1e3) "labile iron pool";
      Modelica.Units.SI.MolarConcentration IRPs_active(
        displayUnit = "mol/L",
        start = 6.889335935e-11 * 1e3) "iron regulatory proteins (active)";
      Modelica.Units.SI.MolarConcentration IRPs_inactive(
        displayUnit = "mol/L",
        start = 7.264345126e-12 * 1e3) "iron regulatory proteins (inactive)";

      //Reactions

      //FT expression ( -> FT-cage;  IRPs_active)
      EnterocyteMucosalBlock.types.AmountRatePerVolume FT_Expression;
      EnterocyteMucosalBlock.types.AmountRatePerVolume k_cat_FT_Expression=
          7.68e-14*1e3;
      parameter Integer n_FT_Expression = 1;
      parameter Modelica.Units.SI.MolarConcentration K_FT_Expression(
        displayUnit = "mol/L") = 1.4e-11 * 1e3;

      //FT degradation (FT-cage -> )
      EnterocyteMucosalBlock.types.AmountRatePerVolume FT_Degradation;
      parameter Modelica.Units.SI.Frequency k_FT_Degradation = 5.461499585e-6;

      //FT degradation core release (core -> LIP; FT-cage core)
      EnterocyteMucosalBlock.types.AmountRatePerVolume
        FT_Degradation_Core_Release;

      //FT Fe oxidation (2 * LIP -> DFP; FT-cage)
      EnterocyteMucosalBlock.types.AmountRatePerVolume FT_Fe_Oxidation;
      parameter Modelica.Units.SI.Frequency k_cat_FT_Fe_Oxidation = 591;
      parameter Modelica.Units.SI.MolarConcentration K_m_FT_Fe_Oxidation(
        displayUnit = "mol/L") = 0.35e-3 * 1e3;
      parameter Real n_FT_Fe_Oxidation = 1.3;

      //FT Fe Reduction (DFP -> 2 * LIP)
      EnterocyteMucosalBlock.types.AmountRatePerVolume FT_Fe_Reduction;
      parameter Modelica.Units.SI.Frequency k_FT_Fe_Reduction = 0.2605;

      //FT Nucleation (2 * DFP -> 4 * core; FT-cage core)
      EnterocyteMucosalBlock.types.AmountRatePerVolume FT_Nucleation;
      parameter Modelica.Units.SI.Frequency k_cat_FT_Nucleation = 5e7;
      parameter Modelica.Units.SI.MolarConcentration K_i_FT_Nucleation(
        displayUnit = "mol/L") = 0.461598e-3 * 1e3;
      parameter Integer n_FT_Nucleation = 4;

      //FT core formation (Mineralization) (DFP -> 2 * core; core)
      EnterocyteMucosalBlock.types.AmountRatePerVolume FT_Core_Formation;
      parameter Modelica.Units.SI.Frequency k_cat_FT_Core_Formation = 0.101564;
      parameter Modelica.Units.SI.MolarConcentration K_m_FT_Core_Formation(
        displayUnit = "mol/L") = 5e-06 * 1e3;
      parameter Modelica.Units.SI.MolarConcentration K_i_FT_Core_Formation(
        displayUnit = "mol/L") = 4.6458e-3 * 1e3;
      parameter Integer n_FT_Core_Formation = 4;
      parameter Integer m_FT_Core_Formation = 8;

      //IRPs degradation (IRPs_active -> IRPs_inactive;  LIP)
      EnterocyteMucosalBlock.types.AmountRatePerVolume IRPs_Degradation;
      parameter EnterocyteMucosalBlock.types.SecondOrderRateConstant
        k_cat_IRPs_Degradation=3.99474*1e-3;

      //IRPs activation (IRPs_inactive -> IRPs_active)
      EnterocyteMucosalBlock.types.AmountRatePerVolume IRPs_Activation;
      parameter Modelica.Units.SI.Frequency k_cat_IRPs_Activation = 4.63671e-6;

      //Global Quantities

      parameter Integer H = 24 "H subunits";
      parameter Integer L = 24 - H "L subunits";
      parameter Real atoms_per_cage_transient = 1550.28325867654 "Transient number of Fe atoms that are stored inside the core of a ferritin cage";
      parameter Integer rN = 50;
      parameter Integer rO = 2;

      //Fe_cell_total

      parameter Real Fe_cell_total_deviation = 0.1;
      Modelica.Units.SI.MolarConcentration Fe_cell_total(
         displayUnit = "mol/L") "Fe in cell (total) = core + LIP + DFP";

      parameter Modelica.Units.SI.MolarConcentration core_start(
        displayUnit = "mol/L") = 3.682217017e-6 * 1e3;
      parameter Modelica.Units.SI.MolarConcentration LIP_start(
        displayUnit = "mol/L") = 1.223884748e-7 * 1e3;
      parameter Modelica.Units.SI.MolarConcentration DFP_start(
        displayUnit = "mol/L") = 1.344769304e-10 * 1e3;

    initial equation

       core = core_start * (1 + Fe_cell_total_deviation);
       LIP = LIP_start * (1 + Fe_cell_total_deviation);
       DFP = DFP_start * (1 + Fe_cell_total_deviation);

    equation

      Fe_cell_total = core + LIP + DFP;

      //10
      FT_Core_Formation = (k_cat_FT_Core_Formation * DFP * core) / (K_m_FT_Core_Formation + DFP)
        * K_i_FT_Core_Formation ^ n_FT_Core_Formation / (K_i_FT_Core_Formation ^ n_FT_Core_Formation + core ^ n_FT_Core_Formation)
        * (4300 ^ m_FT_Core_Formation - atoms_per_cage_transient ^ m_FT_Core_Formation) / 4300 ^ m_FT_Core_Formation;

      //11
      FT_Degradation = k_FT_Degradation * FT_cage;

      //12
      FT_Degradation_Core_Release = k_FT_Degradation * core;

      //13
      FT_Expression = k_cat_FT_Expression * (1 - IRPs_active ^ n_FT_Expression / (K_FT_Expression ^ n_FT_Expression + IRPs_active ^ n_FT_Expression));

      //14
      FT_Fe_Oxidation = (k_cat_FT_Fe_Oxidation * (H + rO) / (24 + rO) * FT_cage * LIP ^ n_FT_Fe_Oxidation)
        / (K_m_FT_Fe_Oxidation ^ n_FT_Fe_Oxidation + LIP ^ n_FT_Fe_Oxidation);

      //15
      FT_Fe_Reduction = k_FT_Fe_Reduction * DFP;

      //16
      FT_Nucleation = k_cat_FT_Nucleation * DFP ^ 2 * FT_cage * (L + rN) / (24 + rN)
        * K_i_FT_Nucleation ^ n_FT_Nucleation / (K_i_FT_Nucleation ^ n_FT_Nucleation + core ^ n_FT_Nucleation);

      //17
      IRPs_Activation = k_cat_IRPs_Activation * IRPs_inactive;

      //18
      IRPs_Degradation = k_cat_IRPs_Degradation * IRPs_active * LIP;

      //cell

      der(FT_cage) * cell = -cell * FT_Degradation + cell * FT_Expression;

      der(core) * cell = 2 * cell * FT_Core_Formation
        + 4 * cell * FT_Nucleation
        - cell * FT_Degradation_Core_Release;

      der(DFP) * cell = -cell * FT_Core_Formation
        - cell * FT_Fe_Reduction
        + cell * FT_Fe_Oxidation
        - 2 * cell * FT_Nucleation;

      der(LIP) * cell = 2 * cell * FT_Fe_Reduction
        - 2 * cell * FT_Fe_Oxidation
        + cell * FT_Degradation_Core_Release;

      der(IRPs_active) * cell = cell * IRPs_Activation
        - cell * IRPs_Degradation;

      der(IRPs_inactive) * cell = -cell * IRPs_Activation
        + cell * IRPs_Degradation;

      annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
            coordinateSystem(preserveAspectRatio=false)));
    end EnterocyteMucosalBlockShortDeviationModel;

    model EnterocyteMucosalBlockShortModel11
      "Enterocyte mucosal block (short)"

      //Compartments

      constant Modelica.Units.SI.Volume cell(
        displayUnit = "L") = 1.4e-12 * 1e-3;

      //Species

      //cell

      Modelica.Units.SI.MolarConcentration FT_cage(
        displayUnit = "mol/L",
        //start = 2.375189822e-9 * 1e3) "FT-cage";
        start = 2.375189822e-9 * 2 * 1e3) "FT-cage";
      Modelica.Units.SI.MolarConcentration core(
        displayUnit = "mol/L",
        start = 3.682217017e-6 * 1e3) "core";
      Modelica.Units.SI.MolarConcentration DFP(
        displayUnit = "mol/L",
        start = 1.344769304e-10 * 1e3) "diferric peroxo complex";
      Modelica.Units.SI.MolarConcentration LIP(
        displayUnit = "mol/L",
        start = 1.223884748e-7 * 1e3) "labile iron pool";
      Modelica.Units.SI.MolarConcentration IRPs_active(
        displayUnit = "mol/L",
        start = 6.889335935e-11 * 1e3) "iron regulatory proteins (active)";
      Modelica.Units.SI.MolarConcentration IRPs_inactive(
        displayUnit = "mol/L",
        start = 7.264345126e-12 * 1e3) "iron regulatory proteins (inactive)";

      //Reactions

      //FT expression ( -> FT-cage;  IRPs_active)
      EnterocyteMucosalBlock.types.AmountRatePerVolume FT_Expression;
      EnterocyteMucosalBlock.types.AmountRatePerVolume k_cat_FT_Expression=
          7.68e-14*1e3;
      parameter Integer n_FT_Expression = 1;
      parameter Modelica.Units.SI.MolarConcentration K_FT_Expression(
        displayUnit = "mol/L") = 1.4e-11 * 1e3;

      //FT degradation (FT-cage -> )
      EnterocyteMucosalBlock.types.AmountRatePerVolume FT_Degradation;
      parameter Modelica.Units.SI.Frequency k_FT_Degradation = 5.461499585e-6;

      //FT degradation core release (core -> LIP; FT-cage core)
      EnterocyteMucosalBlock.types.AmountRatePerVolume
        FT_Degradation_Core_Release;

      //FT Fe oxidation (2 * LIP -> DFP; FT-cage)
      EnterocyteMucosalBlock.types.AmountRatePerVolume FT_Fe_Oxidation;
      parameter Modelica.Units.SI.Frequency k_cat_FT_Fe_Oxidation = 591;
      parameter Modelica.Units.SI.MolarConcentration K_m_FT_Fe_Oxidation(
        displayUnit = "mol/L") = 0.35e-3 * 1e3;
      parameter Real n_FT_Fe_Oxidation = 1.3;

      //FT Fe Reduction (DFP -> 2 * LIP)
      EnterocyteMucosalBlock.types.AmountRatePerVolume FT_Fe_Reduction;
      parameter Modelica.Units.SI.Frequency k_FT_Fe_Reduction = 0.2605;

      //FT Nucleation (2 * DFP -> 4 * core; FT-cage core)
      EnterocyteMucosalBlock.types.AmountRatePerVolume FT_Nucleation;
      parameter Modelica.Units.SI.Frequency k_cat_FT_Nucleation = 5e7;
      parameter Modelica.Units.SI.MolarConcentration K_i_FT_Nucleation(
        displayUnit = "mol/L") = 0.461598e-3 * 1e3;
      parameter Integer n_FT_Nucleation = 4;

      //FT core formation (Mineralization) (DFP -> 2 * core; core)
      EnterocyteMucosalBlock.types.AmountRatePerVolume FT_Core_Formation;
      parameter Modelica.Units.SI.Frequency k_cat_FT_Core_Formation = 0.101564;
      parameter Modelica.Units.SI.MolarConcentration K_m_FT_Core_Formation(
        displayUnit = "mol/L") = 5e-06 * 1e3;
      parameter Modelica.Units.SI.MolarConcentration K_i_FT_Core_Formation(
        displayUnit = "mol/L") = 4.6458e-3 * 1e3;
      parameter Integer n_FT_Core_Formation = 4;
      parameter Integer m_FT_Core_Formation = 8;

      //IRPs degradation (IRPs_active -> IRPs_inactive;  LIP)
      EnterocyteMucosalBlock.types.AmountRatePerVolume IRPs_Degradation;
      parameter EnterocyteMucosalBlock.types.SecondOrderRateConstant
        k_cat_IRPs_Degradation=3.99474*1e-3;

      //IRPs activation (IRPs_inactive -> IRPs_active)
      EnterocyteMucosalBlock.types.AmountRatePerVolume IRPs_Activation;
      parameter Modelica.Units.SI.Frequency k_cat_IRPs_Activation = 4.63671e-6;

      //Global Quantities

      parameter Integer H = 24 "H subunits";
      parameter Integer L = 24 - H "L subunits";
      parameter Real atoms_per_cage_transient = 1550.28325867654 "Transient number of Fe atoms that are stored inside the core of a ferritin cage";
      parameter Integer rN = 50;
      parameter Integer rO = 2;

      Real LIP_rel;

    equation

        LIP_rel = LIP/(LIP + core);

      //10
      FT_Core_Formation = (k_cat_FT_Core_Formation * DFP * core) / (K_m_FT_Core_Formation + DFP)
        * K_i_FT_Core_Formation ^ n_FT_Core_Formation / (K_i_FT_Core_Formation ^ n_FT_Core_Formation + core ^ n_FT_Core_Formation)
        * (4300 ^ m_FT_Core_Formation - atoms_per_cage_transient ^ m_FT_Core_Formation) / 4300 ^ m_FT_Core_Formation;

      //11
      FT_Degradation = k_FT_Degradation * FT_cage;

      //12
      FT_Degradation_Core_Release = k_FT_Degradation * core;

      //13
      FT_Expression = k_cat_FT_Expression * (1 - IRPs_active ^ n_FT_Expression / (K_FT_Expression ^ n_FT_Expression + IRPs_active ^ n_FT_Expression));

      //14
      FT_Fe_Oxidation = (k_cat_FT_Fe_Oxidation * (H + rO) / (24 + rO) * FT_cage * LIP ^ n_FT_Fe_Oxidation)
        / (K_m_FT_Fe_Oxidation ^ n_FT_Fe_Oxidation + LIP ^ n_FT_Fe_Oxidation);

      //15
      FT_Fe_Reduction = k_FT_Fe_Reduction * DFP;

      //16
      FT_Nucleation = k_cat_FT_Nucleation * DFP ^ 2 * FT_cage * (L + rN) / (24 + rN)
        * K_i_FT_Nucleation ^ n_FT_Nucleation / (K_i_FT_Nucleation ^ n_FT_Nucleation + core ^ n_FT_Nucleation);

      //17
      IRPs_Activation = k_cat_IRPs_Activation * IRPs_inactive;

      //18
      IRPs_Degradation = k_cat_IRPs_Degradation * IRPs_active * LIP;

      //cell

      der(FT_cage) * cell = -cell * FT_Degradation + cell * FT_Expression;

      der(core) * cell = 2 * cell * FT_Core_Formation
        + 4 * cell * FT_Nucleation
        - cell * FT_Degradation_Core_Release;

      der(DFP) * cell = -cell * FT_Core_Formation
        - cell * FT_Fe_Reduction
        + cell * FT_Fe_Oxidation
        - 2 * cell * FT_Nucleation;

      der(LIP) * cell = 2 * cell * FT_Fe_Reduction
        - 2 * cell * FT_Fe_Oxidation
        + cell * FT_Degradation_Core_Release;

      der(IRPs_active) * cell = cell * IRPs_Activation
        - cell * IRPs_Degradation;

      der(IRPs_inactive) * cell = -cell * IRPs_Activation
        + cell * IRPs_Degradation;

      annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
            coordinateSystem(preserveAspectRatio=false)));
    end EnterocyteMucosalBlockShortModel11;

    model CellularFerritinIronStorageModel2
      "Cellular ferritin iron storage"

      constant Modelica.Units.SI.Volume V_cell(
        displayUnit = "L") = 1.4e-12 / 1000 "cell volume";

      //Species
      Modelica.Units.SI.MolarConcentration LIP(
        displayUnit = "mol/L",
        start = 1e-05 * 1000) "labile iron pool";
      Modelica.Units.SI.MolarConcentration FT_cage(
        displayUnit = "mol/L",
        start = 5e-09 * 1000) "FT-cage";
      Modelica.Units.SI.MolarConcentration core(
        displayUnit = "mol/L",
        start = 7.5e-06 * 1000) "core";
      Modelica.Units.SI.MolarConcentration DFP(
        displayUnit = "mol/L",
        start = 0) "diferric peroxo complex";

      parameter Integer H = 4 "H subunits";
      parameter Integer L = 24 - H "L subunits";
      parameter Integer atoms_per_cage_transient = 1500 "Transient number of Fe atoms that are stored inside the core of a ferritin cage";
      parameter Integer rN = 50;
      parameter Integer rO = 2;

      parameter Modelica.Units.SI.Frequency k_FTlysis = 1.203e-05;

      parameter Real FT_Expression(
        quantity = "ReactionRate",
        unit = "mol/(m3.s)",
        displayUnit = "mol/(l.s)")
         = 6.015e-14 * 1000;

      //FT degradation
      Real FT_Degradation(
        quantity = "ReactionRate",
        unit = "mol/(m3.s)",
        displayUnit = "mol/(l.s)");

      //FT degradation core release
      Real CoreRelease(
        quantity = "ReactionRate",
        unit = "mol/(m3.s)",
          displayUnit = "mol/(l.s)");

      //Oxidation (2 LIP -> DFP)
      Real Oxidation(
        quantity = "ReactionRate",
        unit = "mol/(m3.s)",
        displayUnit = "mol/(l.s)");

      parameter Modelica.Units.SI.Frequency k_cat_oxidation = 591 "catalytic turnover number";
      parameter Modelica.Units.SI.MolarConcentration K_m_oxidation(
        displayUnit = "mol/L") = 0.35 "Michaelis constant";
      parameter Real n_oxidation = 1.3 "Hill coefficient)";

      //Reduction (DFP -> 2 LIP)
      Real Reduction(
        quantity = "ReactionRate",
        unit = "mol/(m3.s)",
        displayUnit = "mol/(l.s)");

      parameter Modelica.Units.SI.Frequency k_deg = 0.2605 "rate constant";

      //Nucleation (2 DFP -> 4 core)
      Real Nucleation(
        quantity = "ReactionRate",
        unit = "mol/(m3.s)",
        displayUnit = "mol/(l.s)");

      parameter Modelica.Units.SI.Frequency k_cat_nucleation = 5e07 "catalytic turnover number";
      parameter Modelica.Units.SI.MolarConcentration K_i_nucleation(
        displayUnit = "mol/L") = 0.461598 "inhibition constant";
      parameter Integer n_nucleation = 4 "Hill coefficient";

      //Mineralization (DFP -> 2 core)
      Real Mineralization(
        quantity = "ReactionRate",
        unit = "mol/(m3.s)",
        displayUnit = "mol/(l.s)");

      parameter Modelica.Units.SI.Frequency k_cat_mineralization = 0.101564 "catalytic turnover number";
      parameter Modelica.Units.SI.MolarConcentration K_m_mineralization(
        displayUnit = "mol/L") = 5e-03 "Michaelis constant";
      parameter Modelica.Units.SI.MolarConcentration K_i_mineralization(
        displayUnit = "mol/L") = 4.6458 "inhibition constant";
      parameter Integer n_mineralization = 4 "Hill coefficient";
      parameter Integer m_mineralization = 8 "Hill coefficient";

      Modelica.Units.SI.MolarConcentration Fe_total(
        displayUnit = "mol/L") "Fe total";
      Real LIP_rel;


    equation

      Fe_total = core + LIP + DFP;
      LIP_rel = LIP / Fe_total * 100;

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
        coordinateSystem(preserveAspectRatio=false)));

    end CellularFerritinIronStorageModel2;

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
          extent={{100,-64},{120,-44}})));
      Bodylight.Types.Concentration Fe_total;
      //Bodylight.Types.Concentration Fe_total_need;

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

     // BodylightExtension.Types.MolarReactionRate FT_Expression( start = 16.015e-14 * 1000);

    /*    
    parameter Real FT_Expression(
     quantity = "ReactionRate",
     unit = "mol/(m3.s)",
     displayUnit = "mol/(l.s)")
    
    = 16.015e-14 * 1000;
 */

      //FT degradation
      Real FT_Degradation(
        quantity = "ReactionRate",
        unit = "mol/(m3.s)",
        displayUnit = "mol/(l.s)");

      //FT degradation core release
      Real CoreRelease(
        quantity = "ReactionRate",
        unit = "mol/(m3.s)",
          displayUnit = "mol/(l.s)");

      //Oxidation (2 LIP -> DFP)
      Real Oxidation(
        quantity = "ReactionRate",
        unit = "mol/(m3.s)",
        displayUnit = "mol/(l.s)");

      parameter Modelica.Units.SI.Frequency k_cat_oxidation = 591 "catalytic turnover number";
      parameter Modelica.Units.SI.MolarConcentration K_m_oxidation(
        displayUnit = "mol/L") = 0.35 "Michaelis constant";
      parameter Real n_oxidation = 1.3 "Hill coefficient)";

      //Reduction (DFP -> 2 LIP)
      Real Reduction(
        quantity = "ReactionRate",
        unit = "mol/(m3.s)",
        displayUnit = "mol/(l.s)");

      parameter Modelica.Units.SI.Frequency k_deg = 0.2605 "rate constant";

      //Nucleation (2 DFP -> 4 core)
      Real Nucleation(
        quantity = "ReactionRate",
        unit = "mol/(m3.s)",
        displayUnit = "mol/(l.s)");

      parameter Modelica.Units.SI.Frequency k_cat_nucleation = 5e07 "catalytic turnover number";
      parameter Modelica.Units.SI.MolarConcentration K_i_nucleation(
        displayUnit = "mol/L") = 0.461598 "inhibition constant";
      parameter Integer n_nucleation = 4 "Hill coefficient";

      //Mineralization (DFP -> 2 core)
      Real Mineralization(
        quantity = "ReactionRate",
        unit = "mol/(m3.s)",
        displayUnit = "mol/(l.s)");

      parameter Modelica.Units.SI.Frequency k_cat_mineralization = 0.101564 "catalytic turnover number";
      parameter Modelica.Units.SI.MolarConcentration K_m_mineralization(
        displayUnit = "mol/L") = 5e-03 "Michaelis constant";
      parameter Modelica.Units.SI.MolarConcentration K_i_mineralization(
        displayUnit = "mol/L") = 4.6458 "inhibition constant";
      parameter Integer n_mineralization = 4 "Hill coefficient";
      parameter Integer m_mineralization = 8 "Hill coefficient";

    // parameter Modelica.Units.SI.MolarConcentration FT_cage_norm(
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
      LIP = Fe_total_set - (core + DFP);

    equation
      //FT_Expression = 16.015e-14 * 1000;
      // FT_Expression=Ft_expressionIn;

      Fe_in_FT = core + DFP;
      Fe_total = Fe_in_FT + LIP;
      //Fe_total_need = Fe_total_set - Fe_total;

      Fract_Fe_in_Ft = Fe_in_FT / Fe_total;
      Fract_LIP = 1 - Fract_Fe_in_Ft;

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

      //der(FT_cage) = -FT_Degradation + FT_Expression;
      //-FT_Degradation + FT_Expression=0;

      der(LIP) = -2 * Oxidation + 2 * Reduction + CoreRelease;

      der(core) = 2 * Mineralization + 4 * Nucleation - CoreRelease;

      der(DFP) = Oxidation - Mineralization - Reduction - 2 * Nucleation;



      annotation (Diagram(coordinateSystem(extent={{-100,-100},{100,100}})), Icon(
            coordinateSystem(extent={{-100,-100},{100,100}}), graphics={
            Rectangle(
              extent={{102,-100},{-100,100}},
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
              textString="FT_cage")}));
    end FerritinIronStorage;

    model Test_FT_storage
        extends Modelica.Icons.Example;
      Bodylight.Types.Constants.ConcentrationConst Fe_total(k(displayUnit=
              "mol/l") = 0.017499704)
        annotation (Placement(transformation(extent={{-94,70},{-86,78}})));
      Bodylight.Types.Constants.ConcentrationConst FT_cage(k(displayUnit=
              "mmol/l") = 1.33125e-05)
        annotation (Placement(transformation(extent={{-106,20},{-98,28}})));
      FerritinIronStorage ferritinIronStorage
        annotation (Placement(transformation(extent={{-24,-4},{32,52}})));
      FerritinLysis ferritinLysis
        annotation (Placement(transformation(extent={{-84,-20},{-64,0}})));
    equation
      connect(Fe_total.y, ferritinIronStorage.Fe_total_set) annotation (Line(
            points={{-85,74},{-36,74},{-36,43.6},{-25.68,43.6}}, color={0,0,127}));
      connect(FT_cage.y, ferritinLysis.FT_cage) annotation (Line(points={{-97,24},
              {-92,24},{-92,-9.9},{-85.3,-9.9}},           color={0,0,127}));
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
  end models;
  annotation (uses(Modelica(version="4.0.0"), Bodylight(version="1.0")));
end EnterocyteMucosalBlock;
