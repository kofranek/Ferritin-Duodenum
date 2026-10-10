package BodylightExtension

  package Types

    type MolarReactionRate = Real(
      quantity="MolarReactionRate",
      unit = "mol/(m3.s)",
      displayUnit = "mmol/(l.s)"
      );

    type ReactionRateFirstOrder = Real(
      quantity    = "ReactionRateFirstOrder",
      unit        = "1/s"
      );

    type ReactionRateThirdOrder = Real(
      quantity    = "ReactionRateThirdOrder",
      unit        = "m6/(mol2.s)",
      displayUnit = "l2/(mol2.s)"
      );
  end Types;
  annotation (uses(Bodylight(version="1.0")));
end BodylightExtension;

package Modelica "Modelica Standard Library - Version 4.0.0"
extends Modelica.Icons.Package;

  package Icons "Library of icons"
    extends Icons.Package;

    partial model Example "Icon for runnable examples"

      annotation (Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{100,100}}), graphics={
            Ellipse(lineColor = {75,138,73},
                    fillColor={255,255,255},
                    fillPattern = FillPattern.Solid,
                    extent = {{-100,-100},{100,100}}),
            Polygon(lineColor = {0,0,255},
                    fillColor = {75,138,73},
                    pattern = LinePattern.None,
                    fillPattern = FillPattern.Solid,
                    points = {{-36,60},{64,0},{-36,-60},{-36,60}})}), Documentation(info="<html>
<p>This icon indicates an example. The play button suggests that the example can be executed.</p>
</html>"));
    end Example;

    partial package Package "Icon for standard packages"
      annotation (Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{100,100}}), graphics={
            Rectangle(
              lineColor={200,200,200},
              fillColor={248,248,248},
              fillPattern=FillPattern.HorizontalCylinder,
              extent={{-100.0,-100.0},{100.0,100.0}},
              radius=25.0),
            Rectangle(
              lineColor={128,128,128},
              extent={{-100.0,-100.0},{100.0,100.0}},
              radius=25.0)}), Documentation(info="<html>
<p>Standard package icon.</p>
</html>"));
    end Package;

    partial package SourcesPackage "Icon for packages containing sources"
      extends Modelica.Icons.Package;
      annotation (Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,
                -100},{100,100}}), graphics={
            Polygon(origin={23.3333,0.0},
              fillColor={128,128,128},
              pattern=LinePattern.None,
              fillPattern=FillPattern.Solid,
              points={{-23.333,30.0},{46.667,0.0},{-23.333,-30.0}}),
            Rectangle(
              fillColor = {128,128,128},
              pattern = LinePattern.None,
              fillPattern = FillPattern.Solid,
              extent = {{-70,-4.5},{0,4.5}})}),
                                Documentation(info="<html>
<p>This icon indicates a package which contains sources.</p>
</html>"));
    end SourcesPackage;
    annotation (Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,
                -100},{100,100}}), graphics={Polygon(
              origin={-8.167,-17},
              fillColor={128,128,128},
              pattern=LinePattern.None,
              fillPattern=FillPattern.Solid,
              points={{-15.833,20.0},{-15.833,30.0},{14.167,40.0},{24.167,20.0},{
                  4.167,-30.0},{14.167,-30.0},{24.167,-30.0},{24.167,-40.0},{-5.833,
                  -50.0},{-15.833,-30.0},{4.167,20.0},{-5.833,20.0}},
              smooth=Smooth.Bezier), Ellipse(
              origin={-0.5,56.5},
              fillColor={128,128,128},
              pattern=LinePattern.None,
              fillPattern=FillPattern.Solid,
              extent={{-12.5,-12.5},{12.5,12.5}})}), Documentation(info="<html>
<p>This package contains definitions for the graphical layout of components which may be used in different libraries. The icons can be utilized by inheriting them in the desired class using &quot;extends&quot; or by directly copying the &quot;icon&quot; layer.</p>

<h4>Main Authors</h4>

<dl>
<dt><a href=\"http://www.robotic.dlr.de/Martin.Otter/\">Martin Otter</a></dt>
    <dd>Deutsches Zentrum fuer Luft und Raumfahrt e.V. (DLR)</dd>
    <dd>Oberpfaffenhofen</dd>
    <dd>Postfach 1116</dd>
    <dd>D-82230 Wessling</dd>
    <dd>email: <a href=\"mailto:Martin.Otter@dlr.de\">Martin.Otter@dlr.de</a></dd>
<dt>Christian Kral</dt>

    <dd>  <a href=\"https://christiankral.net/\">Electric Machines, Drives and Systems</a><br>
</dd>
    <dd>1060 Vienna, Austria</dd>
    <dd>email: <a href=\"mailto:dr.christian.kral@gmail.com\">dr.christian.kral@gmail.com</a></dd>
<dt>Johan Andreasson</dt>
    <dd><a href=\"https://www.modelon.com/\">Modelon AB</a></dd>
    <dd>Ideon Science Park</dd>
    <dd>22370 Lund, Sweden</dd>
    <dd>email: <a href=\"mailto:johan.andreasson@modelon.se\">johan.andreasson@modelon.se</a></dd>
</dl>

<p>
Copyright &copy; 1998-2020, Modelica Association and contributors
</p>
</html>"));
  end Icons;

  package Units "Library of type and unit definitions"
    extends Modelica.Icons.Package;

    package SI "Library of SI unit definitions"
      extends Modelica.Icons.Package;

      type Frequency = Real (final quantity="Frequency", final unit="Hz");

      type Concentration = Real (final quantity="Concentration", final unit=
              "mol/m3");
      annotation (Icon(graphics={Text(
              extent={{-80,80},{80,-78}},
              textColor={128,128,128},
              fillColor={128,128,128},
              fillPattern=FillPattern.None,
              fontName="serif",
              textString="SI",
              textStyle={TextStyle.Italic})}),
                                       Documentation(info="<html>
<p>This package provides predefined types based on the international standard
on units.
</p>
<p>
For an introduction to the conventions used in this package, have a look at:
<a href=\"modelica://Modelica.Units.UsersGuide.Conventions\">Conventions</a>.
</p>
</html>"));
    end SI;
    annotation (Icon(graphics={
        Polygon(
          fillColor = {128,128,128},
          pattern = LinePattern.None,
          fillPattern = FillPattern.Solid,
          points = {{-80,-40},{-80,-40},{-55,50},{-52.5,62.5},{-65,60},{-65,65},{-35,77.5},{-32.5,60},{-50,0},{-50,0},{-30,15},{-20,27.5},{-32.5,27.5},{-32.5,27.5},{-32.5,32.5},{-32.5,32.5},{2.5,32.5},{2.5,32.5},{2.5,27.5},{2.5,27.5},{-7.5,27.5},{-30,7.5},{-30,7.5},{-25,-25},{-17.5,-28.75},{-10,-25},{-5,-26.25},{-5,-32.5},{-16.25,-41.25},{-31.25,-43.75},{-40,-33.75},{-45,-5},{-45,-5},{-52.5,-10},{-52.5,-10},{-60,-40},{-60,-40}},
          smooth = Smooth.Bezier),
        Polygon(
          fillColor = {128,128,128},
          pattern = LinePattern.None,
          fillPattern = FillPattern.Solid,
          points = {{87.5,30},{62.5,30},{62.5,30},{55,33.75},{36.25,35},{16.25,25},{7.5,6.25},{11.25,-7.5},{22.5,-12.5},{22.5,-12.5},{6.25,-22.5},{6.25,-35},{16.25,-38.75},{16.25,-38.75},{21.25,-41.25},{21.25,-41.25},{45,-48.75},{47.5,-61.25},{32.5,-70},{12.5,-65},{7.5,-51.25},{21.25,-41.25},{21.25,-41.25},{16.25,-38.75},{16.25,-38.75},{6.25,-41.25},{-6.25,-50},{-3.75,-68.75},{30,-76.25},{65,-62.5},{63.75,-35},{27.5,-26.25},{22.5,-20},{27.5,-15},{27.5,-15},{30,-7.5},{30,-7.5},{27.5,-2.5},{28.75,11.25},{36.25,27.5},{47.5,30},{53.75,22.5},{51.25,8.75},{45,-6.25},{35,-11.25},{30,-7.5},{30,-7.5},{27.5,-15},{27.5,-15},{43.75,-16.25},{65,-6.25},{72.5,10},{70,20},{70,20},{80,20}},
          smooth = Smooth.Bezier)}), Documentation(info="<html>
<p>This package provides predefined types, such as <em>Mass</em>,
<em>Angle</em>, <em>Time</em>, based on the international standard
on units, e.g.,
</p>

<blockquote><pre>
<strong>type</strong> Angle = Real(<strong>final</strong> quantity = \"Angle\",
                  <strong>final</strong> unit     = \"rad\",
                  displayUnit   = \"deg\");
</pre></blockquote>

<p>
Some of the types are derived SI units that are utilized in package Modelica
(such as ComplexCurrent, which is a complex number where both the real and imaginary
part have the SI unit Ampere).
</p>

<p>
Furthermore, conversion functions from non SI-units to SI-units and vice versa
are provided in subpackage
<a href=\"modelica://Modelica.Units.Conversions\">Conversions</a>.
</p>

<p>
For an introduction how units are used in the Modelica Standard Library
with package Units, have a look at:
<a href=\"modelica://Modelica.Units.UsersGuide.HowToUseUnits\">How to use Units</a>.
</p>

<p>
Copyright &copy; 1998-2020, Modelica Association and contributors
</p>
</html>",   revisions="<html>
<ul>
<li><em>May 25, 2011</em> by Stefan Wischhusen:<br>Added molar units for energy and enthalpy.</li>
<li><em>Jan. 27, 2010</em> by Christian Kral:<br>Added complex units.</li>
<li><em>Dec. 14, 2005</em> by <a href=\"http://www.robotic.dlr.de/Martin.Otter/\">Martin Otter</a>:<br>Add User&#39;s Guide and removed &quot;min&quot; values for Resistance and Conductance.</li>
<li><em>October 21, 2002</em> by <a href=\"http://www.robotic.dlr.de/Martin.Otter/\">Martin Otter</a> and Christian Schweiger:<br>Added new package <strong>Conversions</strong>. Corrected typo <em>Wavelenght</em>.</li>
<li><em>June 6, 2000</em> by <a href=\"http://www.robotic.dlr.de/Martin.Otter/\">Martin Otter</a>:<br>Introduced the following new types<br>type Temperature = ThermodynamicTemperature;<br>types DerDensityByEnthalpy, DerDensityByPressure, DerDensityByTemperature, DerEnthalpyByPressure, DerEnergyByDensity, DerEnergyByPressure<br>Attribute &quot;final&quot; removed from min and max values in order that these values can still be changed to narrow the allowed range of values.<br>Quantity=&quot;Stress&quot; removed from type &quot;Stress&quot;, in order that a type &quot;Stress&quot; can be connected to a type &quot;Pressure&quot;.</li>
<li><em>Oct. 27, 1999</em> by <a href=\"http://www.robotic.dlr.de/Martin.Otter/\">Martin Otter</a>:<br>New types due to electrical library: Transconductance, InversePotential, Damping.</li>
<li><em>Sept. 18, 1999</em> by <a href=\"http://www.robotic.dlr.de/Martin.Otter/\">Martin Otter</a>:<br>Renamed from SIunit to SIunits. Subpackages expanded, i.e., the SIunits package, does no longer contain subpackages.</li>
<li><em>Aug 12, 1999</em> by <a href=\"http://www.robotic.dlr.de/Martin.Otter/\">Martin Otter</a>:<br>Type &quot;Pressure&quot; renamed to &quot;AbsolutePressure&quot; and introduced a new type &quot;Pressure&quot; which does not contain a minimum of zero in order to allow convenient handling of relative pressure. Redefined BulkModulus as an alias to AbsolutePressure instead of Stress, since needed in hydraulics.</li>
<li><em>June 29, 1999</em> by <a href=\"http://www.robotic.dlr.de/Martin.Otter/\">Martin Otter</a>:<br>Bug-fix: Double definition of &quot;Compressibility&quot; removed and appropriate &quot;extends Heat&quot; clause introduced in package SolidStatePhysics to incorporate ThermodynamicTemperature.</li>
<li><em>April 8, 1998</em> by <a href=\"http://www.robotic.dlr.de/Martin.Otter/\">Martin Otter</a> and Astrid Jaschinski:<br>Complete ISO 31 chapters realized.</li>
<li><em>Nov. 15, 1997</em> by <a href=\"http://www.robotic.dlr.de/Martin.Otter/\">Martin Otter</a> and Hubertus Tummescheit:<br>Some chapters realized.</li>
</ul>
</html>"));
  end Units;
annotation (
preferredView="info",
version="4.0.0",
versionDate="2020-06-04",
dateModified = "2020-06-04 11:00:00Z",
revisionId="6626538a2 2020-06-04 19:56:34 +0200",
uses(Complex(version="4.0.0"), ModelicaServices(version="4.0.0")),
conversion(
 from(version={"3.0", "3.0.1", "3.1", "3.2", "3.2.1", "3.2.2", "3.2.3"}, script="modelica://Modelica/Resources/Scripts/Conversion/ConvertModelica_from_3.2.3_to_4.0.0.mos")),
Icon(coordinateSystem(extent={{-100.0,-100.0},{100.0,100.0}}), graphics={
  Polygon(
    origin={-6.9888,20.048},
    pattern=LinePattern.None,
    fillPattern=FillPattern.Solid,
    points={{-93.0112,10.3188},{-93.0112,10.3188},{-73.011,24.6},{-63.011,31.221},{-51.219,36.777},{-39.842,38.629},{-31.376,36.248},{-25.819,29.369},{-24.232,22.49},{-23.703,17.463},{-15.501,25.135},{-6.24,32.015},{3.02,36.777},{15.191,39.423},{27.097,37.306},{32.653,29.633},{35.035,20.108},{43.501,28.046},{54.085,35.19},{65.991,39.952},{77.897,39.688},{87.422,33.338},{91.126,21.696},{90.068,9.525},{86.099,-1.058},{79.749,-10.054},{71.283,-21.431},{62.816,-33.337},{60.964,-32.808},{70.489,-16.14},{77.368,-2.381},{81.072,10.054},{79.749,19.05},{72.605,24.342},{61.758,23.019},{49.587,14.817},{39.003,4.763},{29.214,-6.085},{21.012,-16.669},{13.339,-26.458},{5.401,-36.777},{-1.213,-46.037},{-6.24,-53.446},{-8.092,-52.387},{-0.684,-40.746},{5.401,-30.692},{12.81,-17.198},{19.424,-3.969},{23.658,7.938},{22.335,18.785},{16.514,23.283},{8.047,23.019},{-1.478,19.05},{-11.267,11.113},{-19.734,2.381},{-29.259,-8.202},{-38.519,-19.579},{-48.044,-31.221},{-56.511,-43.392},{-64.449,-55.298},{-72.386,-66.939},{-77.678,-74.612},{-79.53,-74.083},{-71.857,-61.383},{-62.861,-46.037},{-52.278,-28.046},{-44.869,-15.346},{-38.784,-2.117},{-35.344,8.731},{-36.403,19.844},{-42.488,23.813},{-52.013,22.49},{-60.744,16.933},{-68.947,10.054},{-76.884,2.646},{-93.0112,-12.1707},{-93.0112,-12.1707}},
    smooth=Smooth.Bezier),
  Ellipse(
    origin={40.8208,-37.7602},
    fillColor={161,0,4},
    pattern=LinePattern.None,
    fillPattern=FillPattern.Solid,
    extent={{-17.8562,-17.8563},{17.8563,17.8562}})}),
Documentation(info="<html>
<p>
<img src=\"modelica://Modelica/Resources/Images/Logos/Modelica_Libraries.svg\" width=\"250\">
</p>

<p>
The package <strong>Modelica&reg;</strong> is a <strong>standardized</strong> and <strong>free</strong> package
that is developed by the \"<strong>Modelica Association Project - Libraries</strong>\".</p>
<p>
Its development is coordinated with the Modelica&reg; language from the
Modelica Association, see <a href=\"https://www.Modelica.org\">https://www.Modelica.org</a>.
It is also called <strong>Modelica Standard Library</strong>.
It provides model components in many domains that are based on
standardized interface definitions. Some typical examples are shown
in the next figure:
</p>

<p>
<img src=\"modelica://Modelica/Resources/Images/UsersGuide/ModelicaLibraries.png\">
</p>

<p>
For an introduction, have especially a look at:
</p>
<ul>
<li> <a href=\"modelica://Modelica.UsersGuide.Overview\">Overview</a>
  provides an overview of the Modelica Standard Library
  inside the <a href=\"modelica://Modelica.UsersGuide\">User's Guide</a>.</li>
<li><a href=\"modelica://Modelica.UsersGuide.ReleaseNotes\">Release Notes</a>
 summarizes the changes of new versions of this package.</li>
<li> <a href=\"modelica://Modelica.UsersGuide.Contact\">Contact</a>
  lists the contributors of the Modelica Standard Library.</li>
<li> The <strong>Examples</strong> packages in the various libraries, demonstrate
  how to use the components of the corresponding sublibrary.</li>
</ul>

<p>
This version of the Modelica Standard Library consists of
</p>
<ul>
<li><strong>1417</strong> component models and blocks,</li>
<li><strong>512</strong> example models, and</li>
<li><strong>1219</strong> functions</li>
</ul>
<p>
that are directly usable (= number of public, non-partial, non-internal and non-obsolete classes). It is fully compliant
to <a href=\"https://modelica.org/documents/ModelicaSpec34.pdf\">Modelica Specification version 3.4</a>
and it has been tested with Modelica tools from different vendors.
</p>

<p>
<strong>Licensed by the Modelica Association under the 3-Clause BSD License</strong><br>
Copyright &copy; 1998-2020, Modelica Association and <a href=\"modelica://Modelica.UsersGuide.Contact\">contributors</a>.
</p>

<p>
<em>This Modelica package is <u>free</u> software and the use is completely at <u>your own risk</u>; it can be redistributed and/or modified under the terms of the 3-Clause BSD license. For license conditions (including the disclaimer of warranty) visit <a href=\"https://modelica.org/licenses/modelica-3-clause-bsd\">https://modelica.org/licenses/modelica-3-clause-bsd</a>.</em>
</p>

<p>
<strong>Modelica&reg;</strong> is a registered trademark of the Modelica Association.
</p>
</html>"));
end Modelica;

package Bodylight "Modelica library for Physiology (version 1.0)"
 extends Modelica.Icons.Package;

  package Types "Physiological units with nominals"
    extends Modelica.Icons.Package;

    package Constants
      extends Modelica.Icons.SourcesPackage;

    block ConcentrationConst "Constant signal of type Concentration"
     parameter Types.Concentration k "Constant Concentration output value";
          RealIO.ConcentrationOutput y "Concentration constant"
        annotation (Placement(transformation(extent={{40,-10},{60,10}}),
                    iconTransformation(extent={{40,-10},{60,10}})));
    equation
          y=k;
      annotation (defaultComponentName="concentration",
                 Diagram(coordinateSystem(extent={{-40,-40},{40,40}})), Icon(
            coordinateSystem(extent={{-40,-40},{40,40}}, preserveAspectRatio=false),
                graphics={
            Rectangle(extent={{-40,40},{40,-40}},
              lineColor={0,0,0},
                  radius=10,
              fillColor={236,236,236},
                              fillPattern=FillPattern.Solid),
            Text( extent={{-100,-44},{100,-64}},
              lineColor={0,0,0},
                      fillColor={236,236,236},
              fillPattern=FillPattern.Solid,
                  textString="%name"),
            Text(         extent={{-40,10},{40,-10}},
              lineColor={0,0,0},
                  fillColor={236,236,236},
              fillPattern=FillPattern.Solid,
                      textString="Const")}));
    end ConcentrationConst;
    end Constants;

    package RealIO
      extends Modelica.Icons.Package;

      connector ConcentrationInput = input Concentration
        "input Concentration as connector"
        annotation (defaultComponentName="concentration",
        Icon(graphics={Polygon(
                points={{-100,100},{100,0},{-100,-100},{-100,100}},
                lineColor={0,0,127},
                fillColor={0,0,127},
                fillPattern=FillPattern.Solid)},
             coordinateSystem(extent={{-100,-100},{100,100}}, preserveAspectRatio=true, initialScale=0.2)),
        Diagram(coordinateSystem(
              preserveAspectRatio=true, initialScale=0.2,
              extent={{-100,-100},{100,100}},
              grid={1,1}), graphics={Polygon(
                points={{0,50},{100,0},{0,-50},{0,50}},
                lineColor={0,0,127},
                fillColor={0,0,127},
                fillPattern=FillPattern.Solid), Text(
                extent={{-10,85},{-10,60}},
                lineColor={0,0,127},
                textString="%name")}),
          Documentation(info="<html>
    <p>
    Connector with one input signal of type Concentration.
    </p>
    </html>"));

      connector ConcentrationOutput = output Concentration
        "output Concentration as connector"
        annotation (defaultComponentName="concentration",
        Icon(coordinateSystem(
              preserveAspectRatio=true,
              extent={{-100,-100},{100,100}},
              grid={1,1}), graphics={Polygon(
                points={{-100,100},{100,0},{-100,-100},{-100,100}},
                lineColor={0,0,127},
                fillColor={255,255,255},
                fillPattern=FillPattern.Solid)}),
        Diagram(coordinateSystem(
              preserveAspectRatio=true,
              extent={{-100,-100},{100,100}},
              grid={1,1}), graphics={Polygon(
                points={{-100,50},{0,0},{-100,-50},{-100,50}},
                lineColor={0,0,127},
                fillColor={255,255,255},
                fillPattern=FillPattern.Solid), Text(
                extent={{30,110},{30,60}},
                lineColor={0,0,127},
                textString="%name")}),
          Documentation(info="<html>
  <p>
  Connector with one output signal of type Real.
  </p>
  </html>"));

      connector FractionOutput = output Fraction "output Fraction as connector"
        annotation (defaultComponentName="fraction",
        Icon(coordinateSystem(
              preserveAspectRatio=true,
              extent={{-100,-100},{100,100}},
              grid={1,1}), graphics={Polygon(
                points={{-100,100},{100,0},{-100,-100},{-100,100}},
                lineColor={0,0,127},
                fillColor={255,255,255},
                fillPattern=FillPattern.Solid)}),
        Diagram(coordinateSystem(
              preserveAspectRatio=true,
              extent={{-100,-100},{100,100}},
              grid={1,1}), graphics={Polygon(
                points={{-100,50},{0,0},{-100,-50},{-100,50}},
                lineColor={0,0,127},
                fillColor={255,255,255},
                fillPattern=FillPattern.Solid), Text(
                extent={{30,110},{30,60}},
                lineColor={0,0,127},
                textString="%name")}),
          Documentation(info="<html>
  <p>
  Connector with one output signal of type Real.
  </p>
  </html>"));
    end RealIO;

    type Frequency = Modelica.Units.SI.Frequency (
                                                displayUnit="1/min");

    replaceable type Concentration = Modelica.Units.SI.Concentration(displayUnit="mmol/l", min=0) constrainedby Real;

    type Fraction = Real(final quantity="Fraction",final unit="1", displayUnit="%", nominal=1e-2);
    annotation (Documentation(revisions="<html>
<p>Copyright (c) 2008-2015, Marek Matej&aacute;k, Charles University in Prague </p>
<p>All rights reserved. </p>
<p>Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met: </p>
<ol>
<li>Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer. </li>
<li>Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution. </li>
<li>Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission. </li>
</ol>
<p>THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS &quot;AS IS&quot; AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.</p>
</html>",   info="<html>
<p>The main problem of medical research, articles, and experiments is using obscure units from medicine, pharmacology, biology and non-physics disciplines. The Physiolibrary fulfills the Modelica ideal of using SI units as the main unit for each variable, and the previously described physiological units are also implemented as the displayUnits for each variable. Using these displayUnits the user sets and sees the &quot;physiological&quot; values. The implementation can also be joined to any unit-correct Modelica models and physical equations without crashing due to unit incompatibilities. The unit support of Physiolibrary is so strong that you even can chose the right unit-typed &ldquo;input real&rdquo;/&rdquo;output real&rdquo; from the library package Types.RealIO. As can be expected, only the non-specific packages States and Blocks in the Physiolibrary have variables without units.</p>
</html>"));
  end Types;
  annotation (preferredView="info",
version="1.0",
versionDate="2025-02-20",
dateModified = "2025-02-20 12:49:00Z",
revisionId="",
uses(Modelica(version="4.0.0")),
conversion(
      from(version="Physiolibrary-2.4.1", script="Resources/Scripts/ConvertFromPhysiolibrary_2.4.1.mos")),
  Documentation(revisions="<html>
<p>Copyright (c) 2008-2015, Marek Matej&aacute;k, Charles University in Prague </p>
<p>All rights reserved. </p>
<p>Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met: </p>
<ol>
<li>Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer. </li>
<li>Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution. </li>
<li>Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission. </li>
</ol>
<p>THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS &quot;AS IS&quot; AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.</p>
</html>", info="<html>
<ul>
<li>Web pages: <a href=\"http://www.physiolibrary.org\">www.physiolibrary.org</a></li>
<li><a href=\"modelica://Physiolibrary.UsersGuide.Overview\">Overview</a></li>
<li><a href=\"modelica://Physiolibrary.UsersGuide.Connectors\">Connectors</a></li>
<li><a href=\"modelica://Physiolibrary.UsersGuide.Contact\">Contact</a></li>
</ul>

<p><br/>The origin of this Modelica Physiolibrary was in the first version of our HumMod Golem Edition model implementation, where it was called HumMod.Library. As the successors of Guyton&apos;s Medical Physiology School write, the original HumMod model is &ldquo;The best, most complete, mathematical model of human physiology ever created&rdquo;.</p>
<p>We are also developing many types of smaller physiological models for use in medical education, so it was essential to separate this library from our HumMod Modelica implementation. This separation improves the quality of the next HumMod release and provides a useful Modelica library to modelers in this bioscience.</p>
<p>The library contains only carefully-chosen elementary physiological laws, which are the basis of more complex physiological processes. For example from only two type of blocks (Chemical.ChemicalReaction and Chemical.Substance) it is possible to compose the allosteric transitions or the Michaelis-Menten equation.</p>
<p><br/>Physiology is a very progressive discipline, that examines how the living body works. And it is no surprise that all processes in the human body are driven by physical laws of nature. The great challenge is to marry old empirical experiments with the &ldquo;new&rdquo; physical principles. Many teams and projects in the word deal with this formalization of physiology, for example: Physiome, SBML, EuroPhysiome, VPH, CellML etc. It is our hope that this library helps this unflagging effort of physiologists to exactly describe the processes.</p>
</html>"));
end Bodylight;

package EnterocyteMucosalBlock "Enterocyte mucosal block"

  package models

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
    /*
 der(LIP)=0;
 der(core)=0;
 der(DFP)=0;
*/

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
              textString="FT_cage"),
            Text(
              extent={{22,54},{92,30}},
              textColor={28,108,200},
              textString="LIP")}));
    end FerritinIronStorage;

    model Test_FT_storage
        extends Modelica.Icons.Example;
      Bodylight.Types.Constants.ConcentrationConst Fe_total(k(displayUnit=
              "mmol/l") = 0.06)
        annotation (Placement(transformation(extent={{-88,40},{-80,48}})));
      Bodylight.Types.Constants.ConcentrationConst FT_cage(k(displayUnit=
              "mmol/l") = 3e-05)
        annotation (Placement(transformation(extent={{-88,20},{-80,28}})));
      FerritinIronStorage ferritinIronStorage
        annotation (Placement(transformation(extent={{-24,-4},{32,52}})));
    equation
      connect(Fe_total.y, ferritinIronStorage.Fe_total_set) annotation (Line(
            points={{-79,44},{-52.34,44},{-52.34,43.6},{-25.68,43.6}},
                                                                 color={0,0,127}));
      connect(ferritinIronStorage.FT_cage, FT_cage.y)
        annotation (Line(points={{-25.68,24},{-79,24}}, color={0,0,127}));
      annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
            coordinateSystem(preserveAspectRatio=false)));
    end Test_FT_storage;
  end models;
  annotation (uses(Modelica(version="4.0.0"), Bodylight(version="1.0")));
end EnterocyteMucosalBlock;
model EnterocyteMucosalBlock_models_Test_FT_storage
 extends EnterocyteMucosalBlock.models.Test_FT_storage;
  annotation(experiment(
    StopTime=20000,
    __Dymola_NumberOfIntervals=500,
    Tolerance=0.0001,
    __Dymola_Algorithm="dassl"));
end EnterocyteMucosalBlock_models_Test_FT_storage;
