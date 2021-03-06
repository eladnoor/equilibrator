----------------------------------
Energy Currency and ATP Hydrolysis 
----------------------------------

.. figure:: _static/_images/atp.png
   :alt: Adenosine Triphosphate (ATP)
   :align: center

   Adenosine Triphosphate (ATP)

ATP is considered the energy "currency" of the cell. But how much is that currency worth? How much energy is in an ATP hydrolysis reaction? The hydrolysis of one phosphoanhydride bond to form ADP and phosphate

|atp_hydrolysis|_

.. |atp_hydrolysis| replace:: ATP + H\ :sub:`2`\ O ⇌ ADP + Pi
.. _atp_hydrolysis: http://equilibrator.weizmann.ac.il/search?query=ATP+%2B+Water+%3C%3D%3E+ADP+%2B+Phosphate

has a Δ\ :sub:`r`\ G'° of about -26 kJ / mol. [#atp1]_ That is, if you hydrolyze 1 mol of ATP to ADP in standard conditions, 26 kJ of energy that is usable for work is released. Conversely, it takes 26 kJ of work to form 1 mol of ATP from ADP in standard conditions.

Keep in mind that Δ\ :sub:`r`\ G'° is defined in standard conditions (1 molar concentrations) while `biological concentrations are much closer to 1 millimolar <http://book.bionumbers.org/what-are-the-concentrations-of-free-metabolites-in-cells/>`__ (1 mM). We can calculate the Δ\ :sub:`r`\ G' in these conditions using one of the most basic equations of biochemical thermodynamics

.. math::
	\begin{eqnarray}
	\Delta_r G' &=& \Delta_r G'^{\circ} + RT \ln{Q} \\
	&=& \Delta_r G'^{\circ} + RT \ln{\left( \frac{[ADP][Pi]}{[ATP][H_2O]} \right)}
	\end{eqnarray}

Here Q is the ratio of product concentrations to substrate concentrations for the ATP hydrolysis reaction - concentrations are denoted by square brackets like [ATP]. [#atp2]_ Q is called the "reaction quotient" or sometimes the "mass-action ratio" in textbooks. Because biology takes place in aqueous solution, the water concentration [H\ :sub:`2`\ O] is roughly 55 M. [#atp3]_ Typically, however, the "activity" of water is set to 1 molar (`Alberty et al. 2011 <refs.html>`_). If we set the other three concentrations to 10\ :sup:`-3` molar we get

.. math::
	\begin{eqnarray}
	\Delta_r G' &=& \Delta_r G'^{\circ} + RT \ln{\left( \frac{10^{-3} M \times 10^{-3} M}{10^{-3} M \times 1 M} \right)} \\
	&=& \Delta_r G'^{\circ} + RT \ln{10^{-3}} \\
	&\approx& \Delta_r G'^{\circ} - 17.1 \frac{kJ}{mol}
	\end{eqnarray}

Since R = 8.315 x 10\ :sup:`-3` kJ/mol/K and we assume a temperature of T = 298.15 K (25 °C). [#atp4]_ This calculation gives a Δ\ :sub:`r`\ G' value around -44 kJ/mol. Try using eQuilibrator to check this value. [#atp5]_ We call this value the Δ\ :sub:`r`\ G'\ :sup:`m` - the Δ\ :sub:`r`\ G' value when all reactants have 1 mM concentration (except H\ :sub:`2`\ O). We will often use this value instead of Δ\ :sub:`r`\ G'° in this document and on eQuilibrator because 1 mM concentrations are much more representative of small molecule concentrations inside cells than the usual 1 M standard.

By way of comparison, reactions hydrolyzing peptide and C-C bonds are much less exergonic than ATP hydrolysis. For example, consider the hydrolysis of the dipeptide glyclglycine

|glygly|_

.. |glygly| replace:: Glycylglycine + H\ :sub:`2`\ O = 2 Glycine
.. _glygly: http://equilibrator.weizmann.ac.il/search?query=Glycylglycine+%2B+water+%3D+2+Glycine

This reaction has Δ\ :sub:`r`\ G'° ≈ 0 kJ / mol and Δ\ :sub:`r`\ G'\ :sup:`m` ≈ -16 kJ / mol. Similarly, breaking a C-C bond in a sugar, as in the `fructose bisphosphate aldolase reaction <http://equilibrator.weizmann.ac.il/reaction?reactantsId=C00111&reactantsCoeff=1&reactantsName=Glycerone%20phosphate&reactantsPhase=aqueous&reactantsConcentration=0.001&reactantsId=C00118&reactantsCoeff=1&reactantsName=D-Glyceraldehyde%203-phosphate&reactantsPhase=aqueous&reactantsConcentration=0.001&reactantsId=C00354&reactantsCoeff=-1&reactantsName=D-Fructose-1,6-bisphosphate&reactantsPhase=aqueous&reactantsConcentration=0.001&ph=7.000000&pmg=14.000000&ionic_strength=0.100000&e_reduction_potential=0.000000&max_priority=0&mode=BA&query=D-Fructose-1%2C6-bisphosphate%20%3D%20Glycerone%20phosphate%20%2B%20D-Glyceraldehyde%203-phosphate>`__, is intrinsically unfavorable with a Δ\ :sub:`r`\ G'\ :sup:`m` ≈ +3 kJ / mol.

We can also convert the energy of ATP hydrolysis into other units, for example units of mechanical force relevant on the molecular scale: piconewton nanometers (pN nm). Looking up the `definitions of these units <https://en.wikipedia.org/wiki/KT_(energy)>`__, we find that 1 pN nm = 10\ :sup:`-24` kJ. Hydrolysing 1 mol of ATP yields ~45 kJ of energy (assuming mM concentrations) and so one ATP hydrolysis yields 

.. math::
	\frac{45}{N_A} kJ \times 10^{24} \frac{pN \times nm}{kJ}

Where N\ :sub:`A` is Avogadro's number. Altogether, this gives ~73 pN nm per ATP hydrolysis (`Phillips et al., 2009 <refs.html>`_). But how much energy is that? Well, it’s more than enough to `compact DNA <http://bionumbers.hms.harvard.edu/bionumber.aspx?id=103125>`__ into a viral capsid, for example. In other words, ATP stores an amount of energy that can be converted into a usable amount of work at the molecular scale. ATP acts as a bridge between the chemical and physical realms - between energy metabolism and the mechanical work of living.

Moreover, ATP serves a crucial role in metabolism. As we will see throughout this document, cells often want to perform metabolic transformations that are not intrinsically favorable (reactions that are "uphill" energetically). Through evolution and natural selection, cells "learned" how to make these transformations "go" by coupling them to highly favorable reactions like the hydrolysis of ATP. [#atp6]_ So remember the energy scale of ATP hydrolysis - Δ\ :sub:`r`\ G'\ :sup:`m` ≈ 45 kJ/mol - since we will often use that scale to determine how many ATPs will be needed to make some biochemical reaction go in the direction the cell wants it to.

.. [#atp1] That's about 6 kCal / mol for the chemists in the room.
.. [#atp2] Notice that Q is unitless for this reaction because it has the same number of substrates and products. In fact Q *must* unitless so that we can take its logarithm. 
.. [#atp3] Try to calculate this from the density of pure water - remember that 1 L of water weighs 1 kg at standard temperature and pressure and H\ :sub:`2`\ O has a molecular mass of 18 g/mol. 
.. [#atp4] We assume T = 25 °C = 298.15 K throughout eQuilibrator. We don't curretly have enough information about the temperature dependence of Δ\ :sub:`r`\ G'° to remove this assumption.
.. [#atp5] Hint: click on the ATP hydrolysis reaction above.
.. [#atp6] This coupling between catalysis and ATP hydrolysis can occur in various ways inside protein enzymes. If you want to learn more about how that process works, read about "allostery" in your favorite textbook.
