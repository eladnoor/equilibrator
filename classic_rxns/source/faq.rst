Frequently Asked Questions
==========================================================

How do I search for a compound or enzyme?
----------------------------------------------------------

You can search for a compound or enzyme by entering their common name into the search box. For example, try "`Succinyl-CoA </search?query=Succinyl-CoA>`_" or "`pyruvate formate-lyase </search?query=pyruvate+formate-lyase>`_." As you type, you should see a list of potential matches which should give you a sense of what eQuilibrator thinks you are looking for. To the right of each suggestion you will see they type of the match: "Compound" or "Enzyme". Complete your search by selecting a compound or enzyme from the list or typing the name in full. When you've finished typing your query, click 'Search' or press 'Enter' to see search results.

The search results page shows the top matches for your query. Results are shown in two sections: one for compounds and another for enzymes. If you're interested in a compound, click on it to see more details and to find Δ\ :sub:`f`\ G. If you're interested in an enzyme, click on the enzyme name to see detailed information about it. If you want to explore the thermodynamics of an enzyme-catalyzed reaction, choose a reaction listed under the enzyme name to get a Δ\ :sub:`r`\ G estimation.

How do I search for a reaction?
----------------------------------------------------------

You can search for a reaction by entering an (almost) free-text description of the reaction. For example, you can type "`glucose = 2 ethanol + 2 CO\ :sub:`2` <search?query=glucose+%3D+2+ethanol+%2B+2+CO2>`_".

In general, write your reaction using " = " or " <=> " to separate the substrates from the products. The spaces around " = " are required. Always use spaces before and after compound names, "+" and stoichiometric coefficients. Also note that you can find reactions through the enzymes that catalyze them. For example, instead of typing in the full reaction "`Fructose 1,6-bisphosphate = Glycerone phosphate + Glyceraldehyde 3-phosphate </search?query=glucose+%3D+2+ethanol+%2B+2+CO2>`_" you can instead search for "`aldolase </search?query=aldolase>`_" or "`fructose bisphosphate aldolase </search?query=fructose+bisphosphate+aldolase>`_."

What are "standard conditions?"
----------------------------------------------------------

"Standard conditions" implies that all substrates and products of a reaction are at 1 M concentration. Standard conditions are denoted by a degree symbol "°" as in Δ\ :sub:`r`\ G'°. Standard conditions are simply a reference point - reaction rarely take place in standard conditions, if ever. If concentrations in your reaction system are not 1 M, it will affect the Δ\ :sub:`r`\ G'. You should not assume that the Δ\ :sub:`r`\ G'° value is meaningful in real-world scenarios - realistic biological concentrations can alter the Δ\ :sub:`r`\ G' from Δ\ :sub:`r`\ G'° by tens of kJ/mol. 

In general, Δ\ :sub:`r`\ G'\ :sup:`m` is a better reference than Δ\ :sub:`r`\ G'° for biological systems. You can read more about Δ\ :sub:`r`\ G'\ :sup:`m` below. Moreover, if you know the concentrations of reactants in your system, you can click on the test tube icon to enter them explicitly into eQuilibrator. eQuilibrator will then calculate the appropriate Δ\ :sub:`r`\ G' value for you using the formula 

.. math::
	\begin{eqnarray}
	\Delta_r G' &=& \Delta_r G'^{\circ} + RT \ln{Q} \\
	\end{eqnarray}

Where Q is the "reaction quotient" or "mass-action ratio." Read more about these calculations `here <atp.html>`_.

What are Δ\ :sub:`r`\ G°, Δ\ :sub:`r`\ G'° and Δ\ :sub:`r`\ G'?
----------------------------------------------------------------------------

G is the Gibbs free energy of a physical system. It quantifies the amount of energy available for work in a system that is at constant temperature and pressure. Biological systems do not often experience radical changes in temperature or pressure, so G is appropriate for their analysis. On eQuilibrator we assume that the temperature is 25 °C (298.15 K) and the pressue is 1 bar. As a consequence of the second law of thermodynamics, spontaneous changes in systems at constant temperature and pressure must lead to a reduction in G. Chemical reactions change biochemical systems, so the same must be true of them: in order to occur, a reaction must lead to a reduction in the Gibbs free energy, implying that Δ\ :sub:`r`\ G' < 0.

Δ\ :sub:`r`\ G° is the change in Gibbs free energy due to a chemical reaction in standard conditions and without accounting for pH, ionic strength or any other cellular factors (i.e. they are taken as 0). Δ\ :sub:`r`\ G'° is the change in Gibbs free energy due to a chemical reaction due a reaction at a particular pH and ionic strength. These are convenient reference points, but neither Δ\ :sub:`r`\ G° nor Δ\ :sub:`r`\ G'° account for the effect of reactant concentrations on the free energy.

In cells there are certain prevailing metabolite concentrations. Concentrations affect the entropic component of the free energy: if substrate concentrations are higher than product concentrations then forward flow of the reaction will balance the concentrations and increase the entropy of the system. Δ\ :sub:`r`\ G' accounts for the effect of concentrations. It is very important to use Δ\ :sub:`r`\ G' values derived from physiologically plausible concentrations in thermodynamic analysis because reactant concentrations have a substantial effect on Δ\ :sub:`r`\ G' and because the laws of thermodynamics constrain only Δ\ :sub:`r`\ G', not Δ\ :sub:`r`\ G° or Δ\ :sub:`r`\ G'°.

What are Δ\ :sub:`f`\ G' and Δ\ :sub:`f`\ G'°?
----------------------------------------------------------

The formation energy of a compound is the change in free energy due to the compound's formation reaction - the reaction forming it from elements in their standard form. For example, the formation reaction of carbon dioxide is O2 + C ⇌ CO\ :sub:`2`. Δ\ :sub:`f`\ G° is the formation energy of a compound in standard conditions and without accounting for pH, ionic strength or any other cellular factors (i.e. they are taken as 0).

In contrast, Δ\ :sub:`f`\ G'° is the formation energy of a compound at a particular set of cellular conditions. When we show a Δ\ :sub:`f`\ G'° on eQuilibrator, as we do on compound and reaction pages, it has been transformed to the pH and ionic strength values. The default pH is 7 and the default ionic strength is 0.1 M.

What does the 'm' in Δ\ :sub:`r`\ G'\ :sup:`m`, Δ\ :sub:`f`\ G'\ :sup:`m`, and E'\ :sup:`m` mean?
----------------------------------------------------------------------------------------------------------------

As explained above, the reaction Gibbs energy depends on the concentrations of the reactants (and similarly for the formation energy Δ\ :sub:`f`\ G' and the reduction potential E'). It is useful to have a standard for these concentrations in order to compare reactions without specifying a concentration explicitly. It is standard to use 1 M concentrations to compare reactions in solution. However, in the context of metabolic reactions inside living cells, 1 mM is a much more appropriate standard concentration because metabolite concentrations typically range from 1 nM to 10 mM and never exceed 100 mM. As such, 1 M standard concentrations are entirely non-physiological.

The notation for the reaction Gibbs energy in the standard 1M concentrations is Δ\ :sub:`r`\ G'°, i.e. the degree sign (°) represents 1M. For our physiological standard we use the superscript 'm' to mark the 1mM concentration that is used for all reactants: Δ\ :sub:`r`\ G'm, Δ\ :sub:`f`\ G'm, and E'\ :sup:`m`.

Why is Δ\ :sub:`r`\ G "not available?"
----------------------------------------------------------

There are several reasons that a ΔG value might not be available.

In the case of a reaction, Δ\ :sub:`r`\ G is not available when Δ\ :sub:`f`\ G is missing for one of the reactants. Δ\ :sub:`f`\ G is very difficult to measure experimentally. We use the group contribution method to estimate Δ\ :sub:`f`\ G for compounds that have not been characterized experimentally [2,10]. However, we need certain data about a compound in order to approximate Δ\ :sub:`f`\ G.

Our list of compounds is taken from the Kyoto Encyclopedia of Genomes and Genes (`KEGG <http://www.kegg.jp/>`_) [5]. Occasionally, KEGG compound entries don't contain an InChI identifier. Without an InChI identifier we cannot compute the structure of the compound and without structural information we can't estimate Δ\ :sub:`f`\ G.

In some cases, a compound may contain a chemical group that was not observed in any of the compounds that have experimental Δ\ :sub:`f`\ G measurements. We cannot produce estimates for these compounds as we cannot estimate the contribution of unknown groups. In other cases, our algorithm for analyzing compound structure fails to decompose the compound into groups. If we are unable to decompose the compound then we cannot use the group contribution method. Moreover, we can't estimate Δ\ :sub:`r`\ G for any reaction that contains a compound that we cannot estimate Δ\ :sub:`f`\ G for.

How do you calculate the uncertainty for each estimation?
----------------------------------------------------------

In order to fully understand how to calculate estimation uncertainties, you'll probably need to read our paper on the Component Contribution method [12]. The short answer would be that we ran a cross-validation benchmark using a set of reactions for which the Δ\ :sub:`r`\ G' has been measured. Any reaction that you type in, is decomposed into compounds and these compounds are decomposed into groups. By comparing this decomposition vector with the ones from our database, we can estimate the Δ\ :sub:`r`\ G'. Along the way, we can also evaluate how good our estimation is, by checking how good we were for similar reactions in our benchmark.

How do you deal with gases like O\ :sub:`2` and H\ :sub:`2`?
---------------------------------------------------------------

For gases the standard condition is defined as 1 atmosphere (bar) partial pressure. However, if one knows the soluble concentration of the gas of interest it should be specified by choosing "custom" concentrations. Alternatively, if you know the partial pressure of a reactant in the gas phase of the reaction chamber, and assume there is an equilibrium with the solution phase, then you can set a concentration for the gas (in units of mbar) by clicking on the test tube icon. You can also indicate that you want to use the standard gas phase for the ΔG'° by adding (g) to the end of the compound name. For example, try to search for: 

|pep_carb|_

.. |pep_carb| replace:: CO\ :sub:`2`\ (g) + PEP + H\ :sub:`2`\ O = Oxaloacetate + Pi 
.. _pep_carb: /search?query=CO2%28g%29+%2B+PEP+%2B+H2O+%3D+Oxaloacetate+%2B+Pi

This will work only for compounds for which the formation energy in gas phase is found in our database, namely O\ :sub:`2`\ , N\ :sub:`2`\ , H\ :sub:`2`\ , CO\ :sub:`2`\ , and CO.

Why can't I change the concentration of H\ :sup:`+` ions?
----------------------------------------------------------

eQuilibrator uses the "Alberty method" for biochemical thermodynamics. In the Alberty method, H\ :sup:`+` is defined to have 0 free energy [6,8]. Instead of correcting for H\ :sup:`+` concentration, a pH correction accounts for the abundance of H+. You can use the pH "slider" to see the effect of H\ :sup:`+` concentration on your reaction.

Why can't I change the concentration of water?
----------------------------------------------------------

Biochemical systems are generally assumed to be aqueous environments [6,8]. Therefore, the concentration of water is fixed.

Why can't I change the temperature?
----------------------------------------------------------

The temperature is fixed at 25 °C (298.15 K) for all ΔG values given. The group contribution method enables us to approximate Δ\ :sub:`f`\ G of compounds at a particular temperature (the temperature at which they were measured) [10]. As the change in free energy is defined as ΔG = ΔH - TΔS and we don't know the value of ΔS in most cases, we cannot predict how changes in temperature will affect Δ\ :sub:`f`\ G.

What are CO\ :sub:`2`\ (aq) and CO\ :sub:`2`\ (total)?
----------------------------------------------------------

CO\ :sub:`2` in solution gives rise to several chemical species. It can be quite confusing to think about the equilibrium between these species - doing so requires care. CO\ :sub:`2`\ (aq) is dissolved CO\ :sub:`2`. CO\ :sub:`2`\ (aq) undergoes a spontaneous hydration reaction to form carbonic acid: 

CO\ :sub:`2` + H\ :sub:`2`\ O ⇌ H\ :sub:`2`\ CO\ :sub:`3`

or a similar reaction of: 

CO\ :sub:`2` + OH\ :sup:`-` ⇌ HCO\ :sub:`3`\ :sup:`-`

Several hydrated species form in water through addition or release of protons: carbonic acid (H\ :sub:`2`\ CO\ :sub:`3`), bicarbonate (HCO\ :sub:`3`\ :sup:`-`\ ) and carbonate (CO\ :sub:`3`\ :sup:`2-`\ ). In thermodynamics of biochemical reactions, different ionic states (known as pseudo-isomers) are lumped together. If you search for any of these hydrated species, eQuilibrator will use their lumped form - HCO\ :sub:`3`\ :sup:`-`\ (aq). Sometimes it is easier to measure or analyze the sum of CO\ :sub:`2`\ (aq) and its three hydrated forms (H\ :sub:`2`\ CO\ :sub:`3`, HCO3- and CO32-). This sum of species is termed CO\ :sub:`2`\ (total). Note that in the chemical formula of CO\ :sub:`2`\ (total) there are actually 3 oxygen atoms because it also includes the hydrating water molecule.

.. figure:: _static/_images/co2_hydration.png
   :alt: Hydration of Carbon Dioxide
   :align: center

   Hydration of CO\ :sub:`2`

When only the total concentration is known, it is assumed that there is equilibrium among the four species and one uses CO\ :sub:`2`\ (total) in place of CO\ :sub:`2`\ (aq). If, however, you know or can measure the concentration of CO\ :sub:`2`\ (aq) alone, then it is reasonable to use CO\ :sub:`2`\ (aq) as a reactant. The concentration of CO\ :sub:`2`\ (aq) is usually straightforward to calculate based on Henry’s law dictating, for example, that under standard atmospheric conditions of about 400ppm CO\ :sub:`2`\ (g) the concentration of CO\ :sub:`2`\ (aq) is about 10 uM. The concentrations of bicarbonate and CO\ :sub:`2`\ (total), however, depend strongly on pH. More information is supplied in the figure below and in this link.

The uncatalyzed hydration reaction (CO\ :sub:`2`\ (aq) + H\ :sub:`2`\ O ⇌ H\ :sub:`2`\ CO\ :sub:`3`) takes minutes to equilibrate. In many organisms, however, this reaction is catalyzed by the enzyme carbonic anhydrase, which speeds up the reaction by several orders of magnitude [3]. In cells, therefore, CO\ :sub:`2`\ (aq) is generally considered to be in equilibrium with its hydrated forms (carbonic acid, bicarbonate and carbonate) save in some special cases such as in cyanobacterial carbon concentrating mechanisms where carbonic anhydrase is absent from some parts of the cell.

We note that anaplerotic reactions use bicarbonate as their substrate whereas decarboxylation reactions release CO\ :sub:`2` but in finding the ΔG' they can be written using either CO\ :sub:`2`\ (aq), HCO\ :sub:`3`\ :sup:`-`\ (aq) or CO\ :sub:`2`\ (total) as long as the concentrations used are accurate. This is true because of the equilibrium among these species.

We know that this whole issue is quite confusing. We sincerely hope this explanation helps - please contact us if you have suggestions to explain the topic better.

.. todo:: link to TCA cycle/anaplerotic reactions when we have some content for it. 

What are "half-reactions?""
----------------------------------------------------------------------------

A `half-reaction <http://en.wikipedia.org/wiki/Half-reaction>`_ is the oxidation or reduction component of a `redox reaction <http://en.wikipedia.org/wiki/Redox>`_, without the other component. When you search for such a reaction, eQuilibrator recognizes that the number of electrons is not balanced and automatically switches to 'half-reaction' mode. Without knowing the other half, the change in Gibbs energy is not well defined. The parameter that is used to describe the potential difference (in Volts) between the products and substrates of a half-reaction is called the "`standard redox potential <http://en.wikipedia.org/wiki/Redox_reaction#Standard_electrode_potentials_.28reduction_potentials.29>`_" and is marked by E'°. The redox potential is equal to the voltage at equilibrium under standard conditions of an electrochemical cell in which the cathode reaction is the half-reaction considered and the anode is a standard hydrogen electrode where hydrogen is oxidized: ½ H\ :sub:`2` ⇌ H\ :sup:`+` + e\ :sup:`-`.

Assuming you do want the Gibbs energy of a reaction, you have two options. The first option is to balance the electrons in the half-reaction by supplying the other half. eQuilibrator make this simple providing a link for balancing with the biologically ubiquitous redox donor:acceptor pair `NAD+/NADH <glycolysis.html>`_. Alternatively, you can use the bottom panel of results page to adjust the potential of the electrons in the other half-reaction (i.e. change the value of e- potential in mV). This is useful in cases where eQuilibrator doesn't have a value for the second half-reaction, which is sometimes the case when the donors are complicated or not well-defined. For example, protein-based redox carriers like ferredoxin can vary quite quite a lot in their potential.

What's so complicated about redox reactions involving iron?
----------------------------------------------------------------------------

The reduction or oxidation of the pair Fe(III)/Fe(II) is ubiquitous in biology, for example in the iron-sulfur clusters of ferredoxins. However, the chemical environment of the iron atom can have a large effect on the reduction potential of the Fe(III)/Fe(II) pair with the measured reduction potentials of natural ferredoxins varying by more than 350 mV [4]. That is, variation in the measured reduction potential of ferredoxins equals to reduction potential of NAD/NADH!

Similarly, in dissimilatory iron reduction the specific chemical form of Fe(III) can drastically affect the reduction potential. For example, a half reaction with a well-characterized crystalline form Goethite has a redox potential of about -300 mV while y-FeOOH, (Lepidocrocite), which can be treated as having the same empirical formula, has a redox potential of about -100 mV at pH 7 [9]. As a result we strongly suggest that you enter the iron-free half-reaction of interest (e.g. `reduction of pyruvate to acetyl-CoA </search?query=+pyruvate+%2B+CoA+%3D+acetyl-CoA+%2B+CO2>`_) and use the bottom panel to adjust the potential of the electrons in the reaction to match the iron donor-acceptor pair that interests you.

Why is the value for ATP hydrolysis different than some textbooks?
----------------------------------------------------------------------------

The ΔG'° of the ATP hydrolysis reaction is affected by many factors, notably also by the concentration of free Mg\ :sub:`2`\ :sup:`+`\  ions. The value cited in [1] and used in the original version of eQuilibrator (-36.4 kJ/mol) assumes no magnesium ([Mg\ :sub:`2`\ :sup:`+`\ ] = 0). In the current version of eQuilibrator2.0 we use the component contribution method [12] that uses measurements collected in the NIST thermodynamic database for enzyme-catalyzed reactions [7] that were performed under varying levels of Mg\ :sub:`2`\ :sup:`+`\ . This is also the more relevant situation in vivo. As noted in many studies, when taking into account [Mg\ :sub:`2`\ :sup:`+`\ ], the value changes and is observed to be in the range -26 to -32 kJ/mol depending on the reference. A clear discussion can be found at [8].

What is the total driving force of a pathway?
----------------------------------------------------------

We define the driving force of a reaction or pathway as -ΔG' - i.e. a favorable reaction has a negative ΔG' and a positive driving force. The total driving force for a pathway is the driving force associated with the pathway net reaction at particular metabolite and cofactor concentrations. The minimum total driving force, which we calculate in `analyzing a pathway </pathway>`_, is the smallest driving force associated with that pathway given the limits assumed on metabolite an cofactor concentrations. Similarly, the maximum total driving force is the largest driving force associated with pathway given those same limits.


What is the MDF of a pathway?
----------------------------------------------------------

The MDF of a pathway is a metric of how thermodynamically favorable a pathway can be in physiological conditions. The value of the MDF is smallest -ΔG' obtained by any pathway reaction when metabolite concentrations are chosen to make all pathway reactions as favorable as possible (-ΔG' as positive as possible).

You can read more about the MDF in `this paper <http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1003483>`_ [13].

How can I contact the people behind eQuilibrator?
----------------------------------------------------------

If you have questions about eQuilibrator, please consult the `eQuilibrator Google Group <https://groups.google.com/forum/#!forum/equilibrator-users>`_ to see if your question has been asked and answered before. Please also post your questions to the group so that all eQuilibrator users can benefit from your experience. If you have questions about the data and methods on which eQuilibrator is based, please consult `these references <http://localhost:8000/cite>`_. If you'd like to contact us directly, reach out to the `Milo Lab <http://www.weizmann.ac.il/plants/Milo/>`_, which maintains eQuilibrator.

References
----------------------------------------------------------

#. R. Guynn, R.J. Veech, "The equilibrium constants of the adenosine triphosphate hydrolysis and the adenosine triphosphate-citrate lyase reactions", The Journal of Biological Chemistry (1973) 248:6966-6972
#. M.L. Mavrovouniotis, "Estimation of standard Gibbs energy changes of biotransformations" The Journal of Biological Chemistry (1991) 266(22):14440-14445
#. A. Radzicka, R. Wolfenden, "A proficient enzyme", Science (1995) 267:90-93
#. P. J. Stephens, D. R. Jollie, A. Warshel, "Protein Control of Redox Potentials of Iron−Sulfur Proteins" Chem. Rev. (1996) 96:2491–2514
#. M. Kanehisa, S. Goto, "KEGG: Kyoto Encyclopedia of Genes and Genomes" Nucleic Acids Research (2000) 28(1):27-30
#. R.A. Alberty, "Thermodynamics of biochemical reactions" (Hoboken N.J.: Wiley-Interscience, 2003)
#. R.N. Goldberg, Y.B. Tewari, T.N. Bhat, "Thermodynamics of Enzyme-Catalyzed Reactions - a Database for Quantitative Biochemistry", Bioinformatics (2004) 20(16):2874-2877
#. R.A. Alberty, "Biochemical Thermodynamics" (Hoboken, NJ, USA: John Wiley & Sons, Inc., 2006)
#. A. Navrotsky, L. Mazeina, J. Majzlan, "Size-Driven Structural And Thermodynamic Complexity In Iron-Oxides" Science (2008) 319:1635–1638
#. M.D. Jankowski et al., "Group Contribution Method for Thermodynamic Analysis of Complex Metabolic Networks" Biophysical Journal (2008) 95(3):1487-1499
#. E. Noor, A. Bar-Even, A. Flamholz, Y. Lubling, D. Davidi, R. Milo, "An integrated open framework for thermodynamics of reactions that combines accuracy and coverage" Bioinformatics (2012) 28:2037-2044
#. E. Noor, H.S. Haraldsdóttir, R. Milo, R.M.T. Fleming, "Consistent Estimation of Gibbs Energy Using Component Contributions" PLoS Comput Biol (2013) 9:e1003098
#. E. Noor, A. Bar-Even, A. Flamholz, E. Reznik, W. Liebermeister, R. Milo, "Pathway thermodynamics highlights kinetic obstacles in central metabolism" PLoS Comput Biol (2014) 10:e1003483

