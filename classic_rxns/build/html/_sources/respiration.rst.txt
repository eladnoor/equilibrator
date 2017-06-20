Respiration
==========================================================

When molecular oxygen (O\ :sub:`2`) is available many organisms will opt to respire - to donate electrons withdrawn from glucose to O\ :sub:`2`. Respiration is advantageous (as compared to fermentation) because it allows cells to store much more energy as ATP while metabolizing glucose. Instead of reducing pyruvate to some fermentation product(s) (i.e. donating electrons to pyruvate as in the examples above), the carbon enters the tricarboxylic acid (TCA), which oxidizes it (withdraws its electrons) one step at a time to produce 3 CO\ :sub:`2` and 5 reduced electron carriers (3 NADH, 1 FADH2 and one reduced quinone).

|pyr_resp_carriers|_

.. |pyr_resp_carriers| replace:: Pyruvate + 3 NAD+ + 1 FAD + Ubiquinone + 3 H :sub:`2` O ⇌ 3 CO\ :sub:`2` + 3 NADH + 1 FADH2 + Ubiquinol
.. _pyr_resp_carriers: http://equilibrator.weizmann.ac.il/search?query=Pyruvate+%2B+3+NAD%2B+%2B+1+FAD+%2B+Ubiquinone+%2B+3+H2O+%3C%3D%3E+3+CO2+%2B+3+NADH+%2B+1+FADH2+%2B+Ubiquinol

While the oxidation of pyruvate is quite favorable, one caveat has to be considered in examining the net reaction above: the Ubiquinone electron acceptor used in the TCA cycle is a highly hydrophobic molecule that resides primarily in lipid membranes and not in the cytosol. [1]_ In fact, the enzyme that donates electrons to ubiquinone (`succinate dehydrogenase <http://equilibrator.weizmann.ac.il/enzyme?ec=1.3.5.1>`_) also resides on the membrane. What’s the problem? Well, if some reactants reside in the cytosol (or mitochondrial matrix) and others reside in the membrane, it is difficult to compare their concentrations and calculate Δ\ :sub:`r`\ G'. So we can calculate a Δ\ :sub:`r`\ G'° of about -120 kJ/mol and a Δ\ :sub:`r`\ G'm of about -155 kJ/mol for the reaction above, but it is not clear whethe these values are representative of what happens in the cell. 

Let’s now take a step back and remember that the electrons withdrawn from pyruvate and donated to NAD+, FAD or ubiquinone are ultimately going to be given to the terminal electron acceptor - O\ :sub:`2`. This process is achieved through a remarkable chain of membrane-bound protein complexes called the `electron transport chain <https://en.wikipedia.org/wiki/Electron_transport_chain>`_ (Nelson et al., 2008). The net effect of this chain of redox reactions can be summarized as

|pyruvate_resp_ox|_

.. |pyruvate_resp_ox| replace:: 2 Pyruvate + 5 Oxygen ⇌ 6 CO\ :sub:`2` + 4 H :sub:`2` O
.. _pyruvate_resp_ox: http://equilibrator.weizmann.ac.il/search?query=2+Pyruvate+%2B+5+Oxygen+%3C%3D%3E+6+CO2+%2B+4+H2O

Notice that we’ve written this reaction as consuming two pyruvate because (a) each glucose molecule produces two pyruvate through glycolysis and (b) we can avoid considering fractional numbers of oxygen this way. Aside from noting that this net reaction is extremely favorable, with a Δ\ :sub:`r`\ G'm of about -2300 kJ/mol, we also notice that respiration involves the net production of water and CO\ :sub:`2`. Respiration is then the opposite of photosynthetic carbon fixation, which uses light energy to withdraw electrons from water, donate them to CO\ :sub:`2` and (ultimately) make glucose (Nelson et al., 2008; Buchanan et al., 2015).

While it’s clear that there is a lot of energy released in the oxidation of pyruvate [2]_ it’s not so clear from this reaction scheme how the ATP gets made. The short answer is that some of the chemical energy released during electron transfer to O\ :sub:`2` is stored (by the action of the electron transport chain) in a battery of sorts - a battery powered by the difference in proton (H+) concentrations across the cell (or mitochondrial) membrane. This membrane battery is then tapped by a stunningly beautiful rotary protein motor - the `ATP synthase <https://pdb101.rcsb.org/motm/72>`_ - which uses the energy stored in the gradient of chemical potential across the membrane to catalyze the phosphorylation of ADP to ATP. `This video <https://www.youtube.com/watch?v=GM9buhWJjlA>`_ illustrating the catalytic cycle of the ATP synthase is mesmerizing. You can find a much more detailed explanation of this incredible process in `most biochemistry textbooks <https://www.ncbi.nlm.nih.gov/books/NBK21528/>`_ (Nelson et al., 2008). 

.. [1] Or the `mitochondrial matrix <https://en.wikipedia.org/wiki/Mitochondrion>`_ of eukaryotes.
.. [2] Enough to make more than 50 ATP in theory.

Alternative Respiratory Pathways
----------------------------------------------------------

As humans, we might think that respiring is synonymous with breathing oxygen because for us it is. But respiration is defined as a metabolic process where electrons are withdrawn from food (e.g. glucose) and donated to an external molecule (electron acceptor, e.g. O\ :sub:`2`) that is not a metabolic product of the food. The terminal electron acceptor does not need to be O\ :sub:`2` - it can be any common molecule with a tendency to accept electrons. Oxygen is the most commonly discussed terminal electron acceptor because it has the greatest reduction potential of any acceptor used by biology - the greatest tendency to accept electrons. 

*E. coli*, a common and metabolically flexible bacterium, is capable of using both nitrate (NO\ :sub:`3`\ :sup:`-`) and nitrite (NO\ :sub:`2`\ :sup:`-`) as terminal electron acceptors (Unden and Bongaerts, 1997). We can use eQuilibrator to compare the reduction potential of these three electron acceptors, O\ :sub:`2`, NO3- and NO\ :sub:`2`-

#. `Oxygen + 4 e- ⇌ 2 Water <http://equilibrator.weizmann.ac.il/search?query=oxygen++%3C%3D%3E+2+H2O>`_

#. `Nitrate + 2 e- ⇌ Nitrite + Water <http://equilibrator.weizmann.ac.il/reaction?reactantsId=C00244&reactantsCoeff=-1&reactantsName=Nitrate&reactantsPhase=aqueous&reactantsConcentration=0.001&reactantsId=C00088&reactantsCoeff=1&reactantsName=Nitrite&reactantsPhase=aqueous&reactantsConcentration=0.001&reactantsId=C00001&reactantsCoeff=1&reactantsName=H2O&reactantsPhase=liquid&reactantsConcentration=1&ph=7.000000&pmg=14.000000&ionic_strength=0.100000&e_reduction_potential=0.000000&max_priority=0&mode=BA&query=nitrate%20%20%3C%3D%3E%20nitrite>`_

#. `Nitrite + 6 e- ⇌ Ammonia + 2 Water <http://equilibrator.weizmann.ac.il/reaction?reactantsId=C00088&reactantsCoeff=-1&reactantsName=Nitrite&reactantsPhase=aqueous&reactantsConcentration=0.001&reactantsId=C00014&reactantsCoeff=1&reactantsName=Ammonia&reactantsPhase=aqueous&reactantsConcentration=0.001&reactantsId=C00001&reactantsCoeff=2&reactantsName=H2O&reactantsPhase=liquid&reactantsConcentration=1&ph=7.000000&pmg=14.000000&ionic_strength=0.100000&e_reduction_potential=0.000000&max_priority=0&mode=BA&query=Nitrite%20%3C%3D%3E%20ammonia>`_

The oxygen/water pair has, by far, the greatest tendency to accept electrons with a reduction potential E’\ :sup:`m` of about 800 millivolts per electron transferred. The nitrate/nitrite and nitrite/ammonia pairs both have lower (but still positive) reduction potentials around 350-400 millivolts per electron. 

What do these reduction potentials mean? They represent the amount of energy that is released by donating an electron to one molecule in order to form another. For example, we have to donate 4 e- to O\ :sub:`2` to form 2 H :sub:`2` O or two e- per water molecule (try to verify this for yourself) - the reduction potential denotes the average amount of energy released by donating these electrons to O\ :sub:`2`. The reduction potential is relative, however. In order to know how much energy is released in units of kJ/mol we need to know what potential the electrons started with! 

Let’s say we take electrons from our favorite sugar and make CO\ :sub:`2`

|gluc_ox_half|_

.. |gluc_ox_half| replace:: Glucose + 6 H :sub:`2` O ⇌ 6 CO\ :sub:`2` + 24 e-
.. _gluc_ox_half: http://equilibrator.weizmann.ac.il/reaction?reactantsId=C00031&reactantsCoeff=-1&reactantsName=Glucose&reactantsPhase=aqueous&reactantsConcentration=0.001&reactantsId=C00011&reactantsCoeff=6&reactantsName=CO2&reactantsPhase=aqueous&reactantsConcentration=0.001&reactantsId=C00001&reactantsCoeff=-6&reactantsName=H2O&reactantsPhase=liquid&reactantsConcentration=1&ph=7.000000&pmg=14.000000&ionic_strength=0.100000&e_reduction_potential=0.000000&max_priority=0&mode=BA&query=glucose%20%3D%3E%206co2

This half-reaction has a per-electron potential E’\ :sup:`m` of about -400 millivolts (mV). Since this value is much more negative than the values for O\ :sub:`2`, nitrate or nitrite, we know that withdrawing electrons from glucose and donating them to these molecules will be favorable - remember that electrons tend to flow towards more positive reduction potential. 

We can convert these mV values into more familiar kJ/mol using the Nernst equation Δ\ :sub:`r`\ G' = nFΔE’, where n is the number of electrons transferred and F is the Faraday constant (96.4x103 kJ/mV/mol). Donating all 24 electrons from glucose to O\ :sub:`2` should release [3]_

.. math::
	\Delta_r G'^m \approx 2 \times 96.4 \times 10^{-3} \times (-400 - 800) \approx -2800 \frac{kJ}{mol}

Respiring using nitrate or nitrite as the electron acceptor will release less energy because they have lower (less positive) reduction potentials. If we use 400 mV as a representative value for these two acceptors, then 

.. math::
	\Delta_r G'^m \approx 2 \times 96.4 \times 10^{-3} \times (-400 - 400) \approx -1800 \frac{kJ}{mol}

So it’s not surprising that *E. coli* prefers to use O\ :sub:`2` as a terminal electron acceptor over nitrate and nitrite - respiring glucose using oxygen releases 1000 kJ/mol more energy! That’s more than 20 ATPs worth of energy! In fact, the presence of O\ :sub:`2` represses all of *E. coli*’s alternative respiratory pathways.

.. todo::
	reference for repression statement above

An exercise for the reader: try using eQuilibrator to combine the glucose/CO\ :sub:`2` half-reactions with the electron acceptor half-reactions (O\ :sub:`2`, nitrate and nitrite) to calculate the stoichiometries of these alternate respiratory pathways. 

.. [3] consistent with our earlier of this Δ\ :sub:`r`\ G above using eQuilibrator.
