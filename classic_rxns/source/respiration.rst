Respiration
==========================================================

When molecular oxygen (O\ :sub:`2`) is available many organisms will opt to respire - to donate electrons withdrawn from glucose to O\ :sub:`2` - and store even more energy as ATP in the process. Instead of reducing pyruvate to some fermentation product(s) (i.e. donating electrons to pyruvate as in the examples above), the carbon enters the tricarboxylic acid (TCA), which oxidizes it (withdraws its electrons) one step at a time to produce 3 CO\ :sub:`2` and 5 reduced electron carriers (3 NADH, 1 FADH2 and one reduced quinone).

Pyruvate + 3 NAD+ + 1 FAD + Ubiquinone + 3 H :sub:`2` O ⇌ 3 CO\ :sub:`2` + 3 NADH + 1 FADH2 + Ubiquinol

While the oxidation of pyruvate is quite favorable, one caveat has to be considered in examining the net reaction above: the Ubiquinone electron acceptor used in the TCA cycle is a highly hydrophobic molecule that resides primarily in lipid membranes and not in the cytosol (or the mitochondrial matrix of eukaryotes). In fact, the enzyme that donates electrons to ubiquinone (succinate dehydrogenase) also resides on the membrane. What’s the problem? Well, if some reactants reside in the cytosol (or mitochondrial matrix) and others reside in the membrane, it is difficult to compare their concentrations and calculate Δ\ :sub:`r`\ G'. So we can calculate a Δ\ :sub:`r`\ G'° of about -120 kJ/mol and a Δ\ :sub:`r`\ G'm of about -155 kJ/mol for the reaction above, but it is not clear whethe these values are representative of what happens in the cell. 

Let’s now take a step back and remember that the electrons withdrawn from pyruvate and donated to NAD+, FAD or ubiquinone are ultimately going to be given to the terminal electron acceptor - O\ :sub:`2`. This process is achieved through a remarkable chain of membrane-bound protein complexes called the “electron transport chain” or “respiratory chain” (Nelson et al., 2008). The net effect of this chain of redox reactions can be summarized as

2 Pyruvate + 5 Oxygen ⇌ 6 CO\ :sub:`2` + 4 H :sub:`2` O

Notice that we’ve written this reaction as consuming two pyruvate because (a) each glucose molecule produces two pyruvate through glycolysis and (b) we can avoid considering fractional numbers of oxygen this way. Aside from noting that this net reaction is extremely favorable, with a Δ\ :sub:`r`\ G'm of about -2300 kJ/mol, we also notice that respiration involves the net production of water and CO\ :sub:`2`. Respiration is then the opposite of photosynthetic carbon fixation, which uses light energy to withdraw electrons from water, donate them to CO\ :sub:`2` and (ultimately) make glucose (Nelson et al., 2008; Buchanan et al., 2015).

While it’s clear that there is a lot of energy released in the oxidation of pyruvate (enough for more than 50+ ATP) it’s not so clear from this reaction scheme how the ATP gets made. The short answer is that some of the chemical energy released during electron transfer to O\ :sub:`2` is stored (by the action of the electron transport chain) in a battery of sorts - a battery powered by the difference in proton (H+) concentrations across the cell (or mitochondrial) membrane. This membrane battery is then tapped by a stunningly beautiful rotary protein motor - the ATP synthase - which uses the energy stored in the gradient of chemical potential across the membrane to catalyze the phosphorylation of ADP to ATP. This video illustrating the catalytic cycle of the ATP synthase is mesmerizing. You can find a much more detailed explanation of this incredible process in most biochemistry textbooks (Nelson et al., 2008). 

Alternative Respiratory Pathways
----------------------------------------------------------

As humans, we might think that respiring is synonymous with breathing oxygen because for us it is. But respiration is defined as a metabolic process where electrons are withdrawn from food (e.g. glucose) and donated to an external molecule (electron acceptor, e.g. O\ :sub:`2`) that is not a metabolic product of the food. The terminal electron acceptor does not need to be O\ :sub:`2` - it can be any common molecule with a tendency to accept electrons. Oxygen is the most commonly discussed terminal electron acceptor because it has the greatest reduction potential of any acceptor used by biology - the greatest tendency to accept electrons. 

E. coli, a common and metabolically flexible bacterium, is capable of using both nitrate (NO3-) and nitrite (NO\ :sub:`2`-) as terminal electron acceptors (Unden and Bongaerts, 1997). We can use eQuilibrator to compare the reduction potential of these three electron acceptors, O\ :sub:`2`, NO3- and NO\ :sub:`2`-

Oxygen  ⇌ 2 H :sub:`2` O
Nitrate ⇌ Nitrite + H :sub:`2` O
Nitrite ⇌ Ammonia + 2 H :sub:`2` O

The oxygen/water pair has, by far, the greatest tendency to accept electrons with a reduction potential E’m of about 800 millivolts per electron transferred. The nitrate/nitrite and nitrite/ammonia pairs both have lower (but still positive) reduction potentials around 350-400 millivolts per electron. 

What do these reduction potentials mean? They represent the amount of energy that is released by donating an electron to one molecule in order to form another. For example, we have to donate 4 e- to O\ :sub:`2` to form 2 H :sub:`2` O or two e- per water molecule (try to verify this for yourself) - the reduction potential denotes the average amount of energy released by donating these electrons to O\ :sub:`2`. The reduction potential is relative, however. In order to know how much energy is released in units of kJ/mol we need to know what potential the electrons started with! 

Let’s say we take electrons from our favorite sugar and make CO\ :sub:`2`

Glucose + 6 H :sub:`2` O ⇌ 6 CO\ :sub:`2`

This half-reaction has a per-electron potential E’m of about -400 millivolts (mV). Since this value is much more negative than the values for O\ :sub:`2`, nitrate or nitrite, we know that withdrawing electrons from glucose and donating them to these molecules will be favorable - remember that electrons tend to flow towards more positive reduction potential. 

We can convert these mV values into more familiar kJ/mol using the Nernst equation Δ\ :sub:`r`\ G' = nFΔE’, where n is the number of electrons transferred and F is the Faraday constant (F ≈ 96.4x103 kJ/mV/mol). Donating all 24 electrons from glucose to O\ :sub:`2` should release 24 x 96.4x10-3 x (-400-800) ≈ -2800 kJ/mol (consistent with our calculation of this Δ\ :sub:`r`\ G above using eQuilibrator). Respiring using nitrate or nitrite as the electron acceptor will release less energy because they have lower (less positive) reduction potentials. If we use 400 mV as a representative value for these two acceptors, then Δ\ :sub:`r`\ G’ = 24 x 96.4x10-3 x (-400-400) ≈ -1800 kJ / mol. So it’s not surprising that E. coli prefers to use O\ :sub:`2` as a terminal electron acceptor over nitrate and nitrite - respiring glucose using oxygen releases 1000 kJ/mol more energy! That’s more than 20 ATPs worth of energy! In fact, the presence of O\ :sub:`2` represses all of E. coli’s alternative respiratory pathways.

An exercise for the reader: try using eQuilibrator to combine the glucose/CO\ :sub:`2` half-reactions with the electron acceptor half-reactions (O\ :sub:`2`, nitrate and nitrite) to figure out the stoichiometries of these respiratory pathways. 
