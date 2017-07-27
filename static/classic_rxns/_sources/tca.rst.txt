The TCA Cycle
==========================================================

As we've discussed elsewhere (in the sections on `thioesters <thioester.html>`_ and `fatty acid metabolism <fatty_acid_met.html>`_) acetyl-CoA is a very central molecule in metabolism. 

.. figure:: _static/_images/accoa.png
   :alt: Acetyl-CoA
   :align: center

   Acetyl-CoA is a very central molecule in metabolism and is produced during the metabolism of sugars, alkanes, fatty acids and some amino acids like leucine. 

Indeed, after glucose is metabolized to pyruvate through `glycolysis <glycolysis.html>`_ acetyl-CoA is made from pyruvate and CoA by a very large enzyme called `the pyruvate dehydrogenase complex <pdb101.rcsb.org/motm/153>`_. 

|pdh|_

.. |pdh| replace:: Pyruvate + CoA + NAD\ :sup:`+` ⇌ acetyl-CoA + CO\ :sub:`2` + NADH
.. _pdh: http://equilibrator.weizmann.ac.il/search?query=Pyruvate+%2B+CoA+%2B+NAD%2B+%E2%87%8C+acetyl-CoA+%2B+CO2+%2B+NADH


This reaction is an oxidative decarboxylation - pyruvate is oxidized and decarboxylated (a CO\ :sub:`2` is removed) to form acetate. The free-energy released by this transformation is substantial (Δ\ :sub:`r`\ G'\ :sup:`m` ≈ -80 kJ/mol) and is harnessed to make a thioester bond on acetyl-CoA. You can see directly how favorable this transformation is by omitting the CoA from the reaction and entering it into eQuilibrator:

|pdh_noA|_

.. |pdh_noA| replace:: Pyruvate + H\ :sub:`2`\ O + NAD\ :sup:`+` ⇌ acetate + CO\ :sub:`2` + NADH
.. _pdh_noA: http://equilibrator.weizmann.ac.il/search?query=Pyruvate+%2B+H2O+%2B+NAD%2B+%E2%87%8C+acetate+%2B+CO2+%2B+NADH

.. todo:: mention why we need water in a footnote?

We will not go into detail about this complex enzyme here. Rather, we bring it up to highlight two important concepts. First, that oxidative decarboxylations are very favorable reactions, meaning that NADH can be produced by decarboxylating acids like pyruvate. We will return to this theme as we discuss the TCA cycle below. Second, that acetyl-CoA is a major point of convergence in central metabolism, standing at the crossroads between the catabolism of glucose, amino acids and lipids. But how is acetyl-CoA itself metabolized? 

Metabolism of Acetyl-CoA
----------------------------------------------------------
The tricarboxylic acid (TCA) cycle oxidizes the acetyl group of acetyl-CoA

|TCA_net_acCoA|_

.. |TCA_net_acCoA| replace:: Acetyl-CoA + 3 H\ :sub:`2`\ O ⇌ 2 CO\ :sub:`2` + CoA + 8 e\ :sup:`-` 
.. _TCA_net_acCoA: http://equilibrator.weizmann.ac.il/search?query=Acetyl-CoA+%2B+3+H2O+%E2%87%8C+2+CO2+%2B+CoA+%2B+8+e-

We can see that the net effect of the TCA cycle is to decarboxylate acetyl-CoA twice to make  CO\ :sub:`2` and 8 electrons. In the TCA cycle, electrons are mostly donated to NAD\ :sup:`+` to form `NADH <glycolysis.html#nadh>`_. These electrons are subsequently donated to O\ :sub:`2` in a manner that is coupled to the formation of ATP via the `respiratory electron transport chain <respiration.html#etc>`_.

You might think that the oxidation of acetyl-CoA would proceed directly. However, the direct decarboxylation of acetate would give CO\ :sub:`2` and CH\ :sub:`4` (methane) a reduced (electron-rich) one-carbon molecule that is very stable gas that is quite complex to metabolize (`Lehninger Principles of Biochemistry, 2008 <refs.html>`_). In fact, methanotrophs (organisms that consume methane) use very complex and sensitive enzymes to do so. The direct decarboxylation of acetyl-CoA poses similar problems. Biology solves this problem by using "scaffold" molecules - the acetyl group is first attached to another metabolite and only then does metabolism proceed by decarboxylation and oxidation. These scaffolds aid by "holding onto" the acetyl group and by providing a chemical environment that increases its reactivity (making it more straightforward to metabolize, `Lehninger Principles of Biochemistry, 2008 <refs.html>`_).

We will approach the scaffold-based metabolism of acetyl-CoA through a sort of back door by reminding ourselves about `transamination <transamination.html>`_. Remember that many amino acids are made by transamination of an appropriate alpha-keto acid acceptor. For example, glutamate can be made transaminating a five carbon dicarboxylic acid called α-ketoglutarate (also known as 2-oxoglutarate). 

.. figure:: _static/_images/alphaketoglutarate.png
   :alt: α-ketoglutaric acid or 2-oxoglutarate
   :align: center

   α-ketoglutarate or 2-oxoglutarate is called α-keto acid because it has a carbonyl group (ketone) adjacent to (in the α position relative to) a terminal carboxylic acid.

Glutamate is a proteogenic amino acid in its own right, but also a precursor for the synthesis of several amino acids including glutamine, proline and arginine (`Lehninger Principles of Biochemistry, 2008 <refs.html>`_). Similarly, aspartate can be made by transamination of a four carbon dicarboxylic acid called oxaloacetate. Aspartate is a precursor for the production of asparagine, lysine, arginine, methionine, threonine and isoleucine (`Lehninger Principles of Biochemistry, 2008 <refs.html>`_).

.. figure:: _static/_images/oxaloacetate.png
   :alt: Oxaloacetate
   :align: center

   Oxaloacetate is a four carbon dicarboxylic acid that can be transaminated to form aspartate.

.. figure:: _static/_images/aspartate.png
   :alt: Aspartate
   :align: center

   Aspartate is made by transamination of oxaloacetate - notice how the only difference between the two is the amino group.
   
.. todo:: links for individual AAs and for the biosynthetic pathways. 

.. todo:: footnote about variation in AA biosynthetic pathways. references.

In order to make amino acids, therefore, cells need to make α-ketoglutarate and oxaloacetate continually. One way to think about the tricaboxylic acid cycle, then, is as means for metabolizing acetyl-CoA that uses these very central molecules, α-ketoglutarate and oxaloacetate, as scaffolds. We might imagine that these two very central molecules were already present in the cell at relatively high concentrations because they are necessary for amino acid synthesis, which could explain why they evolved to be so central to acetyl-CoA metabolism.

With this in mind, let's consider again the net reaction of the TCA cycle starting from acetyl-CoA

|TCA_net_acCoA2|_

.. |TCA_net_acCoA2| replace:: Acetyl-CoA + 3 H\ :sub:`2`\ O ⇌ 2 CO\ :sub:`2` + CoA + 8 e\ :sup:`-` 
.. _TCA_net_acCoA2: http://equilibrator.weizmann.ac.il/search?query=Acetyl-CoA+%2B+3+H2O+%E2%87%8C+2+CO2+%2B+CoA+%2B+8+e-

Notice that both α-ketoglutarate and oxaloacetate are absent from the net reaction because they are intermediates of the cycle - the scaffolds are neither created or destroyed in cycle, but rather serve as co-substrates for the metabolism of acetyl-CoA. The net effect of the TCA cycle is to decarboxylate acetyl-CoA twice to make  CO\ :sub:`2` and 8 electrons. These electrons have an average potential of E'\ :sup:`m` ≈ -350 mV, meaning that they can be donated to NAD\ :sup:`+` to form NADH (which has a more positive E'\ :sup:`m` ≈ -330 mV). 

You might notice that these E'\ :sup:`m` values are very close to each other (within 20 mV), meaning small variations in concentrations might make these redox reactions infeasible. Cells solve this problem in two ways. First, as we will see, only 6 of 8  electrons are donated to NADH. The other two e- are donated to a higher (more positive) potential acceptor FAD+ (E'\ :sup:`m` ≈ -240 mV). Second, the NAD\ :sup:`+` and NADH concentrations are not exactly 1 mM. Typically [NAD\ :sup:`+`] is about 30-50 times greater than [NADH] - try using eQuilibrator to see how this concentration ratio affects the E' value for NAD\ :sup:`+` reduction to NADH.

Constructing the TCA Cycle
----------------------------------------------------------
We are left with the task of building a cycle for the metabolism of acetyl-CoA to CO\ :sub:`2` and electrons that uses the central amine acceptors α-ketoglutarate and oxaloacetate as scaffolds. We'll do this by considering the two arms of this cycle - from α-ketoglutarate and oxaloacetate and back again. 

α-Ketoglutarate to Oxaloacetate
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Let's first consider what it would take to make oxaloacetate from α-ketoglutarate.

|AKG_2_OAA|_

.. |AKG_2_OAA| replace:: α-ketoglutarate + 2 H\ :sub:`2`\ O ⇌ oxaloacetate + CO\ :sub:`2` + 6 e\ :sup:`-` 
.. _AKG_2_OAA: http://equilibrator.weizmann.ac.il/search?query=alpha-ketoglutarate+%2B+2+H2O+%E2%87%8C+oxaloacetate+%2B+CO2+%2B+6+e-

As we said, α-ketoglutarate (α-KG for short) contains 5 carbon atoms and oxaloacetate (OAA for short) contains 4. Therefore, production of OAA from α-KG requires the decarboxylation of α-KG - i.e. the removal of a carbon atom. As we foreshadowed above, this proceeds through the mechanism of oxidative decaboxylation of α-KG to make succinyl-CoA

|AKG_2_OAA_NADH|_

.. |AKG_2_OAA_NADH| replace:: α-ketoglutarate + CoA + NAD\ :sup:`+` ⇌ succinyl-CoA + CO\ :sub:`2` + NADH
.. _AKG_2_OAA_NADH: http://equilibrator.weizmann.ac.il/search?query=alpha-ketoglutarate+%2B+CoA+%2B+NAD%2B+%E2%87%8C+succinyl-CoA+%2B+CO2+%2B+NADH

The thioester on succinyl-CoA is `approximately energetically equivalent to ATP <thioester.html>`_, which explains how the next reaction step manages to make ATP while hydrolysing the thioester. [1]_

|ScCoA_synth|_

.. |ScCoA_synth| replace:: succinyl-CoA + ADP + Pi ⇌ succinate + CoA + ATP
.. _ScCoA_synth: http://equilibrator.weizmann.ac.il/search?query=succinyl-CoA+%2B+ADP+%2B+Pi+%E2%87%8C+succinate+%2B+CoA+%2B+ATP

Succinate has 4 carbons, like OAA, but is more reduced - having 4 more electrons. So it must be oxidized twice to make oxaloacetate. 

|SUCC_2_OAA|_

.. |SUCC_2_OAA| replace:: succinate + NAD\ :sup:`+` + FAD + H\ :sub:`2`\ O ⇌ oxaloacetate + NADH + FADH\ :sub:`2`
.. _SUCC_2_OAA: http://equilibrator.weizmann.ac.il/search?query=succinate+%2B+NAD%2B+%2B+FAD+%2B+H2O+%E2%87%8C+oxaloacetate+%2B+NADH+%2B+FADH2

The above reaction is actually a three-step process catalyzed by three different enzymes in the TCA cycle. Notice that two electrons are donated to NAD\ :sup:`+` and two are donated to a similar, but higher-potential donor called FAD (as discussed in above). Altogether, this arm of the TCA cycle has a net reaction of

|AKG_2_OAA_NET|_

.. |AKG_2_OAA_NET| replace:: α-ketoglutarate + 2 H\ :sub:`2`\ O + 2 NAD\ :sup:`+` + FAD ⇌ oxaloacetate + CO\ :sub:`2` + 2 NADH + FADH\ :sub:`2`
.. _AKG_2_OAA_NET: http://equilibrator.weizmann.ac.il/search?query=alpha-ketoglutarate+%2B+2+H2O+%2B+2+NAD%2B+%2B+FAD+%E2%87%8C+oxaloacetate+%2B+CO2+%2B+2+NADH+%2B+FADH2

and a Δ\ :sub:`r`\ G'\ :sup:`m` ≈ 0. As discussed above, this almost infeasible Δ\ :sub:`r`\ G'\ :sup:`m` can be remedied by setting the NAD\ :sup:`+` and NADH concentrations to more physiologically relevant values. For example, `measurements in E. coli <http://book.bionumbers.org/what-are-the-concentrations-of-free-metabolites-in-cells/>`_ give [NAD\ :sup:`+`] ≈ 3 mM and [NADH] ≈ 0.08 mM. Try using these values to calculate Δ\ :sub:`r`\ G' in eQuilibrator - does this help resolve the problem?

Oxaloacetate back to α-Ketoglutarate
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
First of all - it is important to note that the TCA cycle can't possibly follow the same reaction scheme in both directions. If the TCA cycle used one reaction sequence from α-ketoglutarate to oxaloacetate and then the same sequence in the reverse direction to get back to α-ketoglutarate, this would be a closed cycle with nonzero flux [2]_ - a perpetual motion machine in violation of the first law of thermodynamics. Another, perhaps even simpler reason that the TCA cycle can't use the same reaction sequence in the reverse direction is that we haven't yet seen how acetyl-CoA is incorporated into the cycle. Acetyl-CoA needs to find it's way into the cycle in order for the TCA cycle to metabolize acetyl-CoA, after all! 

Acetyl-CoA gets into the TCA cycle is through the citrate synthase reaction

|citrate_synthase|_

.. |citrate_synthase| replace:: acetyl-CoA + oxaloacetate + H\ :sub:`2`\ O ⇌ CoA + citrate
.. _citrate_synthase: http://equilibrator.weizmann.ac.il/search?query=acetyl-CoA+%2B+oxaloacetate+%2B+H2O+%E2%87%8C+CoA+%2B+citrate

which adds acetyl-CoA to oxaloacetate to form the tricarboxylic acid citrate (after which the cycle is named). 

.. figure:: _static/_images/citrate.png
   :alt: Citrate
   :align: center

   Citrate is the a tricarboxylic acid after which the TCA cycle is named.

The citrate synthase reaction is quite favorable (Δ\ :sub:`r`\ G'\ :sup:`m` ≈ -35 kJ/mol) due to the hydrolysis of a thioester bond. If we consider the citrate synthase reaction without thioester hydrolysis, i.e. adding acetate to oxaloacetate directly, we see that the reaction is intrinsically unfavorable

|citrate_synthase_noA|_

.. |citrate_synthase_noA| replace:: acetate + oxaloacetate ⇌ citrate
.. _citrate_synthase_noA: http://equilibrator.weizmann.ac.il/search?query=acetate+%2B+oxaloacetate+%E2%87%8C+citrate

From this we learn that the formation of a thioester on acetyl-CoA in the pyruvate dehydrogenase reaction above essentially "carries forward" the energy output oxidative decarboxylation into the citrate synthase step. Yet another example of an intrinsically unfavorable chemical reaction that biology enables through clever energetic activation.

Since we added two carbons from acetyl-CoA to oxaloacetate to form citrate, citrate must have 6 carbons. This means that it must be decarboxylated one time to form α-ketoglutarate. Due to the mechanism of enzymes carrying out oxidative decarboxylation, `citrate must be isomerized to isocitrate <https://www.ncbi.nlm.nih.gov/books/NBK22427/>`_. It can then undergo oxidative decarboxylation 

|iso_dehy|_

.. |iso_dehy| replace:: isocitrate + NAD\ :sup:`+` ⇌ α-ketoglutarate + CO\ :sub:`2` + NADH
.. _iso_dehy: http://equilibrator.weizmann.ac.il/search?query=isocitrate+%2B+NAD%2B+%E2%87%8C+alpha-ketoglutarate+%2B+CO2+%2B+NADH

to form α-KG and close the cycle. The net reaction of this arm of the cycle is 

|OAA_2_AKG_NET|_

.. |OAA_2_AKG_NET| replace:: oxaloacetate + acetyl-CoA + NAD\ :sup:`+` + H\ :sub:`2`\ O ⇌ α-ketoglutarate + CO\ :sub:`2` + NADH + CoA
.. _OAA_2_AKG_NET: http://equilibrator.weizmann.ac.il/search?query=oxaloacetate+%2B+acetyl-CoA+%2B+NAD%2B+%2B+H2O+%E2%87%8C+alpha-ketoglutarate+%2B+CO2+%2B+NADH+%2B+CoA

and is quite favorable with a Δ\ :sub:`r`\ G'\ :sup:`m` ≈ -40 kJ/mol. Putting together the two arms of the TCA cycle, we see that acetyl-CoA is added to oxaloacetate, oxidatively decarboxylated once to α-ketoglutarate, which is subsequently oxidatively decarboxylated once and oxidized twice to remake oxaloacetate. The two decarboxylations and 4 reduced electron carriers formed (3 NADH and 1 FADH\ :sub:`2`) account for the 2 carbons and 8 e- introduced to the cycle by the acetyl group of acetyl-CoA.

.. figure:: _static/_images/tca_cycle.svg
   :alt: The TCA Cycle
   :align: center

   The TCA cycle and anaplerotic reactions. Red titles denote the abbreviated names of the enzymes catalyzing the labeled reaction in E. coli.

.. [1] The enzyme that catalyzes this reaction, succinyl-CoA synthetase, is unfortunately named for the reverse direction of the reaction.

.. [2] By nonzero flux we mean that the cycle moves in a particular direction. For example the cell presumably "wants" the cycle to move in the direction of acetyl-CoA metabolism and energy production. In equilibrium the forward and reverse fluxes are definitionally equal and the cycle carries no net flux (forward - reverse = 0), meaning (in this case) that acetyl-CoA is neither created or destroyed by the cycle. However if we assume there is some acetyl-CoA degradation happening through the TCA cycle it must therefore carry net flux and not be in equilibrium. Above we implied that using P and P' (a pathway P and its reverse) for both arms of the TCA cycle would imply that we are in equilibrium (because the concentrations of all pathway intermediates are the same in both arms). But if we are in equilibrium there cannot be any net flux.

Anaplerotic Reactions
----------------------------------------------------------
α-Ketoglutarate to oxaloacetate are constantly being consumed by transamination reactions to make various amino acids. However, they are both also intermediates of the TCA cycle, meaning that they are neither created or destroyed by the action of the cycle. So we have a conundrum! If α-KG and OAA are removed from the cycle by transamination reactions but never replenished, their concentrations will eventually dwindle to 0, the TCA cycle would come to a halt, amino acids could no longer be made and the cell would die. [3]_

Cells sidestep this problem by continually replenishing TCA cycle intermediates through "`anaplerotic reactions <https://en.wikipedia.org/wiki/Anaplerotic_reactions>`_" from the `Greek meaning <https://en.wiktionary.org/wiki/anaplerotic>`_ "filling up" or "replenishing" reactions. Several of these reactions are illustrated in the TCA cycle figure above. For example, the pyruvate carboxylase reaction

|pyr_carboxylase|_

.. |pyr_carboxylase| replace:: pyruvate + ATP + CO\ :sub:`2` + H\ :sub:`2`\ O ⇌ oxaloacetate + ADP + Pi
.. _pyr_carboxylase: http://equilibrator.weizmann.ac.il/search?query=Pyruvate+%2B+ATP+%2B+CO2+%2B+H2O+%3C%3D%3E+Oxaloacetate+%2B+ADP+%2B+Pi

replenishes oxaloacetate from pyruvate. Notice that this actually replenishes all the TCA cycle intermediates, including α-KG, because oxaloacetate will be quickly converted into those metabolites through the action of the TCA cycle. In fact, there is no anaplerotic reaction that directly produces α-KG - cells rely on the TCA cycle to do this for them.

.. [3] In technical terms we would say that this configuration - where transamination happens at a constant nonzero rate but cycle intermediates are not replenished - has no nonzero "steady-state". By this we mean that there is no way to arrange this system where the flux through the TCA cycle is greater than 0 for an extended period of time. Try to convince yourself of this. 

pH Dependence
----------------------------------------------------------
The TCA cycle reaction that makes oxaloacetate is called `malate dehydrogenase <http://pdb101.rcsb.org/motm/154>`_ because, well, it dehydrogenates malate

|malate_dehy|_

.. |malate_dehy| replace:: malate + NAD\ :sup:`+` ⇌ oxaloacetate + NADH
.. _malate_dehy: http://equilibrator.weizmann.ac.il/search?query=malate+%2B+NAD%2B+%E2%87%8C+oxaloacetate+%2B+NADH

this reaction is problematic for two reasons. First of all, it is not very favorable, having a Δ\ :sub:`r`\ G'\ :sup:`m` ≈ +30 kJ/mol at pH 7. Using more plausible concentrations for [NAD\ :sup:`+`] ≈ 3 mM and [NADH] ≈ 0.08 mM helps but not enough to make the forward direction favorable (try this for yourself). Moreover, `malate <https://en.wikipedia.org/wiki/Malic_acid>`_ and `oxaloacetate <https://en.wikipedia.org/wiki/Oxaloacetic_acid>`_ have different pKas on their two carboxylic acid groups, meaning that Δ\ :sub:`r`\ G' of the malate dehydrogenase reaction will depend on the pH. Try using the pH slider and pH graphing utility on eQuilibrator to see how Δ\ :sub:`r`\ G' and Δ\ :sub:`r`\ G'° depend on the pH. 

The pKa of an acid is the pH at which that acidic group is 50% protonated (also 50% deprotonated). Overall, the pKas associated with oxaloacetate (pKa = 2.2, 3.9) are lower than those associated with malate (pKa = 3.4, 5.2). This means that as the pH goes down from 7 (i.e. becomes more acidic) and approaches the higher pKa of malate, malate reaches a pH where it has multiple populated protonation states but oxaloacetate does not. [4]_ As a result of this effect, we find that lowering the pH (making the environment more acidic) makes the Δ\ :sub:`r`\ G' more positive, favoring the reverse reaction even more. 

Some organisms maintain a pH of 6 in their cytosol (`Noor et. al, 2014 <refs.html>`_). In these conditions Δ\ :sub:`r`\ G'° = +36 kJ/mol. What would it take to make this reaction flow in the direction of oxaloacetate (i.e. make the Δ\ :sub:`r`\ G' negative) if we assume [NAD\ :sup:`+`] ≈ 3 mM and [NADH] ≈ 0.08? We know that 

.. math::
	\begin{eqnarray}
	\Delta_r G' &=& \Delta_r G'^{\circ} + RT \ln{Q} \\
	&=& 36 \frac{kJ}{mol} + RT \ln{\left( \frac{[NADH][oxaloacetate]}{[NAD^+][malate]} \right)} \\
	&=& 36 \frac{kJ}{mol} + RT \ln{\left( \frac{0.08 mM \times [oxaloacetate]}{3 mM \times [malate]} \right)} \\
	&=& 36 \frac{kJ}{mol} + RT \left( \ln{\left( \frac{0.08 mM}{3 mM} \right)} + \ln{\left(\frac{[oxaloacetate]}{[malate]} \right)} \right) \\
	\end{eqnarray}

Since R = 8.315 x 10\ :sup:`-3` kJ/mol/K and we assume a temperature of T = 298.15 K, RT ln(0.08/3) ≈ -9 kJ/mol. Therefore, in order for Δ\ :sub:`r`\ G' < 0 we need 

.. math::
	\begin{eqnarray}
	RT  \ln{\left(\frac{[oxaloacetate]}{[malate]} \right)} &<& -27 \frac{kJ}{mol} \\
	\frac{[oxaloacetate]}{[malate]} &<& \exp{\left( \frac{-27 \frac{kJ}{mol}}{RT} \right)} \\
	\frac{[oxaloacetate]}{[malate]} &<& 1.8\times10^{-5}\\
	54000 \times [oxaloacetate] &<& [malate]
	\end{eqnarray}

In other words, we need at least ≈54000 times more malate than oxaloacetate to make this reaction flow towards oxaloacetate spontaneously. A 5x10\ :sup:`4`\ -fold difference is `biologically plausibile <http://book.bionumbers.org/what-are-the-concentrations-of-free-metabolites-in-cells/>`_, but we must remember that the previous reaction in the TCA cycle, the one that makes malate, must also be thermodynamically feasible for the TCA cycle to work. Forcing a very high malate concentration will strongly constrain the operation of the TCA cycle! 

One solution employed by some organisms is to use a different, higher potential `quinone <https://en.wikipedia.org/wiki/Coenzyme_Q10>`_ electron carrier (`Noor et. al, 2014 <refs.html>`_). This has the effect of increasing the intrinsic favorability of the reaction so that such extreme malate concentrations are not required. Since the quinone has a higher potential, however, less energy is released and less ATP can be formed when electrons carried by the quinone are ultimately donated to O\ :sub:`2`. 

.. [4] Remember that the pH is a log base 10 scale, meaning that the 1.3 pH point difference between the higher pKa of malate and that of oxaloacetate indicates a very large difference more than 10-fold difference in the abundance of protonated carboxylic acid at pH 5.2. 
