The TCA Cycle
==========================================================

As we've discussed elsewhere (in the sections on `thioesters <thioester.html>`_ and `fatty acid metabolism <fatty_acid_met.html>`_) acetyl-CoA is a very central molecule in metabolism. 

.. figure:: _static/_images/accoa.png
   :alt: Acetyl-CoA
   :align: center

   Acetyl-CoA is a very central molecule in metabolism and is produced during the metabolism of sugars, alkanes, fatty acids and some amino acids like leucine. 

Indeed, after glucose is metabolized to pyruvate through `glycolysis <glycolysis.html>`_ acetyl-CoA is made from pyruvate by a very complex enzyme called `the pyruvate dehydrogenase complex <pdb101.rcsb.org/motm/153>`_. 

Pyruvate + CoA + NAD+ ⇌ acetyl-CoA + CO\ :sub:`2` + NADH

This reaction is an oxidative decarboxylation - pyruvate is oxidized and decarboxylated (a CO\ :sub:`2` is removed) to form acetate. The free-energy released by this transformation is substantial (Δ\ :sub:`r`\ G'\ :sup:`m` ≈ -80 kJ/mol) and is harnessed to make a thioester bond on acetyl-CoA. You can see directly how favorable this transformation is by omitting the CoA from the reaction and entering it into eQuilibrator:

Pyruvate + H\ :sub:`2`\ O + NAD+ ⇌ acetate + CO\ :sub:`2` + NADH

.. todo:: mention why we need water in a footnote?

We will not go into detail about this complex enzyme here. Rather, we bring it up to highlight two important concepts. First, that oxidative decarboxylations are very favorable reactions, meaning that NADH can be produced by decarboxylating acids like pyruvate. We will return to this theme as we discuss the TCA cycle below. Second, that acetyl-CoA is a major point of convergence in central metabolism, bringing together the catabolism of glucose, amino acids and lipids. But how is acetyl-CoA itself metabolized? 

Metabolism of Acetyl-CoA
----------------------------------------------------------
We will approach the metabolism of acetyl-CoA through a sort of back door by reminding ourselves about `transamination <transamination.html>`_. Remember that many amino acids are made by transamination of an appropriate alpha-keto acid acceptor. For example, glutamate can be made transaminating a five carbon dicarboxylic acid called α-ketoglutarate (also known as 2-oxoglutarate). 

.. figure:: _static/_images/alphaketoglutarate.png
   :alt: α-ketoglutaric acid or 2-oxoglutarate
   :align: center

   α-ketoglutarate or 2-oxoglutarate is called α-keto acid because it has a carbonyl group (ketone) adjacent to (in the α position relative to) a terminal carboxylic acid.

Glutamate is a proteogenic amino acid in its own right, but also a precursor for the synthesis of several amino acids including glutamine, proline and arginine (ref). Similarly, aspartate can be made by transamination of a four carbon dicarboxylic acid called oxaloacetate. Aspartate is a precursor for the production of asparagine, lysine, arginine, methionine, threonine and isoleucine (ref).


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

In order to make amino acids, therefore, cells need to make α-ketoglutarate and oxaloacetate continually. One way to think about the tricaboxylic acid cycle (TCA cycle) is as a cycle for the production of ATP through metabolism of acetyl-CoA that also continually produces both α-ketoglutarate and oxaloacetate. Let's consider the net reaction of the TCA cycle starting from acetyl-CoA

Acetyl-CoA + 3 H2O ⇌ 2 CO\ :sub:`2` + CoA + 8 e- 

Notice that both α-ketoglutarate and oxaloacetate are absent from the net reaction because they are intermediates of the cycle - they are neither created or destroyed in cycle, but rather serve as co-substrates for the metabolism of acetyl-CoA. We can see that the net effect of the TCA cycle is to decarboxylate acetyl-CoA twice to make  CO\ :sub:`2` and electrons. These electrons have an average potential of E'\ :sup:`m` ≈ -350 mV, meaning that they can be donated to NAD+ to form NADH (which has a more positive E'\ :sup:`m` ≈ -330 mV). In the TCA cycle, electrons are mostly donated to NAD+ to form NADH. They are subsequently donated to O\ :sub:`2` in a manner that is coupled to the formation of ATP via the `respiratory electron transport chain <respiration.html>`_. 

You might notice that the E'\ :sup:`m` values above are very close to each other, meaning small variations in concentrations might make these redox reactions infeasible. In practice this problem is "solved" in two ways. First, as we will see, not all the electrons are donated to NADH. Two of the 8 electrons produced by oxidizing acetyl-CoA are donated to a higher (more positive) potential acceptor FAD+ (E'\ :sup:`m` ≈ -240 mV). Second, the NAD+ and NADH concentrations are not exactly 1 mM. Typically [NAD+] is about 30-50 times greater than [NADH] - try using eQuilibrator to see how this concentration ratio affects the E' value for NAD+ reduction to NADH.

Constructing the TCA Cycle
----------------------------------------------------------
We are left with the task of building a cycle for the metabolism of acetyl-CoA to CO2 and electrons that passes produces the central amine acceptors α-ketoglutarate and oxaloacetate "on the way." 

α-Ketoglutarate to Oxaloacetate
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Let's first consider what it would take to make oxaloacetate from alpha-ketoglutarate.

alpha-ketoglutarate + 2 H2O ⇌ oxaloacetate + CO\ :sub:`2` + 6 e-

As we said, α-ketoglutarate (α-KG for short) contains 5 carbon atoms and oxaloacetate (OAA for short) contains 4. Therefore, production of OAA from α-KG requires the decarboxylation of α-KG - i.e. the removal of a carbon atom. As we foreshadowed above, this proceeds through the mechanism of oxidative decaboxylation of α-KG to make succinyl-CoA

alpha-ketoglutarate + CoA + NAD+ ⇌ succinyl-CoA + CO\ :sub:`2` + NADH

The thioester on succinyl-CoA is `approximately energetically equivalent to ATP <thioester.html>`_, which explains how the next reaction step manages to make ATP while hydrolysing the thioester. [1]_

succinyl-CoA + ADP + Pi ⇌ succinate + CoA + ATP

Succinate has 4 carbons, like OAA, but is more reduced - having 4 more electrons. So it must be oxidized twice to make oxaloacetate. 

succinate + NAD+ + FAD + H2O ⇌ oxaloacetate + NADH + FADH2

Notice that two electrons are donated to NAD+ and two are donated to a similar, but higher-potential donor called FAD (as discussed in above). Altogether, this arm of the TCA cycle has a net reaction of

alpha-ketoglutarate + 2 H2O + 2 NAD+ + FAD ⇌ oxaloacetate + CO2 + 2 NADH + FADH2

Notice also that this net reaction has a Δ\ :sub:`r`\ G'\ :sup:`m` ≈ 0. As discussed above, this can be remedied by setting the NAD+ and NADH concentrations to more physiologically relevant values. For example, `measurements in E. coli <http://book.bionumbers.org/what-are-the-concentrations-of-free-metabolites-in-cells/>`_ give [NAD+] ≈ 3 mM and [NADH] ≈ 0.08 mM. Try using these values to calculate Δ\ :sub:`r`\ G' in eQuilibrator - does this help resolve the problem?

Oxaloacetate back to α-Ketoglutarate
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
First of all - it is important to note that the TCA cycle can't possibly follow the same reaction scheme in both directions. If the TCA cycle used one reaction sequence from α-ketoglutarate to oxaloacetate and then the same sequence in the reverse direction to get back to α-ketoglutarate, this would be a closed cycle with nonzero flux [2]_ - a perpetual motion machine in violation of the first law of thermodynamics. Another, perhaps even simpler reason that the TCA cycle can't use the same reaction sequence in the reverse direction is that we haven't yet seen how acetyl-CoA is incorporated into the cycle. Acetyl-CoA needs to find it's way into the cycle in order for the TCA cycle to metabolize acetyl-CoA, after all! 

Acetyl-CoA gets into the TCA cycle is through the citrate synthase reaction

acetyl-CoA + oxaloacetate + H2O ⇌ CoA + citrate

which adds acetyl-CoA to oxaloacetate to form the tricarboxylic acid citrate (after which the cycle is named). 

.. figure:: _static/_images/citrate.png
   :alt: Citrate
   :align: center

   Citrate is the a tricarboxylic acid after which the TCA cycle is named.

The citrate synthase reaction is quite favorable (Δ\ :sub:`r`\ G'\ :sup:`m` ≈ -35 kJ/mol) due to the hydrolysis of a thioester bond. If we consider the citrate synthase reaction without thioester hydrolysis, i.e. adding acetate to oxaloacetate directly, we see that the reaction is intrinsically unfavorable

acetate + oxaloacetate ⇌ citrate

From this we learn that the formation of a thioester on acetyl-CoA in the pyruvate dehydrogenase reaction above essentially "carries forward" the energy output oxidative decarboxylation into the citrate synthase step. Yet another example of an intrinsically unfavorable chemical reaction that biology enables through clever energetic activation.

Since we added two carbons from acetyl-CoA to oxaloacetate to form citrate, citrate must have 6 carbons. This means that it must be decarboxylated one time to form α-ketoglutarate. Due to the mechanism of enzymes carrying out oxidative decarboxylation, `citrate must be isomerized to isocitrate <https://www.ncbi.nlm.nih.gov/books/NBK22427/>`_. It can then undergo oxidative decarboxylation 

isocitrate + NAD+ ⇌ α-ketoglutarate + CO2 + NADH

to form αKG and close the cycle. The net reaction of this arm of the cycle is 

oxaloacetate + acetyl-CoA + NAD+ + H2O ⇌ α-ketoglutarate + CO2 + NADH + CoA

and is quite favorable with a Δ\ :sub:`r`\ G'\ :sup:`m` ≈ -40 kJ/mol. Putting together the two arms of the TCA cycle, we see that acetyl-CoA is added to oxaloacetate, oxidatively decarboxylated once to α-ketoglutarate, which is subsequently oxidatively decarboxylated once and oxidized twice to remake oxaloacetate. The two decarboxylations and 4 reduced electron carriers formed (3 NADH and 1 FADH2) account for the 2 carbons and 8 e- introduced to the cycle by the acetyl group of acetyl-CoA.

.. [1] The enzyme that catalyzes this reaction, succinyl-CoA synthetase, is unfortunately named for the reverse direction of the reaction.

.. [2] By nonzero flux we mean that the cycle moves in a particular direction. For example the cell presumably "wants" the cycle to move in the direction of acetyl-CoA metabolism and energy production. In equilibrium the forward and reverse fluxes are definitionally equal and the cycle carries no net flux (forward - reverse = 0), meaning (in this case) that acetyl-CoA is neither created or destroyed by the cycle. However if we assume there is some acetyl-CoA degradation happening through the TCA cycle it must therefore carry net flux and not be in equilibrium. Above we implied that using P and P' (a pathway P and its reverse) for both arms of the TCA cycle would imply that we are in equilibrium (because the concentrations of all pathway intermediates are the same in both arms). But if we are in equilibrium there cannot be any net flux.

Anaplerotic Reactions
----------------------------------------------------------
As discussed above and in the section on `transamination <transamination.html>`_ α-ketoglutarate to oxaloacetate are constantly being consumed by transamination reactions to make various amino acids. However, we also discussed how they are both intermediates of the TCA cycle meaning that they are neither created or destroyed by the action of the cycle. So we have a problem! If α-KG and OAA are removed from the cycle by transamination reactions but never replenished then their concentrations will eventually dwindle to 0 - the TCA cycle would stop, amino acids could no longer be made and the cell would die. [3]_

Cells sidestep this problem by continually replenishing TCA cycle intermediates through "`anaplerotic reactions <https://en.wikipedia.org/wiki/Anaplerotic_reactions>`_" from the `Greek meaning <https://en.wiktionary.org/wiki/anaplerotic>`_ "filling up" or "replenishing" reactions. For example, the pyruvate carboxylase reaction

pyruvate + ATP + CO2 ⇌ oxaloacetate + ADP + Pi

replenishes oxaloacetate from pyruvate. Notice that this actually replenishes all the TCA cycle intermediates, including α-KG, because oxaloacetate will be quickly converted into those metabolites through the action of the TCA cycle. In fact, there is no anaplerotic reaction that directly produces α-KG - cells rely on the TCA cycle to do this for them.

.. [3] In technical terms we would say that this configuration - where transamination happens at a constant nonzero rate but cycle intermediates are not replenished - has no nonzero "steady-state". By this we mean that there is no way to arrange this system where the flux through the TCA cycle is greater than 0. 

pH Dependence
----------------------------------------------------------
Add content here

