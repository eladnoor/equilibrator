Thioester Bonds are Energetic Too
==========================================================

Phosphoanhydride bonds like those in ATP are not the only energy-carrying bonds used by biology. In fact, phosphate was likely scarce during the early history of Earth and some researchers hypothesize that `thioester bonds <https://en.wikipedia.org/wiki/Thioester>`_ - bonds between a sulfur and a carbonyl carbon (R–S–CO–R') - were used as energy currency during the evolution of early life (Goldford et al.).

.. figure:: _static/_images/accoa.png
   :alt: Acetyl-CoA
   :align: center

   Acetyl-CoA is the protypical example of a biological molecule containing a thioester bond.

So how much energy is there in a thioester bond? We can use eQuilibrator to examine the hydrolysis of the common biological thioester, acetyl-CoA

|thio_hydrolysis|_

.. |thio_hydrolysis| replace:: Acetyl-CoA + H\ :sub:`2`\ O ⇌ Acetate + CoA
.. _thio_hydrolysis: http://equilibrator.weizmann.ac.il/reaction?reactantsId=C00024&reactantsCoeff=-1&reactantsName=Acetyl-CoA&reactantsPhase=aqueous&reactantsConcentration=0.001&reactantsId=C00033&reactantsCoeff=1&reactantsName=Acetate&reactantsPhase=aqueous&reactantsConcentration=0.001&reactantsId=C00010&reactantsCoeff=1&reactantsName=CoA&reactantsPhase=aqueous&reactantsConcentration=0.001&reactantsId=C00001&reactantsCoeff=-1&reactantsName=H2O&reactantsPhase=liquid&reactantsConcentration=1&ph=7.000000&pmg=14.000000&ionic_strength=0.100000&e_reduction_potential=0.000000&max_priority=0&mode=BA&query=acetyl-CoA%20%3D%20acetate%20%2B%20CoA%20%2B%20h2o

We find that this hydrolysis reaction has a Δ\ :sub:`r`\ G'\ :sup:`m` around -40 kJ/mol - very similar to the energetic scale of ATP hydrolysis. In fact, the thioester bond is sometimes exchanged with a phosphoester bond (like the one in ATP) during metabolism of acetyl-CoA

`Acetyl-CoA + Pi ⇌ Acetyl-phosphate + CoA <http://equilibrator.weizmann.ac.il/reaction?query=acetyl-CoA+%2B+pi+%3D%3E+acetyl-phosphate+%2B+CoA&ph=7.0&ionic_strength=0.1&reactantsCoeff=-1&reactantsId=C00024&reactantsName=Acetyl-CoA&reactantsConcentration=1&reactantsConcentrationPrefactor=0.001&reactantsPhase=aqueous&reactantsCoeff=-1&reactantsId=C00009&reactantsName=Pi&reactantsConcentration=10&reactantsConcentrationPrefactor=0.001&reactantsPhase=aqueous&reactantsCoeff=1&reactantsId=C00227&reactantsName=Acetyl+phosphate&reactantsConcentration=1&reactantsConcentrationPrefactor=0.001&reactantsPhase=aqueous&reactantsCoeff=1&reactantsId=C00010&reactantsName=CoA&reactantsConcentration=1&reactantsConcentrationPrefactor=0.001&reactantsPhase=aqueous&max_priority=0&submit=Update>`_

Notice that this reaction has a positive Δ\ :sub:`r`\ G'\ :sup:`m`  ≈ 10 kJ/mol, meaning that it would flow in the reverse direction (forming acetyl-CoA) if all the reactants had 1 mM concentrations. However, concentrations are of course not exactly 1 mM in cells. Cellular concentrations of inorganic phosphate are, for example, typically closer to 10 mM [1]_.

Try using eQuilibrator to set the phosphate concentration to 10 mM, leaving everything else at 1 mM. Notice that the Δ\ :sub:`r`\ G'\ :sup:`m` value changed by about 6 kJ/mol. As an exercise: show with some simple math that a 10-fold change in one reactant concentration will always alter the Δ\ :sub:`r`\ G' by about 6 kJ/mol [2]_

.. [1] For typical concentrations in *E. coli* see `BioNumbers ID 105540 <http://bionumbers.hms.harvard.edu/bionumber.aspx?&id=105540&ver=3&trm=inorganic%20phosphate%20concentration>`_
.. [2] Hint: think about the formula for Δ\ :sub:`r`\ G' = Δ\ :sub:`r`\ G'° + RT ln Q.
