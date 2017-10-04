-------------------------------------
Transamination for Making Amino Acids
-------------------------------------

Proteins are made of amino acids and amino acids have amine groups (-NH\ :sub:`2`). Here we consider two related questions: (a) where do the amines come from and (b) how are different amino acids made? The most common way that cells make new amines is to incorporate ammonia into glutamate, forming glutamine

|glutamine_synth|_

.. |glutamine_synth| replace:: Glutamate + NH\ :sub:`3` ⇌ Glutamine + H\ :sub:`2`\ O
.. _glutamine_synth: http://equilibrator.weizmann.ac.il/search?query=Glutamate+%2B+NH3+%3C%3D%3E+Glutamine+%2B+H2O

Using eQuilibrator we can see that this reaction is highly unfavorable with a positive Δ\ :sub:`r`\ G'm around 28 kJ/mol. In fact, the enzyme that catalyzes this reaction (`glutamine synthetase <http://equilibrator.weizmann.ac.il/enzyme?ec=6.3.1.2>`__) `couples the reaction to ATP hydrolysis <http://equilibrator.weizmann.ac.il/reaction?reactantsId=C00002&reactantsCoeff=-1&reactantsName=ATP&reactantsPhase=aqueous&reactantsConcentration=0.001&reactantsId=C00008&reactantsCoeff=1&reactantsName=ADP&reactantsPhase=aqueous&reactantsConcentration=0.001&reactantsId=C00009&reactantsCoeff=1&reactantsName=Orthophosphate&reactantsPhase=aqueous&reactantsConcentration=0.001&reactantsId=C00014&reactantsCoeff=-1&reactantsName=NH3&reactantsPhase=aqueous&reactantsConcentration=0.001&reactantsId=C00025&reactantsCoeff=-1&reactantsName=L-Glutamate&reactantsPhase=aqueous&reactantsConcentration=0.001&reactantsId=C00064&reactantsCoeff=1&reactantsName=L-Glutamine&reactantsPhase=aqueous&reactantsConcentration=0.001&ph=7.000000&pmg=14.000000&ionic_strength=0.100000&e_reduction_potential=0.000000&max_priority=0&mode=BA&query=ATP%20%2B%20NH3%20%2B%20L-Glutamate%20%3D%20ADP%20%2B%20Orthophosphate%20%2B%20L-Glutamine>`__ in order to make it favorable (Δ\ :sub:`r`\ G'\ :sup:`m` = -15 kJ / mol). Glutamine is an amino acid with two nitrogen atoms - one derived from glutamate and another from ammonia. Once glutamine is formed, cells use an enzyme (`glutamine:oxoglutarate aminotransferase <http://equilibrator.weizmann.ac.il/enzyme?ec=1.4.1.13>`__) to transfer the ammonia-derived nitrogen to another molecule (an alpha-keto acid called 2-oxoglutarate) to end up with two glutamates.

`NADPH + 2-Oxoglutarate + Glutamine ⇌ NADP+ + 2 Glutamate <http://equilibrator.weizmann.ac.il/reaction?query=NADP++%2B+2+L-Glutamate+%3C%3D%3E+NADPH+%2B+2-Oxoglutarate+%2B+L-Glutamine&ph=7.0&ionic_strength=0.1&reactantsCoeff=1.0&reactantsId=C00005&reactantsName=NADPH&reactantsConcentration=1&reactantsConcentrationPrefactor=0.001&reactantsPhase=aqueous&reactantsCoeff=-1.0&reactantsId=C00006&reactantsName=NADP+&reactantsConcentration=1&reactantsConcentrationPrefactor=0.001&reactantsPhase=aqueous&reactantsCoeff=-2.0&reactantsId=C00025&reactantsName=L-Glutamate&reactantsConcentration=1&reactantsConcentrationPrefactor=0.001&reactantsPhase=aqueous&reactantsCoeff=1.0&reactantsId=C00026&reactantsName=2-Oxoglutarate&reactantsConcentration=1&reactantsConcentrationPrefactor=0.001&reactantsPhase=aqueous&reactantsCoeff=1.0&reactantsId=C00064&reactantsName=L-Glutamine&reactantsConcentration=1&reactantsConcentrationPrefactor=0.001&reactantsPhase=aqueous&max_priority=0&submit=Reverse>`__

This reaction is highly favorable, in part because it couples to the oxidation of NADPH to the reduction of a carbonyl group on 2-oxoglutarate (see the images below and `Bar-Even et al., 2012 <refs.html>`_). So by tracing the pathway you can see that it has the net effect of converting ammonia to amines stored on glutamate. [#tra1]_ 

The most common way that amino acids are made is by transferring amines from glutamate to the right alpha-keto acid acceptor (as there are several such alpha keto acid compounds). 

.. figure:: _static/_images/alphaketoglutarate.png
   :alt: α-ketoglutaric acid or 2-oxoglutarate
   :align: center

   α-ketoglutaric acid or 2-oxoglutarate is called α-keto acid because it has a carbonyl group (ketone) adjacent to (in the α position relative to) a terminal carboxylic acid.

.. todo::
    charges on carboxylic acids in figures???

.. figure:: _static/_images/glutamate.png
   :alt: L-glutamate
   :align: center

   L-glutamate: notice that glutamate is identical to α-ketoglutaric acid except the the α carbonyl group of has been replaced by an amine. 

Phenylalanine, for example, can be made by transferring an amine to phenylpyruvate

`Glutamate + Phenylpyruvate ⇌ 2-Oxoglutarate + Phenylalanine <http://equilibrator.weizmann.ac.il/reaction?query=2-Oxoglutarate+%2B+L-Phenylalanine+%3C%3D%3E+L-Glutamate+%2B+Phenylpyruvate&ph=7.0&ionic_strength=0.1&reactantsCoeff=1.0&reactantsId=C00025&reactantsName=L-Glutamate&reactantsConcentration=1&reactantsConcentrationPrefactor=0.001&reactantsPhase=aqueous&reactantsCoeff=-1.0&reactantsId=C00026&reactantsName=2-Oxoglutarate&reactantsConcentration=1&reactantsConcentrationPrefactor=0.001&reactantsPhase=aqueous&reactantsCoeff=-1.0&reactantsId=C00079&reactantsName=L-Phenylalanine&reactantsConcentration=1&reactantsConcentrationPrefactor=0.001&reactantsPhase=aqueous&reactantsCoeff=1.0&reactantsId=C00166&reactantsName=Phenylpyruvate&reactantsConcentration=1&reactantsConcentrationPrefactor=0.001&reactantsPhase=aqueous&max_priority=0&submit=Reverse>`__

Similarly, aspartate can be made by transferring an amine to oxaloacetate 

`Glutamate + Oxaloacetate ⇌ 2-Oxoglutarate + Aspartate <http://equilibrator.weizmann.ac.il/reaction?query=2-Oxoglutarate+%2B+L-Aspartate+%3C%3D%3E+L-Glutamate+%2B+Oxaloacetate&ph=7.0&ionic_strength=0.1&reactantsCoeff=1.0&reactantsId=C00025&reactantsName=L-Glutamate&reactantsConcentration=1&reactantsConcentrationPrefactor=0.001&reactantsPhase=aqueous&reactantsCoeff=-1.0&reactantsId=C00026&reactantsName=2-Oxoglutarate&reactantsConcentration=1&reactantsConcentrationPrefactor=0.001&reactantsPhase=aqueous&reactantsCoeff=1.0&reactantsId=C00036&reactantsName=Oxaloacetate&reactantsConcentration=1&reactantsConcentrationPrefactor=0.001&reactantsPhase=aqueous&reactantsCoeff=-1.0&reactantsId=C00049&reactantsName=L-Aspartate&reactantsConcentration=1&reactantsConcentrationPrefactor=0.001&reactantsPhase=aqueous&max_priority=0&submit=Reverse>`__

If we inspect these reactions in eQuilibrator, we find that they have Δ\ :sub:`r`\ G'\ :sup:`m` values that are very close to 0 kJ / mol. This makes intuitive sense because both reactions involve the transfer of a group (an amine) from one molecule to the same position on another molecule (both the substrate and product are amino acids). In other words, these reactions are not intrinsically thermodynamically driven in either direction. Rather, it is the extremely high concentration of glutamate [#tra2]_  that drives these reactions forward (in the direction of making various amino acids). Try setting the glutamate concentration to 50 mM using eQuilibrator to see what happens to the Δ\ :sub:`r`\ G' of these transamination reactions. 

.. [#tra1] Exercise for the reader: explain why this reaction is a reduction reaction, i.e. why it requires electrons from NADPH? 

.. [#tra2] `Upwards of 50 mM in many cell types <http://book.bionumbers.org/what-are-the-concentrations-of-free-metabolites-in-cells/>`__

