Glycolysis and Fermentation
==========================================================

The reason that glucose oxidation produces so much energy is that molecular oxygen has a strong "preference" to accept electrons (very high reduction potential). This preference of electrons to flow to oxygen is so powerful that it can be used to drive the formation of several molecules of ATP. This process is called oxidative phosphorylation. In anaerobic conditions, however, oxygen is not available (by definition) and it's impossible to drive oxidative phosphorylation by donating electrons to O2.

The anaerobic breakdown of glucose for energy production is called fermentation. Several approaches to fermenting glucose occur in nature, but here we'll discuss only one in this section: fermentation of glucose to lactate. Fermentation of glucose begins in the same way that complete oxidation of glucose does: by breaking glucose down into two pyruvate molecules through a process known as glycolysis:

Glucose ⇌ 2 Pyruvate

If you inspect this reaction in eQuilibrator, you'll notice that it is not "electron-balanced." That is, there are 4 fewer electrons in two pyruvate molecules than there are in one glucose molecule. We should write this reaction as a “half-reaction” 

Glucose ⇌ 2 Pyruvate + 4 e-

In order track the electrons explicitly. These electrons don't just float around the cell: they are carried by specific compounds called electron carriers. In glycolysis, the electron carrier of choice is nicotinamide adenine dinucleotide (NAD). To fix the electron imbalance, try clicking on the "Balance with NAD+/NADH" link.

The full glycolysis reaction, including electron carriers and ATP production is:

Glucose + 2 NAD+ + 2 ADP + 2 Pi ⇌ 2 Pyruvate + 2 NADH + 2 ATP + 2 H\ :sub:`2`\ O

Clearly glycolysis produces ATP without any oxygen. However, if you ran glycolysis over and over again to provide ATP, the amount of NADH in the cell would grow continually. The ratio of NAD+ concentration to NADH concentration - sometimes called the "redox state" - is very important to living cells because it affects the energetics of redox reactions, in particular in glycolysis. In order to keep this ratio constant (homeostasis), the electrons stored in NADH need to move on to a different molecule. Fermentation to lactate solves this problem by taking the electrons carried by NADH and giving them to pyruvate to form lactate:

2 Pyruvate + 2 NADH ⇌ 2 Lactate + 2 NAD+

Since lactate build-up is problematic (in much the same way that NADH build up is) the lactate is then excreted from the cell. Fermentation to lactate, then, produces 2 ATP from one glucose while forming 2 lactate:

Glucose + 2 ADP + 2 Pi ⇌ 2 Lactate + 2 ATP + 2 H\ :sub:`2`\ O

This neat biochemical trick enables organisms as varied as E. coli and humans to produce energy from sugars even when oxygen isn't available.

NADH as an Electron Carrier
----------------------------------------------------------

Many reactions in metabolism redox reactions (reduction/oxidation reactions) involving the transfer of electrons between molecules. One classic such reaction is the lactate dehydrogenase reaction described above. 

NADH + Pyruvate <=> NAD+ + Lactate

We can identify that this is a redox reaction from the presence of the two-electron carrier NAD+/NADH. One way to understand the energetics of redox reactions is to split them into “half-reactions.”

Pyruvate + 2e- <=> Lactate E'm ≈ -190 mV
NADH <=> NAD+ + 2e- E'm ≈ -330 mV

If we add these two half-reactions together we see that the electrons balance and cancel: the two half-cells produce a balanced chemical reaction. 

The energetics of these half-reactions can be characterized by a number called the “reduction potential,” given above as E'm. What do these reduction potentials mean? The -190 mV E'm given for converting pyruvate to lactate represents the amount of energy that each electron needs to carry in order to convert pyruvate into lactate assuming, as above, 1 mM concentrations for pyruvate and lactate (denoted by the m superscript). Here, for example, we have to donate 2 e- to pyruvate to form lactate and each electron must carry roughly -190 mV of potential. 

The reduction potential is relative value - in order to know whether the pyruvate/lactate half-reaction will be favorable, we need to know where electrons came from. We can tell that NADH is a suitable donor for this reaction because the electrons in NADH have a potential of -330 mV. Electrons flow towards more positive potentials, and so we can tell that donating electrons from NADH to pyruvate to form NAD+ and lactate will be favorable - each electron gains about 140 mV of potential. 

We can convert these mV values into more familiar kJ/mol using the Nernst equation ΔrG' = nFΔE’, where n is the number of electrons transferred and F is the Faraday constant (F ≈ 96.4x103 kJ/mV/mol). Reducing pyruvate using NADH (i.e. donating electrons from NADH to pyruvate) will therefore have a ΔrGm ≈ 2 x 96.4x10-3 (-330-190) ≈ -27 kJ/mol (consistent with our calculation above using eQuilibrator). 

There are many common biological electron carriers aside from NADH and its close cousin NADPH. These include ferredoxin, glutathione, quinones and various flavins. But NAD(P)H is the most common biological electron carrier for a number of reasons. Most importantly, molecular oxygen (O2) is abundant in our atmosphere and has a very high reduction potential (Nelson et al., 2008). Many biological electron carriers can spontaneously donate electrons to O2 and so must be “protected” from oxygen in various ways (e.g. by being buried deep inside proteins). As compared to other biological electron carriers, NAD(P)H is relatively insensitive to O2 (Nelson et al., 2008). NADH is also a two electron carrier, which is suited to many metabolic reactions where electrons are often transferred in pairs. Finally, the reduction potential of NAD(P)H around -320 mV is well-suited to many common biological transformations, which mainly alter the “redox state” (number of electrons) associated with carbon atoms (Bar-Even et al., 2012). 

Ethanol Fermentation
----------------------------------------------------------

Fermentation isn’t limited to making lactate - cells can produce anything that is “redox neutral” relative to their growth substrate (i.e. has the same total number of electrons as glucose). So long as the product has the same number of electrons as the substrate, the process doesn’t force accumulation of reduced electron carriers like NADH and is called fermentation. Ethanol is a famous and important fermentation product of glucose - one that is imbibed and used as a fuel additive worldwide. 

We can use eQuilibrator to check that ethanol production is in fact redox neutral (relative to glucose) by searching for the reaction producing ethanol from pyruvate

2 pyruvate + 4 e- ⇌ 2 ethanol + 2 CO2

This reaction is categorized as a “half-reaction” by eQuilibrator, meaning that there are excess electrons on one side of the reaction. The 4 electrons required can be withdrawn from NADH as discussed above (eQuilibrator does this automatically if you click the “Balance with NAD+/NADH” link).

2 Pyruvate + 2 NADH ⇌ 2 Ethanol + 2 CO2 + 2 NAD+

As with lactate, we see that production of ethanol from two pyruvate molecules perfectly balances the production of two pyruvate from glucose -- glucose breakdown produces two NADH and ethanol production consumes them both. As a result, the NADH cancel and the net reaction for ethanol fermentation contains no electron carriers.

Glucose ⇌ 2 Ethanol + 2 CO2

Indeed, production of ethanol and CO2 from pyruvate is quite favorable (ΔrG'm = -114 kJ / mol) which helps explain why yeast are so content to make large quantities of ethanol for our enjoyment. 

Mixtures of Fermentation Products
----------------------------------------------------------

Bacteria produce many different kinds of fermentation products. The only requirement is that the products of a fermentation pathway contain the same number of electrons as the substrates. Indeed, as we saw with ethanol fermentation, it’s possible to make multiple fermentation products (ethanol + CO2) so long as redox balance is preserved.

In mixed fermentation, cells produce (surprise!) a precise mixture of products that together are redox balanced with their substrate. E. coli, for example, can use this strategy to make an extra ATP by converting pyruvate to a 1:1:1 mixture of acetate, ethanol and formate. 

2 Pyruvate + 2 NADH + H\ :sub:`2`\ O <=> Acetate + Ethanol + 2 Formate + 2 NAD+

Formate is fairly toxic, so it is subsequently converted to CO2 and molecular hydrogen by an enzyme called formate hydrogenlyase

Formate <=> CO2 + H2

Considering the net reaction of this process in eQuilibrator 

2 Pyruvate + 2 NADH + H :sub:`2` O <=> Acetate + Ethanol + 2 CO2 + 2 NAD+ + 2 H2

We see that it has a ΔrG'm around -100 kJ/mol, more than sufficient to make an additional ATP. Altogether, the net reaction from glucose makes 3 ATP and is still quite favorable

Glucose + 3 ADP + 3 Pi <=> Acetate + Ethanol + 2 CO2 + 2 H2 + 3 ATP + 2 H\ :sub:`2`\ O

So how does E. coli get an extra ATP out of this transformation from pyruvate to acetate, ethanol and formate? To see how, consider the conversion of pyruvate to acetate, which takes place in the following 3 steps 

#. CoA + Pyruvate <=> Acetyl-CoA + Formate

#. Acetyl-CoA + Pi <=> Acetyl-Phosphate + CoA

#. ADP + Acetyl phosphate <=> ATP + Acetate

Acetate is formed via the intermediate of acetyl-CoA, which allows for the production of ATP (as discussed above). Notice, however, that this 3-step pathway is redox neutral - it does not consume any of the NADH that would have been produced in glycolysis to make pyruvate. So E. coli can’t ferment glucose entirely to acetate and formate because that would not be a redox neutral transformation, as you can see by balancing the net reaction on eQuilibrator.

Glucose + 2 NAD+ + H\ :sub:`2`\ O <=> 2 Acetate + 2 Formate + 2 NADH

Producing one ethanol for every acetate ensures that the overall fermentation pathway is redox balanced. This can be seen by considering how ethanol is made from pyruvate in this case 

#. CoA + Pyruvate <=> Acetyl-CoA + Formate

#. NADH + Acetyl-CoA <=> NAD+ + CoA + Acetaldehyde

#. NADH + Acetaldehyde <=> NAD+ + Ethanol

We notice that the production of ethanol from pyruvate in these three steps involves the consumption of 2 NADH. So if the cell makes sure to make one ethanol molecule for every acetate, then redox balance will be preserved and one extra ATP will be made for every two pyruvates metabolized - a balancing act indeed!

There are several other pathways of this sort that produce a defined mixture of fermentation products that are collectively redox neutral compared to glucose. The pathways are termed “mixed acid fermentation pathways” because they usually produce a mixture of acids (Kim and Gadd, 2008). But this name can be confusing, as in the case of 1:1:1 production of acetate, ethanol and formate because not all of the products are acids (e.g. ethanol is an alcohol). You can learn more about the variety of mixed fermentation pathways on Wikipedia and Biocyc. 
