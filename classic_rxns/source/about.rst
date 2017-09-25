About Us
==========================================================

.. _about-equilibrator:

About eQuilibrator
----------------------------------------------------------

eQuilibrator is a simple web interface designed to enable easy thermodynamic analysis of biochemical systems. eQuilibrator enables free-text search for biochemical compounds and reactions and provides thermodynamic estimates for both in a variety of conditions. Estimation of thermodynamic parameters (Δ\ :sub:`r`\ G and Δ\ :sub:`f`\ G) elucidates how much energy is required to drive a particular biochemical reaction and in which direction the reaction will flow in particular cellular conditions [#Alberty2003]_ [#Alberty2006]_.

Because experimental measurement of the free energy of formation (Δ\ :sub:`f`\ G°) of compounds is technically challenging, the vast majority of known metabolites have not been thermodynamically characterized. eQuilibrator uses a well-studied approximation of Δ\ :sub:`f`\ G called *group contribution* [#Mavrovouniotis1991]_ [#Jankowski2008]_ [#Noor2012]_, enabling thermodynamic analysis of many biochemical reactions and pathways.

Currently, eQuilibrator can provide estimates for many compounds in the KEGG database [#Kanehisa2000]_ (about 4500). Individual compounds and enzymes can be searched for by their common names ("water", "glucosamine", "hexokinase"), and reactions can be entered in a simple, free-text format ("ribulose bisphosphate + CO2 + water => 2 3-phosphoglycerate") that eQuilibrator parses automatically. eQuilibrator also allows manipulation of the conditions of a reaction - pH, ionic strength, and reactant and product concentrations - to help explore the thermodynamic landscape of a biochemical reaction.

eQuilibrator is a project of the `Milo Lab <http://www.weizmann.ac.il/plants/Milo/>`_ at the Weizmann Institute in Rehovot, Israel. If you have any thoughts or questions feel free to write us on the `eQuilibrator Users Google Group <https://groups.google.com/forum/#!forum/equilibrator-users>`_

.. _implementation:

Implementation of eQuilibrator
----------------------------------------------------------

eQuilibrator makes heavy use of open source libraries and frameworks and is open-source itself. You can find the eQuilibrator source code on
our GitHub repository called `equilibrator <https://github.com/eladnoor/equilibrator/>`_.

The eQuilibrator back-end is implemented in pure `Python <http://www.python.org>`_ using the `Django web framework <http://www.djangoproject.com/>`_ and several other open-source Python libraries, including `NumPy <http://numpy.scipy.org/>`_ and `Pyparsing <http://pyparsing.wikispaces.com/>`_.
The eQuilibrator user-interface is implemented using HTML, CSS and the `jQuery Javascript framework <http://jquery.com/>`_.

All thermodynamic data presented by eQuilibrator was generated using the Component Contribution method [`DOI: 10.1371/journal.pcbi.1003098 <http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1003098>`_] [#Noor2013]_ which was implemented also in Python and is open-source as well. You can find the source code on our GitHub repository called `component-contribution <https://github.com/eladnoor/component-contribution/>`_.

.. _acknowledgements:

Acknowledgements
----------------------------------------------------------

*   Avi Flamholz, Elad Noor, Ron Milo and Arren Bar-Even designed the eQuilibrator interface
*   Avi Flamholz and Elad Noor implemented eQuilibrator
*   Elad Noor, Avi Flamholz, Hulda Haraldsdóttir, Arren Bar-Even and Ron Milo designed and tuned eQuilibrator's version of the component contribution

References
----------------------------------------------------------

.. [#Alberty2003] R.A. Alberty, "Thermodynamics of biochemical reactions" (Hoboken N.J.: Wiley-Interscience, 2003)
.. [#Alberty2006] R.A. Alberty, "Biochemical Thermodynamics" (Hoboken, NJ, USA: John Wiley & Sons, Inc., 2006)
.. [#Mavrovouniotis1991] M.L. Mavrovouniotis, "Estimation of standard Gibbs energy changes of biotransformations" The Journal of Biological Chemistry (1991) 266(22):14440-14445
.. [#Jankowski2008] M.D. Jankowski et al., "Group Contribution Method for Thermodynamic Analysis of Complex Metabolic Networks" Biophysical Journal (2008) 95(3):1487-1499
.. [#Noor2012] E. Noor, A. Bar-Even, A. Flamholz, Y. Lubling, D. Davidi, R. Milo, "An integrated open framework for thermodynamics of reactions that combines accuracy and coverage" Bioinformatics (2012) 28:2037-2044
.. [#Noor2013] E. Noor, H.S. Haraldsdóttir, R. Milo, R.M.T. Fleming, "Consistent Estimation of Gibbs Energy Using Component Contributions" PLoS Comput Biol (2013) 9:e1003098
.. [#Kanehisa2000] M. Kanehisa, S. Goto, "KEGG: Kyoto Encyclopedia of Genes and Genomes" Nucleic Acids Research (2000) 28(1):27-30


