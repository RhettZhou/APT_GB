# APT_GB
Help identify grain boundaries from atom probe datasets for mapping interfacial excesses and concentrations.

This app is used for atom probe tomography grain boundary analysis. The 3D atom probe tip is projected to a 2D concentration map. After applying filter and machine learning identification, grain boundary could be automatically detected and its interface excess and composition could be quantified. 

Grain boundaries are planar lattice defects that govern the properties of many types of polycrystalline materials. Hence, their structures have been investigated in great detail. However, much less is known about their chemical features, owing to the experimental difficulties to probe these features at the atomic length scale inside bulk material specimens. Atom probe tomography is a tool capable of accomplishing this task, with an ability to quantify chemical characteristics at a near-atomic scale.  Using APT data sets, we present here a machine-learning-based approach for the automated quantification of chemical features of grain boundaries. We trained a convolutional neural network using twenty thousand synthesized images of grain interiors, grain boundarie, or triple junctions. Such a trained convolutional neural network automatically detected the location of the grain boundaries from atom probe data set. Those grain boundaies are then subjected to compositional mapping and analysis, including revealing in-plane chemical patterns. 

The whole process has been interpreted in this software package. We provide graphical interfaces for all analysis procedures. The details for each step are documented in the manual slides with an synthetic data set as an example. 

- Please use the defaut folder for the installation.
- Users must also install Blender 2.92.0. Please also use the default folder.

Acknowledgement: Peter Felfer, Anna Ceguerra, Varvara Efremova, Ye Wei, Markus Kühbach, Huan Zhao, Florian Vogel, Reza Darvishi Kamachali, Baptise Gault, Gregory B Thompson, Dierk Raabe.

