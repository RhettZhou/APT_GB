# APT_GB
Help identify grain boundaries from atom probe datasets for mapping interfacial excesses and compositions.

This application is used for atom probe tomography grain boundary analysis. The 3D atom probe tip is projected to a 2D composition map. After applying filter and machine learning identification, grain boundary could be automatically detected and its interfacial excess and composition could be quantified. 

Grain boundaries are planar lattice defects that govern the properties of many types of polycrystalline materials. Hence, their structures have been investigated in great detail. However, much less is known about their chemical features, owing to the experimental difficulties to probe these features at the atomic length scale inside bulk material specimens. Atom probe tomography is a tool capable of accomplishing this task, with an ability to quantify chemical characteristics at a near-atomic scale.  Using APT data sets, we present here a machine-learning-based approach for the automated quantification of chemical features of grain boundaries. We trained a convolutional neural network using twenty thousand synthesized images of grain interiors, grain boundarie, or triple junctions. Such a trained convolutional neural network automatically detected the location of the grain boundaries from atom probe data set. Those grain boundaies are then subjected to compositional mapping and analysis, including revealing in-plane chemical patterns. 

Flow chart of working principle

![image](https://user-images.githubusercontent.com/51905661/122538918-58613680-d027-11eb-8fe4-46f8dc899667.png)

Flow chart summarizes the steps of machine learning-enhanced mapping of grain boundary (GB) composition and interfacial excess from atom probe tomography (APT) data sets. ‘AP Suite’ refers to the atom prober's toolkit for data analysis workstations. Matlab is a numerical computing environment developed by MathWorks, Inc. Spyder is an open-source integrated development environment for scientific programming in Python. Blender is an open-source computer graphics software. POS is an APT file format for the atomic positions and their respective mass-to-charge-state ratios. RRNG is a range file format identifying the chemical information of each ion species in APT data by associating an element with its mass-to-charge-state ratios. OBJ is a geometry definition file format. 2D and 3D refer to two- and three- dimensions, respectively.

The whole process has been interpreted in this software package. We provide graphical interfaces for all analysis procedures. The details for each step are documented in the manual slides with an synthetic data set as an example. 

Installation:

- Please use the defaut folder for the installation.
- Users must also install Blender 2.92.0. Please also use the default folder.
- Details on the installation process can be found in the manual slides. 

Please use, edit and distribute it freely!  Please do not use it for commercial purposes!

Xuyang Zhou, Ye Wei, Markus Kühbach, Huan Zhao, Florian Vogel, Reza Darvishi Kamachali, Gregory B. Thompson, Dierk Raabe, Baptiste Gault, Revealing in-plane grain boundary composition features through machine learning from atom probe tomography data, Preprint arXiv (2021)

Acknowledgement: Peter Felfer, Anna Ceguerra, Varvara Efremova, Ye Wei, Markus Kühbach, Huan Zhao, Florian Vogel, Reza Darvishi Kamachali, Baptise Gault, Gregory B Thompson, Dierk Raabe.

