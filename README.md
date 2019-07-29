# polybuild

### Dispersed polymer configurations for Molecular Dynamics

Polybuild, corruption of 'polymer builder', is a software developed in Shell Script 
language, interpreted by GNU bash (version 4.3.30), which facilitates building initial 
configurations for molecular dynamics simulations of dispersed polymers.


### What does polybuild do?

From a concise specification of the microstructure of a linear architecture dispersed 
polymeric system, the quantities of different molecules in the simulated sample and the 
parameters referring to the simulation box, the polymer builder (polybuild) creates the 
configuration files, called xyz files, containing the positions of all particles for all 
molecular types in the sample. After that, it creates a file in a format recognized by 
the software Playmol[1] that packs the molecules in a simulation box with all the required 
description.


### Why polybuild is useful?

Particle-based simulations require an initial configuration, which consists in the 
initial positions of all particles of the simulated system. The task of creating such 
initial conditions may not be trivial. It consists of an optimization problem with several 
constraints, such as box size and shape restrictions, reproduction of the system density, 
and the most difficult one, which is that short-range repulsions are not excessively large, 
thus resulting in numerical instability during the course of the simulation.[2,3]
It is a common practice to use computer programs that solve this optimization problem by 
packing the particles in simulation boxes subject to the constraints that force fields 
impose.[4] 
Generally, these programs create an initial configuration by replicating atoms and molecules 
that make up the system in a defined region of space, respecting a series of spatial constraints, 
thus creating complex initial configurations. 
However, the number of different species in a dispersed polymer system can be in the order of 
thousands. This requires, therefore, the description of the topology of these thousands of 
chains. Even if the simulation reproduces only a small portion of space or if it involves a 
coarse-grained model, obtaining such information and feeding it to computer programs that 
generate initial configurations can be costly. For example, some simulations performed in 
the present study have about 300 oligomeric chains that differ in size or composition, and 
the definitions of the particles, chemical bonds and topology necessary to generate initial 
configurations occupy, for each system, about of 13000 lines of code, even without considering 
the information about angles and dihedrals.



### Criteria

The concise writing of microstructures of all chains demands a code that can represent single 
molecules or polymer chains with any architecture. Here a coding system is proposed for linear 
chains, since this is the scope of the simulations of the present article. Each block is 
represented by two pieces of information: the first one specifies the type of particle that 
constitutes the block and the second represents the block length, which must necessarily be a 
cardinal number. For example, A5B5 represents a diblock-type oligomer, consisting of a block with 
five beads of type A attached to a block with five beads of type B.
This coding is a representation of the structural formula of the molecule, so that it encloses all 
compositional information and the connectivity between the component particles of the molecule. 
However, it does not carry any conformational information, which is also necessary for the 
specification of the molecule to be replicated. Therefore, in addition to the sequence that 
translates the structural information, it is necessary to specify conformational parameters for 
writing the xyz files. The simplest choice is adopted, with the particles inserted to form a straight 
line in a given cartesian direction. Thus, each particle is inserted at a distance lb from the previous 
one, forming a connecting angle of 180º.

### Requirements
- Playmol - see reference [1].
- a bash interpreter (the version provided here is tested for a GNU bash version 4.3.30).

### How to use polybuild? 
The use is quite straightforward. The *'polybuild.sh'* script receives an input file that contains 
the microstructure specifications, system parameters, and Playmol parameters to be used for packaging - 
see the template file called *'input_file_template'*.

The command that processes the script has the following syntax and is capable of receiving multiple input files at once:

```
> ./polybuild.sh input_file1 [input_file2 … input_fileN]
```

### If this script is used in the preparation of scientific publications, the following paper should be cited.

LEMOS, T. S. M., ABREU, C. R. A., PINTO, J. C. 
Mesoscopic Simulation of Dispersed Copolymers - Effects of Chain Length, Chemical
    Composition, and Block Length Distributions on Self-Assembly, 
*Macromolecular Theory and Simulations*, **2019**

### References


&nbsp;
[1]	Playmol, http://atoms.peq.coppe.ufrj.br/playmol.


&nbsp;
[2]	D. Frenkel, B. Smit, *Understanding Molecular Simulation*, Academic Press, London, **2001**. 


&nbsp;
[3]	M. P. Allen, D. J. Tildesley, *Computer Simulation of Liquids*, Oxford University Press, New York, **1991**.


&nbsp;
[4]	L. Martínez, R. Andrade, E. G. Birgin, J. M. Martínez, *Journal of Computational Chemistry* **2009**, 30, 2157-2164.
