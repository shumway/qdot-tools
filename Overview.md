# Introduction #

I still need to post documentation. Below is a demo I gave for a simulation lab class at Harvard a couple years ago.
This write up should be updated with python scripts.


# InGaAs/GaAs Self-Assembled Dot Demo #

The steps to modeling a self-assembled dot are illustrated in this figure.

<img src='http://physics2.asu.edu/shumway/codes/demo/InGaAsSAD/method.png' />

<ol>
<li>Position atoms to define the quantum dot. Based on experimental knowledge, design a model for your stucture. This includes the shape, size, and composition of the dot and any quantum wells or wetting layers. Also specify the size of the periodic simulation cell. We define the structure in an XML file, struct.xml, and use a utitily program, geom, to create the atom positions and other model information in struct.h5.</li>
<li>Relax the atomic positions. Because different chemical species have different bond lengths, the heterostructure is strained. Use a ball-and-spring model, such as VFF or Stillinger-Weber, to model the strain energy. The utility relax finds the equilibrium atom positions that minimize the strain energy, and stores the relaxed positions and residual stress fields in the file relaxed.h5.</li>
<li>Calculate the strained band offsets. The path integral method uses effective mass models. The utility getEMA creates a file emagrids.h5 that have the strain-modified confinement potentials and other information about the nanostructure.</li>
<li>Simulate electrons and holes with path integral Monte Carlo. The program paths is used to run path integral simulations of electrons and holes in the dots.</li>
</ol>

In the following sections we give detailed instructions for how to perform each of these four steps to simulate an InGaAs/GaAs dot.


## 1. Position atoms to define the quantum dot ##
The file struct.xml defines the geometry of the nanostructure. Lengths are in units of the GaAs lattice constant, aGaAs = 10.683 a0 = 0.565 nm.
```
<Nanostructure>
  <!-- This is a 200 A (160) x 35 A cone, as in PRB v. 64, p. 125302 (2001) -->
  <!-- For speed of demo, we make supercell way too small (should be 80x80x80). --> 
  <Materials>
    <!--A pure GaAs binary and two random InGaAs alloys -->
    <Binary name="GaAs" anion="As" cation="Ga"/>
    <CommonAnionTernary name="InGaAs" anion="As" cation1="In" cation2="Ga" x1="0.5"/>
    <CommonAnionTernary name="In03GaAs" anion="As" cation1="In" cation2="Ga" x1="0.3"/>
  </Materials>
  <SuperCell nx="48" ny="48" nz="48" a="10.683" material="GaAs"/>
  <Structures>
    <!-- 200 A (160) x 35 A cone -->
    <Cone x="24" y="24" z="18.48" height="6.0" base="35.3" top="28.3"
       material="InGaAs"/>
    <!-- A 6 monolayer (3 a) quantum well as the wetting layer -->
    <Well z="16.98" thickness="3.0" material="In03GaAs"/>
  </Structures>
</Nanostructure>
```
To construct the atomic structure, run geom. This will produce an HDF5 file, struct.h5 (10.9 MB), that has all the atomic positions and other information. Use the command h5dump --headers struct.h5 to see a list of all the information contained in this compressed binary file. This structure has 884,736 atoms in a 27.1 x 27.1 x 27.1 nm supercell.

Here's an image of the unrelaxed atoms in the struct.h5 file (slice taken through center, 16 pixels/nm). The dark bands on the left and right side are the edge of the periodic supercell, which should really be about twice as large as the small cell used for this demo.

<img src='http://physics2.asu.edu/shumway/codes/demo/InGaAsSAD/unrelaxed.png' />

(This image was generated with struct2png.cc, available as executable on nnin-cluster in /home/jshumway/bin/struct2png.)

## 2. Relax atomic positions ##
The atoms in struct.h5 are located on the ideal GaAs zinc-blende lattice. The larger indium atoms in the dot should distort the lattice to minimize strain energy.
The program relax moves the atoms to reduce strain energy. It needs an input file, relax.xml, that describes the relaxation algorithm and system parameters.
```
<AtomisticStructure>
  <!-- Input file to relax an InGaAs heterostructure.-->
  <Structure infile="struct.h5"/>
<Relax maxiter="300" ftol="1e-6" outfile="relaxed.h5">
  <Checkpoint interval="20" file="checkpoint.h5"/>
  <CellRelax interval="10" mode="zonly"/>
  </Relax>
  <CalculateStresses outfile="relaxed.h5"/>
  <ForceModel>
    <VFF/> 
  </ForceModel>
</AtomisticStructure>
```
The relaxation routine will run for 300 line minimization or until the maximum force is less than 1x10-6 (atomic units). This may take several hours. At the conclusion you will have an output file relaxed.h5. This HDF5 file has the same format as struct.h5, but with strain-relaxed atomic coordinates. If the structure does not fully relax, you can copy the file relaxed.h5 onto struct.h5 and rerun.

Note that we allow the supercell to relax in the z-direction. This is necessary because the InGaAs wetting layer is thicker than the GaAs. In the limit of low dot density (large supercell), only the wetting layer (not the dot layer) should stretch the supercell in the z-direction. For this demo, the supercell z-dimension expands 0.55%, or 0.15 nm, from 27.13 to 27.28 nm.

Here's an image of the relaxed structure. The atom positions only shift a few percent, so structural differences from the unrelaxed structure are difficult to see. Mostly the strain fields and band-offsets change, which are crucial to electron and hole properties.

<img src='http://physics2.asu.edu/shumway/codes/demo/InGaAsSAD/relaxed.png' />

(This image was generated with struct2png.cc, available as executable on nnin-cluster in /home/jshumway/bin/struct2png.)

## 3. Calculate strained band offsets ##
Next, use getEMA to calculate the strained band offsets. This program does not require any input parameters at this time because it is hard-coded for InGaAs. The program reads the atomic positions and strain fields from relaxed.h5 and produces the output file emagrids.h5 (8.9 MB), which contains a set of rectangular grids for the composition, strain, and strain-modified band offsets. To see a list of the contents of the emagrids.h5 file, use the command h5dump --headers emagrids.h5
Use can use Matlab to view the grids. Transfer emagrids.h5 to your computer's desktop, then start up Matlab.

```
>> ve = 27.211*double(hdf5read('emagrids.h5','boffset/ve'));
>> [C,h]=contourf(ve(:,:,24)');
>> colorbar(C);
```
The factor 27.211 (2\*13.6 eV) converts to electron volts. The zero of energy is the GaAs valence band edge, and the GaAs conduction band edge is at 15.2 eV. You should find an image like this:

<img src='http://physics2.asu.edu/shumway/codes/demo/InGaAsSAD/ve.png' />

You can repeat the calculation for the valence band, which is stored in vh:

<img src='http://physics2.asu.edu/shumway/codes/demo/InGaAsSAD/vh.png' />

Other datasets in emagrids.h5 you can plot include:
  * strain/trace
  * strain/biaxial
  * composition/Ga
  * composition/In
  * composition/As (note, this is uniform in InGaAs)

## 4. Simulate exciton with path integral Monte Carlo ##
The program
<a href='http://code.google.com/p/pi-qmc'><b><em>pi</em></b></a> is used to sample thermal properties of electrons and holes using path integrals. To run the demos, download the XML file and run with pi filename.xml.

(Use "View Source" if you have viewing or copying raw XML in your browser.)

You can look at the energy under "thermo\_energy" in the log file. Output is in pimc.dat and pimc.h5.
<ol>
<li>ONE ELECTRON<br />
Use this input file: 1e.xml</li>
<li>ONE HOLE<br />
Use this input file: 1h.xml</li>
<li>EXCITON WITHOUT COULOMB INTERACTION<br />
Use this input file: 1h.xml</li>
<li>EXCITON WITH COULOMB INTERACTION<br />
Use this input file: 1ex.xml</li>
<li>BIEXCITON<br />
Use this input file: 2ex.xml</li>
</ol>


Work supported by NSF Grant DMR 02-39819.