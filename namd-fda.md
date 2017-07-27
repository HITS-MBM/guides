## Using GROMACS force distribution analysis (FDA) tool with NAMD trajectories

Install the [FDA tool for GROMACS](https://github.com/HITS-MBM/gromacs-fda).

Then install the python tool [PyTopol](https://github.com/resal81/PyTopol/wiki/PyTopol-Installation). After the installation you need to download the [psf2top.py](https://github.com/resal81/PyTopol/tree/master/scripts) script.

Now run the script: `psf2top.py -p psf_file.psf -c parameter_file1.prm parameter_file2.inp` One can add as many parameter files as needed. This script will generate a `top.top` file and maybe some itp files. This are now the topology files in GROMACS format.

The PyTopol tool recommends to check if the topology conversion was succesful. This can be done by using the python tool [MDAnalysis](https://github.com/MDAnalysis/mdanalysis/wiki/TopologyDataStructures) which can read binary files like trr, tpr, psf and dcd.

Next the `trajectory.dcd` file needs to be converted to a `trajectory.trr` file. Open VMD: `vmd psf_file.psf trajectory.dcd`. After VMD was opened select the molecule. Then click on `File` and select `Save Coordinates`. Now choose the trr format and save it.
Now create the file named `mdp_file.mdp` which contains for example the following code:

	; Run parameters
	integrator	= md		; leap-frog integrator
	nsteps          = 500000	; 2 * 500000fs = 1ns
	dt		= 0.002		; 2 fs
	; Output control
	nstxout         = 10000		; save coordinates every 0.02 ns
	nstvout         = 10000		; save velocities every 0.02 ns
	nstfout		= 0		; do not save forces
	nstxtcout	= 1000		; write xtc every 2 ps
	;nstxout-compressed = 1000	; write xtc every 2 ps
	nstenergy	= 1000		; write energies every 2 ps
	nstlog		= 10000		; update log file every 0.02 ns
	; Constraints
	constraint_algorithm = lincs	; holonomic constraints
	constraints	= hbonds	; hbonds (heavy atom-H bonds) constrained
	lincs_iter	= 1		; accuracy of LINCS
	lincs_order	= 4		; also related to accuracy
	continuation	= yes		; continuing after NPT
	; Neighborsearching
	cutoff-scheme   = Verlet        ; for avx and cuda
	verlet-buffer-tolerance = 0.005	; 5.0 syntax
	ns_type         = grid          ; search neighboring grid cells
	nstlist         = 20            ; steps
	rlist		= 1.0		; short-range neighborlist cutoff (in nm)
	rcoulomb	= 1.0		; short-range electrostatic cutoff (in nm)
	rvdw		= 1.0		; short-range van der Waals cutoff (in nm)
	; Electrostatics
	coulombtype	= PME		; Particle Mesh Ewald for long-range electrostatics
	pme_order	= 4		; cubic interpolation
	fourierspacing	= 0.16		; grid spacing for FFT
	; Temperature coupling is on
	tcoupl		= V-rescale	; modified Berendsen thermostat
	tc-grps		= System	; two coupling groups - more accurate
	tau_t		= 0.1		; time constant, in ps
	ref_t		= 300		; reference temperature, one for each group, in K
	; Pressure coupling is on
	pcoupl		= Parrinello-Rahman	; Pressure coupling on in NPT
	pcoupltype	= isotropic	; uniform scaling of box vectors
	tau_p		= 2.0		; time constant, in ps
	ref_p		= 1.0		; reference pressure, in bar
	compressibility = 4.5e-5	; isothermal compressibility of water, bar^-1
	; Periodic boundary conditions
	pbc		= xyz		; 3-D PBC
	; Dispersion correction
	DispCorr	= EnerPres	; account for cut-off vdW scheme
	; Velocity generation
	gen_vel		= no		; Velocity generation is off
	; COM PULLING
	; Pull type: no, umbrella, constraint or constant_force
	pull		= no

The mdp file will not be used but is needed to creat the tpr file.

Having created the mdp file run: `gmx editconf -f pdb_file.pdb -o gro_file.gro -d 1` which creates the gro file from the pdb file. The -d 1 adds a box around the System with a distance of 1nm between the wall of the box and the System. Be careful the box might need to be bigger if the system rotates during the simulation.
Now one can create the tpr file: `gmx grompp -f mdp_file.mdp -c gro_file.gro -p top.top -o tpr_file.tpr`
Now use `gmx make_ndx -f tpr_file.tpr` to create the `index.ndx` file. In this file one can specify different atom groups in the system, for example by atom type.
After this one has to create an `input.pfi` which looks for example like the following:

	onepair = summed
	group1 =  System
	group2 =  System
	atombased = pairwise_forces_scalar
	residuebased = no
	type =  all

`group1` and `group2` have to be specified in the index file and are the two groups between which the forces are calculated. In the index file every group starts with its name in brackets [ name ] followed by a list of atom ids which are in this group. A default group is the System group which contains all atoms of the whole system. Force calculations can be either residual based or atom based. type specifies which types of interactions should be considered for the force distribution analysis. The output of the FDA can be either vectors, signed scalars or unsigned scalars. For more information regarding the index and input file check the [FDA manual](https://github.com/HITS-MBM/gromacs-fda/blob/master-fda/fda-manual/fda-manual.pdf).

Now one has all needed files and can start the actual FDA by running: `gmx_fda mdrun -nt 1 -rerun trr_file.trr -pfi input.pfi -pfn index.ndx -s tpr_file.tpr -pfa output.pfa` For more information check the [FDA manual](https://github.com/HITS-MBM/gromacs-fda/blob/master-fda/fda-manual/fda-manual.pdf).
