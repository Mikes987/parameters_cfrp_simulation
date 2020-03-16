# Parameters for CFRP Simulation
Python scripts in order to successfully simulate progressive damage of carbon fiber reinforced polymers.

In my master thesis, I simulated temperature dependent progressive damage in Abaqus 2017. The foundation was the Nasa [CompDam_DGD](https://github.com/nasa/CompDam_DGD) code and I expanded the script by adding a module that recalculated mechanical properties based on the operating temperature.

For a successful simulation, two processes had to be accomplished.

1. Transform temperature dependent results of real mechanical tests into a sufficient mathematical equation and implement it into the CompDam Code.

2. Create a CFRP-model in Abaqus that fulfills the requirements of the CompDam Code.


## Temperature dependent mechanical Properties
While mechanical properties in longitudinal direction are assumed to be temperature independent within a certain temperature range under tensile load and results of mechanical testing have led to the assumption of nearly linear correlations between temperature and mechanical properties in longidutinal direction under compressive load as well transverse direction, correlations under 12-shear load were not that obvious.

Mechanical testing in shear direction do not show linear but non linear mechanical behavior, i.e. it is not only necessary to find a mathematical equation for temperaturedependent shear strength and shear modulus but obtain non linear behavior during simulation.
CompDam_DGD handles this nonlinear behavior in quasi static simulation by using the Ramberg Osgood equation. This equation uses two parameter. For further information e.g. visit [wikipedia](https://en.wikipedia.org/wiki/Ramberg%E2%80%93Osgood_relationship).

In order to evaluate the Ramberg Osgood Parameters, two scripts were written in Python.
- In ```ROfit.py``` the txt files, which represent the results of mechanical shear tests, are loaded. The script summarizes them according to their operating temperature, creates a mean curve and transforms it into the Ramberg Osgood equation. The CompDam Code does not use the equation as shown in Wikipedia but an adjusted one. Furthermore, diagrams are drawn, a file that includes all RO parameters according to the temperature and a logfile will be printed.
- In ```RO_T_Fit.py```the file ```rop.txt``` will be loaded and the columns transferred into np.arrays. By checking the datapoints and varying the shape of the y-axis, it was clear that the RO-parameter n can be predicted by using a arctan function whereas a can be predicted by using a exp(arctan) function within that temperature range. These function were essentially implemented into CompDam.

