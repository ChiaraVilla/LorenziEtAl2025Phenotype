# LorenziEtAl2025Phenotype
Simulating phenotype-structured PDE models of collective cell migration

## Generalities

**Public gitlab repository LorenziEtAl2025Phenotype** <br />
This repository provides Python code to simulate a phenotype-structured PDE (PS-PDE) models of collective cell migration. <br />
It is supporting material for the tutorial-style review article: <br />

Tommaso Lorenzi, Kevin J. Painter, Chiara Villa (2025) <br />
<i>Phenotype structuring in collective cell migration:
a tutorial of mathematical models and methods</i> <br />
Available of ArXiv [arXiv:2410.13629] and HAL [hal-04851615] <br />

**Authors** <br />
Kevin J. Painter (Politecnico di Torino) and Chiara Villa (Centre Inria de Saclay)

**Citation** <br />
Painter, K.J. and Villa, C. (2025). Python code to simulate phenotype-structured PDE models of collective cell migration. <br />
If you use this software in your work then please cite the above named paper.

**Copyright notice** <br />
Python code to simulate phenotype-structured PDE models of collective cell migration. <br />
Copyright (C) 2025 K.J. Painter & C. Villa

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.
This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
You should have received a copy of the GNU General Public License
along with this program.  If not, see https://www.gnu.org/licenses/.

## Repository content and how to use

The code is set up to simulate the PS-PDE (7) in Lorenzi et al. (2025), using the numerical scheme outlined in Section 4.3 of the paper.
This corresponds to simulating diffusion-driven movement of a cell population with phenotype-dependent diffusion coefficient. <br />

**Numerical scheme** <br />
The numerical scheme is based on the method of lines, employing a finite volume approximation in physical and phenotype space. See Section 4.3 in Lorenzi et al. (2025) for more details. <br />

**Files content** <br />
- 'PhenMotion_Diff.py' : this is the main python script to run in order to simulate the model. You can change here: model parameters, grid definition, time-integration method. The script calls the function 'odeRHSeps' in which the ODE system (obtained upon finite volume approximation of the spatial derivatives) to solve is stored, and the function 'plot_TW' to plot the solution over time: make sure that the files 'odeRHSeps_Diff.py' and 'plot_Diff.py' are in the same folder as 'PhenMotion_Diff.py' when running the file. <br />
- 'odeRHSeps_Diff.py' : this file stores the function 'odeRHSeps' (called from the main file 'PhenMotion_Diff.py') in which the ODE system obtained upon finite volume approximation of the spatial derivatives is stored. You can change here: the definition of the diffusion coefficient D(y) and the cell proliferation and death term (defined in "kinetics"), as well as details of the finite volume scheme. <br />
- 'plot_Diff.py' : this file stores the function 'plot_TW' (called from the main file 'PhenMotion_Diff.py') to plot the solution at a given time t. You can change here: the plot.<br />
- 'PhenMotion_Diff_notebook.ipynb' : this file is a Jubiter notebook encoding all parts of the python script found in the files 'PhenMotion_Diff.py', 'odeRHSeps_Diff.py' and 'plot_Diff.py', should you prefer to work with this instead. Just make sure to run the different sections of the notebook in the correct order.

## For students

<i>Are you a student?</i><br />
Try extending the code to simulate the PS-PDE model (17), capturing pressure-based cell motion, and the PS-PDE model (20), capturing taxis-based motion!<br />
Just follow the scheme outlined in Section 4.3 ("Discretisation in phenotype and physical space") and apply first-order upwind for the advection term.<br />

<i>Want a challenge?</i><br />
Try speeding up the code, by experimenting with time-splitting schemes or employing the WKB ansatz, as hinted at at the end of Section 4.3 ("Time integration") of the manuscript.<br />

Buon lavoro!<br />
Tommaso, Kevin & Chiara <br />

<i>You can find our contact information on our websites:</i> <br />
https://staff.polito.it/tommaso.lorenzi/contacts.html <br />
https://sites.google.com/view/kevinjpainter/home <br />
https://chiaravilla.github.io/website/contact.html







