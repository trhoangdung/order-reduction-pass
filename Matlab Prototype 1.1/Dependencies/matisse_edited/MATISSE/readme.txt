% MATISSE
% Version 1.0 July-2005
% Metrics for Approximate TransItion Systems Simulation and Equivalence.
%
% Author: Antoine Girard
% Department of Electrical and Systems Engineering
% University of Pennsylvania
% agirard@seas.upenn.edu



MATISSE Toolbox installation notes

Section 1
=========

In order to use MATISSE, set a Matlab path to the "MATISSE/" directory and to all it's subdirectories. If you are using Matlab for Windows, go to the "File - Set Path..." menu, choose "Add with Subfolders..." and pick up the MATISSE directory. Click on the "Save" button to store the updated path setting. Under Unix, you can either manually edit the file "startup.m", or to use the same procedure described above.

To explore functionality of MATISSE, try one of the following:

help matisse
matisse_demo1
matisse_demo2
matisse_demo3

If you have any questions or comments, or you observe buggy behaviour of the
toolbox, please send your reports to

  agirard@seas.upenn.edu

Section 2 (Additional software requirements)
============================================

The MATISSE toolbox uses the following packages:

- The Multi Parametric Toolbox (MPT) by M. Kvasnica and P. Grieder and M. Baotic,
  http://control.ee.ethz.ch/~mpt/

- The YALMIP interface by Johan Lofberg
  included in the MPT, also at http://control.ee.ethz.ch/~joloef/

- SEDUMI: the Matlab toolbox for solving optimization problems over symmetric cones by Jos F. Sturm
  included in the MPT, also at http://fewcal.kub.nl/sturm/software/sedumi.html

- The MATLAB optimization toolbox by MathWorks
