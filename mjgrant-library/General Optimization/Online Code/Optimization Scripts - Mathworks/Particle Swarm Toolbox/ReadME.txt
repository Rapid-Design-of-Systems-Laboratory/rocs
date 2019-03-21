PSOt, particle swarm optimization toolbox for matlab. 

May be distributed freely as long as none of the files are modified. 

Send suggestions to bkbirge@unity.ncsu.edu. 

Updates will be posted periodically at www4.ncsu.edu/~bkbirge

To install:
Extract into any directory you want but make sure the matlab path points to that directory and the subdirectories 'hiddenutils' and 'testfunctions'. Enjoy! - Brian Birge


PSO toolbox implementing Common, Clerc 1", and Trelea types along with an alpha version of tracking changing environments. Can search for min, max, or 'distance' of user developed cost function. Very easy to use and hack with reasonably good documentation and will take advantage of vectorized cost functions. It uses similar syntax to Matlab's optimization toolbox. Includes a suite of static and dynamic test functions. Run 'DemoPSOBehavior' to get a quick example of maximizing the static F6 Shaffer function. It is in constant development and I welcome suggestions.

Usage ideas: to find a global min/max, to optimize training of neural nets, error topology change tracking, teaching PSO, investigate Emergence, tune control systems/filters, paradigm for multi-agent interaction, etc.

The download also includes a pdf of a powerpoint presentation that gives an overview of PSO.

Files included:

in main directory:

0) ReadMe.txt - this file, duh
1) A Particle Swarm Optimization (PSO) Primer.pdf  -  powerpoint converted to pdf presentation explaining the very basics of PSO
2) DemoPSOBehavior.m - demo script, useful to see how the pso main function is called
3) goplotpso4demo.m - plotting routine called by the demo script, useful to see how custom plotting can be developed though this routine slows down the PSO a lot
4) goplotpso.m - default plotting routine used by pso algorithm
5) pso_Trelea_vectorized.m - main PSO algorithm function, implements Common, Trelea 1&2, Clerc 1", and an alpha version of tracking environmental changes.
6) pso_Trelea.m - PSO algo that doesn't require cost function to have a single input, only implements Common and Trelea types, will be phased out in the future but useful if you have some of your own functions already written with several inputs
in 'hiddenutils'

1) forcerow, forcecol.m - utils to force a vector to be a row or column, superseded by Matlab 7 functions I believe but I think they are still called in the main algo
2) normalize.m - takes a matrix and reformats the data to fit between a new range, very flexible

in 'testfunctions'

1) 20 different testfunction scripts, see help for each one for specifics