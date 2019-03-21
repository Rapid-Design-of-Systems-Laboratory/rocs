                       DISTRIBUTION DISKETTE

       FOR THE SHUFFLED COMPLEX EVOLUTION (SCE-UA) METHOD


                           Version 2.2





                          Qingyun Duan - 1
                        Hoshin V. Gupta - 2
                       Soroosh Sorooshian - 3
                       
 1 - Lawrence Livermore National Laboratory, 7000 East Avenue, 
     Livermore, CA 94568, ph. 925-422-7704, email: qduan@llnl.gov
 
 2 - Dept of HWR, University of Arizona, Tucson, AZ 85721
     Ph: (520) 626-6974, email: hoshin.gupta@hwr.arizona.edu
 
 3 - Dept. of Civil & Env. Engr., University of California,
     Irvine, CA 92697, ph. 949-824-8825, email: soroosh@uci.edu



                           May, 1994


Table of Contents:
==================

I.   Statement and Disclaimer by the Authors

II.  Introduction

III. Files Contained in This Diskette

     a.  Description Files

     b.  FORTRAN Files

     c.  Input Data Files

     d.  Sample Output Files

IV.  How to Use This Diskette

     a.  How to Use SCE-UA Method on the Standard Test Problems

     b.  How to Use SCE-UA Method on the SIXPAR Model Problem

     c.  How to Use SCE-UA Method on Your Own Problem

V.   How to Prepare Input Data File SCEIN.DAT

VI.  How to Prepare Input Data File SXPINP.DAT

VII. More Information and Questions

I.   Statement and Disclaimer by the Authors:
=============================================

The users of the programs contained in this diskette can copy and
use these programs freely, without seeking authors' permission.
The authors request all users make appropriate references to the
use of these programs.  The authors disclaim any responsibility
resulting from use of these programs.

II.   Introduction:
===================

The SCE-UA method is a general purpose global optimization
program.  It was originally developed by Dr. Qingyun Duan as part
of his doctoral dissertation work at the Department of Hydrology
and Water Resources, University of Arizona, Tucson, AZ 85721, USA. 
The dissertation is entitled "A Global Optimization Strategy for
Efficient and Effective Calibration of Hydrologic Models".  The
program has since been modified to make it easier for use on
problems of users' interests.  The algorithm has been described
in detail in an article entitled "Effective and Efficient Global
Optimization for Conceptual Rainfall-Runoff Models", Water Resources
Research, Vol 28(4), pp.1015-1031, 1992; and in an article entitled
"A Shuffled Complex Evolution Approach for Effective and Efficient
Global Minimization" by Q. Duan, V.K. Gupta and S. Sorooshian,
Journal of Optimization Theory and its Applications, Vol 76(3), 
pp 501-521, 1993.  A paper entitled "Optimal Use of the SCE-UA Global
Optimization Method for Calibrating Watershed Models", by Q. Duan,
S. Sorooshian, & V.K. Gupta, Journal of Hydrology, Vol.158, 265-284,
1994, discussed how to use the SCE-UA Method in an efficient and 
effective manner.  The SCE-UA algorithm is also briefly described
in the accompanying ASCii file 'CONCEPT_SCE.txt' contained in this
package (Figures referred in this file are in separate ".gif" files).


III.  Files Contained in This Diskette:
=======================================

This diskette contains the FORTRAN codes of the main programs,
the input subroutines, the subroutine implementing the Shuffled
Complex Evolution (SCE-UA) global optimization algorithm, and
some standard test functions from the optimization literature.
The SIXPAR model calibration problem is also included as a test
problem.  The SIXPAR model is a simplified example of a hydrologic
rainfall-runoff model.  More complex versions of such a model have
been used operationally in the United States and in many other
countries for hydrologic forecasting purposes.  The problems
presented here are the ones tested in the above mentioned papers. 
This diskette also contains two sample input data files and two
sample output files.

Specifically, this diskette should contain the following files:


a.   Description Files:
-----------------------

CONCEPT.SCE - ASCii file briefly describing the SCE-UA algorithm
README      - This file


b.   FORTRAN Files:
-------------------

FUNCTN1.FOR - Subroutine for computing function value for Goldstein-
              Price function and for checking constraints
FUNCTN2.FOR - Subroutine for computing function value for Rosenbrock
              function and for checking constraints
FUNCTN3.FOR - Subroutine for computing function value for Six-Hump
              Camel-Back function and for checking constraints
FUNCTN4.FOR - Subroutine for computing function value for Rastrigin
              function and for checking constraints
FUNCTN5.FOR - Subroutine for computing function value for Griewank
              function and for checking constraints
FUNCTN6.FOR - Subroutine for computing function value for Shekel
              function and for checking constraints
FUNCTN7.FOR - Subroutine for computing function value for Hartman
              function and for checking constraints
SCEIN.FOR   - Subroutine for reading the input data needed by
              subroutine SCEUA.FOR
SCEMAIN.FOR - Main program calling subroutines SCEIN and SCEUA
SCEUA.FOR   - Main subroutine implementing the SCE-UA algorithm
SXPINP.FOR  - Subroutine for reading the input data needed to run
              SIXPAR model
SXPMAIN.FOR - Main program to run SCE-UA method on the SIXPAR model


c.   Input Data Files:
----------------------

SAMPLE1.IN  - Input data file for the SCE-UA method for Goldstein-Price
              that generates output file SAMPLE1.OUT
SAMPLE2.IN  - Input data file for the SCE-UA method for the SIXPAR model
              that generates output file SAMPLE2.OUT
SCEIN.DAT   - Input data file for the SCE-UA method
SXPINP.DAT  - Input data file for the SIXPAR model


d.   Sample Output Files:
-------------------------
SAMPLE1.OUT  - Sample output file for the Goldstein-Price function
SAMPLE2.OUT  - Sample output file for the SIXPAR model problem

IV. How to Use This Diskette:
=============================

a.   How to Use SCE-UA Method on the Standard Test Problems:
------------------------------------------------------------

For illustration purpose, the Goldstein-Price function is used to
show you how to use this diskette.  The procedure is outlined as
follows:

i)   Compile and link FORTRAN programs SCEMAIN.FOR, SCEIN.FOR,
     SCEUA.FOR, and FUNCTN1.FOR.  You should get an executable
     file named SCEMAIN.EXE (or by other name given by you).

ii)  Prepare input data file SCEIN.DAT (Section V explains how
     to prepare this input data file).

iii) Run the executable file SCEMAIN.EXE you just created.  When
     the run is completed, a file named as SCEOUT.DAT should be
     created.  If input data file SCEIN.DAT you used has the same
     content as file SAMPLE1.IN, SCEOUT.DAT should be the same
     as SAMPLE1.OUT.


b.   How to Use SCE-UA Method on the SIXPAR Model Problem:
----------------------------------------------------------

The procedure to use SCE-UA method on the SIXPAR model calibration
problem is outlined as follows:

i)   Compile and link FORTRAN programs SXPMAIN.FOR, SXPINP.FOR,
     SCEIN.FOR, SCEUA.FOR, and FUNCTN8.FOR.  You should get an
     executable file named SXPMAIN.EXE (or by other name given by
     you).

ii)  Prepare input data file SCEIN.DAT.

iii) Prepare input data file SXPINP.DAT (Section VI explains how
     to prepare the input data file SXPINP.DAT).

iv)  Run the executable file SXPMAIN.EXE you just created.  When
     the run is completed, a file named as SCEOUT.DAT should be
     created.  If input data file SCEIN.DAT you used has the same
     content as file SAMPLE2.IN, SCEOUT.DAT should be the same as
     SAMPLE2.OUT.


c.   How to Use SCE-UA Method on Your Own Problem:
--------------------------------------------------

To use the SCE-UA method on your problem, you have to write one or
two subroutines of your own, depending on the nature of your
problem.  If you have a pure mathematical function optimization
problem, you will have to write a subroutine that has the same
structure as in files FUNCTN1.FOR through FUNCTN7.FOR.  If you
face a model calibration problem that is similar to the SIXPAR
model, then you have to write a subroutine that has the same
structure as FUNCTN8.FOR which computes the appropriate objective
function.  You also have to write a separate input subroutine
similar to SXPINP.FOR to read the input data information specific
to your model.  Once you complete your own subroutines, follow the
appropriate procedures as outlined above to run SCE-UA method.


V.  Input Summary for SCEIN.DAT:
================================

Card  Format  Cols                    Contents
----  ------  ----                    --------

1     i5      1-5   MAXN
                    Maximum number of trials allowed before
                    optimization is terminated.  The purpose of
                    MAXN is to stop an optimization search before
                    too much computer time is expended.  MAXN
                    should be set large enough so that optimization
                    is generally completed before MAXN trials are
                    performed.
                    Recommended value is 10,000 (increase or
                    decrease as necessary).

      i5      6-10  KSTOP
                    Number of shuffling loops in which the 
                    criterion must improve by the specified
                    percentage or else optimization will be
                    terminated.
                    Recommended value: 5.

      f5.2    11-15 PCENTO
                    Percentage by which the criterion value must
                    change in the specified number of shuffling 
                    loops or else optimization is terminated
                    (Use decimal equivalent: Percentage/100).
                    Recommended value: 0.01.

      i5      16-20 NGS
                    Number of complexes used for optimization
                    search.  Minimum value is 1.
                    Recommended value is between 2 and 20 depending
                    on the number of parameters to be optimized and
                    on the degree of difficulty of the problem.

      i5      21-25 ISEED
                    Random seed used in optimization search.  Enter
                    any integer number.  Default value (=1969) is
                    assumed if this field is left blank.
                    Recommended value: any large integer.

      i5      26-30 IDEFLT
                    Flag for setting the control variables of the
                    SCE-UA algorithm.  Enter '0' or leave the field
                    blank for default setting.  Enter '1' for user
                    specified setting.
                    If option '1' is chosen, type the following
                    input card.  Else, leave the next card blank.

2     i5      1-5   NPG
                    Number of points in each complex.  NPG should
                    be greater than or equal to 2.  The default
                    value is equal to (2 * number of optimized
                    parameters + 1).

      i5      6-10  NPS
                    Number of points in each sub-complex.  NPS
                    should be greater than or equal to 2 and less
                    than NPG.  The default value is equal to 
                    (number of optimized parameters + 1).

      i5      11-15 NSPL
                    Number of evolution steps taken by each complex
                    before next shuffling.  Default value is equal
                    to NPG.

      i5      16-20 MINGS
                    Minimum number of complexes required for
                    optimization search, if the number of complexes
                    is allowed to reduce as the optimization search
                    proceeds.  The default value is equal to NGS.

      i5      21-25 INIFLG
                    Flag on whether to include the initial point in
                    the starting population (see Card 3, Field 1
                    below).  Enter '1' if the initial point is to be
                    included.  The default value is equal to '0'.

      i5      26-30 IPRINT
                    Print-out control flag.  Enter '0' to print out
                    the best estimate of the global optimum at the
                    end of each shuffling loop.  Enter '1' to print
                    out every point in the entire sample population
                    at the end of each shuffling loop.  The default
                    value is equal to 0.

3*    f10.3   1-10  Initial estimate of the parameter to be 
                    optimized (used if INIFLG = 1 in Card 2 above).

      f10.3   11-20 Lower bound of the parameter to be optimized.

      f10.3   21-30 Upper bound of the parameter to be optimized.

 *    Repeat Card 3 as many times as necessary until all parameters
      to be optimized are defined.

      For all test problems, the lower and the upper bounds are given
      at the beginning of the function subroutines, along with the
      location of the global optima.

VI.   Input Summary for SXPINP.DAT:
===================================

Card  Format  Cols                   Contents
----  ------  ----                   --------

1     a80     1-80  Input heading (It can be left blank).

2     i5      1-5   NPAR
                    Total number of parameters in the SIXPAR model
                    (= 6).

      i5      6-10  NDATA
                    Number of data points used for calibration.

      i5      11-15 NS
                    Number of state variables in the SIXPAR model
                    (= 2).

      i5      16-20 IOBJ
                    Flag on which objective function to be used.
                    Enter '1' to choose Simple Least Square (SLS).
                    Enter '0' to choose Heteroscedastic Maximum
                    Likelihood Estimator (HMLE).

      i5      21-25 IDATA
                    Flag on what data set to be used.
                    Enter '0' to use the existing hydrologic time
                    series data given in this input file.  Enter
                    '1' to generate hydrologic time series using
                    "true" parameter set and precipitation time 
                    series data given in this input file.

3     a80     1-80  Input heading

4     f10.4   1-10  UM
                    "True"* parameter value for upper zone capacity.

      f10.4   11-20 UK
                    "True" parameter value for upper zone depletion
                    rate.

      f10.4   21-30 BM
                    "True" parameter value for lower zone capacity.

      f10.4   31-40 BK
                    "True" parameter value for lower zone depletion
                    rate.

      f10.4   41-50 A
                    "True" parameter value for the reparameterized
                    percolation equation

      f10.4   51-60 X
                    "True" parameter value for the exponent in
                    the percolation equation.

*     "True" parameter values are used for generating synthetic
      streamflow data (if IDATA = 1) and are used for defining values
      for parameters that are not optimized (see Card 6 below).

5     a80     1-80  Input heading

6     i5      1-5   Index* of the first optimized parameter.

      i5      6-10  Index of the second optimized parameter.

      i5      11-15 Index of the third optimized parameter.

      i5      16-20 Index of the fourth optimized parameter.

      i5      21-25 Index of the fifth optimized parameter.

      i5      26-30 Index of the sixth optimized parameter.

*     Each of the 6 parameters in the SIXPAR model is assigned to
      a unique index, e.g., the index for UM is equal to 1, UK=2,
      BM=3, BK=4, A=5, and X=6.  These indexes allow users to 
      optimize on a subset of model parameters.

      For example, if we decide to optimize two parameters BM and
      BK, enter '3' in Field 1 of Card 6; enter '4' in Field 2;
      leave the rest of the fields blank.  Properly prepare the
      input file SCEIN.DAT.  The optimization search defined as
      such will try to search the optimal values for BM and BK
      while other parameters take on values specified in Card 4
      of this input file.

7     a80     1-80  Input heading

8     f10.4   1-10  Initial upper zone content.

      f10.4   11-20 Initial lower zone content.

9     a80     1-80  Input heading

10    i5      1-5   Data point number

      f10.4   6-15  Precipitation input value.

      f10.4   16-25 Observed streamflow value*.

*     Observed streamflow values are used only when IDATA is
      equal to '1' (see Field 5 of Card 1).

      Repeat Card 10 as many times as necessary until all data
      points are defined.

VII.   More Information and Questions:
======================================

For more information and questions regarding the SCE-UA method
and the programs contained in this diskette, please contact the
authors through the address, telephone, or email listed on the
front page.
