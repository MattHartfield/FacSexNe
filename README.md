README FOR FAC SEX NE

Supplemental material for manuscript "Approximating the coalescent under facultative sex". Contains both code used to create simulation data, and Mathematica notebook of derivations.

*** C Code ***

The C program is used to simulate spread of neutral alleles under facultative sex. Code based on that used in Agrawal and Hartfield 2016 (https://github.com/MattHartfield/BalSelSims).

Simulation uses routines found with the GNU Scientific Library (GSL) (http://www.gnu.org/software/gsl/). Since GSL is distributed under the GNU General Public License (http://www.gnu.org/copyleft/gpl.html), you must download it separately from this file.

This program can be compiled in e.g. GCC using a command like:

gcc FacSexNe.c -I/usr/local/include -L/usr/local/lib -lm -lgsl -lgslcblas -o FacSexNe

(Replace the -I and -L with locations of the GSL install on your own computer.)

Then run by executing:

./FacSexNe N sex gc reps

Where:

- N is the population size
- sex is the frequency of sex (a value between 0 = obligate asex, and 1 = obligate sex)
- gc is the frequency of mitotic gene conversion
- reps is how many times to introduce the linked neutral allele

*** Mathematica notebook ***

"File S1.zip" is a Mathematica notebook (and PDF printout) outlining the mathematical derivations used in the manuscript. The Wolfram Player (https://www.wolfram.com/player/) can be used to view the file if you do not have the Mathematica software.

Comments to m.hartfield@ed.ac.uk.
