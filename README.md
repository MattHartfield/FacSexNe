README FOR FAC SEX NE

Scripts used to simulate spread of neutral alleles under facultative sex, as used in the study "Approximating the coalescent under facultative sex". Code based on that used in Agrawal and Hartfield 2016 (https://github.com/MattHartfield/BalSelSims).

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

Comments to m.hartfield@ed.ac.uk.
