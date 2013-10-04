bedToEnrichments
================

A genomics tool to calculate if a set of regions are enriched for being near a certain set of genes, or any
other genomic locations that can be assigned attributes.

Installation
============

This code uses the Kent libraries from UCSC as well as the GNU Scientific Library, so those must be installed
to compile bedToEnrichments

<ol>
<li> Download and compile the Kent Libraries

<ol>
<li> Check the setting of your machtype variable with:<br />
echo $MACHTYPE<br />
it should be something in this list: i386 i686 sparc alpha x86_64 ppc.  If it is not, set your machtype variable.
<li> Go to a folder on your computer where you want the kent source tree to reside and type:<br />
git clone git://genome-source.cse.ucsc.edu/kent.git<br />
to download the repository onto your own computer.
<li> go to the src/lib directory within the kent source repo that you just cloned:<br />
cd kent/src/lib<br />
<li> Compile the libraries<br />
make
<li> If this was successful, you should have a file here:<br />
kent/src/lib/x86_64/jkweb.a<br />
the x86_64 will be the machtype of your machine.</br />
</ol>

If this was not successful then you should look at the build instructions in the kent repo itself
by looking at this file:<br />
kent/src/product/README.building.source

<li> Download and compile the GNU Scientific Library

<ol>
<li> Get GSL:<br />
wget ftp://ftp.gnu.org/gnu/gsl/gsl-1.16.tar.gz
<li> Configure the make files and you may want to use the prefix option to change the install path:<br />
./configure --prefix=/home/lowec/src/gsl/gsl-1.16_install
<li> Compile the source:<br />
make
<li> Move the library to the install location:<br />
make install
</ol>

<li> Compile bedToEnrichments
<ol>
<li> Edit the bedToEnrichments makefile so that it points to the kent and gsl libraries on your system.  These
are the four lines you will have to modify:<br />
HG_INC += -I/home/lowec/kent/src/hg/inc -I/home/lowec/kent/src/inc<br />
L += /home/lowec/kent/src/lib/${MACHTYPE}/jkweb.a<br />
HG_INC += -I/home/lowec/src/gsl/gsl-1.16_install/include<br />
L += /home/lowec/src/gsl/gsl-1.16_install/lib/libgsl.a /home/lowec/src/gsl/gsl-1.16_install/lib/libgslcblas.a

<li> Compile bedToEnrichments:<br />
make

<li> Running bedToEnrichments with no parameters should give a brief help message.
</ol>
</ol>
