# proteins-rigidity-flexibility
Calculates protein rigidity/ flexibility from protein crystal structures
Rigidity is represented by packing factor of individual non-H atoms of a protein.
Flexibility is represented by atomic fluctuations obtained from the B-factors of a crystal structure
Error correction in atomic fluctuations based on a new method has been developed that considers 
contributions of errors from different factors in a consolidated manner.
Temperature coorections have been done and validated.
Linear correlation coefficient between two data series has been calculated
Residuewise packing factors and fluctuations have been calculated.

Individual README files for each code file is given below.

______________________________
Filename: Pack.f
Description:
Codes for computing the atom-wise packing fraction following our definition and also for computing the atomic fluctuations using the B-factors from the respective pdb of the crystal structure of the protein.
-----------------------------------

Input data file contents
------------------------
file0 (file0 = param.data the parameter file)
nset, rsph (nset=number of datasets, rsph=radius of sampling sphere).
npro, nall (npro=protein atoms in crystal, nall=all atoms in protein.
-------------------------------------------------------------

The name of file0 used is param.data
The content of the file ‘param.data’ is
9 
C 1.90
N 1.50
O 1.40
S 1.85
H 0.50
E 0.50 
G 0.50
D 0.50
P 1.90


Example of input file content
-----------------------------
param.data    !file0
3 6.0         !nset, rsph 
1bqc.pdb      !file1 = The PDB file of the crystal structure
2336 2336     !npro, nall 
2man.pdb   
2302 2324
1obr.pdb   
2581 2581

In this example, for each of the pdb files two output files containing the packing value and fluctuation for each atom of the protein will be generated. Example: For the pdb file ‘1bqc.pdb’ two output files ‘1bqcpack.out’ and ‘1bqcfluk.out’ will be generated containing the packing value and fluctuation respectively for each atom of the protein. 


___________________________________________________
File name: errorcorr.f
Description:
File for the error correction in the atomic fluctuation (directly obtained values from the B-factors) in a data set following our method. 
-----------------------------------------------------

Input data file contents
------------------------
nset, gmin 
xmin, N, file1, file2

Definitions of each parameter.
-------------------------------  
nset  = no. of dataset in a data series.
gmin  = global minimum fluctuation from a separate pdb set
xmin  = minimum fluctuation in a particular series from file1
N     = Number of data points in data series. 
file1 = File containing fluctuations directly from cryst. struc. 
file2 = File containing corrected atomic fluctuation values.   
----------------------------------------------------------------

In our work we have used gmin=0.05


Example:
Input data file for contents for 2 proteins.

2 0.05
0.23 2336 1bqcff.pdb  1bqcffc.out
0.24 2302 2manlf.pdb  2manlfc.out


_______________________________________
Filename: tempcorr.f

Description:
File for computing the atomic fluctuations at a target temperature based on the atomic fluctuations at the crystal temperature. It also computes the average fluctuation and rms fluctuations at the target temperature.
It also identifies the minimum fluctuation value and position in the array.
!-----------------------------------------------------

Input data file contents
------------------------
nset
N, file1, file2, tcry, ttar

Definitions of each parameter.
-------------------------------  
Nset  = no. of data set in the imput file.
N     = No. data points in the data file file1
file1 = Data file containing the atomic fluctuation data at crystal temp.
file2 = Output file storing the atomic fluctuation data at target temp.
tcry  = crystal temp
ttar  = target temp


Example of input file:
2
2302 2manlfc.out  2manlftc.out  298.0 100.0
2581 1obrffc.out  1obrfftc.out  297.0 100.0


___________________________________________
File name: corrcoef.f
Description:
Codes for computing the correlation coefficient between two data series containing N data points each.
----------------------------------------------

Input data file contents
------------------------
lim, file3
file1, file2, N

Definitions of each parameter.
-------------------------------  
lim   = number of data sets
N     = Total number of data points in each dataset.
file1 = first data array file
file2 = second data array file
file3 = output file containing results.

Example of input data file:
---------------------------
1bqcpack.pdb     1BQCfluk.pdb             2336
2manpack.pdb     2manfluk.pdb             2302



_____________________________________________
Filename: resipac.f

Description:
Computation of packing factor of individual residues by taking the average of the packing values of all the non-H atoms of the residue.
!-----------------------------------------------------------------

Input data file contents
------------------------
nset
file1, file2, file3, nall

Definitions of each parameter.
-------------------------------  
nset  = no. of dataset in a data series.
nall  = Number of data points in data series. 
file1 = PDB file of the protein
file2 = File containing corrected atomic packing factor values.   
file3 = Output file containing residue-wise packing factor values.   


Example of the input data file.
-------------------------------

2
1azm.pdb 1azmpack.out 1azmres.out
1kfw.pdb 1kfwpack.out 1kfwres.out


