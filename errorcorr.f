!File name: errorcorr.f
!File for  making error correction in the atomic fluctuation
!(obtained values from the B-factors)in a data set following 
!our method. 
!-----------------------------------------------------
!Data format.
!nset, gmin 
!xmin, N, file1, file2 
!nset = no. of dataset in a data series.
!gmin = global minimum fluctuation from a separate pdb set
!xmin = minimum fluctuation in a particular series from file1
!N = Number of data points in data series. 
!file1 = File containing fluctuations directly from cryst. struc. 
!file2 = File containing corrected atomic fluctuation values.   
!----------------------------------------------------------------
       character*25 file1
       character*25 file2
       read(*,*)nset,gmin
       lev=1
       open (unit=10, file='avecorr.out', status='new')
10     read(*,*) xmin, N, file1, file2 
       write(*,*) lev, nset 
       write(*,65) lev, xmin, N, file1, file2 
       open (unit=11, file=file1, status='old')
       open (unit=12, file=file2, status='new')
       del=xmin-gmin
       sum=0.0
       nneg=0
       do 20, I=1,N
          read(11,*)na,x
25        xnew=x-del
          write(12,*)na, xnew
          sum=sum + xnew
20     continue
       close (unit=11)
       close (unit=12)
       ave=sum/(N-nneg)
!!       write(*,*) lev, nneg, ave
!Reopen the file to compute RMSF of the values.
!-----------------------------------------------
       open (unit=12, file=file2, status='old')
       sumds=0.0
       do 30, I=1,N-nneg
          read(12,*) na,xnew
          diff=xnew - ave
          sumds=sumds + diff*diff
30     continue
       close (unit=12)
       rmsf=sqrt(sumds/(N-nneg))
       write(10,56) lev, file2, ave, rmsf
       lev=lev + 1
       if (lev .le. nset) go to 10
55     format (A7,1x,F4.2)
56     format (I3,1x,A12,1x,F5.2,1x,F5.2)
61     format (F4.2,1x,I5,1x,A14,1x,A14,I3)
64     format (F4.2,1x,I4,1x,A12,1x,A12,3x,I4)
65     format (I4,1x,F4.2,1x,I4,1x,A14,1x,A14)  
       stop
       end

