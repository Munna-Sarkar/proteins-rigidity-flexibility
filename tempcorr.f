!Filename: tempcorr.f
!File for computing the atomic fluctuations at a target temperature 
!based on the atomic fluctuations at the crystal temperature. 
!It also computes the average fluctuation and rms fluctuations of at 
!the target temperature.
!It also identify the minimum fluctuation value and position in the !array.
!-----------------------------------------------------
!Data structure.
!nset
!N, file1, file2, tcry, ttar
!nset = no. of data set in the imput file.
!N = No. dfata points in the data file file1
!file1=Data file containing the atomic fluctuation data at crystal temp.
!file2 = Output file storing the atomic fluctuation data at target temp.
!tcry = crystal temp
!ttar = target temp
!-----------------------------------------------------
       character*25 file1
       character*25 file2
       dimension x(20000),xnew(20000),na(20000),xi(20000)
       open (unit=12, file='flucTcorr.out', status='new')
       read(*,*) nset
       write(*,*) nset
       lev=1
10     read(*,*) N, file1, file2, tcry, ttar
!       write(*,*) N, file1, file2, tcry, ttar
       open (unit=10, file=file1, status='old')
       open (unit=11, file=file2, status='new')
!Read the fluctuation data.
!---------------------------
       do 15, I=1,N
          read(10,*) na(I),x(I)
          xi(I)=(x(I)*x(I))/tcry
          xnew(I)=sqrt(xi(I)*ttar)
          write(11,*) na(I), xnew(I)
15     continue
!Calculation of average and rms values.
!--------------------------------------
       sum=0.0
       do 20, I=1,N
          sum=sum + xnew(I)
20     continue
       aver=sum/N
!Calculation of RMS fluctuations.
!--------------------------------
       ssum=0.0
       do 30, I=1,N
          ssum=ssum + (xnew(I)-aver)*(xnew(I)-aver)
30     continue
       rms=sqrt(ssum/N)
!-----------------------------
       write(*,51) 'aver=',aver,' rms=',rms
!-----------------------------------------
!Finding min, max values.
!------------------------
       xmax=-999
       xmin=999
       do 40, I=1,N
          if (xnew(I) .lt. xmax) go to 50
          xmax=xnew(I)
          Imax=I
          go to 40
50        if (xnew(I) .gt. xmin) go to 40
          xmin=xnew(I)
          Imin=I
40     continue
       write(12,56)xmin, N, file1, file2, Imin    
!---------------------------------------
       close (unit=10)
       close (unit=11)
!Going for the next dataset.
!---------------------------
       lev=lev + 1
       if (lev .le. nset) go to 10
51     format(A5,1x,F4.2,1x,A5,1x,F4.2)
55     format (I5,1x,F7.4)
56     format(f4.2,1x,I5,1x,A12,1x,A12,1x,I4)
       write(*,*) 'Job successfully completed.'
       stop
       end


