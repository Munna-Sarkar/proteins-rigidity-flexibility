!Computing the correlation coefficient between two data series containing !N datapoints each.
!N== Total number of datapoints.
!file1== first data array file
!file2== secons array data file
!file3== output file containing results.
!-------------------------------------------------------------
       character*25 file1
       character*25 file2
       character*25 file3
       dimension x(15000), y(15000)
       read(*,*) lim, file3
!       write(*,*) lim, file3
       lev=1
       open (unit=12, file=file3, status='new')
222    read(*,*) file1, file2, N
!       write(*,*) file1, file2, N
       open (unit=10, file=file1, status='old')
       open (unit=11, file=file2, status='old')

       do 10, I=1,N
          read(10,*) nsa, x(I)
          read(11,*) nsb, y(I)
10     continue
       close (unit=10)
       close (unit=11)
!-------------------------
15     sumx=0.0
       sumy=0.0
       sds1=0.0
       sds2=0.0
       sds3=0.0
       do 20, I=1,N
          sumx=sumx + x(I)
          sumy=sumy + y(I)
20     continue
       avex=sumx/N
       avey=sumy/N
!Now compute the rms fluctuations of the data.
!---------------------------------------------
       sumx2=0.0
       sumy2=0.0
       do 30, I=1,N
          sumx2=sumx2 + (x(I) - avex)**2.0
          sumy2=sumy2 + (y(I) - avey)**2.0
30     continue
       xfluc=sqrt((sumx2)/N)
       yfluc=sqrt((sumy2)/N)
!Now compute the correlation between the two data series.
!--------------------------------------------------------
       do 40, I=1,N 
          sds1=sds1 + (x(I)- avex)*(y(I)- avey)
          sds2=sds2 + (x(I)- avex)*(x(I)- avex)
          sds3=sds3 + (y(I)- avey)*(y(I)- avey)
40     continue
       cnum=sds1
       cden=sqrt(sds2*sds3)
       corr=cnum/cden
       write(*,99) file1, '  corr-coef = ', corr
       write(12,75) lev, file1, file2, corr
       lev=lev + 1
       if (lev .le. lim) go to 222
75     format (I5,1x,A15,1x,A15,1x,F5.2)
90     format (A14,1x,F4.2,1x,F4.2,4x,F4.2,1x,F4.2,5x,F4.2)
99     format (A14,1x,A14,1x,F5.2)
       stop
       end

