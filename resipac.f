!Filename: resipac.f
!Computation of packing factor of individual residues by taking the !average of the packing values of all the non-H atoms of the residue.
!------------------------------------------------------------------
      character*25 file1   !pdb file
      character*25 file2   !pack file
      character*25 file3   !respac file (output file).
      dimension xa(50000),ya(50000),za(50000),
     & atom(50000),numa(50000),bnam(50000),res(50000),ct(50000),
     & anam(50000),atyp(50000), btyp(50000),
     & nseg(50000),rf(50000)
!---------------------------------------------------------------------
      read(*,*) nset
      levset=1
10    read(*,*) file1, file2, file3, nall

      open (unit=10, file=file1, status='old')
      open (unit=11, file=file2, status='old')
      open (unit=12, file=file3, status='new')
!Now open the PDB file.
!------------------------------
      open (unit=10, file=file1, status='old')
!Now read all the atoms in the PDB file including water O and ligand.
!---------------------------------------------------------------------
      sum=0.0
      num=0
      do 20, I=1,nall
       read(10,99)atom(I),numa(I),ct(I),atyp(I),bnam(I),res(I),nseg(I),
     & xa(I),ya(I),za(I),cf,rf(I)
20     continue
       close(unit=10)
       I=1  
22     if (I .ne. 1) go to 30 
       refres=res(I) 
       write(*,56) refres, nseg(I)
30     if (res(I) .ne. refres) go to 35
       read(11,*)x,y
       sum=sum + y
       num=num+1
       write(*,54)refres, y, num, I
       go to 25
35     ave=sum/num
       nref=nseg(I-1)
       go to 50
45     nref=nseg(I)
50     write(12,*)nref, ave
       write(*,57)nref, ave, num
       write(*,*)'-----------------------------'
       write(*,*)'Going to next refres'
       refres=res(I) 
       sum=0.0
       num=0
      write(*,56) refres, nseg(I)
      go to 60
25    I=I+1 
60    If (I .le. nall) go to 22
21    format(A4,A5)    !format(A5,A4)
54    format (A3,1x,F5.2,1x,I3,1x,I5)
55    format (A3,1x,F5.2,1x,F5.2,I3,1x,A5)
56    format (A3,1x,I3)
57    format (I5,2x,F5.2,2x,I2)
       close (unit=10)
       close (unit=11)
       close (unit=12)
99    format(A4,1x,I6,1x,A1,A1,A2,1x,A3,2x,I4,4x,3(F8.3),     
     & 2x,F4.2,F6.2)
      levset=levset + 1
      if (levset .le. nset) go to 10
      write(*,*) 'Job completed successfully'
      stop
      end





