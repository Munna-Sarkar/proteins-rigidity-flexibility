!Filename: Pack.f
!File for computing the atom-wise packing fraction following our
!definition and also for computing the atomic fluctuations using
!the B-factors from the respective pdb of the crystal structure
!of the protein. 
!-------------------------------------------------------------
!Data format.
!-------------
!file0
!nset, rsph !number of datasets, radius of sampling sphere.
!-------------------------------------------------------------
      character*25 file0
      character*15 file1
      character*4  naam 
      character*5  ext
      character*4  ext2
      character*4  pack
      character*4  fluk
      character*15 file2
      character*15 file3

      dimension xa(50000),ya(50000),za(50000),radf(50000),
     & atom(50000),numa(50000),bnam(50000),res(50000),ct(50000),
     & anam(50000),atyp(50000), btyp(50000), ctyp(2),r(15),
     & nseg(50000),rf(50000),aq(50000),q(50000),e(12),s(12),
     & ep(50000),si(50000),anm(50000),resi(50000),xnem(50000)

!---------------------------------------------------------------
!!!!Reading the atom type and vdw radius.
!========================================
      read(*,*) file0
      open (unit=15, file=file0, status='old')
      read(15,*)npar !No. of parameters.
         do 7, k=1,npar
           read(15,95) btyp(k), r(k)
!!           write(*,95) btyp(k), r(k)
7        continue
95    format(A1,1x,F4.2)  
96    format(A1,1x,A1)  
      close (unit=15)

!-------------------------------------------------
      ext2='.out'
      pack='pack'
      fluk='fluk'
      read(*,*) nset, rsph 
      open (unit=34, file='result.out', status='new')
      open (unit=35, file='corrcoef.data', status='new')
      open (unit=36, file='list.out', status='old')      
      open (unit=37, file='allpac.out', status='new')
      open (unit=38, file='allfluc.out', status='new')

      levset=1
555   read(36,21)naam,ext
      write(*,21)naam,ext
      write(*,*)levset

      open (unit=45, file='tm.out', status='old')
      write(45,*)naam,ext         !file1
      write(45,*)naam,pack,ext2   !file2
      write(45,*)naam,fluk,ext2   !file3
      close (unit=45)

      open(unit=45, file='tm.out', status='old')
      read(45,*) file1
      write(*,*) file1
        read(45,*) file2
        write(*,*) file2
          read(45,*) file3
          write(*,*) file3
      close (unit=45)
      write(34,*)file1

!Now open the NEW files for Storing data.
!-----------------------------------------
      open (unit=16, file=file2, status='new')
      open (unit=32, file=file3, status='new')
      open (unit=10, file=file1, status='old')
      open (unit=39, file='respac.data', status='old')
      write(39,*)file1, file2, file2

!Now read all the atoms in the PDB file.
!nall=all atoms in crystal. npro=all atoms in protein.
!---------------------------------------
      read(10,*)npro,nall 
      do 510, I=1,nall
       read(10,99)atom(I),numa(I),ct(I),atyp(I),bnam(I),res(I),nseg(I),
     & xa(I),ya(I),za(I),cf,rf(I)

!!Now assign atomic radius to each atom.
!---------------------------------------
515      do 520, k=1,npar
            if (btyp(k) .ne. atyp(I)) go to 520 
              radf(I)=r(k)
            go to 510
520      continue
       write(*,*) 'Terminated. Unknown atom type at ', 'I=', I
       go to 175
510   continue
21    format(A4,A5) 
       close (unit=8)
       close (unit=10)
       close (unit=20)
99    format(A4,1x,I6,1x,A1,A1,A2,1x,A3,2x,I4,4x,3(F8.3),     
     & 2x,F4.2,F6.2)
!==========================================================
        k1=1
        k2=2
        k3=3
        k4=4
        k5=5
!--------------------------------------------------------------
!Computation of the OCCUPIED VOLUME and atomic PACKINNG FRACTION 
!using the non-H atoms from the PDB file.
!=======================================
!Computing the packing fraction pacf.
!-------------------------------------
      nha=1 !Initial value of no. of non-H atoms.
540   cont=0.0
!---------------------------------------------------------
!Now generate random points uniformly within the big sphere 
!of radius (rsph)around the selected non-H atom. 
!-----------------------------------------------------------
       IA=16807
       IM=2147483647
       AM=1./IM
       IQ=127773
       IR=2836
       MASK=123459876
       idum=254739  !352743 used as a seed for random numbers.
       lev=1
!----------------------------------------
       vsph=(4.0*3.141/3.0)*(rsph**3.0)
       np=2000 
       do 55, I=1,np
            idum=ieor(idum,MASK)
            k=idum/IQ
            idum=IA*(idum-k*IQ) - IR*k
            if (idum .lt. 0) idum=idum + IM
            rand1=AM*idum
!
              idum=ieor(idum,MASK)
              k=idum/IQ
              idum=IA*(idum-k*IQ) - IR*k
              if (idum .lt. 0) idum=idum + IM
              rand2=AM*idum
!
                idum=ieor(idum,MASK)
                k=idum/IQ
                idum=IA*(idum-k*IQ) - IR*k
                if (idum .lt. 0) idum=idum + IM
                rand3=AM*idum
!
                    idum=ieor(idum,MASK)
                    k=idum/IQ
                    idum=IA*(idum-k*IQ) - IR*k
                    if (idum .lt. 0) idum=idum + IM
                    rand4=AM*idum

                 arad=rand1*rsph
                 thet=rand2*180.0
                 phi=rand3*360.0

              xpt=xa(nha) + arad*sin(thet)*cos(phi)
              ypt=ya(nha) + arad*sin(thet)*sin(phi)
              zpt=za(nha) + arad*cos(thet)
!!!!       write(*,*) 'Points coordinates ', xpt, ypt, zpt

!Check if the generated point is on any protein atoms within the sphere.
!---------------------------------------------------------
       do 60, Im=1,nall   !nall=the total atoms of the protein.
           dxp=xa(Im) - xpt
           dyp=ya(Im) - ypt
           dzp=za(Im) - zpt
            disp=sqrt(dxp*dxp + dyp*dyp + dzp*dzp)
            if (disp .le. radf(Im)) go to 65
60      continue
        go to 55
65      cont=cont + 1.0
55      continue
        occuvol=vsph*cont/np
        pacf=cont/np           !pacf=packing fraction
        write(16,72) nha, pacf
!Calculate average pacf and its rms fluctuation.
!-----------------------------------------------
        sumpac=sumpac + pacf
!-----------------------------------------
!Computing the atomic fluctuations from the B-factors in the pdb.
!------------------------------------------------------------------
        w=8.0*3.14*3.14
        fluc=sqrt((rf(nha)/w))
        write(32,72) nha,fluc
!----------------------------
170     nha=nha + 1
        if (nha .le. npro) go to 540
        nha=nha - 1
        close (unit=16)
        close (unit=32)
!--------------------------------------------------------
!Calculate average pacf, fluc and their rms fluctuations.
      sumpac=0.0
      sumfluc=0.0
      open (unit=16, file=file2, status='old')
      open (unit=32, file=file3, status='old')
      do 80, I=1,nha
          read(16,*)is,pac
          read(32,*)is,fluc
            sumpac=sumpac + pac
            sumfluc=sumfluc + fluc
80    continue
        avepac=sumpac/nha
        avefluc=sumfluc/nha
        close (unit=16)
        close (unit=32)
!----------------------------------------------
      open (unit=16, file=file2, status='old')
      open (unit=32, file=file3, status='old')
      sumpacs=0.0
      sumflu=0.0
      do 82, I=1,nha 
          read(16,72)is,pac
          dpac=pac - avepac
          sumpacs=sumpacs + dpac*dpac
          read(32,72) is,fluc 
          dflu=fluc - avefluc
          sumflu=sumflu + dflu*dflu
82    continue
          pacflu=sqrt(sumpacs/nha)
          fluflu=sqrt(sumflu/nha)
      write(*,777)'Average pac =',avepac, '  RMS pacf fluc =',pacflu
      write(*,777)'Average fluc =',avefluc,' RMS fluc fluc =',fluflu 
      write(*,*) 'Calculation of RMS fluc for pac done'
      write(*,*) 'NHA = ', nha
!----------------------------------------------
      write(35,*) file2, file3, nha
!========================================
      write(34,777)'Average pac =',avepac, '  RMS pacf fluc =',pacflu
      write(34,777)'Average fluc =',avefluc,' RMS fluc fluc =',fluflu 
      write(34,*)'--------------------------------------------------'
      write(*,*)'-------------------------------'
      write(37,*) levset,avepac
      write(38,*) levset,avefluc
      levset=levset + 1
      if (levset .le. nset) go to 555
71     format(I6,2x,F5.3)
72     format(I4,1x,F5.2)
73     format(I4,1x,F7.2)
777    format(A14,1x,F4.2,1x,A18,1x,F4.2)
778    format(A14,1x,F4.2)
100   format(A4,1x,I6,1x,A1,A1,A2,1x,A3,2x,I4,4x,3(F8.3))
       write(*,*)'       '
!-------------------------------------------
      write(*,*) '**  JOB SUCCESSFULLY COMPLETED **'
175    stop
       end
