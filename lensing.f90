Program lensing
  !---------------
  use entree_lensing
  use my_rng

  Implicit none
  
  integer::i=0,j=0,l=0,m=0,ios
  integer::nlNGVS=0,nlSDSS=0
  
  integer::dim1=9 !dispersion interval
  integer::dim2=9  !magnitude interval 
  
  character(len=80)::NGVS_maille(116)
  character(len=80)::SDSS_maille(116)
  double precision,allocatable::vdinf(:), vdsup(:)
  double precision,allocatable::maginf(:), magsup(:)
  character(len=10),allocatable::fich_sortie(:)

  double precision,dimension(:),allocatable:: RA,DEC,e1,e2,de1,de2,z,err_z,snr_win,mag_y
  double precision,dimension(:),allocatable:: objid,z2,zErr,zPhoto,zPhotoErr,type,ra2,dec2 &
       ,veldisp,ved,mabs_r
  
  double precision,parameter::rmin=0.00333,rmax=0.1667
  double precision::lrmin=0.,lrmax=0.,dlr=0.,dr=0.,lr=0.
  double precision::nb=11.
  integer::dim=10
  
  double precision::x=0., y=0., r=0.
  double precision::alpha=0.
  double precision::et=0., er=0.
  integer::k=0.

  double precision,dimension(:),allocatable::shearT, shearR, N, r2
  double precision::sigma 
  double precision::moy=0.
  double precision::var=0.
!------------------------------------------------------------
  double precision::moysig
  double precision,parameter::c=299792458.
  double precision,parameter::Ho=70.
  double precision::Betha=0.53,B, B2
  double precision::cmpt=0.
  double precision::ds=0., dl=0. 

  double precision,dimension(:),allocatable::sumB
  double precision::moydisp=0. ,w=0., wiso=0.
  double precision::moysumB=0.
  double precision::err_vstar=0.

  double precision::mag=0.
  integer::Nlens=0

  Call set_cosmology(0.3d0,0.7d0,0.3d0,Ho)
  call init_random_seed(1)   

  allocate( vdinf(dim1), vdsup(dim1), maginf(dim2), magsup(dim2), fich_sortie(dim1) ) 

!------------------------------------------------------------
  open(10,file='/home/amarchal/file_NGVS.txt',status='old',action='read',iostat=ios)
    if(ios .ne. 0)stop"Pb !"
    do i=1,116
       read(10,*)NGVS_maille(i)
    enddo
  close(10)

  open(11,file='/home/amarchal/file_SDSS.txt',status='old',action='read',iostat=ios)
    if(ios .ne. 0)stop"Pb !"
    do i=1,116
       read(11,*)SDSS_maille(i)
    enddo
  close(11)

  open(12,file='/home/amarchal/file_disp.txt',status='old',action='read',iostat=ios)
     if(ios .ne. 0)stop"Pb !"
     read(12,*)
     do i=1,dim1
        read(12,*)vdinf(i),vdsup(i)
     enddo
  close(12)

  open(13,file='/home/amarchal/file_mag.txt',status='old',action='read',iostat=ios)
     if(ios .ne. 0)stop"Pb !"
     read(13,*)
     do i=1,dim2
        read(13,*)maginf(i), magsup(i)
     enddo
  close(13)

  open(14,file='/home/amarchal/file_sortie.txt',status='old',action='read',iostat=ios)
     if(ios .ne. 0)stop"Pb !"
     read(14,*)
     do i=1,dim1
        read(14,*)fich_sortie(i)
     enddo
  close(14)
  !------------------------------------------------------------
  lrmin=log(rmin)
  lrmax=log(rmax)
  dlr=(lrmax-lrmin)/nb

  do m=1,dim2 
  
     allocate(shearT(dim),shearR(dim),N(dim),r2(dim),sumB(dim))

     do i=1,dim
        shearT(i)=0.; N(i)=0.; r2=0.; sumB(i)=0.; shearR(i)=0.
     enddo

!-------------------------------------------------------------------------------------
     do l=1,116
        call lecture_NGVS('/home/amarchal/NGVS_data/'//NGVS_maille(l),nlNGVS,RA,DEC,e1,e2,de1,de2,z,err_z,snr_win,mag_y)
        call lecture_SDSS('/home/amarchal/new_SDSS_data/'//SDSS_maille(l),nlSDSS,objid,z2,zErr,zPhoto,zPhotoErr,type,ra2,dec2,veldisp,ved,mabs_r)
!-----------------------calcul de sigmaS----------------------
        do i=1,nlNGVS-1  
           if (e1(i) .ne. 0 .and. e2(i) .ne. 0 .and. e1(i) .gt. -1 &
           .and. e1(i) .lt. 1) then                    
              moy=moy+(1./nlNGVS)*e1(i)
           endif
        enddo
  
        do i=1,nlNGVS-1 
           if (e1(i) .ne. 0 .and. e2(i) .ne. 0 .and. e1(i) .gt. -1 &
           .and. e1(i) .lt. 1) then 
              var=var+(1./nlNGVS)*((e1(i)-moy)**2.)
           endif
        enddo
        sigma=sqrt(var)
        moy=0.; var=0.
!------------------------------------------------------------
!   wiso = 1. / sigma/sigma 
        do j=1,nlSDSS-1   
!        if (z2(j) .gt. 0 .and. z2(j) .lt. 0.5 .and. veldisp(j) .gt. 0. .and. veldisp(j) .lt. 180. ) then 
           !if (z2(j) .gt. 0.05 .and. z2(j) .lt. 0.6 .and. veldisp(j) .gt. vdinf(m)  &
           !.and. veldisp(j) .lt. vdsup(m)) then 
            
           if (z2(j) .gt. 0.05 .and. z2(j) .lt. 0.6 .and. mabs_r(j) .gt. magsup(m) &
           .and. mabs_r(j) .lt. maginf(m) .and. veldisp(j) .gt. 10. .and. veldisp(j) .lt. 500.) then 

              do i=1,nlNGVS
                 if (z(i)+err_z(i) .gt. 0.7 .and. e1(i) .ne. 0 .and. e2(i) .ne. 0 .and. e1(i) .gt. -1 &
                 .and. e1(i) .lt. 1) then 

                    wiso = 1./( sigma*sigma + 0.5*(de1(i)*de1(i) + de2(i)*de2(i)))

                    x=-(RA(i)-ra2(j))/cos(dec2(j)*pi/180)
                    y=DEC(i)-dec2(j)
                    r=sqrt(x**2+y**2)
                    lr=log(r)
                    k=(lr-lrmin)/dlr
                 
                    if (k > 0 .and. k < nb) then
                  
!                   call sample_dlsdos( z2(j), z(i), err_z(i), B, B2, 0 )
                       call sample_dlsdos( z2(j), z(i), err_z(i), B, B2, 100 )
!                    ds=dda(0.d0,z(i))
!                    dl=dda(z2(j),z(i))
!                    B = max(dl/ds,0.) ; B2=B*B

!                    dl=(c/Ho)*(z2(j)+0.5*(1-qo)*z2(j)**2.) !Ordre 2 loi de Hubble
!                    ds=(c/Ho)*(z(i)+0.5*(1-qo)*z(i)**2.)
!                    B = max((ds-dl),0.)/ds

                       alpha=atan2(y,x)
                       et=-(e1(i)*cos(2.*alpha) + e2(i)*sin(2.*alpha))
                       er= e2(i)*cos(2.*alpha) - e1(i)*sin(2.*alpha)
                       shearT(k)=shearT(k)+(et)*B*wiso
                       shearR(k)=shearR(k)+(er)*B*wiso
                       N(k)=N(k)+1
                    
                       sumB(k)=sumB(k)+B2*wiso
           
 ! write(14,*)x*60.,y*60.
                    endif
                    x=0.; y=0.; r=0.; k=0.; et=0.; er=0.; alpha=0.; lr=0.; dl=0.; ds=0.; B=0.; B2=0; wiso=0.
                 endif
              enddo
              moydisp=moydisp+(veldisp(j)*veldisp(j)*(1/ved(j)**2.))
              w=w+(1/ved(j)**2.)
              mag=mag+mabs_r(j)
              Nlens=Nlens+1
           endif
        enddo

        do i=1,dim
           e1(i)=0.; e2(i)=0.; RA(i)=0.; ra2(i)=0.; DEC(i)=0.; dec2(i)=0.
        enddo
        deallocate(objid,z2,zErr,zPhoto,zPhotoErr,type,ra2,dec2 &
          ,veldisp,ved,mabs_r,RA,DEC,e1,e2,de1,de2,z,err_z,snr_win,mag_y)
        sigma=0.
        nlNGVS=0.; nlSDSS=0.
        write(*,*)'done',l
     enddo

     do i=1,dim
        r2(i)=exp(lrmin+(i)*dlr)
!     sumB(i)=sumB(i)/N(i)
!     shearT(i)=shearT(i)/(N(i)*sumB(i))
!     ErrT(i)=0.26331349542305638/sqrt(N(i))!*sumB(i))   !A vérifier
        shearT(i)=shearT(i)/sumB(i)
        shearR(i)=shearR(i)/sumB(i)
        write(*,'(30(x,F15.6))') r2(i)*60,shearT(i),shearR(i),1./sqrt(sumB(i)),N(i)
     enddo
!-----------------------Resultats----------------------
     open(14,file=fich_sortie(m),status='replace', action='write',iostat=ios)
     if(ios .ne. 0)stop"Pb !"
     write(14,*)"r(')    shearT   shearR   Err_ShearT    N"
     do i=1,dim
           write(14,'(30(x,F15.6))') r2(i)*60,shearT(i),shearR(i),1./sqrt(sumB(i)),N(i)
     enddo
     close(14) 
!-------------------------------------------------------------
!---Vel disp moy (moyenne pondérée) star + mag abs (moyenne)---
     write(*,*)'--------------------------------'
     moydisp=moydisp/w
     moydisp=sqrt(moydisp)    
     err_vstar=1./(sqrt(w))
     mag=mag/Nlens
!-------------------------------------------------------------
!-------------Ecriture resultat vitesses de dispersions-------
     open(15,file='test/resultat.txt',status='old',position='append', action='write',iostat=ios)
     if(ios .ne. 0)stop"Pb !"
        write(15,'(30(x,F15.6))')moydisp,err_vstar,mag
     close(15) 
  
     moydisp=0.; err_vstar=0.; w=0.; mag=0.; Nlens=0
     deallocate(r2,shearT,shearR,sumB,N)
  enddo
end program lensing
