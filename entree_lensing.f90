Module entree_lensing
  !---------------
  Implicit none
Contains
  !-----------------------------------------------------------------------------------------
  
  Subroutine lecture_NGVS(NGVS_maille,nlBG,RA,DEC,e1,e2,de1,de2,z,err_z,snr_win,mag_y)
    
    Implicit none
    integer::ios=0,i
    integer::nlBG
    character(len=60)::NGVS_maille
    character(len=400) :: line
    double precision,dimension(:),allocatable:: RA,DEC,e1,e2,de1,de2,z,err_z,snr_win,mag_y
    
    open(10,file=NGVS_maille,status='old',action='read',iostat=ios)
    if(ios .ne. 0)stop"Pb !"
    nlBG=0
    do while(ios .eq. 0)
       read(10,*,iostat=ios) line

       if(index(trim(line),'#')==1) cycle
       nlBG=nlBG+1
    enddo
    nlBG=nlBG-1
    close(10)
    
    allocate(RA(nlBG),DEC(nlBG),e1(nlBG),e2(nlBG),de1(nlBG),de2(nlBG),z(nlBG),err_z(nlBG),snr_win(nlBG),mag_y(nlBG))
    
    open(11,file=NGVS_maille,status='old',action='read',iostat=ios)
    if(ios .ne. 0)stop"Pb !"
    read(11,*) 
    do i=1,nlBG-1
       read(11,'(A)',iostat=ios) line
       if(index(trim(line),'#')==1) cycle
       read(line,*) RA(i),DEC(i),e1(i),e2(i),de1(i),de2(i),z(i),err_z(i),snr_win(i),mag_y(i)
    enddo
    
    close(11)
  End Subroutine lecture_NGVS
  !---------------------------------------------------------------------------------------
  
  Subroutine lecture_SDSS(SDSS_maille,nlFG,objid,z2,zErr,zPhoto,zPhotoErr,type,ra2,dec2,veldisp,ved,mabs_r)
    
    Implicit none
    integer::ios=0,i
    integer::nlFG
    character(len=60)::SDSS_maille
    double precision,dimension(:),allocatable:: objid,z2,zErr,zPhoto,zPhotoErr,type &
               ,ra2,dec2,veldisp,ved,mabs_r
    
    open(12,file=SDSS_maille,status='old',action='read',iostat=ios)
    if(ios .ne. 0)stop"Pb !"
    nlFG=0
    do while(ios .eq. 0)
       read(12,*,iostat=ios)
       nlFG=nlFG+1
    enddo
    close(12)
    
    allocate(objid(nlFG),z2(nlFG),zErr(nlFG),zPhoto(nlFG),zPhotoErr(nlFG),type(nlFG),ra2(nlFG),dec2(nlFG),veldisp(nlFG),ved(nlFG),mabs_r(nlFG))
    
    open(13,file=SDSS_maille,status='old',action='read',iostat=ios)
    if(ios .ne. 0)stop"Pb !"
    
    do i=1,nlFG-1
       read(13,*)objid(i),z2(i),zErr(i),zPhoto(i),zPhotoErr(i),type(i) &
               ,ra2(i),dec2(i),veldisp(i),ved(i),mabs_r(i)
    enddo
    
    close(13)
  End Subroutine lecture_SDSS
  
  !-------------------------------------------------------------------------------------------
End Module entree_lensing
