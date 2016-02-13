module my_rng
#ifdef __INTEL_COMPILER
  Use IFPORT, Only : getpid
#endif
  Use par_zig_mod
  use cosmology
  implicit none 

  Integer, Allocatable :: seedlist(:)

contains


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Subroutine init_random_seed(nval)
  Integer, Intent(in) :: nval
    Integer, Allocatable :: seed(:)
    Integer :: i, n, un, istat, dt(8), pid, t(2), s
    Integer(8) :: count, tms
    Real(8) :: r
    
    Call Random_seed(size = n)
    Allocate(seed(n))
    ! First try if the OS provides a random number generator
    Open(newunit=un,file="/dev/urandom", access="stream", &
         form="unformatted", action="read", status="old", iostat=istat)
    If (istat == 0) Then
       Read(un) seed
       Close(un)
    Else
       ! Fallback to XOR:ing the current time and pid. The PID is
       ! useful in case one launches multiple instances of the same
       ! program in parallel.
       Call System_clock(count)
       If (count /= 0) Then
          t = Transfer(count, t)
       Else
          Call Date_and_time(values=dt)
          tms = ((((((dt(1)-1970)*365_8 + dt(2) )*31_8 + dt(3) )*24 + &
               dt(5) )*60 + dt(6) )*60 + dt(7) )*1000 + dt(8)
          t = Transfer(tms, t)
       End If
       s = Ieor(t(1), t(2))
       pid = getpid() + 1099279 ! Add a prime
       s = Ieor(s, pid)
       If (n >= 3) Then
          seed(1) = t(1) + 36269
          seed(2) = t(2) + 72551
          seed(3) = pid
          If (n > 3) Then
             seed(4:) = s + 37 * (/ (i, i = 0, n - 4) /)
          End If
       Else
          seed = s + 37 * (/ (i, i = 0, n - 1 ) /)
       End If
    End If
    Call Random_seed(put=seed)
    ! stores 
    Allocate(seedlist(0:nval-1))
    Do i=0,nval-1
       Call Random_number(r)
       seedlist(i) = 123456789*r
    Enddo
    Call par_zigset( nval, seedlist, 4 )

  End Subroutine init_random_seed


  subroutine sample_dlsdos( zl, zs, errzs, B, B2, N )
     real(8), intent(in) :: zl,zs, errzs
     real(8), intent(out) :: B,B2
     integer, optional, intent(in) :: N
     integer :: nin,i
     real(8) :: zs_tmp, b_tmp

     nin = merge(N,0, present(N))
     if(nin<1) then
         B=dlsdos(zl,zs)
         B2=B*B         
     else
!        allocate( zsvec(nin), dlsvec(nin), dsvec(nin) )
        B=0. ; B2=0.
        do i=1,nin
            zs_tmp = zs + par_rnor(0)*errzs
            b_tmp = dlsdos(zl,zs_tmp)
            B = B + b_tmp
            B2 = B2 + b_tmp*b_tmp
        enddo
        B=B/nin  ; B2=B2/nin
     endif
 end subroutine


end module my_rng
