module scatter_module
  use, intrinsic:: iso_fortran_env, only : d8 => real64
  implicit none
  real(kind=d8),parameter ::pi = dacos(-1.0_d8)

  abstract interface
  subroutine potential_ab(x,V)
      import :: d8
      real(kind=d8), intent(in)::x 
      real (kind=d8), intent(inout)::V 
  end subroutine
 end interface

 type onedscatter 
 integer :: n_dvr
 real(kind=d8)::r_start,r_end,mass,r0,esk,delta
 real(kind=d8)::t_total,dt,r_flux,e_start,e_end,de 
 real(kind=d8)::pabs,cabs1,cabs2,rabs1,rabs2
 procedure(potential_ab),pointer,nopass::potential=>null()
 contains 
 procedure,public ::propagation_scatter_so
 procedure,public ::propagation_scatter_crwp
 end type onedscatter
 
 contains 

 subroutine propagation_scatter_so(this,nz,rs1,rs2,rms,r_0,esk_ev,delt,totalt,deltt,posflux,e1,e2,delte,p_abs,c_abs1,c_abs2,r_abs1,r_abs2,potential)
 implicit none
 class(onedscatter) ::this
 integer,intent(in)::nz
 real(kind=d8),intent(in)::rs1,rs2,rms,r_0,esk_ev,delt,totalt,deltt
 real(kind=d8),intent(in)::posflux,e1,e2,delte,p_abs,c_abs1,c_abs2,r_abs1,r_abs2
 procedure(potential_ab)::potential
 
 integer::ne,i,j,k,nabsl,nabsr,nstep
 real(kind=d8),allocatable::esk0(:),rq0(:),trans(:,:),t_eng(:)
 real(kind=d8),allocatable::t_dvr(:,:),fsin(:),fcos(:),dvec(:,:)
 real(kind=d8),allocatable::fabs1(:),fabs2(:),vpot(:)
 complex(kind=8),allocatable::w0(:),w1(:),w2(:),aef(:),ae0(:)
 complex(8),allocatable::phie0(:),phie0p(:),phie(:)
 real(kind=d8)::pi,ev2hartree,length,amu2au,dr,time,V,probs
 complex(kind=8):: zi,s
 pi=dacos(-1.0_d8)
 ev2hartree=27.2114_d8
 zi=(0.0_d8,1.0_d8)
 amu2au=1822.888_d8
 
 this%n_dvr=nz
 this%r_start=rs1
 this%r_end=rs2
 this%r0=r_0
 this%mass=rms
 this%esk=esk_ev
 this%delta=delt
 this%t_total=totalt
 this%dt=deltt
 this%r_flux=posflux
 this%e_start=e1
 this%e_end=e2
 this%de=delte
 this%pabs=p_abs
 this%cabs1=c_abs1
 this%cabs2=c_abs2
 this%rabs1=r_abs1
 this%rabs2=r_abs2
 this%potential=>potential

 this%esk=this%esk/ev2hartree
 this%mass=this%mass*amu2au
 length=this%r_end-this%r_start
 !energy_grid
 ne=(this%e_end-this%e_start)/this%de+1
 allocate(esk0(ne))
 this%e_start=this%e_start/ev2hartree
 this%e_end=this%e_end/ev2hartree
 this%de=this%de/ev2hartree
 do i=1,ne
  esk0(i)=this%e_start+(i-1)*this%de
 end do
 !spatial_grid
 allocate(trans(this%n_dvr,this%n_dvr),t_eng(this%n_dvr),rq0(this%n_dvr))
 dr=(this%r_end-this%r_start)/(this%n_dvr+1)
 do i=1,this%n_dvr
  rq0(i)=this%r_start+i*dr
  t_eng(i)=(i*pi/length)**2/2.0_d8/this%mass
  do j=1,i
  trans(i,j)=dsqrt(2.0/length*dr)*dsin(i*j*pi/(this%n_dvr+1))
  trans(j,i)=trans(i,j)
  end do
 end do
 !T in dvr
 allocate(t_dvr(this%n_dvr,this%n_dvr))
 do i=1,this%n_dvr
  t_dvr(i,i)=pi*pi/4.0_d8/this%mass/length/length*((2*(this%n_dvr+1)**2+1.0_d8) &
   /3-1.0_d8/dsin(pi*i/(this%n_dvr+1))**2)
  do j=1,i-1
  t_dvr(i,j)=pi*pi/4.0_d8/this%mass/length/length*(-1)**(i-j)*(1.0_d8/dsin(pi*(i-j) &
   /2.0_d8/(this%n_dvr+1))**2-1.0_d8/dsin(pi*(i+j)/2.0_d8/(this%n_dvr+1))**2)
  t_dvr(j,i)=t_dvr(i,j)
  end do
 end do
 !flux sin and cos 
 allocate(fsin(this%n_dvr),fcos(this%n_dvr))
 do i=1,this%n_dvr
 fsin(i)=dsqrt(2.0_d8/length)*dsin(i*pi*(this%r_flux-this%r_start)/length)
 fcos(i)=dsqrt(2.0_d8/length)*(i*pi/length)*dcos(i*pi*(this%r_flux-this%r_start)/length)
 end do
 !init
 allocate(w0(this%n_dvr),w1(this%n_dvr),w2(this%n_dvr),ae0(ne),aef(ne))
 call init(this%n_dvr,rq0,this%mass,this%r_flux,this%esk,this%delta,w2,ne,esk0,aef)
 call transf(this%n_dvr,trans,w2,w0)
 call init(this%n_dvr,rq0,this%mass,this%r0,this%esk,this%delta,w2,ne,esk0,ae0)
 call transf(this%n_dvr,trans,w2,w1)
 allocate(dvec(this%n_dvr,this%n_dvr))
 call sincdivi(this%n_dvr,rq0,dvec)
 !init v
 allocate(vpot(this%n_dvr))
 do i=1,this%n_dvr
  call this%potential(rq0(i), V)
  vpot(i)=V
 end do
 !absorb pot
 allocate(fabs1(this%n_dvr),fabs2(this%n_dvr))
 call wfabs(this%pabs,this%n_dvr,rq0,this%cabs1,this%rabs1,this%cabs2,this%rabs2, &
       nabsl,fabs1,nabsr,fabs2,this%dt)
  
 allocate(phie(ne),phie0(ne),phie0p(ne))
 phie0=(0.0_d8,0.0_d8)
 phie=(0.0_d8,0.0_d8)
 phie0p=(0.0_d8,0.0_d8)
 
 time=0.0_d8
 nstep=this%t_total/this%dt



 !propagation
 do i=1,nstep
  call prog_so(this%n_dvr,w1,w2,vpot,t_eng,trans,this%dt,nabsL,fabs1,nabsR,fabs2)
  time=time+this%dt
  call flux(this%n_dvr,w1,fsin,fcos,ne,esk0,phie0,phie0p,time)
  call flux1(this%n_dvr,w1,fsin,fcos,ne,esk0,phie,w0,time)
  
 end do
 open(111,file='probability_sincdvr_so',status='unknown') 
 do i=1,ne 
 s = phie(i) / (pi*2 * ae0(i) * conjg(aef(i))) * this%dt 
 probs=conjg(s)*s 
 write(111,*) esk0(i)*ev2hartree,probs 
 end do 
 close(111)
  
 end subroutine

 subroutine propagation_scatter_crwp(this,nz,rs1,rs2,rms,r_0,esk_ev,delt,totalt,posflux,e1,e2,delte,p_abs,c_abs1,c_abs2,r_abs1,r_abs2,potential)
 implicit none
 class(onedscatter) ::this
 integer,intent(in)::nz
 real(kind=d8),intent(in)::rs1,rs2,rms,r_0,esk_ev,delt,totalt
 real(kind=d8),intent(in)::posflux,e1,e2,delte,p_abs,c_abs1,c_abs2,r_abs1,r_abs2
 procedure(potential_ab)::potential
 
 integer::ne,i,j,k,nabsl,nabsr,nstep,idabsr,anansr
 real(kind=d8),allocatable::esk0(:),rq0(:),trans(:,:),t_eng(:),vpot(:)
 real(kind=d8),allocatable::esk0s(:),fsin(:),fcos(:),fabs1(:),fabs2(:)
 real(kind=d8)::pi,ev2hartree,length,amu2au,Vmin,Vmax,EKmin,EKmax
 real(kind=d8)::Hmax,Hmin,ascale,bscale,dt,V,dr,en,factid,probs
 complex(kind=8),allocatable::w0(:),w1(:),w2(:),aef(:),ae0(:),cwk1(:),cwk2(:)
 complex(8),allocatable::phie0(:),phie0p(:),phie(:),cwk3(:)
 complex(8)::s

 pi=dacos(-1.0_d8)
 ev2hartree=27.2114_d8
 amu2au=1822.888_d8


 this%n_dvr=nz
 this%r_start=rs1
 this%r_end=rs2
 this%r0=r_0
 this%mass=rms
 this%esk=esk_ev
 this%delta=delt
 this%t_total=totalt
 this%r_flux=posflux
 this%e_start=e1
 this%e_end=e2
 this%de=delte
 this%pabs=p_abs
 this%cabs1=c_abs1
 this%cabs2=c_abs2
 this%rabs1=r_abs1
 this%rabs2=r_abs2
 this%potential=>potential
 
 this%esk=this%esk/ev2hartree
 this%mass=this%mass*amu2au
 length=this%r_end-this%r_start
 
 !spatial_grid
 allocate(trans(this%n_dvr,this%n_dvr),t_eng(this%n_dvr),rq0(this%n_dvr))
 dr=(this%r_end-this%r_start)/(this%n_dvr+1)
 do i=1,this%n_dvr
  rq0(i)=this%r_start+i*dr
  t_eng(i)=(i*pi/length)**2/2.0_d8/this%mass
  do j=1,i
  trans(i,j)=dsqrt(2.0/length*dr)*dsin(i*j*pi/(this%n_dvr+1))
  trans(j,i)=trans(i,j)
  end do
 end do
 Vmin=-1.1_d8/ev2hartree
 !init V
 allocate(vpot(this%n_dvr))
 do i=1,this%n_dvr
  call this%potential(rq0(i), V)
  vpot(i)=V+Vmin
 end do
 !scale
 Vmax=0.6_d8/ev2hartree+Vmin
 EKmax=1.0_d8/dr**2*pi*pi/2.0_d8/this%mass
 EKmin=0.0_d8
 Hmax=Vmax+EKmax
 Hmin=Vmin+EKmin
 ascale=2.0_d8/(Hmax-Hmin)
 bscale=-1.0_d8-ascale*Hmin
 dt=ascale/dsqrt(1.0_d8-((Vmin+this%esk)*Ascale-Bscale)**2)
 nstep=ceiling(this%t_total/dt)
 do i=1,this%n_dvr
  vpot(i)=vpot(i)*ascale+bscale
  t_eng(i)=t_eng(i)*ascale
 end do 
 !init energy
 ne=(this%e_end-this%e_start)/this%de+1
 allocate(esk0(ne),esk0s(ne))
 this%e_start=this%e_start/ev2hartree
 this%e_end=this%e_end/ev2hartree
 this%de=this%de/ev2hartree
 do i=1,ne
  esk0(i)=this%e_start+(i-1)*this%de
  esk0s(i)=-dacos((esk0(i)+Vmin)*ascale+bscale)
 end do 
 !flux sin and cos 
 allocate(fsin(this%n_dvr),fcos(this%n_dvr))
 do i=1,this%n_dvr
 fsin(i)=dsqrt(2.0_d8/length)*dsin(i*pi*(this%r_flux-this%r_start)/length)
 fcos(i)=dsqrt(2.0_d8/length)*(i*pi/length)*dcos(i*pi*(this%r_flux-this%r_start)/length)
 end do
 !absorb pot
 allocate(fabs1(this%n_dvr),fabs2(this%n_dvr))
 !absorb pot
 call wfabs(this%pabs,this%n_dvr,rq0,this%cabs1,this%rabs1,this%cabs2,this%rabs2, &
       nabsl,fabs1,nabsr,fabs2,dt)
 !init
 allocate(w0(this%n_dvr),w1(this%n_dvr),w2(this%n_dvr),ae0(ne),aef(ne))
 w2=0.0_d8
 call init2(this%n_dvr,rq0,this%mass,this%r_flux,this%esk,this%delta,w2,ne,esk0,aef)
 w0=w2
 w2=0.0_d8
 call init2(this%n_dvr,rq0,this%mass,this%r0,this%esk,this%delta,w2,ne,esk0,ae0)
 w1=w2

 allocate(phie(ne),phie0(ne),phie0p(ne))
 phie0=(0.0_d8,0.0_d8)
 phie=(0.0_d8,0.0_d8)
 phie0p=(0.0_d8,0.0_d8)

 allocate(cwk1(this%n_dvr),cwk2(this%n_dvr))
 call flux_b(this%n_dvr,w1,cwk1,fsin,fcos,ne,esk0s,phie0,phie0p,trans,0)
 call flux1_b(this%n_dvr,w1,fsin,fcos,ne,esk0s,phie,w0,0)
 
 do i=1,this%n_dvr
  cwk1(i)=w1(i)*vpot(i)
 end do
 call fttrans(this%n_dvr,trans,w1,cwk2)
 do i=1,this%n_dvr
 cwk2(i)=cwk2(i)*t_eng(i)
 end do
 call fttrans(this%n_dvr,trans,cwk2,w2)
 do i=1,this%n_dvr
 w2(i)=w2(i)+cwk1(i)
 end do
 
 en=0
 do i=1,this%n_dvr
 en=en+w2(i)*w1(i)
 end do

 call flux_b(this%n_dvr,w2,cwk1,fsin,fcos,ne,esk0s,phie0,phie0p,trans,1)
 call flux1_b(this%n_dvr,w2,fsin,fcos,ne,esk0s,phie,w0,1)

 factid=1.0_d8/this%mass
 idabsr=0
 anansr=0
 
 allocate(cwk3(this%n_dvr))

 do i=2,nstep
  call prog_crwp(this%n_dvr,w1,w2,cwk1,cwk2,cwk3,vpot,t_eng,trans,dt,nabsl,fabs1,&
                idabsr,anansr,nabsr,fabs2,ascale,bscale)
  call flux_b(this%n_dvr,w2,cwk1,fsin,fcos,ne,esk0s,phie0,phie0p,trans,i)
  call flux1_b(this%n_dvr,w2,fsin,fcos,ne,esk0s,phie,w0,i)
 end do
 
 open(222,file='probability_sincdvr_crwp',status='unknown') 
 do i=1,ne
 s=phie(i)/(2*pi*ae0(i)*conjg(aef(i)))
 s=s*ascale/sqrt(1.0_d8-dcos(esk0s(i))**2)
 probs=conjg(s)*s
 write(222,*) esk0(i)*ev2hartree,probs 
 end do
 close(222)


 end subroutine


 subroutine init(nq1,rq1,mm,rr0,et,del,w_2,n_e,et0,ae)
  implicit none
  integer, intent(in) :: nq1, n_e
  real(8), intent(in) :: rq1(nq1), mm, rr0, et, del
  real(8), intent(in) :: et0(n_e)
  complex(8), intent(out) :: w_2(nq1)
  complex(8), intent(out) :: ae(n_e)
  complex(8) :: w1(n_e),wi
  real(8)::ak,dr,pi
  integer :: i,j
  real(8)::aki,fac0
  pi=dacos(-1.0_d8)
  wi=(0.0_d8,1.0_d8)
  w1=(0.0_d8,0.0_d8)
  ak=dsqrt(2.0_d8*mm*et)
  dr=rq1(2)-rq1(1)
  fac0=sqrt(2.0*mm/pi)
  do i=1,nq1
   w_2(i)=(1.0_d8/(pi*del**2)**(0.25))*dsqrt(dr)*dexp(-((rq1(i)-rr0)/del)**2/2.0_d8) &
         *cdexp(-wi*ak*rq1(i))
   do j=1,n_e
   aki=sqrt(2.0*mm*et0(j))
  
   w1(j)=w1(j)+fac0/dsqrt(aki)*w_2(i) &
        *cdexp(wi*aki*rq1(i))/2.0_d8*dsqrt(dr)
   end do
  end do
  
  do i=1,n_e
   ae(i)=w1(i)
  end do
end subroutine 

subroutine init2(nq1, rq1, rms, r0, esk, delta, w2, ne0, esk0, ae0)
  implicit none
  integer, intent(in) :: nq1, ne0
  real(8), intent(in) :: rq1(nq1), rms, r0, esk, delta
  real(8), intent(in) :: esk0(ne0)
  complex(8), intent(out) :: w2(nq1), ae0(ne0)
  integer :: i, ne
  real(8) :: ak, aki, dr, fac, fac0
  complex(8) :: wi
  real(8), parameter :: pi = 3.141592653589793238d0

  wi = (0.d0, 1.d0)
  ak = sqrt(2.d0 * rms * abs(esk))

  ae0 = (0.d0, 0.d0)
  dr = rq1(2) - rq1(1)
  fac = 1.d0 / (pi * delta**2)**0.25 * sqrt(2.d0 * dr)
  fac0 = sqrt(2.d0 * rms / pi)

  do i = 1, nq1
     w2(i) = fac * exp(-((rq1(i)-r0)/delta)**2 / 2.d0) * cos(ak*rq1(i))
     do ne = 1, ne0
        aki = sqrt(2.d0 * rms * esk0(ne))
        ae0(ne) = ae0(ne) + fac0 / sqrt(aki) * w2(i) * exp(wi*aki*rq1(i)) / 2.d0 * sqrt(dr)
     end do
  end do

  do i = 1, ne0
  end do

end subroutine init2


subroutine transf(nd,tr,w_1,w_2)
implicit none
integer, intent(in)::nd
complex(8),intent(out)::w_2(nd)
complex(8),intent(in)::w_1(nd)
real(8),intent(in)::tr(nd,nd)
integer::i,j
  w_2 = (0.0d0,0.0d0)
  do i = 1, nd
     do j = 1, nd
        w_2(i) = w_2(i) + tr(j, i) * w_1(j)
     end do
  end do
end subroutine

subroutine sincdivi(nd,rq,tdiv)
  implicit none
  integer, intent(in) :: nd
  real(8), intent(in) :: rq(nd)
  real(8), intent(out) :: tdiv(nd, nd)
  integer :: i, j
  real(8) :: dr
  dr = rq(2) - rq(1)
  tdiv = 0.0_d8
  do i = 1, nd
     do j = 1, i-1
        tdiv(i, j) = (-1.0_d8)**(i-j) / ((i-j) * dr)
        tdiv(j, i) = -tdiv(i, j)
     end do
  end do
end subroutine

subroutine fttrans(nz, tmat, w1, w2)
    use, intrinsic :: iso_fortran_env, only: real64
    implicit none
    integer, intent(in) :: nz
    real(real64), intent(in)  :: tmat(nz, nz)
    complex(real64), intent(in)  :: w1(nz)
    complex(real64), intent(out) :: w2(nz)
    integer :: i, j

    w2 = (0.0_real64, 0.0_real64)

    do i = 1, nz
        do j = 1, nz
            w2(i) = w2(i) + tmat(j, i) * w1(j)
        end do
    end do
end subroutine 

subroutine wfabs(pabs, nq1, rq1, cabs1, rabsl, cabs2, rabsr, &
                 nabsl, fabs1, nabsr, fabs2, dt)
  implicit none
  integer, intent(in) :: nq1
  real(8), intent(in) :: pabs, cabs1, rabsl, cabs2, rabsr, dt
  real(8), intent(in) :: rq1(nq1)
  integer, intent(out) :: nabsl, nabsr
  real(8), intent(out) :: fabs1(nq1), fabs2(nq1)
  integer :: i
  real(8) :: r1, r2, rabs1, rabs

  fabs1 = 1.0d0
  fabs2 = 1.0d0

  r1 = rq1(1) - (rq1(2) - rq1(1))
  rabs1 = rabsl + r1

  do i = 1, nq1
     if (rq1(i) > rabs1) exit
     fabs1(i) = dexp(-cabs1 * ((rabs1 - rq1(i)) / rabsl)**pabs * dt)
  end do
  nabsl = i

  r2 = rq1(nq1) + (rq1(2) - rq1(1))
  rabs = r2 - rabsr

  do i = nabsl, nq1
     if (rq1(i) >= rabs) then
        fabs2(i) = dexp(-cabs2 * ((rq1(i) - rabs)/(r2 - rabs))**pabs * dt)
     else
        nabsr = i
     end if
  end do

end subroutine

subroutine prog_so(nq1, w1, w2, vpot, eint, tmat, dt, nabsL, fabsL, nabsR, fabsR)
    implicit none
    integer, intent(in) :: nq1, nabsL, nabsR
    real(8), intent(in) :: vpot(nq1), eint(nq1), dt
    real(8), intent(in) :: tmat(nq1, nq1)
    real(8), intent(in) :: fabsL(nq1), fabsR(nq1)
    complex(8), intent(inout) :: w1(nq1), w2(nq1)

    complex(8) :: wi
    integer :: i,j

    wi = (0.0_8, 1.0_8)

    do i = 1,nq1
        w1(i) = w1(i) * cdexp(-wi * dt / 2.0_8 * eint(i))
    end do

    do i = 1, nq1
     w2(i) = (0.0_8, 0.0_8)
    end do
    do i=1,nq1
     do j = 1, nq1
      w2(i) = w2(i) + tmat(j, i) * w1(j)
     end do
    end do

    w2(1:nq1) = w2(1:nq1) * cdexp(-wi * dt * vpot(1:nq1))
    w2(1:nabsL) = w2(1:nabsL) * fabsL(1:nabsL)
    w2(nabsR:nq1) = w2(nabsR:nq1) * fabsR(nabsR:nq1)

    do i = 1, nq1
     w1(i) = (0.0_8, 0.0_8)
    end do
    do i=1,nq1
     do j = 1, nq1
      w1(i) = w1(i) + tmat(j, i) * w2(j)
     end do
    end do

    do i = 1,nq1
        w1(i) = w1(i) * exp(-wi * dt / 2.0_8 * eint(i))
    end do

end subroutine 

subroutine prog_crwp(nq1, w1, w2, cwk1, cwk2, cwk3, vpot, eint, tmat, dt, &
                nabsL, fabsL, idabsR, anabsR, nabsR, fabsR, ascale, bscale)
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none

    integer, intent(in) :: nq1, nabsL, idabsR, anabsR, nabsR
    real(dp), intent(in) :: dt, ascale, bscale
    real(dp), intent(in) :: vpot(nq1), eint(nq1), fabsL(nabsL), fabsR(nq1), tmat(nq1,*)
    complex(dp), intent(inout) :: w1(nq1), w2(nq1), cwk1(nq1), cwk2(nq1), cwk3(nq1)

    integer :: i
    real(dp) :: tmp

    do i = 1, nq1
        cwk2(i) = w2(i)
    end do

    do i = 1, nq1
        cwk1(i) = w2(i) * vpot(i)
    end do

    call fttrans(nq1, tmat, w2,cwk2)

    do i = 1, nq1
        cwk2(i) = cwk2(i) * eint(i)
    end do

    call fttrans(nq1, tmat, cwk2,cwk3)

    do i = 1, nq1
        cwk3(i) = cwk3(i) + cwk1(i)
    end do

    do i = 1, nabsL
        cwk3(i) = cwk3(i) * fabsL(i)
        w1(i)   = w1(i) * fabsL(i) * fabsL(i)
    end do

    do i = nabsR, nq1
        cwk3(i) = cwk3(i) * fabsR(i)
        w1(i)   = w1(i) * fabsR(i) * fabsR(i)
    end do

    do i = 1, nq1
        tmp   = w2(i)
        w2(i) = 2.0_dp * cwk3(i) - w1(i)
        w1(i) = tmp
    end do

end subroutine 

subroutine flux1(nz, wa1, fsin, fcos, ne0, esk0, phie0, w0, ttot)
    implicit none
    integer, intent(in) :: nz, ne0
    real(kind=8), intent(in) :: esk0(ne0), ttot
    complex(kind=8), intent(in) :: wa1(nz), w0(nz)
    complex(kind=8), intent(inout) :: phie0(ne0)
    real(kind=8), intent(in) :: fsin(nz), fcos(nz)
    complex(kind=8) :: c1, ctmp
    integer :: i, ne

    c1 = (0.0d0, 0.0d0)
    do i = 1, nz
        c1 = c1 + wa1(i) * conjg(w0(i))
    end do

    do ne = 1, ne0
        ctmp= exp((0.0d0,1.0d0)*esk0(ne)*ttot)
        phie0(ne) = phie0(ne) + c1 * ctmp
    end do

end subroutine 

subroutine flux1_b(nz, wa1, fsin, fcos, ne0, esk0, phie0, w0, nn)
    implicit none
    integer, intent(in) :: nz, ne0,nn
    real(kind=8), intent(in) :: esk0(ne0)
    complex(kind=8), intent(in) :: wa1(nz), w0(nz)
    complex(kind=8), intent(inout) :: phie0(ne0)
    real(kind=8), intent(in) :: fsin(nz), fcos(nz)
    complex(kind=8) :: c1, ctmp
    integer :: i, ne

    c1 = (0.0d0, 0.0d0)
    do i = 1, nz
        c1 = c1 + wa1(i) * conjg(w0(i))
    end do

    if(nn/=0) then
    c1=c1*2.0d0
    end if

    do ne = 1, ne0
        ctmp= exp((0.0d0,1.0d0)*esk0(ne)*nn)
        phie0(ne) = phie0(ne) + c1 * ctmp
    end do

end subroutine 

subroutine flux(nz, wa1, fsin, fcos, ne0, esk0, phie0, phie0p, ttot)
    implicit none
    integer, intent(in) :: nz, ne0
    real(kind=8), intent(in) :: esk0(ne0), ttot
    complex(kind=8), intent(in) :: wa1(nz)
    complex(kind=8), intent(inout) :: phie0(ne0), phie0p(ne0)
    real(kind=8), intent(in) :: fsin(nz), fcos(nz)
    complex(kind=8) :: c1, cd1, ctmp
    integer :: i, ne

    c1 = (0.0d0, 0.0d0)
    cd1 = (0.0d0, 0.0d0)
    do i = 1, nz
        c1 = c1 + wa1(i) * fsin(i)
        cd1 = cd1 + wa1(i) * fcos(i)
    end do

    do ne = 1, ne0
        ctmp =  exp((0.0d0,1.0d0)*esk0(ne)*ttot)
        phie0(ne) = phie0(ne) + c1 * ctmp
        phie0p(ne) = phie0p(ne) + cd1 * ctmp
    end do

end subroutine 

subroutine flux_b(nz, wa1, cwk1,fsin, fcos, ne0, esk0, phie0, phie0p, trans,nn)
    implicit none
    integer, intent(in) :: nz, ne0,nn
    real(kind=8), intent(in) :: esk0(ne0),trans(nz,nz)
    complex(kind=8), intent(in) :: wa1(nz)
    complex(kind=8), intent(inout) :: phie0(ne0), phie0p(ne0),cwk1(nz)
    real(kind=8), intent(in) :: fsin(nz), fcos(nz)
    complex(kind=8) :: c1, cd1, ctmp
    integer :: i, ne

    c1 = (0.0d0, 0.0d0)
    cd1 = (0.0d0, 0.0d0)
    cwk1=wa1
    do i = 1, nz
        c1 = c1 + wa1(i) * fsin(i)
        cd1 = cd1 + wa1(i) * fcos(i)
    end do
     
    if(nn/=0) then
    c1=c1*2.0d0
    cd1=cd1*2.0d0
    end if

    do ne = 1, ne0
        ctmp =  exp((0.0d0,1.0d0)*esk0(ne)*nn)
        phie0(ne) = phie0(ne) + c1 * ctmp
        phie0p(ne) = phie0p(ne) + cd1 * ctmp
    end do

end subroutine 

end module 

