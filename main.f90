  program oned  
    use scatter_module
    use, intrinsic:: iso_fortran_env, only : d8 => real64
    implicit none
     
    type(onedscatter):: H2scatter
     integer :: n_dvr
     real(kind=d8)::r_start,r_end,mass,r0,esk,delta
     real(kind=d8)::t_total,dt,r_flux,e_start,e_end,de 
     real(kind=d8)::pabs,cabs1,cabs2,rabs1,rabs2
     n_dvr=512
     r_start=-30
     r_end=30
     mass=0.5
     r0=15
     esk=1.5
     delta=0.3
     t_total=8000
     dt=10
     r_flux=-15
     e_start=0.1
     e_end=3.0
     de=0.01
     pabs=2
     cabs1=0.04
     cabs2=0.02
     rabs1=12
     rabs2=12   
     
    call H2scatter%propagation_scatter_so(n_dvr,r_start,r_end,mass,r0,esk,delta,t_total,dt,r_flux,e_start,e_end,de,pabs,cabs1,cabs2,rabs1,rabs2,potential_eckart)
 
    call H2scatter%propagation_scatter_crwp(n_dvr,r_start,r_end,mass,r0,esk,delta,t_total,r_flux,e_start,e_end,de,pabs,cabs1,cabs2,rabs1,rabs2,potential_eckart)

contains
    
    subroutine potential_eckart(x, V)
        use, intrinsic :: iso_fortran_env, only : d8 => real64
        real(kind=d8), intent(in)  :: x
        real(kind=d8), intent(inout) :: V
        real(kind=d8) :: V0,a0,R0
        V0 = 0.5_d8/27.2116_d8    
        a0 = 0.5_d8    
        R0 = 0.0_d8      
        V  = V0 * 4.0_d8*exp(-(x-R0)/2/a0)/(1.0_d8+exp(-(x-R0)/a0))**2
    end subroutine 
 end program
