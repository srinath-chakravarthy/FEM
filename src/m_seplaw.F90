#define ALP_PC

module seplaw

  contains
  
  ! Get traction for a single integration point that
  ! forms gap and vgap
    subroutine Seplaw_1_Tract(props, gap, vgap, dtime, tract)
      ! ---- Xu - Needleman Cohesive law
      ! ---- Look at Gao and Bower -- MSMSE 2004 for definition or Xu and Needleman -- JMPS -2003
      ! ---- Parameters to be input as props are
      ! ---- Props(1) ---> sigma_max --- Maximum normal stress of cohesive law
      ! ---- Props(2) ---> delta_n --- Maximum normal separation (max gap(2)) of cohesive law
      ! ---- Props(3) ----> delta_t --- Maximum tangential separation (max gap(1))
      ! ---- Props(4) -----> q value --- relates max normal to tangential energy
      ! ---- Props(5) -----> r value --- dictates motion of cohesive surfaces
      ! ---- Props(6) ----> zeta = viscosity --- from Gao and Bower paper 
    implicit none
    ! Constants
    integer, parameter  :: PDIM=2
    
    ! Params
    real(8), intent(in) :: props(:), gap(PDIM), vgap(PDIM), dtime
    real(8), intent(out) :: tract(PDIM)
    
    ! Conveniences
    real(8) :: sepwrk, sigma_max, dn, dt, q, r, zeta, c1, c2
    
    sepwrk = dexp(1.d0)*props(1)*props(2)
    dn = props(2)
    dt = props(3)
    q = props(4)
    r = props(5)
    zeta = props(6)
    
    C1 = ( 1.D0-DEXP(-GAP(2)*GAP(2)/(DT*DT)) )
    C1 = C1 * (1.D0-Q)/(R-1.D0) * (R - GAP(1)/DN)
    C2 = (GAP(1)/DN)*DEXP(-GAP(2)*GAP(2)/(DT*DT))
    TRACT(1) = (SEPWRK/DN)*DEXP(-GAP(1)/DN)*(C2+C1)

    C1 = Q + (R-Q)/(R-1.D0)*(GAP(1)/DN)
    C1 = C1 * DEXP(-GAP(1)/DN) * DEXP(-GAP(2)*GAP(2)/(DT*DT))
    C1 = C1 * 2.D0*(DN/DT)*(SEPWRK/DN)
    TRACT(2) = C1*GAP(2)/DT
    tract(1)=tract(1)+zeta*props(1)*vgap(1)/dn
  end subroutine Seplaw_1_Tract
  
  ! Get stiffness for a single integration point that
  ! forms gap and vgap
  subroutine Seplaw_1_Stiff(props, gap, vgap, dtime, stiff)
    implicit none
    ! Constants
    integer, parameter  :: PDIM=2
    
    ! Params
    real(8), intent(in) :: props(:), gap(PDIM), vgap(PDIM), dtime
    real(8), intent(out) :: stiff(PDIM,PDIM)
    
    ! Conveniences
    real(8) :: c1,sep_work, sigma_max, dn, dt, q, r, zeta
    
    sep_work = dexp(1.d0)*props(1)*props(2)
    dn = props(2)
    dt = props(3)
    q = props(4)
    r = props(5)
    zeta = props(6)

    C1 = (1.D0-Q)/(R-1)*(1.D0-DEXP(-GAP(2)*GAP(2)/(DT*DT)))
    C1 = C1 * (R+1.D0-GAP(1)/DN)
    C1 = (1.D0-GAP(1)/DN)*DEXP(-GAP(2)*GAP(2)/(DT*DT)) - C1
    stiff(1,1) = (sep_work/(DN*DN))*DEXP(-GAP(1)/DN) * C1

    C1 = Q + (GAP(1)/DN) * (R-Q)/(R-1.D0)
    C1 = C1 * DEXP(-GAP(1)/DN)*DEXP(-GAP(2)*GAP(2)/(DT*DT))
    C1 = 2.D0*(sep_work/(DT*DT)) * C1
    stiff(2,2) = C1*(1.D0-2.D0*(GAP(2)*GAP(2)/(DT*DT)))

    C1 = -GAP(1)/DN + (1.D0-Q)/(R-1.D0)*(R-GAP(1)/DN)
    C1 = C1 * DEXP(-GAP(1)/DN)*DEXP(-GAP(2)*GAP(2)/(DT*DT))
    C1 = 2.D0*(sep_work/(DT*DN))*C1
    stiff(1,2) = (GAP(2)/DT)*C1
    stiff(2,1) = stiff(1,2)
     
    stiff(1,1) = stiff(1,1) + zeta*props(1)/dn/dtime
  
  end subroutine Seplaw_1_Stiff

end module seplaw
