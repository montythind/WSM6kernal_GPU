
MODULE module_mp_wsm6
!
   IMPLICIT NONE
!
   REAL, PARAMETER, PUBLIC :: dtcldcr     = 120. ! maximum time step for minor loops
   REAL, PARAMETER, PUBLIC :: n0r = 8.e6         ! intercept parameter rain
   REAL, PARAMETER, PUBLIC :: n0g = 4.e6         ! intercept parameter graupel
   REAL, PARAMETER, PUBLIC :: avtr = 841.9       ! a constant for terminal velocity of rain
   REAL, PARAMETER, PUBLIC :: bvtr = 0.8         ! a constant for terminal velocity of rain
   REAL, PARAMETER, PUBLIC :: r0 = .8e-5         ! 8 microm  in contrast to 10 micro m
   REAL, PARAMETER, PUBLIC :: peaut = .55        ! collection efficiency
   REAL, PARAMETER, PUBLIC :: xncr = 3.e8        ! maritime cloud in contrast to 3.e8 in tc80
   REAL, PARAMETER, PUBLIC :: xmyu = 1.718e-5    ! the dynamic viscosity kgm-1s-1
   REAL, PARAMETER, PUBLIC :: avts = 11.72       ! a constant for terminal velocity of snow
   REAL, PARAMETER, PUBLIC :: bvts = .41         ! a constant for terminal velocity of snow
   REAL, PARAMETER, PUBLIC :: avtg = 330.        ! a constant for terminal velocity of graupel
   REAL, PARAMETER, PUBLIC :: bvtg = 0.8         ! a constant for terminal velocity of graupel
   REAL, PARAMETER, PUBLIC :: deng = 500.        ! density of graupel
   REAL, PARAMETER, PUBLIC :: n0smax =  1.e11    ! maximum n0s (t=-90C unlimited)
   REAL, PARAMETER, PUBLIC :: lamdarmax = 8.e4   ! limited maximum value for slope parameter of rain
   REAL, PARAMETER, PUBLIC :: lamdasmax = 1.e5   ! limited maximum value for slope parameter of snow
   REAL, PARAMETER, PUBLIC :: lamdagmax = 6.e4   ! limited maximum value for slope parameter of graupel
   REAL, PARAMETER, PUBLIC :: dicon = 11.9       ! constant for the cloud-ice diamter
   REAL, PARAMETER, PUBLIC :: dimax = 500.e-6    ! limited maximum value for the cloud-ice diamter
   REAL, PARAMETER, PUBLIC :: n0s = 2.e6         ! temperature dependent intercept parameter snow
   REAL, PARAMETER, PUBLIC :: alpha = .12        ! .122 exponen factor for n0s
   REAL, PARAMETER, PUBLIC :: pfrz1 = 100.       ! constant in Biggs freezing
   REAL, PARAMETER, PUBLIC :: pfrz2 = 0.66       ! constant in Biggs freezing
   REAL, PARAMETER, PUBLIC :: qcrmin = 1.e-9     ! minimun values for qr, qs, and qg
   REAL, PARAMETER, PUBLIC :: eacrc = 1.0        ! Snow/cloud-water collection efficiency
   REAL, PARAMETER, PUBLIC :: dens  =  100.0     ! Density of snow
   REAL, PARAMETER, PUBLIC :: qs0   =  6.e-4     ! threshold amount for aggretion to occur

   REAL, SAVE ::                                      &
             qc0, qck1,bvtr1,bvtr2,bvtr3,bvtr4,g1pbr, &
             g3pbr,g4pbr,g5pbro2,pvtr,eacrr,pacrr,    &
             bvtr6,g6pbr,                             &
             precr1,precr2,roqimax,bvts1,             &
             bvts2,bvts3,bvts4,g1pbs,g3pbs,g4pbs,     &
             g5pbso2,pvts,pacrs,precs1,precs2,pidn0r, &
             pidn0s,xlv1,pacrc,pi,                    &
             bvtg1,bvtg2,bvtg3,bvtg4,g1pbg,           &
             g3pbg,g4pbg,g5pbgo2,pvtg,pacrg,          &
             precg1,precg2,pidn0g,                    &
             rslopermax,rslopesmax,rslopegmax,        &
             rsloperbmax,rslopesbmax,rslopegbmax,     &
             rsloper2max,rslopes2max,rslopeg2max,     &
             rsloper3max,rslopes3max,rslopeg3max

CONTAINS
!===================================================================
!
!===================================================================
!
!ACC$ROUTINE(FORTRAN)
  SUBROUTINE copyConstantsToGPU(pvtr,pvts,pvtg, &
             rslopermax,rslopesmax,rslopegmax, &
             rsloperbmax,rslopesbmax,rslopegbmax,pidn0r,     &
             rsloper2max,rslopes2max,rslopeg2max,pidn0s,     &
             rsloper3max,rslopes3max,rslopeg3max,pidn0g,    &
             qc0,qck1,pacrr,g6pbr,precr1,precr2,roqimax,     &
             precs1,precs2,xlv1,pacrc,pi,pacrg,precg1,precg2)

real,intent(IN) :: pvtr,pvts,pvtg, &
             rslopermax,rslopesmax,rslopegmax, &
             rsloperbmax,rslopesbmax,rslopegbmax,pidn0r,     &
             rsloper2max,rslopes2max,rslopeg2max,pidn0s,     &
             rsloper3max,rslopes3max,rslopeg3max,pidn0g,     &
             qc0,qck1,pacrr,g6pbr,precr1,precr2,roqimax,     &
             precs1,precs2,xlv1,pacrc,pi,pacrg,precg1,precg2

!ACC$DATA(<alpha,qcrmin,n0s,n0smax,pvtr,pvts,pvtg:in,constant>, &
!ACC$> <rslopermax,rsloperbmax,rsloper2max,rsloper3max,bvtr,pidn0r, &
!ACC$>  rslopesmax,rslopesbmax,rslopes2max,rslopes3max,bvts,pidn0s, &
!ACC$>  rslopegmax,rslopegbmax,rslopeg2max,rslopeg3max,bvtg,pidn0g:in,constant>)

!ACC$DATA(<qc0,qck1,pacrr,g6pbr,pvtr,precr1,precr2,roqimax,pvts, &
!ACC$> precs1,precs2,xlv1,pacrc,pi,pvtg,pacrg,precg1,precg2:in,global>)
  return
  END SUBROUTINE copyConstantsToGPU

!ACC$ROUTINE(FORTRAN)
  SUBROUTINE wsm62D(t, q                                          &   
                   ,qci, qrs, den, p, delz                        &
                   ,delt,g, cpd, cpv, rd, rv, t0c                 &
                   ,ep1, ep2, qmin                                &
                   ,XLS, XLV0, XLF0, den0, denr                   &
                   ,cliq,cice,psat                                &
                   ,lat                                           &
                   ,rain,rainncv                                  &
                   ,sr                                            &
                   ,ids,ide, jds,jde, kds,kde                     &
                   ,ims,ime, jms,jme, kms,kme                     &
                   ,its,ite, jts,jte, kts,kte                     &
                   ,snow,snowncv                                  &
                   ,graupel,graupelncv                            &
                                                                  )
!-------------------------------------------------------------------
  IMPLICIT NONE
!-------------------------------------------------------------------
!
!  This code is a 6-class GRAUPEL phase microphyiscs scheme (WSM6) of the 
!  Single-Moment MicroPhyiscs (WSMMP). The WSMMP assumes that ice nuclei
!  number concentration is a function of temperature, and seperate assumption
!  is developed, in which ice crystal number concentration is a function
!  of ice amount. A theoretical background of the ice-microphysics and related
!  processes in the WSMMPs are described in Hong et al. (2004).
!  All production terms in the WSM6 scheme are described in Hong and Lim (2006).
!  All units are in m.k.s. and source/sink terms in kgkg-1s-1.
!
!  WSM6 cloud scheme
!
!  Coded by Song-You Hong and Jeong-Ock Jade Lim (Yonsei Univ.)
!           Summer 2003
!
!  Implemented by Song-You Hong (Yonsei Univ.) and Jimy Dudhia (NCAR)
!           Summer 2004
!
!  History :  semi-lagrangian scheme sedimentation(JH), and clean up
!             Hong, August 2009
!
!  Reference) Hong, Dudhia, Chen (HDC, 2004) Mon. Wea. Rev.
!             Hong and Lim (HL, 2006) J. Korean Meteor. Soc.
!             Dudhia, Hong and Lim (DHL, 2008) J. Meteor. Soc. Japan
!             Lin, Farley, Orville (LFO, 1983) J. Appl. Meteor.
!             Rutledge, Hobbs (RH83, 1983) J. Atmos. Sci.
!             Rutledge, Hobbs (RH84, 1984) J. Atmos. Sci.
!             Juang and Hong (JH, 2010) Mon. Wea. Rev.
!
  INTEGER,      INTENT(IN   )    ::   ids,ide, jds,jde, kds,kde , &
                                      ims,ime, jms,jme, kms,kme , &
                                      its,ite, jts,jte, kts,kte,  &
                                      lat
  REAL, DIMENSION( its:ite , kts:kte ),                           &
        INTENT(INOUT) ::                                          &
                                                               t
  REAL, DIMENSION( its:ite , kts:kte, 2 ),                        &
        INTENT(INOUT) ::                                          &
                                                             qci
  REAL, DIMENSION( its:ite , kts:kte, 3 ),                        &
        INTENT(INOUT) ::                                          &
                                                             qrs
  REAL, DIMENSION( ims:ime , kms:kme ),                           &
        INTENT(INOUT) ::                                          &
                                                               q
  REAL, DIMENSION( ims:ime , kms:kme ),                           &
        INTENT(IN   ) ::                                          &
                                                             den, &
                                                               p, &
                                                            delz
  REAL, INTENT(IN   ) ::                                    delt, &
                                                               g, &
                                                             cpd, &
                                                             cpv, &
                                                             t0c, &
                                                            den0, &
                                                              rd, &
                                                              rv, &
                                                             ep1, &
                                                             ep2, &
                                                            qmin, &
                                                             XLS, &
                                                            XLV0, &
                                                            XLF0, &
                                                            cliq, &
                                                            cice, &
                                                            psat, &
                                                            denr
  REAL, DIMENSION( ims:ime ),                                     &
        INTENT(INOUT) ::                                    rain, &
                                                         rainncv, &
                                                              sr
  REAL, DIMENSION( ims:ime, jms:jme ), OPTIONAL,                  &
        INTENT(INOUT) ::                                    snow, &
                                                         snowncv
  REAL, DIMENSION( ims:ime, jms:jme ), OPTIONAL,                  &
        INTENT(INOUT) ::                                 graupel, &
                                                      graupelncv
! LOCAL VAR
  REAL, DIMENSION( its:ite , kts:kte , 3) ::                      &
                                                              rh, &
                                                              qs, &
                                                          rslope, &
                                                         rslope2, &
                                                         rslope3, &
                                                         rslopeb, &
                                                         qrs_tmp, & 
                                                            falk, &
                                                            fall, &
                                                           work1
  REAL, DIMENSION( its:ite , kts:kte ) ::                         &
                                                           fallc, &
                                                           falkc, &
                                                          work1c, &
                                                          work2c, &
                                                           workr, &
                                                           worka
  REAL, DIMENSION( its:ite , kts:kte ) ::                         &
                                                         den_tmp, &
                                                        delz_tmp
  REAL, DIMENSION( its:ite , kts:kte ) ::                         &
                                                           pigen, &
                                                           pidep, &
                                                           pcond, &
                                                           prevp, &
                                                           psevp, &
                                                           pgevp, &
                                                           psdep, &
                                                           pgdep, &
                                                           praut, &
                                                           psaut, &
                                                           pgaut, &
                                                           piacr, &
                                                           pracw, &
                                                           praci, &
                                                           pracs, &
                                                           psacw, &
                                                           psaci, &
                                                           psacr, &
                                                           pgacw, &
                                                           pgaci, &
                                                           pgacr, &
                                                           pgacs, &
                                                           paacw, &
                                                           psmlt, &
                                                           pgmlt, &
                                                           pseml, &
                                                           pgeml
  REAL, DIMENSION( its:ite , kts:kte ) ::                         &
                                                            qsum, &
                                                              xl, &
                                                             cpm, &
                                                           work2, &
                                                          denfac, &
                                                             xni, &
                                                         denqrs1, &
                                                         denqrs2, &
                                                         denqrs3, &
                                                          denqci, & 
                                                          delta2, &
                                                          delta3, &
                                                           n0sfac
  REAL, DIMENSION( its:ite ) ::                          delqrs1, &
                                                         delqrs2, &
                                                         delqrs3, &
                                                           delqi  
  REAL, DIMENSION( its:ite ) ::                        tstepsnow, &
                                                      tstepgraup
  INTEGER, DIMENSION( its:ite ) ::                         mstep, &
                                                           numdt
  LOGICAL, DIMENSION( its:ite ) ::                        flgcld
  REAL  ::                                                        &
            cpmcal, xlcal, diffus,                                &
            viscos, xka, venfac, conden, diffac,                  &
            x, y, z, a, b, c, d, e,                               &
            qdt, holdrr, holdrs, holdrg, supcol, supcolt, pvt,    &
            coeres, supsat, dtcld, xmi, eacrs, satdt,             &
            qimax, diameter, xni0, roqi0,                         &
            fallsum, fallsum_qsi, fallsum_qg,                     &
            vt2i,vt2r,vt2s,vt2g,acrfac,egs,egi,                   &
            xlwork2, factor, source, value,                       &
            xlf, pfrzdtc, pfrzdtr, supice, alpha2
  REAL  :: vt2ave
  REAL  :: holdc, holdci
  INTEGER :: i, j, k, mstepmax,                                   &
            iprt, latd, lond, loop, loops, ifsat, n, idim, kdim
  INTEGER :: itest,ktest
! Temporaries used for inlining fpvs function
  REAL  :: dldti, xb, xai, tr, xbi, xa, hvap, cvap, hsub, dldt, ttp
! variables for optimization
  REAL, DIMENSION( its:ite ) ::                             tvec1
  REAL                       ::                              temp
! MWG: iintermediate variables for F2C translation
  REAL                       :: diffus1,viscos1,xka1

!ACC$DATA(<t,qci,qrs,den,delz,rain,rainncv,snow,snowncv,graupel,graupelncv,sr,p,q:in,global>)
!ACC$DATA(<rh,qs,mstep,rslope,rslope2,rslope3,rslopeb,tstepsnow,tstepgraup:in,global>)
!ACC$DATA(<delqrs1,delqrs2,delqrs3,delqi,denqrs1,denqrs2,denqrs3,denqci:in,global>)

!=================================================================
!   compute internal functions
!
#ifndef INLINE
      cpmcal(x) = cpd*(1.-max(x,qmin))+max(x,qmin)*cpv
      xlcal(x) = xlv0-xlv1*(x-t0c)
!----------------------------------------------------------------
!     diffus: diffusion coefficient of the water vapor
!     viscos: kinematic viscosity(m2s-1)
!     Optimizatin : A**B => exp(log(A)*(B))
!
      diffus(x,y) = 8.794e-5 * exp(log(x)*(1.81)) / y        ! 8.794e-5*x**1.81/y
      viscos(x,y) = 1.496e-6 * (x*sqrt(x)) /(x+120.)/y  ! 1.496e-6*x**1.5/(x+120.)/y
      xka(x,y) = 1.414e3*viscos(x,y)*y
      diffac(a,b,c,d,e) = d*a*a/(xka(c,d)*rv*c*c)+1./(e*diffus(c,b))
      venfac(a,b,c) = exp(log((viscos(b,c)/diffus(b,a)))*((.3333333)))         &
                     /sqrt(viscos(b,c))*sqrt(sqrt(den0/c))
      conden(a,b,c,d,e) = (max(b,qmin)-c)/(1.+d*d/(rv*e)*c/(a*a))
#endif
!
#ifdef ALIGN_OK
!DIR$ ASSUME_ALIGNED t:64,qci:64,qrs:64,q:64,den:64,p:64,delz:64,rain:64,rainncv:64,sr:64
!DIR$ ASSUME_ALIGNED snow:64,snowncv:64,graupel:64,graupelncv:64,mstep:64,flgcld:64
!DIR$ ASSUME_ALIGNED worka:64,pgacs:64,n0sfac:64,work2:64,psmlt:64,rslope:64
!DIR$ ASSUME_ALIGNED rslope2:64,rslope3:64,rslopeb:64,cpm:64,pgmlt:64
!DIR$ ASSUME_ALIGNED tstepsnow:64,tstepgraup:64
#endif
!
      idim = ite-its+1
      kdim = kte-kts+1
itest=979
ktest=1
!
!----------------------------------------------------------------
!     padding 0 for negative values generated by dynamics
!
!ACC$REGION(<128:chunk>,<(ite-its+1)/128+1,(kte-kts+1)>,<lat:in>, &
!ACC$> <xlv1:none,extern>, &
!ACC$> <t,qci,qrs,den,delz,rainncv,snowncv,graupelncv,sr,q:none,extern>, &
!ACC$> <tstepgraup,tstepsnow:none,extern>) BEGIN
!$acc parallel vector_length(128)
!ACC$DO PARALLEL(2,kts:kte)
      do k = kts, kte
!$acc loop gang vector
!ACC$DO PARALLEL(1,its:ite)
        do i = its, ite
          qci(i,k,1) = max(qci(i,k,1),0.0)
          qrs(i,k,1) = max(qrs(i,k,1),0.0)
          qci(i,k,2) = max(qci(i,k,2),0.0)
          qrs(i,k,2) = max(qrs(i,k,2),0.0)
          qrs(i,k,3) = max(qrs(i,k,3),0.0)
        enddo
      enddo
!
!----------------------------------------------------------------
!     latent heat for phase changes and heat capacity. neglect the
!     changes during microphysical process calculation
!     emanuel(1994)
!
!ACC$DO PARALLEL(2,kts:kte)
      do k = kts, kte
!ACC$DO PARALLEL(1,its:ite)
!$acc loop gang vector
        do i = its, ite
#ifdef INLINE
          cpm(i,k) = cpd*(1.-max(q(i,k),qmin))+max(q(i,k),qmin)*cpv
          xl(i,k) = xlv0-xlv1*(t(i,k)-t0c)
#else
          cpm(i,k) = cpmcal(q(i,k))
          xl(i,k) = xlcal(t(i,k))
#endif
        enddo
      enddo
!ACC$DO PARALLEL(2,kts:kte)
      do k = kts, kte
!ACC$DO PARALLEL(1,its:ite)
!$acc loop gang vector
        do i = its, ite
          delz_tmp(i,k) = delz(i,k)
          den_tmp(i,k) = den(i,k)
        enddo
      enddo
!
!----------------------------------------------------------------
!    initialize the surface rain, snow, graupel
!
! MWG "k" loop added for F2C parallelization to avoid multi-block execution on
! the GPU
!ACC$INSERTC if (blockIdx.y == 0) {
!ACC$DO PARALLEL(1,its:ite)
!$acc loop gang vector
      do i = its, ite
        rainncv(i) = 0.
#ifndef F2C
        if(PRESENT (snowncv) .AND. PRESENT (snow)) snowncv(i,lat) = 0.
        if(PRESENT (graupelncv) .AND. PRESENT (graupel)) graupelncv(i,lat) = 0.
#else
        snowncv(i,lat) = 0.
        graupelncv(i,lat) = 0.
#endif
        sr(i) = 0.
! new local array to catch step snow and graupel
        tstepsnow(i) = 0.
        tstepgraup(i) = 0.
      enddo
!ACC$INSERTC }
!$acc end parallel 
!ACC$REGION END
!
!----------------------------------------------------------------
!     compute the minor time steps.
!
      loops = max(nint(delt/dtcldcr),1)
      dtcld = delt/loops
      if(delt.le.dtcldcr) dtcld = delt
!
      do loop = 1,loops
!
!----------------------------------------------------------------
!     initialize the large scale variables
!
#ifdef ALIGN_OK
!DIR$ VECTOR ALIGNED
#endif
!ACC$REGION(<160:chunk>,<(ite-its+1)/160+1,kte-kts+1>,<its,ite,kts,kte:in>, &
!ACC$> <t,qci,qrs,den,p,rh,qs,q,mstep,rslope,rslopeb,rslope2,rslope3:none,extern>, &
!ACC$> <qrs,den,denqrs1,denqrs2,denqrs3:none,extern>) BEGIN
!$acc parallel vector_length(128)
!ACC$DO PARALLEL(1,its:ite)
!ACC$INSERTC if (blockIdx.y == 0) {
!$acc loop gang vector
      do i = its, ite
        mstep(i) = 1
        flgcld(i) = .true.
      enddo
!ACC$INSERTC }

!ACC$DO PARALLEL(2)
      do k = kts, kte
!ACC$DO PARALLEL(1,its:ite)
!$acc loop gang vector
        do i = its, ite
          denfac(i,k) = sqrt(den0/den(i,k))
        enddo
      enddo
!
! Inline expansion for fpvs
!         qs(i,k,1) = fpvs(t(i,k),0,rd,rv,cpv,cliq,cice,xlv0,xls,psat,t0c)
!         qs(i,k,2) = fpvs(t(i,k),1,rd,rv,cpv,cliq,cice,xlv0,xls,psat,t0c)
      hsub = xls
      hvap = xlv0
      cvap = cpv
      ttp=t0c+0.01
      dldt=cvap-cliq
      xa=-dldt/rv
      xb=xa+hvap/(rv*ttp)
      dldti=cvap-cice
      xai=-dldti/rv
      xbi=xai+hsub/(rv*ttp)
!ACC$DO PARALLEL(2)
      do k = kts, kte
!ACC$DO PARALLEL(1,its:ite)
!$acc loop gang vector
        do i = its, ite
          tr=ttp/t(i,k)
          qs(i,k,1)=psat*exp(log(tr)*(xa))*exp(xb*(1.-tr))
          qs(i,k,1) = min(qs(i,k,1),0.99*p(i,k))
          qs(i,k,1) = ep2 * qs(i,k,1) / (p(i,k) - qs(i,k,1))
          qs(i,k,1) = max(qs(i,k,1),qmin)
          rh(i,k,1) = max(q(i,k) / qs(i,k,1),qmin)
! MWG: when run on the CPU, this division operation leads to different results
! between Intel Fortran CPU and the NVCC compiler run on the CPU
!          rh(i,k,1) = q(i,k)/qs(i,k,1)
          tr=ttp/t(i,k)
          if(t(i,k).lt.ttp) then
            qs(i,k,2)=psat*exp(log(tr)*(xai))*exp(xbi*(1.-tr))
          else
            qs(i,k,2)=psat*exp(log(tr)*(xa))*exp(xb*(1.-tr))
          endif
          qs(i,k,2) = min(qs(i,k,2),0.99*p(i,k))
          qs(i,k,2) = ep2 * qs(i,k,2) / (p(i,k) - qs(i,k,2))
          qs(i,k,2) = max(qs(i,k,2),qmin)
          rh(i,k,2) = max(q(i,k) / qs(i,k,2),qmin)
! MWG: when run on the CPU, this division operation leads to different results
! between Intel Fortran CPU and the NVCC compiler run on the CPU
!          rh(i,k,2) = q(i,k)/qs(i,k,2)
        enddo
      enddo
!
!----------------------------------------------------------------
!     initialize the variables for microphysical physics
!
!
!ACC$DO PARALLEL(2)
      do k = kts, kte
!ACC$DO PARALLEL(1,its:ite)
!$acc loop gang vector
        do i = its, ite
          prevp(i,k) = 0.
          psdep(i,k) = 0.
          pgdep(i,k) = 0.
          praut(i,k) = 0.
          psaut(i,k) = 0.
          pgaut(i,k) = 0.
          pracw(i,k) = 0.
          praci(i,k) = 0.
          piacr(i,k) = 0.
          psaci(i,k) = 0.
          psacw(i,k) = 0.
          pracs(i,k) = 0.
          psacr(i,k) = 0.
          pgacw(i,k) = 0.
          paacw(i,k) = 0.
          pgaci(i,k) = 0.
          pgacr(i,k) = 0.
          pgacs(i,k) = 0.
          pigen(i,k) = 0.
          pidep(i,k) = 0.
          pcond(i,k) = 0.
          psmlt(i,k) = 0.
          pgmlt(i,k) = 0.
          pseml(i,k) = 0.
          pgeml(i,k) = 0.
          psevp(i,k) = 0.
          pgevp(i,k) = 0.
          falk(i,k,1) = 0.
          falk(i,k,2) = 0.
          falk(i,k,3) = 0.
          fall(i,k,1) = 0.
          fall(i,k,2) = 0.
          fall(i,k,3) = 0.
          fallc(i,k) = 0.
          falkc(i,k) = 0.
          xni(i,k) = 1.e3
        enddo
      enddo
!-------------------------------------------------------------
! Ni: ice crystal number concentraiton   [HDC 5c]
!-------------------------------------------------------------
!ACC$DO PARALLEL(2)
      do k = kts, kte
!ACC$DO PARALLEL(1,its:ite)
!$acc loop gang vector
        do i = its, ite
          temp = (den(i,k)*max(qci(i,k,2),qmin))
          temp = sqrt(sqrt(temp*temp*temp))
          xni(i,k) = min(max(5.38e7*temp,1.e3),1.e6)
        enddo
      enddo
!
!----------------------------------------------------------------
!     compute the fallout term:
!     first, vertical terminal velosity for minor loops
!----------------------------------------------------------------
!ACC$DO PARALLEL(2)
      do k = kts, kte
!ACC$DO PARALLEL(1,its:ite)
!$acc loop gang vector
        do i = its, ite
          qrs_tmp(i,k,1) = qrs(i,k,1)
          qrs_tmp(i,k,2) = qrs(i,k,2)
          qrs_tmp(i,k,3) = qrs(i,k,3)
        enddo
      enddo
!$acc end parallel
      call slope_wsm6(qrs_tmp,den_tmp,denfac,t,rslope,rslopeb,rslope2,rslope3, & 
                     work1,its,ite,kts,kte)
!
!$acc parallel vector_length(128)
!ACC$DO PARALLEL(2)
      do k = kts, kte
!MWG      do k = kte, kts, -1  |F2C translation error
#ifdef ALIGN_OK
!DIR$ VECTOR ALIGNED
#endif
!ACC$DO PARALLEL(1,its:ite)
!$acc loop gang vector
        do i = its, ite
          workr(i,k) = work1(i,k,1)
          qsum(i,k) = max( (qrs(i,k,2)+qrs(i,k,3)), 1.E-15)
          IF ( qsum(i,k) .gt. 1.e-15 ) THEN
            worka(i,k) = (work1(i,k,2)*qrs(i,k,2) + work1(i,k,3)*qrs(i,k,3)) &
                      /qsum(i,k)
          ELSE
            worka(i,k) = 0.
          ENDIF
          denqrs1(i,k) = den(i,k)*qrs(i,k,1)
          denqrs2(i,k) = den(i,k)*qrs(i,k,2)
          denqrs3(i,k) = den(i,k)*qrs(i,k,3)
          if(qrs(i,k,1).le.0.0) workr(i,k) = 0.0
        enddo
      enddo
!$acc end parallel
!ACC$REGION END

!ACC$DATA(<t:out,extern>)
      call nislfv_rain_plm(its,ite,kts,kte,den_tmp,denfac,t,delz_tmp,workr,denqrs1,  &
                           delqrs1,dtcld,1,1)
      call nislfv_rain_plm6(its,ite,kts,kte,den_tmp,denfac,t,delz_tmp,worka,         & 
                           denqrs2,denqrs3,delqrs2,delqrs3,dtcld,1,1)
!ACC$DATA(<t:in,extern>)

!ACC$REGION(<128:chunk>,<(ite-its+1)/128+1,(kte-kts+1)>,<its,ite,kts,kte:in>, &
!ACC$>      <precs1,precs2,pi,precg1,precg2:none,extern>, &
!ACC$>      <t,qci,qrs,den,delz,delqrs1,delqrs2,delqrs3,denqrs1,denqrs2,denqrs3:none,extern>, &
!ACC$>      <p,mstep,denqci,rslope,rslopeb,rslope2,rslope3:none,extern>) BEGIN
!$acc parallel vector_length(128)
!ACC$DO PARALLEL(2)
      do k = kts, kte
!ACC$DO PARALLEL(1,its:ite)
!$acc loop gang vector
        do i = its, ite
          qrs(i,k,1) = max(denqrs1(i,k)/den(i,k),0.)
          qrs(i,k,2) = max(denqrs2(i,k)/den(i,k),0.)
          qrs(i,k,3) = max(denqrs3(i,k)/den(i,k),0.)
          fall(i,k,1) = denqrs1(i,k)*workr(i,k)/delz(i,k)
          fall(i,k,2) = denqrs2(i,k)*worka(i,k)/delz(i,k)
          fall(i,k,3) = denqrs3(i,k)*worka(i,k)/delz(i,k)
        enddo
      enddo
!ACC$INSERTC if (blockIdx.y == 0) {
!ACC$DO PARALLEL(1,its:ite)
!$acc loop gang vector
      do i = its, ite
        fall(i,1,1) = delqrs1(i)/delz(i,1)/dtcld
        fall(i,1,2) = delqrs2(i)/delz(i,1)/dtcld
        fall(i,1,3) = delqrs3(i)/delz(i,1)/dtcld
      enddo
!ACC$INSERTC }
!ACC$DO PARALLEL(2)
      do k = kts, kte
!ACC$DO PARALLEL(1,its:ite)
!$acc loop gang vector
        do i = its, ite
          qrs_tmp(i,k,1) = qrs(i,k,1)
          qrs_tmp(i,k,2) = qrs(i,k,2)
          qrs_tmp(i,k,3) = qrs(i,k,3)
        enddo
      enddo
!$acc end parallel
      call slope_wsm6(qrs_tmp,den_tmp,denfac,t,rslope,rslopeb,rslope2,rslope3, &
                     work1,its,ite,kts,kte)
!
!$acc parallel vector_length(128)
!MWG      do k = kte, kts, -1  |F2C translation error
!ACC$DO PARALLEL(2)
      do k = kts, kte
#ifdef ALIGN_OK
!DIR$ VECTOR ALIGNED
#endif
!ACC$DO PARALLEL(1,its:ite)
!$acc loop gang vector
        do i = its, ite
          supcol = t0c-t(i,k)
          n0sfac(i,k) = max(min(exp(alpha*supcol),n0smax/n0s),1.)
          if(t(i,k).gt.t0c) then
!---------------------------------------------------------------
! psmlt: melting of snow [HL A33] [RH83 A25]
!       (T>T0: S->R)
!---------------------------------------------------------------
            xlf = xlf0
#ifdef INLINE
            diffus1 = 8.794e-5 * exp(log(t(i,k))*(1.81))/p(i,k) 
            viscos1 = 1.496e-6 * (t(i,k)*sqrt(t(i,k))) /(t(i,k)+120.)/den(i,k)
            xka1 = 1.414e3*1.496e-6 * (t(i,k)*sqrt(t(i,k))) /(t(i,k)+120.)/den(i,k)*den(i,k)

            work2(i,k) = exp(log((viscos1/diffus1))*((.3333333)))         &
                     /sqrt(viscos1)*sqrt(sqrt(den0/den(i,k)))
#else
            work2(i,k) = venfac(p(i,k),t(i,k),den(i,k))
#endif
            if(qrs(i,k,2).gt.0.) then
              coeres = rslope2(i,k,2)*sqrt(rslope(i,k,2)*rslopeb(i,k,2))
#ifdef INLINE
              psmlt(i,k) = 1.414e3*1.496e-6 * (t(i,k)*sqrt(t(i,k)))          &
               /(t(i,k)+120.)/den(i,k)*den(i,k)/xlf*(t0c-t(i,k))*pi/2.       &
                         *n0sfac(i,k)*(precs1*rslope2(i,k,2)                 &
                         +precs2*work2(i,k)*coeres)
#else
              psmlt(i,k) = xka(t(i,k),den(i,k))/xlf*(t0c-t(i,k))*pi/2.       &
                         *n0sfac(i,k)*(precs1*rslope2(i,k,2)                 &
                         +precs2*work2(i,k)*coeres)
#endif
              psmlt(i,k) = min(max(psmlt(i,k)*dtcld/mstep(i),                &
                          -qrs(i,k,2)/mstep(i)),0.)
              qrs(i,k,2) = qrs(i,k,2) + psmlt(i,k)
              qrs(i,k,1) = qrs(i,k,1) - psmlt(i,k)
              t(i,k) = t(i,k) + xlf/cpm(i,k)*psmlt(i,k)
            endif
!---------------------------------------------------------------
! pgmlt: melting of graupel [HL A23]  [LFO 47]
!       (T>T0: G->R)
!---------------------------------------------------------------
            if(qrs(i,k,3).gt.0.) then
              coeres = rslope2(i,k,3)*sqrt(rslope(i,k,3)*rslopeb(i,k,3))
#ifdef INLINE
              pgmlt(i,k) = 1.414e3*1.496e-6 * (t(i,k)*sqrt(t(i,k)))          &
                         /(t(i,k)+120.)/den(i,k)*den(i,k)/xlf                &
                         *(t0c-t(i,k))*(precg1*rslope2(i,k,3)                &
                         +precg2*work2(i,k)*coeres)
#else
              pgmlt(i,k) = xka(t(i,k),den(i,k))/xlf                          &
                         *(t0c-t(i,k))*(precg1*rslope2(i,k,3)                &
                         +precg2*work2(i,k)*coeres)
#endif
              pgmlt(i,k) = min(max(pgmlt(i,k)*dtcld/mstep(i),                &
                          -qrs(i,k,3)/mstep(i)),0.)                          
              qrs(i,k,3) = qrs(i,k,3) + pgmlt(i,k)
              qrs(i,k,1) = qrs(i,k,1) - pgmlt(i,k)
              t(i,k) = t(i,k) + xlf/cpm(i,k)*pgmlt(i,k)
            endif
          endif
        enddo
      enddo
!---------------------------------------------------------------
! Vice [ms-1] : fallout of ice crystal [HDC 5a]
!---------------------------------------------------------------
!MWG      do k = kte, kts, -1  |F2C translation error
!ACC$DO PARALLEL(2)
      do k = kts, kte
!ACC$DO PARALLEL(1,its:ite)
!$acc loop gang vector
        do i = its, ite
          if(qci(i,k,2).le.0.) then
            work1c(i,k) = 0.
          else
            xmi = den(i,k)*qci(i,k,2)/xni(i,k)
            diameter  = max(min(dicon * sqrt(xmi),dimax), 1.e-25)
            work1c(i,k) = 1.49e4*exp(log(diameter)*(1.31))
          endif
        enddo
      enddo
!
!  forward semi-laglangian scheme (JH), PCM (piecewise constant),  (linear)
!
!MWG      do k = kte, kts, -1  |F2C translation error
!ACC$DO PARALLEL(2)
      do k = kts, kte
!ACC$DO PARALLEL(1,its:ite)
!$acc loop gang vector
        do i = its, ite
          denqci(i,k) = den(i,k)*qci(i,k,2)
        enddo
      enddo
!$acc end parallel
!ACC$REGION END

!ACC$DATA(<t:out,extern>)
      call nislfv_rain_plm(its,ite,kts,kte,den_tmp,denfac,t,delz_tmp,work1c,denqci,  &
                           delqi,dtcld,1,0)
!ACC$DATA(<t:in,extern>)

!ACC$REGION(<128:chunk>,<(ite-its+1)/128+1,(kte-kts+1)>,<its,ite,kts,kte:in>, &
!ACC$> <qc0,qck1,pacrr,g6pbr,pvtr,precr1,precr2,roqimax,pvts,precs1,precs2,pacrc,pi,pvtg,pacrg,precg1,precg2:none,extern>, &
!ACC$> <t,qci,qrs,den,delz,p,q,rh,qs,rain,rainncv,snow,snowncv,graupel,graupelncv,sr:none,extern>, &
!ACC$> <rslope,rslopeb,rslope2,rslope3,tstepgraup,tstepsnow,delqi,denqci:none,extern>) BEGIN
!$acc parallel vector_length(128)
!ACC$DO PARALLEL(2)
      do k = kts, kte
!ACC$DO PARALLEL(1,its:ite)
!$acc loop gang vector
        do i = its, ite
          qci(i,k,2) = max(denqci(i,k)/den(i,k),0.)
        enddo
      enddo
!ACC$INSERTC if (blockIdx.y == 0) {
!ACC$DO PARALLEL(1,its:ite)
!$acc loop gang vector
      do i = its, ite
        fallc(i,1) = delqi(i)/delz(i,1)/dtcld
      enddo
!
!----------------------------------------------------------------
!      rain (unit is mm/sec;kgm-2s-1: /1000*delt ===> m)==> mm for wrf
!
#ifdef ALIGN_OK
!DIR$ VECTOR ALIGNED
#endif
!ACC$DO PARALLEL(1,its:ite)
!$acc loop gang vector
      do i = its, ite
        fallsum = fall(i,kts,1)+fall(i,kts,2)+fall(i,kts,3)+fallc(i,kts)
        fallsum_qsi = fall(i,kts,2)+fallc(i,kts)
        fallsum_qg = fall(i,kts,3)
        if(fallsum.gt.0.) then
          rainncv(i) = fallsum*delz(i,kts)/denr*dtcld*1000. + rainncv(i)
          rain(i) = fallsum*delz(i,kts)/denr*dtcld*1000. + rain(i)
        endif
        if(fallsum_qsi.gt.0.) then
          tstepsnow(i)   = fallsum_qsi*delz(i,kts)/denr*dtcld*1000.            &
                           +tstepsnow(i)
#ifndef F2C
        IF ( PRESENT (snowncv) .AND. PRESENT (snow)) THEN
#endif
          snowncv(i,lat) = fallsum_qsi*delz(i,kts)/denr*dtcld*1000.            & 
                           +snowncv(i,lat)
          snow(i,lat) = fallsum_qsi*delz(i,kts)/denr*dtcld*1000. + snow(i,lat)
#ifndef F2C
        ENDIF
#endif
        endif
        if(fallsum_qg.gt.0.) then
          tstepgraup(i)  = fallsum_qsi*delz(i,kts)/denr*dtcld*1000.            &
                           +tstepgraup(i)
#ifndef F2C
        IF ( PRESENT (graupelncv) .AND. PRESENT (graupel)) THEN
#endif
          graupelncv(i,lat) = fallsum_qg*delz(i,kts)/denr*dtcld*1000.          &   
                              + graupelncv(i,lat)
          graupel(i,lat) = fallsum_qg*delz(i,kts)/denr*dtcld*1000. + graupel(i,lat)
#ifndef F2C
        ENDIF
#endif
        endif
!       if(fallsum.gt.0.)sr(i)=(snowncv(i,lat) + graupelncv(i,lat))/(rainncv(i)+1.e-12)
        if(fallsum.gt.0.)sr(i)=(tstepsnow(i) + tstepgraup(i))/(rainncv(i)+1.e-12)
      enddo
!ACC$INSERTC }
!
!---------------------------------------------------------------
! pimlt: instantaneous melting of cloud ice [HL A47] [RH83 A28]
!       (T>T0: I->C)
!---------------------------------------------------------------
!ACC$DO PARALLEL(2)
      do k = kts, kte
!ACC$DO PARALLEL(1,its:ite)
!$acc loop gang vector
        do i = its, ite
          supcol = t0c-t(i,k)
          xlf = xls-xl(i,k)
          if(supcol.lt.0.) xlf = xlf0
          if(supcol.lt.0.and.qci(i,k,2).gt.0.) then
            qci(i,k,1) = qci(i,k,1) + qci(i,k,2)
            t(i,k) = t(i,k) - xlf/cpm(i,k)*qci(i,k,2)
            qci(i,k,2) = 0.
          endif
!---------------------------------------------------------------
! pihmf: homogeneous freezing of cloud water below -40c [HL A45]
!        (T<-40C: C->I)
!---------------------------------------------------------------
          if(supcol.gt.40..and.qci(i,k,1).gt.0.) then
            qci(i,k,2) = qci(i,k,2) + qci(i,k,1)
            t(i,k) = t(i,k) + xlf/cpm(i,k)*qci(i,k,1)
            qci(i,k,1) = 0.
          endif
!---------------------------------------------------------------
! pihtf: heterogeneous freezing of cloud water [HL A44]
!        (T0>T>-40C: C->I)
!---------------------------------------------------------------
          if(supcol.gt.0..and.qci(i,k,1).gt.qmin) then
!           pfrzdtc = min(pfrz1*(exp(pfrz2*supcol)-1.)                         &
!              *den(i,k)/denr/xncr*qci(i,k,1)**2*dtcld,qci(i,k,1))
            supcolt=min(supcol,50.)
            pfrzdtc = min(pfrz1*(exp(pfrz2*supcolt)-1.)                        &
            *den(i,k)/denr/xncr*qci(i,k,1)*qci(i,k,1)*dtcld,qci(i,k,1))
            qci(i,k,2) = qci(i,k,2) + pfrzdtc
            t(i,k) = t(i,k) + xlf/cpm(i,k)*pfrzdtc
            qci(i,k,1) = qci(i,k,1)-pfrzdtc
          endif
!---------------------------------------------------------------
! pgfrz: freezing of rain water [HL A20] [LFO 45]
!        (T<T0, R->G)
!---------------------------------------------------------------
          if(supcol.gt.0..and.qrs(i,k,1).gt.0.) then
!           pfrzdtr = min(20.*pi**2*pfrz1*n0r*denr/den(i,k)                    &
!                 *(exp(pfrz2*supcol)-1.)*rslope3(i,k,1)**2                    &
!                 *rslope(i,k,1)*dtcld,qrs(i,k,1))
            temp = rslope3(i,k,1)
            temp = temp*temp*rslope(i,k,1)
            supcolt=min(supcol,50.)
            pfrzdtr = min(20.*(pi*pi)*pfrz1*n0r*denr/den(i,k)                  &
                  *(exp(pfrz2*supcolt)-1.)*temp*dtcld,                         &
                  qrs(i,k,1))
            qrs(i,k,3) = qrs(i,k,3) + pfrzdtr
            t(i,k) = t(i,k) + xlf/cpm(i,k)*pfrzdtr
            qrs(i,k,1) = qrs(i,k,1)-pfrzdtr
          endif
        enddo
      enddo
!
!
!----------------------------------------------------------------
!     update the slope parameters for microphysics computation
!
!ACC$DO PARALLEL(2)
      do k = kts, kte
!ACC$DO PARALLEL(1,its:ite)
!$acc loop gang vector
        do i = its, ite
          qrs_tmp(i,k,1) = qrs(i,k,1)
          qrs_tmp(i,k,2) = qrs(i,k,2)
          qrs_tmp(i,k,3) = qrs(i,k,3)
        enddo
      enddo
!$acc end parallel
      call slope_wsm6(qrs_tmp,den_tmp,denfac,t,rslope,rslopeb,rslope2,rslope3, &
                     work1,its,ite,kts,kte)
!------------------------------------------------------------------
!     work1:  the thermodynamic term in the denominator associated with
!             heat conduction and vapor diffusion
!             (ry88, y93, h85)
!     work2: parameter associated with the ventilation effects(y93)
!
!$acc parallel vector_length(128)
!ACC$DO PARALLEL(2)
      do k = kts, kte
!ACC$DO PARALLEL(1,its:ite)
!$acc loop gang vector
        do i = its, ite
#ifdef INLINE
          diffus1 = 8.794e-5 * exp(log(t(i,k))*(1.81))/p(i,k) 
          viscos1 = 1.496e-6 * (t(i,k)*sqrt(t(i,k))) /(t(i,k)+120.)/den(i,k)
          work1(i,k,1) = den(i,k)*xl(i,k)*xl(i,k)/ (1.414e3*1.496e-6 *                &
                             (t(i,k)*sqrt(t(i,k))) /(t(i,k)+120.)/den(i,k)*den(i,k) &
                                *rv*t(i,k)*t(i,k))+1./(qs(i,k,1)*diffus1)
          work1(i,k,2) = den(i,k)*xls*xls        / (1.414e3*1.496e-6 *                &
                             (t(i,k)*sqrt(t(i,k))) /(t(i,k)+120.)/den(i,k)*den(i,k) &
                                *rv*t(i,k)*t(i,k))+1./(qs(i,k,2)*diffus1)
          work2(i,k) = exp(log((viscos1/diffus1))*((.3333333)))         &
                     /sqrt(viscos1)*sqrt(sqrt(den0/den(i,k)))
#else
          work1(i,k,1) = diffac(xl(i,k),p(i,k),t(i,k),den(i,k),qs(i,k,1))
          work1(i,k,2) = diffac(xls,p(i,k),t(i,k),den(i,k),qs(i,k,2))
          work2(i,k) = venfac(p(i,k),t(i,k),den(i,k))
#endif
        enddo
      enddo
!
!===============================================================
!
! warm rain processes
!
! - follows the processes in RH83 and LFO except for autoconcersion
!
!===============================================================
!
!ACC$DO PARALLEL(2)
      do k = kts, kte
#ifdef ALIGN_OK
!DIR$ VECTOR ALIGNED
#endif
!ACC$DO PARALLEL(1,its:ite)
!$acc loop gang vector
        do i = its, ite
          supsat = max(q(i,k),qmin)-qs(i,k,1)
          satdt = supsat/dtcld
!---------------------------------------------------------------
! praut: auto conversion rate from cloud to rain [HDC 16]
!        (C->R)
!---------------------------------------------------------------
          if(qci(i,k,1).gt.qc0) then
            praut(i,k) = qck1*qci(i,k,1)**(7./3.)
            praut(i,k) = min(praut(i,k),qci(i,k,1)/dtcld)
          endif
!---------------------------------------------------------------
! pracw: accretion of cloud water by rain [HL A40] [LFO 51]
!        (C->R)
!---------------------------------------------------------------
          if(qrs(i,k,1).gt.qcrmin.and.qci(i,k,1).gt.qmin) then
            pracw(i,k) = min(pacrr*rslope3(i,k,1)*rslopeb(i,k,1)               &
                        *qci(i,k,1)*denfac(i,k),qci(i,k,1)/dtcld)
          endif
!---------------------------------------------------------------
! prevp: evaporation/condensation rate of rain [HDC 14]
!        (V->R or R->V)
!---------------------------------------------------------------
          if(qrs(i,k,1).gt.0.) then
            coeres = rslope2(i,k,1)*sqrt(rslope(i,k,1)*rslopeb(i,k,1))
            prevp(i,k) = (rh(i,k,1)-1.)*(precr1*rslope2(i,k,1)                 &
                         +precr2*work2(i,k)*coeres)/work1(i,k,1)
            if(prevp(i,k).lt.0.) then
              prevp(i,k) = max(prevp(i,k),-qrs(i,k,1)/dtcld)
              prevp(i,k) = max(prevp(i,k),satdt/2)
            else
              prevp(i,k) = min(prevp(i,k),satdt/2)
            endif
          endif
        enddo
      enddo
!
!===============================================================
!
! cold rain processes
!
! - follows the revised ice microphysics processes in HDC
! - the processes same as in RH83 and RH84  and LFO behave
!   following ice crystal hapits defined in HDC, inclduing
!   intercept parameter for snow (n0s), ice crystal number
!   concentration (ni), ice nuclei number concentration
!   (n0i), ice diameter (d)
!
!===============================================================
!
!ACC$DO PARALLEL(2)
      do k = kts, kte
#ifdef ALIGN_OK
!DIR$ VECTOR ALIGNED
#endif
!ACC$DO PARALLEL(1,its:ite)
!$acc loop gang vector
        do i = its, ite
          supcol = t0c-t(i,k)
          n0sfac(i,k) = max(min(exp(alpha*supcol),n0smax/n0s),1.)
          supsat = max(q(i,k),qmin)-qs(i,k,2)
          satdt = supsat/dtcld
          ifsat = 0
!-------------------------------------------------------------
! Ni: ice crystal number concentraiton   [HDC 5c]
!-------------------------------------------------------------
!         xni(i,k) = min(max(5.38e7*(den(i,k)                                  &
!                      *max(qci(i,k,2),qmin))**0.75,1.e3),1.e6)
          temp = (den(i,k)*max(qci(i,k,2),qmin))
          temp = sqrt(sqrt(temp*temp*temp))
          xni(i,k) = min(max(5.38e7*temp,1.e3),1.e6)
          eacrs = exp(0.07*(-supcol))
!
          xmi = den(i,k)*qci(i,k,2)/xni(i,k)
          diameter  = min(dicon * sqrt(xmi),dimax)
          vt2i = 1.49e4*diameter**1.31
          vt2r=pvtr*rslopeb(i,k,1)*denfac(i,k)
          vt2s=pvts*rslopeb(i,k,2)*denfac(i,k)
          vt2g=pvtg*rslopeb(i,k,3)*denfac(i,k)
          qsum(i,k) = max( (qrs(i,k,2)+qrs(i,k,3)), 1.E-15)
          if(qsum(i,k) .gt. 1.e-15) then
          vt2ave=(vt2s*qrs(i,k,2)+vt2g*qrs(i,k,3))/(qsum(i,k))
          else
          vt2ave=0.
          endif
          if(supcol.gt.0.and.qci(i,k,2).gt.qmin) then
            if(qrs(i,k,1).gt.qcrmin) then
!-------------------------------------------------------------
! praci: Accretion of cloud ice by rain [HL A15] [LFO 25]
!        (T<T0: I->R)
!-------------------------------------------------------------
              acrfac = 2.*rslope3(i,k,1)+2.*diameter*rslope2(i,k,1)            &
                      +diameter**2*rslope(i,k,1)
              praci(i,k) = pi*qci(i,k,2)*n0r*abs(vt2r-vt2i)*acrfac/4.
              praci(i,k) = min(praci(i,k),qci(i,k,2)/dtcld)
!-------------------------------------------------------------
! piacr: Accretion of rain by cloud ice [HL A19] [LFO 26]
!        (T<T0: R->S or R->G)
!-------------------------------------------------------------
              piacr(i,k) = pi**2*avtr*n0r*denr*xni(i,k)*denfac(i,k)            &
                          *g6pbr*rslope3(i,k,1)*rslope3(i,k,1)                 &
                          *rslopeb(i,k,1)/24./den(i,k)
              piacr(i,k) = min(piacr(i,k),qrs(i,k,1)/dtcld)
            endif
!-------------------------------------------------------------
! psaci: Accretion of cloud ice by snow [HDC 10]
!        (T<T0: I->S)
!-------------------------------------------------------------
            if(qrs(i,k,2).gt.qcrmin) then
              acrfac = 2.*rslope3(i,k,2)+2.*diameter*rslope2(i,k,2)            &
                      +diameter**2*rslope(i,k,2)
              psaci(i,k) = pi*qci(i,k,2)*eacrs*n0s*n0sfac(i,k)                 &
                          *abs(vt2ave-vt2i)*acrfac/4.
              psaci(i,k) = min(psaci(i,k),qci(i,k,2)/dtcld)
            endif
!-------------------------------------------------------------
! pgaci: Accretion of cloud ice by graupel [HL A17] [LFO 41]
!        (T<T0: I->G)
!-------------------------------------------------------------
            if(qrs(i,k,3).gt.qcrmin) then
              egi = exp(0.07*(-supcol))
              acrfac = 2.*rslope3(i,k,3)+2.*diameter*rslope2(i,k,3)            &
                      +diameter**2*rslope(i,k,3)
              pgaci(i,k) = pi*egi*qci(i,k,2)*n0g*abs(vt2ave-vt2i)*acrfac/4.
              pgaci(i,k) = min(pgaci(i,k),qci(i,k,2)/dtcld)
            endif
          endif
!-------------------------------------------------------------
! psacw: Accretion of cloud water by snow  [HL A7] [LFO 24]
!        (T<T0: C->S, and T>=T0: C->R)
!-------------------------------------------------------------
          if(qrs(i,k,2).gt.qcrmin.and.qci(i,k,1).gt.qmin) then
            psacw(i,k) = min(pacrc*n0sfac(i,k)*rslope3(i,k,2)*rslopeb(i,k,2)   &    
                        *qci(i,k,1)*denfac(i,k),qci(i,k,1)/dtcld)
          endif
!-------------------------------------------------------------
! pgacw: Accretion of cloud water by graupel [HL A6] [LFO 40]
!        (T<T0: C->G, and T>=T0: C->R)
!-------------------------------------------------------------
          if(qrs(i,k,3).gt.qcrmin.and.qci(i,k,1).gt.qmin) then
            pgacw(i,k) = min(pacrg*rslope3(i,k,3)*rslopeb(i,k,3)               &
                        *qci(i,k,1)*denfac(i,k),qci(i,k,1)/dtcld)
          endif
!-------------------------------------------------------------
! paacw: Accretion of cloud water by averaged snow/graupel 
!        (T<T0: C->G or S, and T>=T0: C->R) 
!-------------------------------------------------------------
          if(qrs(i,k,2).gt.qcrmin.and.qrs(i,k,3).gt.qcrmin) then
            paacw(i,k) = (qrs(i,k,2)*psacw(i,k)+qrs(i,k,3)*pgacw(i,k))         & 
                        /(qsum(i,k))
           endif      
!-------------------------------------------------------------
! pracs: Accretion of snow by rain [HL A11] [LFO 27]
!         (T<T0: S->G)
!-------------------------------------------------------------
          if(qrs(i,k,2).gt.qcrmin.and.qrs(i,k,1).gt.qcrmin) then
            if(supcol.gt.0) then
              acrfac = 5.*rslope3(i,k,2)*rslope3(i,k,2)*rslope(i,k,1)          &
                      +2.*rslope3(i,k,2)*rslope2(i,k,2)*rslope2(i,k,1)         &
                      +.5*rslope2(i,k,2)*rslope2(i,k,2)*rslope3(i,k,1)
              pracs(i,k) = pi**2*n0r*n0s*n0sfac(i,k)*abs(vt2r-vt2ave)          &
                          *(dens/den(i,k))*acrfac
              pracs(i,k) = min(pracs(i,k),qrs(i,k,2)/dtcld)
            endif
!-------------------------------------------------------------
! psacr: Accretion of rain by snow [HL A10] [LFO 28]
!         (T<T0:R->S or R->G) (T>=T0: enhance melting of snow)
!-------------------------------------------------------------
            acrfac = 5.*rslope3(i,k,1)*rslope3(i,k,1)*rslope(i,k,2)            &
                    +2.*rslope3(i,k,1)*rslope2(i,k,1)*rslope2(i,k,2)           &
                    +.5*rslope2(i,k,1)*rslope2(i,k,1)*rslope3(i,k,2)
            psacr(i,k) = pi**2*n0r*n0s*n0sfac(i,k)*abs(vt2ave-vt2r)            &
                        *(denr/den(i,k))*acrfac
            psacr(i,k) = min(psacr(i,k),qrs(i,k,1)/dtcld)
          endif
!-------------------------------------------------------------
! pgacr: Accretion of rain by graupel [HL A12] [LFO 42]
!         (T<T0: R->G) (T>=T0: enhance melting of graupel)
!-------------------------------------------------------------
          if(qrs(i,k,3).gt.qcrmin.and.qrs(i,k,1).gt.qcrmin) then
            acrfac = 5.*rslope3(i,k,1)*rslope3(i,k,1)*rslope(i,k,3)            &
                    +2.*rslope3(i,k,1)*rslope2(i,k,1)*rslope2(i,k,3)           &
                    +.5*rslope2(i,k,1)*rslope2(i,k,1)*rslope3(i,k,3)
            pgacr(i,k) = pi**2*n0r*n0g*abs(vt2ave-vt2r)*(denr/den(i,k))        &
                        *acrfac
            pgacr(i,k) = min(pgacr(i,k),qrs(i,k,1)/dtcld)
          endif
!
!-------------------------------------------------------------
! pgacs: Accretion of snow by graupel [HL A13] [LFO 29]
!        (S->G): This process is eliminated in V3.0 with the 
!        new combined snow/graupel fall speeds
!-------------------------------------------------------------
          if(qrs(i,k,3).gt.qcrmin.and.qrs(i,k,2).gt.qcrmin) then
            pgacs(i,k) = 0.
          endif
          if(supcol.le.0) then
            xlf = xlf0
!-------------------------------------------------------------
! pseml: Enhanced melting of snow by accretion of water [HL A34]
!        (T>=T0: S->R)
!-------------------------------------------------------------
            if(qrs(i,k,2).gt.0.)                                               &
              pseml(i,k) = min(max(cliq*supcol*(paacw(i,k)+psacr(i,k))         &
                          /xlf,-qrs(i,k,2)/dtcld),0.)
!-------------------------------------------------------------
! pgeml: Enhanced melting of graupel by accretion of water [HL A24] [RH84 A21-A22]
!        (T>=T0: G->R)
!-------------------------------------------------------------
            if(qrs(i,k,3).gt.0.)                                               &
              pgeml(i,k) = min(max(cliq*supcol*(paacw(i,k)+pgacr(i,k))         &
                          /xlf,-qrs(i,k,3)/dtcld),0.)
          endif
          if(supcol.gt.0) then
!-------------------------------------------------------------
! pidep: Deposition/Sublimation rate of ice [HDC 9]
!       (T<T0: V->I or I->V)
!-------------------------------------------------------------
            if(qci(i,k,2).gt.0.and.ifsat.ne.1) then
              pidep(i,k) = 4.*diameter*xni(i,k)*(rh(i,k,2)-1.)/work1(i,k,2)
              supice = satdt-prevp(i,k)
              if(pidep(i,k).lt.0.) then
                pidep(i,k) = max(max(pidep(i,k),satdt/2),supice)
                pidep(i,k) = max(pidep(i,k),-qci(i,k,2)/dtcld)
              else
                pidep(i,k) = min(min(pidep(i,k),satdt/2),supice)
              endif
              if(abs(prevp(i,k)+pidep(i,k)).ge.abs(satdt)) ifsat = 1
            endif
!-------------------------------------------------------------
! psdep: deposition/sublimation rate of snow [HDC 14]
!        (T<T0: V->S or S->V)
!-------------------------------------------------------------
            if(qrs(i,k,2).gt.0..and.ifsat.ne.1) then
              coeres = rslope2(i,k,2)*sqrt(rslope(i,k,2)*rslopeb(i,k,2))
              psdep(i,k) = (rh(i,k,2)-1.)*n0sfac(i,k)*(precs1*rslope2(i,k,2)   &    
                           + precs2*work2(i,k)*coeres)/work1(i,k,2)
              supice = satdt-prevp(i,k)-pidep(i,k)
              if(psdep(i,k).lt.0.) then
                psdep(i,k) = max(psdep(i,k),-qrs(i,k,2)/dtcld)
                psdep(i,k) = max(max(psdep(i,k),satdt/2),supice)
              else
                psdep(i,k) = min(min(psdep(i,k),satdt/2),supice)
              endif
              if(abs(prevp(i,k)+pidep(i,k)+psdep(i,k)).ge.abs(satdt))          &
                ifsat = 1
            endif
!-------------------------------------------------------------
! pgdep: deposition/sublimation rate of graupel [HL A21] [LFO 46]
!        (T<T0: V->G or G->V)
!-------------------------------------------------------------
            if(qrs(i,k,3).gt.0..and.ifsat.ne.1) then
              coeres = rslope2(i,k,3)*sqrt(rslope(i,k,3)*rslopeb(i,k,3))
              pgdep(i,k) = (rh(i,k,2)-1.)*(precg1*rslope2(i,k,3)               &
                              +precg2*work2(i,k)*coeres)/work1(i,k,2)
              supice = satdt-prevp(i,k)-pidep(i,k)-psdep(i,k)
              if(pgdep(i,k).lt.0.) then
                pgdep(i,k) = max(pgdep(i,k),-qrs(i,k,3)/dtcld)
                pgdep(i,k) = max(max(pgdep(i,k),satdt/2),supice)
              else
                pgdep(i,k) = min(min(pgdep(i,k),satdt/2),supice)
              endif
              if(abs(prevp(i,k)+pidep(i,k)+psdep(i,k)+pgdep(i,k)).ge.          &
                abs(satdt)) ifsat = 1
            endif
!-------------------------------------------------------------
! pigen: generation(nucleation) of ice from vapor [HL 50] [HDC 7-8]
!       (T<T0: V->I)
!-------------------------------------------------------------
            if(supsat.gt.0.and.ifsat.ne.1) then
              supice = satdt-prevp(i,k)-pidep(i,k)-psdep(i,k)-pgdep(i,k)
              xni0 = 1.e3*exp(0.1*supcol)
              roqi0 = 4.92e-11*xni0**1.33
              pigen(i,k) = max(0.,(roqi0/den(i,k)-max(qci(i,k,2),0.))/dtcld)
              pigen(i,k) = min(min(pigen(i,k),satdt),supice)
            endif
!
!-------------------------------------------------------------
! psaut: conversion(aggregation) of ice to snow [HDC 12]
!        (T<T0: I->S)
!-------------------------------------------------------------
            if(qci(i,k,2).gt.0.) then
              qimax = roqimax/den(i,k)
              psaut(i,k) = max(0.,(qci(i,k,2)-qimax)/dtcld)
            endif
!
!-------------------------------------------------------------
! pgaut: conversion(aggregation) of snow to graupel [HL A4] [LFO 37]
!        (T<T0: S->G)
!-------------------------------------------------------------
            if(qrs(i,k,2).gt.0.) then
              alpha2 = 1.e-3*exp(0.09*(-supcol))
              pgaut(i,k) = min(max(0.,alpha2*(qrs(i,k,2)-qs0)),qrs(i,k,2)/dtcld)
            endif
          endif
!
!-------------------------------------------------------------
! psevp: Evaporation of melting snow [HL A35] [RH83 A27]
!       (T>=T0: S->V)
!-------------------------------------------------------------
          if(supcol.lt.0.) then
            if(qrs(i,k,2).gt.0..and.rh(i,k,1).lt.1.) then
              coeres = rslope2(i,k,2)*sqrt(rslope(i,k,2)*rslopeb(i,k,2))
              psevp(i,k) = (rh(i,k,1)-1.)*n0sfac(i,k)*(precs1                  &
                           *rslope2(i,k,2)+precs2*work2(i,k)                   &
                           *coeres)/work1(i,k,1)
              psevp(i,k) = min(max(psevp(i,k),-qrs(i,k,2)/dtcld),0.)
            endif
!-------------------------------------------------------------
! pgevp: Evaporation of melting graupel [HL A25] [RH84 A19]
!       (T>=T0: G->V)
!-------------------------------------------------------------
            if(qrs(i,k,3).gt.0..and.rh(i,k,1).lt.1.) then
              coeres = rslope2(i,k,3)*sqrt(rslope(i,k,3)*rslopeb(i,k,3))
              pgevp(i,k) = (rh(i,k,1)-1.)*(precg1*rslope2(i,k,3)               &
                         +precg2*work2(i,k)*coeres)/work1(i,k,1)
              pgevp(i,k) = min(max(pgevp(i,k),-qrs(i,k,3)/dtcld),0.)
            endif
          endif
        enddo
      enddo
!
!
!----------------------------------------------------------------
!     check mass conservation of generation terms and feedback to the
!     large scale
!
!ACC$DO PARALLEL(2)
      do k = kts, kte
#ifdef ALIGN_OK
!DIR$ VECTOR ALIGNED
#endif
!ACC$DO PARALLEL(1,its:ite)
!$acc loop gang vector
        do i = its, ite
          delta2(i,k)=0.
          delta3(i,k)=0.
          if(qrs(i,k,1).lt.1.e-4.and.qrs(i,k,2).lt.1.e-4) delta2(i,k)=1.
          if(qrs(i,k,1).lt.1.e-4) delta3(i,k)=1.
        enddo
      enddo
!ACC$DO PARALLEL(2)
      do k = kts, kte
#ifdef ALIGN_OK
!DIR$ VECTOR ALIGNED
#endif
!ACC$DO PARALLEL(1,its:ite)
!$acc loop gang vector
        do i = its, ite
!
!     cloud water
!
          value = max(qmin,qci(i,k,1))
          source = (praut(i,k)+pracw(i,k)+paacw(i,k)+paacw(i,k))*dtcld
          if(t(i,k).le.t0c) then
            if (source.gt.value) then
              factor = value/source
              praut(i,k) = praut(i,k)*factor
              pracw(i,k) = pracw(i,k)*factor
              paacw(i,k) = paacw(i,k)*factor
            endif
          else
            if (source.gt.value) then
              factor = value/source
              praut(i,k) = praut(i,k)*factor
              pracw(i,k) = pracw(i,k)*factor
              paacw(i,k) = paacw(i,k)*factor
            endif
          endif
        enddo
      enddo
!ACC$DO PARALLEL(2)
      do k = kts, kte
#ifdef ALIGN_OK
!DIR$ VECTOR ALIGNED
#endif
!ACC$DO PARALLEL(1,its:ite)
!$acc loop gang vector
        do i = its, ite
!
!     cloud ice
!
          if(t(i,k).le.t0c) then
            value = max(qmin,qci(i,k,2))
            source = (psaut(i,k)-pigen(i,k)-pidep(i,k)+praci(i,k)+psaci(i,k)   &
                     +pgaci(i,k))*dtcld
            if (source.gt.value) then
              factor = value/source
              psaut(i,k) = psaut(i,k)*factor
              pigen(i,k) = pigen(i,k)*factor
              pidep(i,k) = pidep(i,k)*factor
              praci(i,k) = praci(i,k)*factor
              psaci(i,k) = psaci(i,k)*factor
              pgaci(i,k) = pgaci(i,k)*factor
            endif
          endif
        enddo
      enddo
!ACC$DO PARALLEL(2)
      do k = kts, kte
#ifdef ALIGN_OK
!DIR$ VECTOR ALIGNED
#endif
!ACC$DO PARALLEL(1,its:ite)
!$acc loop gang vector
        do i = its, ite
!
!     rain
!
          value = max(qmin,qrs(i,k,1))
          if(t(i,k).le.t0c) then
            source = (-praut(i,k)-prevp(i,k)-pracw(i,k)+piacr(i,k)+psacr(i,k)  &
                      +pgacr(i,k))*dtcld
            if (source.gt.value) then
              factor = value/source
              praut(i,k) = praut(i,k)*factor
              prevp(i,k) = prevp(i,k)*factor
              pracw(i,k) = pracw(i,k)*factor
              piacr(i,k) = piacr(i,k)*factor
              psacr(i,k) = psacr(i,k)*factor
              pgacr(i,k) = pgacr(i,k)*factor
            endif
          else
            source = (-paacw(i,k)-praut(i,k)+pseml(i,k)+pgeml(i,k)-pracw(i,k)  &
                      -paacw(i,k)-prevp(i,k))*dtcld
            if (source.gt.value) then
              factor = value/source
              praut(i,k) = praut(i,k)*factor
              prevp(i,k) = prevp(i,k)*factor
              pracw(i,k) = pracw(i,k)*factor
              paacw(i,k) = paacw(i,k)*factor
              pseml(i,k) = pseml(i,k)*factor
              pgeml(i,k) = pgeml(i,k)*factor
            endif
          endif
        enddo
      enddo
!ACC$DO PARALLEL(2)
      do k = kts, kte
#ifdef ALIGN_OK
!DIR$ VECTOR ALIGNED
#endif
!ACC$DO PARALLEL(1,its:ite)
!$acc loop gang vector
        do i = its, ite
!
!     snow
!
          if(t(i,k).le.t0c) then
            value = max(qmin,qrs(i,k,2))
            source = -(psdep(i,k)+psaut(i,k)-pgaut(i,k)+paacw(i,k)+piacr(i,k)  &
                     *delta3(i,k)+praci(i,k)*delta3(i,k)                       &
                     -pracs(i,k)*(1.-delta2(i,k))                              &
                     +psacr(i,k)*delta2(i,k)+psaci(i,k)-pgacs(i,k) )*dtcld
            if (source.gt.value) then
              factor = value/source
              psdep(i,k) = psdep(i,k)*factor
              psaut(i,k) = psaut(i,k)*factor
              pgaut(i,k) = pgaut(i,k)*factor
              paacw(i,k) = paacw(i,k)*factor
              piacr(i,k) = piacr(i,k)*factor
              praci(i,k) = praci(i,k)*factor
              psaci(i,k) = psaci(i,k)*factor
              pracs(i,k) = pracs(i,k)*factor
              psacr(i,k) = psacr(i,k)*factor
              pgacs(i,k) = pgacs(i,k)*factor
            endif
          else
            value = max(qcrmin,qrs(i,k,2))
            source=(pgacs(i,k)-pseml(i,k)-psevp(i,k))*dtcld
            if (source.gt.value) then
              factor = value/source
              pgacs(i,k) = pgacs(i,k)*factor
              psevp(i,k) = psevp(i,k)*factor
              pseml(i,k) = pseml(i,k)*factor
            endif
          endif
        enddo
      enddo
!ACC$DO PARALLEL(2)
      do k = kts, kte
#ifdef ALIGN_OK
!DIR$ VECTOR ALIGNED
#endif
!ACC$DO PARALLEL(1,its:ite)
!$acc loop gang vector
        do i = its, ite
!
!     graupel
!
          if(t(i,k).le.t0c) then
            value = max(qmin,qrs(i,k,3))
            source = -(pgdep(i,k)+pgaut(i,k)                                   &
                     +piacr(i,k)*(1.-delta3(i,k))+praci(i,k)*(1.-delta3(i,k))  &
                     +psacr(i,k)*(1.-delta2(i,k))+pracs(i,k)*(1.-delta2(i,k))  &
                     +pgaci(i,k)+paacw(i,k)+pgacr(i,k)+pgacs(i,k))*dtcld
            if (source.gt.value) then
              factor = value/source
              pgdep(i,k) = pgdep(i,k)*factor
              pgaut(i,k) = pgaut(i,k)*factor
              piacr(i,k) = piacr(i,k)*factor
              praci(i,k) = praci(i,k)*factor
              psacr(i,k) = psacr(i,k)*factor
              pracs(i,k) = pracs(i,k)*factor
              paacw(i,k) = paacw(i,k)*factor
              pgaci(i,k) = pgaci(i,k)*factor
              pgacr(i,k) = pgacr(i,k)*factor
              pgacs(i,k) = pgacs(i,k)*factor
            endif
          else
            value = max(qcrmin,qrs(i,k,3))
            source=-(pgacs(i,k)+pgevp(i,k)+pgeml(i,k))*dtcld
            if (source.gt.value) then
              factor = value/source
              pgacs(i,k) = pgacs(i,k)*factor
              pgevp(i,k) = pgevp(i,k)*factor
              pgeml(i,k) = pgeml(i,k)*factor
            endif
          endif
        enddo
      enddo
!ACC$DO PARALLEL(2)
      do k = kts, kte
#ifdef ALIGN_OK
!DIR$ VECTOR ALIGNED
#endif
!ACC$DO PARALLEL(1,its:ite)
!$acc loop gang vector
        do i = its, ite
!
!     update
!
          if(t(i,k).le.t0c) then
            work2(i,k)=-(prevp(i,k)+psdep(i,k)+pgdep(i,k)+pigen(i,k)+pidep(i,k))
            q(i,k) = q(i,k)+work2(i,k)*dtcld
            qci(i,k,1) = max(qci(i,k,1)-(praut(i,k)+pracw(i,k)                 &
                           +paacw(i,k)+paacw(i,k))*dtcld,0.)
            qrs(i,k,1) = max(qrs(i,k,1)+(praut(i,k)+pracw(i,k)                 &
                           +prevp(i,k)-piacr(i,k)-pgacr(i,k)                   &
                           -psacr(i,k))*dtcld,0.)
            qci(i,k,2) = max(qci(i,k,2)-(psaut(i,k)+praci(i,k)                 &
                           +psaci(i,k)+pgaci(i,k)-pigen(i,k)-pidep(i,k))       &
                           *dtcld,0.)
            qrs(i,k,2) = max(qrs(i,k,2)+(psdep(i,k)+psaut(i,k)+paacw(i,k)      &
                           -pgaut(i,k)+piacr(i,k)*delta3(i,k)                  &
                           +praci(i,k)*delta3(i,k)+psaci(i,k)-pgacs(i,k)       &
                           -pracs(i,k)*(1.-delta2(i,k))+psacr(i,k)*delta2(i,k))&
                           *dtcld,0.)
            qrs(i,k,3) = max(qrs(i,k,3)+(pgdep(i,k)+pgaut(i,k)                 &
                           +piacr(i,k)*(1.-delta3(i,k))                        &
                           +praci(i,k)*(1.-delta3(i,k))                        &
                           +psacr(i,k)*(1.-delta2(i,k))                        &
                           +pracs(i,k)*(1.-delta2(i,k))+pgaci(i,k)+paacw(i,k)  &
                           +pgacr(i,k)+pgacs(i,k))*dtcld,0.)
            xlf = xls-xl(i,k)
            xlwork2 = -xls*(psdep(i,k)+pgdep(i,k)+pidep(i,k)+pigen(i,k))       &
                      -xl(i,k)*prevp(i,k)-xlf*(piacr(i,k)+paacw(i,k)           &
                      +paacw(i,k)+pgacr(i,k)+psacr(i,k))
            t(i,k) = t(i,k)-xlwork2/cpm(i,k)*dtcld
          else
            work2(i,k)=-(prevp(i,k)+psevp(i,k)+pgevp(i,k))
            q(i,k) = q(i,k)+work2(i,k)*dtcld
            qci(i,k,1) = max(qci(i,k,1)-(praut(i,k)+pracw(i,k)                 &
                    +paacw(i,k)+paacw(i,k))*dtcld,0.)
            qrs(i,k,1) = max(qrs(i,k,1)+(praut(i,k)+pracw(i,k)                 &
                    +prevp(i,k)+paacw(i,k)+paacw(i,k)-pseml(i,k)               &
                    -pgeml(i,k))*dtcld,0.)
            qrs(i,k,2) = max(qrs(i,k,2)+(psevp(i,k)-pgacs(i,k)                 &
                    +pseml(i,k))*dtcld,0.)
            qrs(i,k,3) = max(qrs(i,k,3)+(pgacs(i,k)+pgevp(i,k)                 &
                    +pgeml(i,k))*dtcld,0.)
            xlf = xls-xl(i,k)
            xlwork2 = -xl(i,k)*(prevp(i,k)+psevp(i,k)+pgevp(i,k))              &
                      -xlf*(pseml(i,k)+pgeml(i,k))
            t(i,k) = t(i,k)-xlwork2/cpm(i,k)*dtcld
          endif
        enddo
      enddo
!
! Inline expansion for fpvs
!         qs(i,k,1) = fpvs(t(i,k),0,rd,rv,cpv,cliq,cice,xlv0,xls,psat,t0c)
!         qs(i,k,2) = fpvs(t(i,k),1,rd,rv,cpv,cliq,cice,xlv0,xls,psat,t0c)
      hsub = xls
      hvap = xlv0
      cvap = cpv
      ttp=t0c+0.01
      dldt=cvap-cliq
      xa=-dldt/rv
      xb=xa+hvap/(rv*ttp)
      dldti=cvap-cice
      xai=-dldti/rv
      xbi=xai+hsub/(rv*ttp)
!ACC$DO PARALLEL(2)
      do k = kts, kte
!ACC$DO PARALLEL(1,its:ite)
!$acc loop gang vector
        do i = its, ite
          tr=ttp/t(i,k)
          qs(i,k,1)=psat*exp(log(tr)*(xa))*exp(xb*(1.-tr))
          qs(i,k,1) = min(qs(i,k,1),0.99*p(i,k))
          qs(i,k,1) = ep2 * qs(i,k,1) / (p(i,k) - qs(i,k,1))
          qs(i,k,1) = max(qs(i,k,1),qmin)
          tr=ttp/t(i,k)
          if(t(i,k).lt.ttp) then
            qs(i,k,2)=psat*exp(log(tr)*(xai))*exp(xbi*(1.-tr))
          else
            qs(i,k,2)=psat*exp(log(tr)*(xa))*exp(xb*(1.-tr))
          endif
          qs(i,k,2) = min(qs(i,k,2),0.99*p(i,k))
          qs(i,k,2) = ep2 * qs(i,k,2) / (p(i,k) - qs(i,k,2))
          qs(i,k,2) = max(qs(i,k,2),qmin)
        enddo
      enddo
!
!----------------------------------------------------------------
!  pcond: condensational/evaporational rate of cloud water [HL A46] [RH83 A6]
!     if there exists additional water vapor condensated/if
!     evaporation of cloud water is not enough to remove subsaturation
!
!ACC$DO PARALLEL(2)
      do k = kts, kte
!ACC$DO PARALLEL(1,its:ite)
!$acc loop gang vector
        do i = its, ite
#ifdef INLINE
          !conden(a,b,c,d,e) = (max(b,qmin)-c)/(1.+d*d/(rv*e)*c/(a*a))
          !work1(i,k,1) = conden(t(i,k),q(i,k),qs(i,k,1),xl(i,k),cpm(i,k))
          work1(i,k,1) = (max(q(i,k),qmin)-qs(i,k,1))/(1.+xl(i,k)*xl(i,k)   &
                       / (rv*cpm(i,k))*qs(i,k,1)/(t(i,k)*t(i,k)))
#else
          work1(i,k,1) = conden(t(i,k),q(i,k),qs(i,k,1),xl(i,k),cpm(i,k))
#endif
          work2(i,k) = qci(i,k,1)+work1(i,k,1)
          pcond(i,k) = min(max(work1(i,k,1)/dtcld,0.),max(q(i,k),0.)/dtcld)
          if(qci(i,k,1).gt.0..and.work1(i,k,1).lt.0.)                          &
            pcond(i,k) = max(work1(i,k,1),-qci(i,k,1))/dtcld
          q(i,k) = q(i,k)-pcond(i,k)*dtcld
          qci(i,k,1) = max(qci(i,k,1)+pcond(i,k)*dtcld,0.)
          t(i,k) = t(i,k)+pcond(i,k)*xl(i,k)/cpm(i,k)*dtcld
        enddo
      enddo
!
!
!----------------------------------------------------------------
!     padding for small values
!
!ACC$DO PARALLEL(2)
      do k = kts, kte
#ifdef ALIGN_OK
!DIR$ VECTOR ALIGNED
#endif
!ACC$DO PARALLEL(1,its:ite)
!$acc loop gang vector
        do i = its, ite
          if(qci(i,k,1).le.qmin) qci(i,k,1) = 0.0
          if(qci(i,k,2).le.qmin) qci(i,k,2) = 0.0
        enddo
      enddo
!$acc end parallel
!ACC$REGION END
      enddo                  ! big loops
!ACC$DATA(<t,qci,qrs,rain,rainncv,snow,snowncv,graupel,graupelncv,sr,q:out,extern>)

  END SUBROUTINE wsm62d
!--------------------------------------------------------------------------
!ACC$ROUTINE(GPU:<96:chunk>,<(ite-its+1)/96+1>,<n0sfac:none,local,demote(1,2)>, &
!ACC$> <alpha,qcrmin,n0s,n0smax,pvtr,pvts,pvtg:none,constant>, &
!ACC$> <rslopermax,rsloperbmax,rsloper2max,rsloper3max,bvtr,pidn0r, &
!ACC$>  rslopesmax,rslopesbmax,rslopes2max,rslopes3max,bvts,pidn0s, &
!ACC$>  rslopegmax,rslopegbmax,rslopeg2max,rslopeg3max,bvtg,pidn0g:none,constant>)
      subroutine slope_wsm6(qrs,den,denfac,t,rslope,rslopeb,rslope2,rslope3,   &
                            vt,its,ite,kts,kte)
  IMPLICIT NONE
  INTEGER       ::               its,ite, jts,jte, kts,kte
  REAL, DIMENSION( its:ite , kts:kte,3) ::                                     &
                                                                          qrs, &
                                                                       rslope, &
                                                                      rslopeb, &                                                 
                                                                      rslope2, &                                                 
                                                                      rslope3, &                                                 
                                                                           vt
  REAL, DIMENSION( its:ite , kts:kte) ::                                       &
                                                                          den, &
                                                                       denfac, &
                                                                            t
  REAL, PARAMETER  :: t0c = 273.15
  REAL, DIMENSION( its:ite , kts:kte ) ::                                      &
                                                                       n0sfac
  REAL       ::  lamdar, lamdas, lamdag, x, y, z, supcol
  integer :: i, j, k
!----------------------------------------------------------------
!     size distributions: (x=mixing ratio, y=air density):
!     valid for mixing ratio > 1.e-9 kg/kg.
#ifndef INLINE
      lamdar(x,y)=   sqrt(sqrt(pidn0r/(x*y)))      ! (pidn0r/(x*y))**.25
      lamdas(x,y,z)= sqrt(sqrt(pidn0s*z/(x*y)))    ! (pidn0s*z/(x*y))**.25
      lamdag(x,y)=   sqrt(sqrt(pidn0g/(x*y)))      ! (pidn0g/(x*y))**.25
#endif
!
#ifdef ALIGN_OK
!DIR$ ASSUME_ALIGNED qrs:64,den:64,denfac:64,t:64,rslope:64,rslopeb:64,rslope2:64,rslope3:64,vt:64
#endif
!$acc parallel vector_length(128)
!ACC$DO PARALLEL(2)
      do k = kts, kte
#ifdef ALIGN_OK
!DIR$ VECTOR ALIGNED
#endif
!ACC$DO PARALLEL(1,its:ite)
!$acc loop gang vector
        do i = its, ite
          supcol = t0c-t(i,k)
!---------------------------------------------------------------
! n0s: Intercept parameter for snow [m-4] [HDC 6]
!---------------------------------------------------------------
          n0sfac(i,k) = max(min(exp(alpha*supcol),n0smax/n0s),1.)
          if(qrs(i,k,1).le.qcrmin)then
            rslope(i,k,1) = rslopermax
            rslopeb(i,k,1) = rsloperbmax
            rslope2(i,k,1) = rsloper2max
            rslope3(i,k,1) = rsloper3max
          else
#ifdef INLINE
            rslope(i,k,1) = 1./sqrt(sqrt(pidn0r/(qrs(i,k,1)*den(i,k))))
#else
            rslope(i,k,1) = 1./lamdar(qrs(i,k,1),den(i,k))
#endif
            rslopeb(i,k,1) = rslope(i,k,1)**bvtr
            rslope2(i,k,1) = rslope(i,k,1)*rslope(i,k,1)
            rslope3(i,k,1) = rslope2(i,k,1)*rslope(i,k,1)
          endif
          if(qrs(i,k,2).le.qcrmin)then
            rslope(i,k,2) = rslopesmax
            rslopeb(i,k,2) = rslopesbmax
            rslope2(i,k,2) = rslopes2max
            rslope3(i,k,2) = rslopes3max
          else
#ifdef INLINE
            rslope(i,k,2) = 1./sqrt(sqrt(pidn0s*n0sfac(i,k)/(qrs(i,k,2)*den(i,k))))
#else
            rslope(i,k,2) = 1./lamdas(qrs(i,k,2),den(i,k),n0sfac(i,k))
#endif
            rslopeb(i,k,2) = rslope(i,k,2)**bvts
            rslope2(i,k,2) = rslope(i,k,2)*rslope(i,k,2)
            rslope3(i,k,2) = rslope2(i,k,2)*rslope(i,k,2)
          endif
          if(qrs(i,k,3).le.qcrmin)then
            rslope(i,k,3) = rslopegmax
            rslopeb(i,k,3) = rslopegbmax
            rslope2(i,k,3) = rslopeg2max
            rslope3(i,k,3) = rslopeg3max
          else
#ifdef INLINE
            rslope(i,k,3) = 1./sqrt(sqrt(pidn0g/(qrs(i,k,3)*den(i,k))))
#else
            rslope(i,k,3) = 1./lamdag(qrs(i,k,3),den(i,k))
#endif
            rslopeb(i,k,3) = rslope(i,k,3)**bvtg
            rslope2(i,k,3) = rslope(i,k,3)*rslope(i,k,3)
            rslope3(i,k,3) = rslope2(i,k,3)*rslope(i,k,3)
          endif
          vt(i,k,1) = pvtr*rslopeb(i,k,1)*denfac(i,k)
          vt(i,k,2) = pvts*rslopeb(i,k,2)*denfac(i,k)
          vt(i,k,3) = pvtg*rslopeb(i,k,3)*denfac(i,k)
          if(qrs(i,k,1).le.0.0) vt(i,k,1) = 0.0
          if(qrs(i,k,2).le.0.0) vt(i,k,2) = 0.0
          if(qrs(i,k,3).le.0.0) vt(i,k,3) = 0.0
        enddo
      enddo
!$acc end parallel

  END subroutine slope_wsm6
!-----------------------------------------------------------------------------
!MWG: F2C-ACC ignores the thread and block arguments here
!ACC$ROUTINE(GPU:<1>,<1>,<rslopermax,rsloperbmax,rsloper2max,rsloper3max,    &
!ACC$>                      bvtr,pidn0r,qcrmin,pvtr:none,constant>)
      subroutine slope_rain(qrs,den,denfac,t,rslope,rslopeb,rslope2,rslope3,   & 
                            vt,kts,kte)
  IMPLICIT NONE
  INTEGER       ::               kts,kte
  REAL, DIMENSION( kts:kte) ::                                                 &
                                                                          qrs, &
                                                                       rslope, &
                                                                      rslopeb, &
                                                                      rslope2, &
                                                                      rslope3, &
                                                                           vt, &      
                                                                          den, &
                                                                       denfac, &
                                                                            t
  REAL, PARAMETER  :: t0c = 273.15
  REAL, DIMENSION( kts:kte ) ::                                                &
                                                                       n0sfac
  REAL       ::  lamdar, x, y, z, supcol
  integer :: k
!----------------------------------------------------------------
!     size distributions: (x=mixing ratio, y=air density):
!     valid for mixing ratio > 1.e-9 kg/kg.
#ifndef INLINE
      lamdar(x,y)=   sqrt(sqrt(pidn0r/(x*y)))      ! (pidn0r/(x*y))**.25
#endif
!
#ifdef ALIGN_OK
!DIR$ ASSUME_ALIGNED qrs:64,den:64,denfac:64,t:64,rslope:64,rslopeb:64,rslope2:64,rslope3:64,vt:64
!DIR$ VECTOR ALIGNED
#endif
      do k = kts, kte
          if(qrs(k).le.qcrmin)then
            rslope(k) = rslopermax
            rslopeb(k) = rsloperbmax
            rslope2(k) = rsloper2max
            rslope3(k) = rsloper3max
          else
#ifdef INLINE
            rslope(k) = 1./sqrt(sqrt(pidn0r/(qrs(k)*den(k)))) 
#else
            rslope(k) = 1./lamdar(qrs(k),den(k))
#endif
            rslopeb(k) = rslope(k)**bvtr
            rslope2(k) = rslope(k)*rslope(k)
            rslope3(k) = rslope2(k)*rslope(k)
          endif
          vt(k) = pvtr*rslopeb(k)*denfac(k)
          if(qrs(k).le.0.0) vt(k) = 0.0
      enddo
  END subroutine slope_rain
!------------------------------------------------------------------------------
!MWG: F2C-ACC ignores the thread and block arguments here
!ACC$ROUTINE(GPU:<1>,<1>,<rslopesmax,rslopesbmax,rslopes2max,rslopes3max, &
!ACC$> n0smax,n0s,pidn0s,bvts,pvts,alpha,qcrmin:none,constant>)
      subroutine slope_snow(qrs,den,denfac,t,rslope,rslopeb,rslope2,rslope3,   &
                            vt,kts,kte)

  IMPLICIT NONE
  INTEGER       ::               kts,kte
  REAL, DIMENSION( kts:kte) ::                                                 &
                                                                          qrs, &
                                                                       rslope, &
                                                                      rslopeb, &
                                                                      rslope2, &
                                                                      rslope3, &
                                                                           vt, &  
                                                                          den, &
                                                                       denfac, &
                                                                            t
  REAL, PARAMETER  :: t0c = 273.15
  REAL, DIMENSION( kts:kte ) ::                                                &
                                                                       n0sfac
  REAL       ::  lamdas, x, y, z, supcol
  integer :: k
!----------------------------------------------------------------
!     size distributions: (x=mixing ratio, y=air density):
!     valid for mixing ratio > 1.e-9 kg/kg.
#ifndef INLINE
      lamdas(x,y,z)= sqrt(sqrt(pidn0s*z/(x*y)))    ! (pidn0s*z/(x*y))**.25
#endif
!
#ifdef ALIGN_OK
!DIR$ ASSUME_ALIGNED qrs:64,den:64,denfac:64,t:64,rslope:64,rslopeb:64,rslope2:64,rslope3:64,vt:64
!DIR$ VECTOR ALIGNED
#endif
      do k = kts, kte
          supcol = t0c-t(k)
!---------------------------------------------------------------
! n0s: Intercept parameter for snow [m-4] [HDC 6]
!---------------------------------------------------------------
          n0sfac(k) = max(min(exp(alpha*supcol),n0smax/n0s),1.)
          if(qrs(k).le.qcrmin)then
            rslope(k) = rslopesmax
            rslopeb(k) = rslopesbmax
            rslope2(k) = rslopes2max
            rslope3(k) = rslopes3max
          else
#ifdef INLINE
            rslope(k) = 1./sqrt(sqrt(pidn0s*n0sfac(k)/(qrs(k)*den(k))))
#else
            rslope(k) = 1./lamdas(qrs(k),den(k),n0sfac(k))
#endif
            rslopeb(k) = rslope(k)**bvts
            rslope2(k) = rslope(k)*rslope(k)
            rslope3(k) = rslope2(k)*rslope(k)
          endif
          vt(k) = pvts*rslopeb(k)*denfac(k)
          if(qrs(k).le.0.0) vt(k) = 0.0
      enddo
  END subroutine slope_snow
!------------------------------------------------------------------------------
!MWG: F2C-ACC ignores the thread and block arguments here
!ACC$ROUTINE(GPU:<1>,<1>,<rslopegmax,rslopegbmax,rslopeg2max,rslopeg3max,    &
!ACC$>                   qcrmin,pidn0g,bvtg,pvtg:none,constant>)
      subroutine slope_graup(qrs,den,denfac,t,rslope,rslopeb,rslope2,rslope3,  &
                            vt,kts,kte) 
  IMPLICIT NONE
  INTEGER       :: kts,kte
  REAL, DIMENSION( kts:kte) ::                                                 &
                                                                          qrs, &
                                                                       rslope, &
                                                                      rslopeb, &
                                                                      rslope2, &
                                                                      rslope3, &
                                                                           vt, &  
                                                                          den, &
                                                                       denfac, &
                                                                            t
  REAL, PARAMETER  :: t0c = 273.15
  REAL, DIMENSION( kts:kte ) ::                                                &
                                                                       n0sfac
  REAL       ::  lamdag, x, y, z, supcol
  integer :: j, k
!----------------------------------------------------------------
!     size distributions: (x=mixing ratio, y=air density):
!     valid for mixing ratio > 1.e-9 kg/kg.
#ifndef INLINE
      lamdag(x,y)=   sqrt(sqrt(pidn0g/(x*y)))      ! (pidn0g/(x*y))**.25
#endif
!
#ifdef ALIGN_OK
!DIR$ ASSUME_ALIGNED qrs:64,den:64,denfac:64,t:64,rslope:64,rslopeb:64,rslope2:64,rslope3:64,vt:64
!DIR$ VECTOR ALIGNED
#endif
      do k = kts, kte
!---------------------------------------------------------------
! n0s: Intercept parameter for snow [m-4] [HDC 6]
!---------------------------------------------------------------
          if(qrs(k).le.qcrmin)then
            rslope(k) = rslopegmax
            rslopeb(k) = rslopegbmax
            rslope2(k) = rslopeg2max
            rslope3(k) = rslopeg3max
          else
#ifdef INLINE
            rslope(k) = 1./sqrt(sqrt(pidn0g/(qrs(k)*den(k))))
#else
            rslope(k) = 1./lamdag(qrs(k),den(k))
#endif
            rslopeb(k) = rslope(k)**bvtg
            rslope2(k) = rslope(k)*rslope(k)
            rslope3(k) = rslope2(k)*rslope(k)
          endif
          vt(k) = pvtg*rslopeb(k)*denfac(k)
          if(qrs(k).le.0.0) vt(k) = 0.0
      enddo
  END subroutine slope_graup

!MWG: F2C-ACC ignores the thread and block arguments here
!ACC$ROUTINE(CPU:<1>,<1>,<pidn0r,bvtr,pvtr,rslopermax,rsloperbmax,rsloper2max,rsloper3max:none,constant>)
      SUBROUTINE nislfv_rain_plm(its,ite,kts,kte,denl,denfacl,tkl,dzl,wwl,rql,precip,dt,id,iter)
!-------------------------------------------------------------------
!
! for non-iteration semi-Lagrangain forward advection for cloud
! with mass conservation and positive definite advection
! 2nd order interpolation with monotonic piecewise linear method
! this routine is under assumption of decfl < 1 for semi_Lagrangian
!
! dzl    depth of model layer in meter
! wwl    terminal velocity at model layer m/s
! rql    cloud density*mixing ration
! precip precipitation
! dt     time step
! id     kind of precip: 0 test case; 1 raindrop
! iter   how many time to guess mean terminal velocity: 0 pure forward.
!        0 : use departure wind for advection
!        1 : use mean wind for advection
!        > 1 : use mean wind after iter-1 iterations
!
! author: hann-ming henry juang <henry.juang@noaa.gov>
!         implemented by song-you hong
!
      implicit none
      integer,intent(IN) :: its,ite,kts,kte,id
      real  dt
      real  dzl(its:ite,kts:kte),wwl(its:ite,kts:kte),rql(its:ite,kts:kte),precip(its:ite)
      real  denl(its:ite,kts:kte),denfacl(its:ite,kts:kte),tkl(its:ite,kts:kte)
!
      integer  i,k,n,m,kk,kb,kt,iter
      real  tl,tl2,qql,dql,qqd
      real  th,th2,qqh,dqh
      real  zsum,qsum,dim,dip,c1,con1,fa1,fa2
      real  allold, allnew, zz, dzamin, cflmax, decfl
      real  dz(kts:kte), ww(kts:kte), qq(kts:kte), wd(kts:kte), wa(kts:kte), was(kts:kte)
      real  den(kts:kte), denfac(kts:kte), tk(kts:kte)
      real  wi(kts:kte+1), zi(kts:kte+1), za(kts:kte+1)
      real  qn(kts:kte), qr(kts:kte),tmp(kts:kte),tmp1(kts:kte),tmp2(kts:kte),tmp3(kts:kte)
      real  dza(kts:kte+1), qa(kts:kte+1), qmi(kts:kte+1), qpi(kts:kte+1)
!
#ifdef ALIGN_OK
!DIR$ ASSUME_ALIGNED denl:64,denfacl:64,tkl:64,dzl:64,wwl:64,rql:64,precip:64
#endif
!ACC$REGION(<64:chunk>,<(ite-its+1)/64+1>,<dz,ww,qq,wd,wa,was,den,denfac,tk,wi,zi,za,qn,qr,tmp,tmp1,tmp2,tmp3,dza,qa,qmi,qpi:none,local>, &
!ACC$> <pidn0r,bvtr,pvtr,rslopermax,rsloperbmax,rsloper2max,rsloper3max:none,constant>, &
!ACC$> <rql=>denqci[iter==0],rql=>denqrs1[iter==1]:none,extern>, &
!ACC$> <tkl=>t,precip=>delqi[iter==0],precip=>delqrs1[iter==1]:none,extern> ) BEGIN
#ifdef CRAYOPENACC
!$acc parallel vector_length(128) 
#endif
!ACC$DO PARALLEL(1,its:ite)
#ifdef CRAYOPENACC
!$acc loop gang vector
#endif
      do i=its,ite
        precip(i) = 0.0
      enddo
!
!ACC$DO PARALLEL(1,its:ite)
#ifdef CRAYOPENACC
!$acc loop gang vector private(dz,ww,qq,wd,wa,was,den,denfac,tk,wi,zi,za,qn,qr,tmp,tmp1,tmp2,tmp3,dza,qa,qmi,qpi)
#endif

      i_loop : do i=its,ite
! -----------------------------------
      do k=kts,kte
        dz(k) = dzl(i,k)
        qq(k) = rql(i,k)
        ww(k) = wwl(i,k)
        den(k) = denl(i,k)
        denfac(k) = denfacl(i,k)
        tk(k) = tkl(i,k)
      enddo
! skip for no precipitation for all layers
      allold = 0.0
      do k=kts,kte
        allold = allold + qq(k)
      enddo
      if(allold.le.0.0) then
        cycle i_loop
      endif
!
! compute interface values
      zi(kts)=0.0
      do k=kts,kte
        zi(k+1) = zi(k)+dz(k)
      enddo
!
! save departure wind
      do k=kts,kte
        wd(k) = ww(k)
      enddo
      n=1
 100  continue
! plm is 2nd order, we can use 2nd order wi or 3rd order wi
! 2nd order interpolation to get wi
      wi(kts) = ww(kts)
      wi(kte+1) = ww(kte)
      do k=kts+1,kte
        wi(k) = (ww(k)*dz(k-1)+ww(k-1)*dz(k))/(dz(k-1)+dz(k))
      enddo
! 3rd order interpolation to get wi
      fa1 = 9./16.
      fa2 = 1./16.
      wi(kts) = ww(kts)
      wi(kts+1) = 0.5*(ww(kts+1)+ww(kts))
      do k=kts+2,kte-1
        wi(k) = fa1*(ww(k)+ww(k-1))-fa2*(ww(k+1)+ww(k-2))
      enddo
      wi(kte) = 0.5*(ww(kte)+ww(kte-1))
      wi(kte+1) = ww(kte)
!
! terminate of top of raingroup
      do k=kts+1,kte
        if( ww(k).eq.0.0 ) wi(k)=ww(k-1)
      enddo
!
! diffusivity of wi
      con1 = 0.05
      do k=kte,kts,-1
        decfl = (wi(k+1)-wi(k))*dt/dz(k)
        if( decfl .gt. con1 ) then
          wi(k) = wi(k+1) - con1*dz(k)/dt
        endif
      enddo
! compute arrival point
      do k=kts,kte+1
        za(k) = zi(k) - wi(k)*dt
      enddo
!
      do k=kts,kte
        dza(k) = za(k+1)-za(k)
      enddo
      dza(kte+1) = zi(kte+1) - za(kte+1)
!
! computer deformation at arrival point
      do k=kts,kte
        qa(k) = qq(k)*dz(k)/dza(k)
        qr(k) = qa(k)/den(k)
      enddo
      qa(kte+1) = 0.0
!     call maxmin(kte-kts+1,1,qa,' arrival points ')
!
! compute arrival terminal velocity, and estimate mean terminal velocity
! then back to use mean terminal velocity
      if( n.le.iter ) then
        call slope_rain(qr,den,denfac,tk,tmp,tmp1,tmp2,tmp3,wa,kts,kte)
        if( n.ge.2 ) then
          do k=kte,kte
            wa(k)=0.5*(wa(k)+was(k))
          enddo
        endif
        do k=kts,kte
!#ifdef DEBUG
!        print*,' slope_wsm3 ',qr(k)*1000.,den(k),denfac(k),tk(k),tmp(k),tmp1(k),tmp2(k),ww(k),wa(k)
!#endif
! mean wind is average of departure and new arrival winds
          ww(k) = 0.5* ( wd(k)+wa(k) )
        enddo
        do k=kte,kte
        was(k) = wa(k)
        enddo
        n=n+1
        go to 100
      endif
!
! estimate values at arrival cell interface with monotone
      do k=kts+1,kte
        dip=(qa(k+1)-qa(k))/(dza(k+1)+dza(k))
        dim=(qa(k)-qa(k-1))/(dza(k-1)+dza(k))
        if( dip*dim.le.0.0 ) then
          qmi(k)=qa(k)
          qpi(k)=qa(k)
        else
          qpi(k)=qa(k)+0.5*(dip+dim)*dza(k)
          qmi(k)=2.0*qa(k)-qpi(k)
          if( qpi(k).lt.0.0 .or. qmi(k).lt.0.0 ) then
            qpi(k) = qa(k)
            qmi(k) = qa(k)
          endif
        endif
      enddo
      qpi(kts)=qa(kts)
      qmi(kts)=qa(kts)
      qmi(kte+1)=qa(kte+1)
      qpi(kte+1)=qa(kte+1)
!
! interpolation to regular point
      do k = kts,kte
        qn(k) = 0.0
      enddo
      kb=kts
      kt=kts
      intp : do k=kts,kte
             kb=max(kb-1,kts)
             kt=max(kt-1,kts)
! find kb and kt
             if( zi(k).ge.za(kte+1) ) then
               exit intp
             else
               find_kb : do kk=kb,kte
                         if( zi(k).le.za(kk+1) ) then
                           kb = kk
                           exit find_kb
                         else
                           cycle find_kb
                         endif
               enddo find_kb
               find_kt : do kk=kt,kte
                         if( zi(k+1).le.za(kk) ) then
                           kt = kk
                           exit find_kt
                         else
                           cycle find_kt
                         endif
               enddo find_kt
               kt = kt - 1
! compute q with piecewise constant method
               if( kt.eq.kb ) then
                 tl=(zi(k)-za(kb))/dza(kb)
                 th=(zi(k+1)-za(kb))/dza(kb)
                 tl2=tl*tl
                 th2=th*th
                 qqd=0.5*(qpi(kb)-qmi(kb))
                 qqh=qqd*th2+qmi(kb)*th
                 qql=qqd*tl2+qmi(kb)*tl
                 qn(k) = (qqh-qql)/(th-tl)
               else if( kt.gt.kb ) then
                 tl=(zi(k)-za(kb))/dza(kb)
                 tl2=tl*tl
                 qqd=0.5*(qpi(kb)-qmi(kb))
                 qql=qqd*tl2+qmi(kb)*tl
                 dql = qa(kb)-qql
                 zsum  = (1.-tl)*dza(kb)
                 qsum  = dql*dza(kb)
                 if( kt-kb.gt.1 ) then
                 do m=kb+1,kt-1
                   zsum = zsum + dza(m)
                   qsum = qsum + qa(m) * dza(m)
                 enddo
                 endif
                 th=(zi(k+1)-za(kt))/dza(kt)
                 th2=th*th
                 qqd=0.5*(qpi(kt)-qmi(kt))
                 dqh=qqd*th2+qmi(kt)*th
                 zsum  = zsum + th*dza(kt)
                 qsum  = qsum + dqh*dza(kt)
                 qn(k) = qsum/zsum
               endif
               cycle intp
             endif
!
       enddo intp
!
! rain out
      sum_precip: do k=kts,kte
                    if( za(k).lt.0.0 .and. za(k+1).lt.0.0 ) then
                      precip(i) = precip(i) + qa(k)*dza(k)
                      cycle sum_precip
                    else if ( za(k).lt.0.0 .and. za(k+1).ge.0.0 ) then
                      precip(i) = precip(i) + qa(k)*(0.0-za(k))
                      exit sum_precip
                    endif
                    exit sum_precip
      enddo sum_precip
!
! replace the new values
#ifdef ALIGN_OK
!DIR$ VECTOR ALIGNED
#endif
      do k=kts,kte
        rql(i,k) = qn(k)
      enddo
!
! ----------------------------------
      enddo i_loop
#ifdef CRAYOPENACC
!$acc end parallel
#endif
!ACC$REGION END
!
  END SUBROUTINE nislfv_rain_plm

!-----------------------------------------------------------------------------
!MWG: F2C-ACC ignores the thread and block arguments here
!ACC$ROUTINE(CPU:<1>,<1>,<pidn0s,pidn0g,pvts,pvtg,rslopesmax,rslopesbmax,rslopes2max,rslopes3max,rslopegmax,rslopegbmax,rslopeg2max,rslopeg3max:none,constant>)
      SUBROUTINE nislfv_rain_plm6(its,ite,kts,kte,denl,denfacl,tkl,dzl,wwl,rql,rql2, precip1, precip2,dt,id,iter)
!-------------------------------------------------------------------
!
! for non-iteration semi-Lagrangain forward advection for cloud
! with mass conservation and positive definite advection
! 2nd order interpolation with monotonic piecewise linear method
! this routine is under assumption of decfl < 1 for semi_Lagrangian
!
! dzl    depth of model layer in meter
! wwl    terminal velocity at model layer m/s
! rql    cloud density*mixing ration
! precip precipitation
! dt     time step
! id     kind of precip: 0 test case; 1 raindrop
! iter   how many time to guess mean terminal velocity: 0 pure forward.
!        0 : use departure wind for advection
!        1 : use mean wind for advection
!        > 1 : use mean wind after iter-1 iterations
!
! author: hann-ming henry juang <henry.juang@noaa.gov>
!         implemented by song-you hong
!
      implicit none
      integer  its,ite,kts,kte,id
      real  dt
      real  dzl(its:ite,kts:kte),wwl(its:ite,kts:kte),rql(its:ite,kts:kte),rql2(its:ite,kts:kte),precip(its:ite),precip1(its:ite),precip2(its:ite)
      real  denl(its:ite,kts:kte),denfacl(its:ite,kts:kte),tkl(its:ite,kts:kte)
!
      integer  i,k,n,m,kk,kb,kt,iter,ist
      real  tl,tl2,qql,dql,qqd
      real  th,th2,qqh,dqh
      real  zsum,qsum,dim,dip,c1,con1,fa1,fa2
      real  allold, allnew, zz, dzamin, cflmax, decfl
      real  dz(kts:kte), ww(kts:kte), qq(kts:kte), qq2(kts:kte), wd(kts:kte), wa(kts:kte), wa2(kts:kte), was(kts:kte)
      real  den(kts:kte), denfac(kts:kte), tk(kts:kte)
      real  wi(kts:kte+1), zi(kts:kte+1), za(kts:kte+1)
      real  qn(kts:kte), qr(kts:kte),qr2(kts:kte),tmp(kts:kte),tmp1(kts:kte),tmp2(kts:kte),tmp3(kts:kte)
      real  dza(kts:kte+1), qa(kts:kte+1), qa2(kts:kte+1),qmi(kts:kte+1), qpi(kts:kte+1)
!
!ACC$REGION(<96:chunk>,<(ite-its+1)/96+1>,<dz,ww,qq,qq2,wd,wa,wa2,was,den,denfac,tk,wi,zi,za,qn,qr,qr2,tmp,tmp1,tmp2,tmp3,dza,qa,qa2,qmi,qpi:none,local>, &
!ACC$> <tkl=>t,precip1=>delqrs2,precip2=>delqrs3,rql=>denqrs2,rql2=>denqrs3:none,extern>) BEGIN
#ifdef CRAYOPENACC
!$acc parallel vector_length(128) 
#endif
!ACC$INSERTC double dim;
#ifdef ALIGN_OK
!DIR$ ASSUME_ALIGNED denl:64,denfacl:64,tkl:64,dzl:64,wwl:64,rql:64,rql2:64,precip1:64,precip2:64
#endif
!ACC$DO PARALLEL(1,its:ite)
!$acc loop gang vector
    do i=its,ite
      precip(i) = 0.0
      precip1(i) = 0.0
      precip2(i) = 0.0
    enddo
!
!ACC$DO PARALLEL(1,its:ite)
#ifdef CRAYOPENACC
!$acc loop gang vector private(dz,ww,qq,qq2,wd,wa,wa2,was,den,denfac,tk,wi,zi,za,qn,qr,qr2,tmp,tmp1,tmp2,tmp3,dza,qa,qa2,qmi,qpi)
#endif
    i_loop : do i=its,ite
! -----------------------------------
      do k=kts,kte
        dz(k) = dzl(i,k)
        qq(k) = rql(i,k)
        qq2(k) = rql2(i,k)
        ww(k) = wwl(i,k)
        den(k) = denl(i,k)
        denfac(k) = denfacl(i,k)
        tk(k) = tkl(i,k)
      enddo
! skip for no precipitation for all layers
      allold = 0.0
      do k=kts,kte
        allold = allold + qq(k)
      enddo
      if(allold.le.0.0) then
        cycle i_loop
      endif
!
! compute interface values
      zi(kts)=0.0
      do k=kts,kte
        zi(k+1) = zi(k)+dz(k)
      enddo
!
! save departure wind
      do k=kts,kte
        wd(k) = ww(k)
      enddo
      n=1
 100  continue
! plm is 2nd order, we can use 2nd order wi or 3rd order wi
! 2nd order interpolation to get wi
      wi(kts) = ww(kts)
      wi(kte+1) = ww(kte)
      do k=kts+1,kte
        wi(k) = (ww(k)*dz(k-1)+ww(k-1)*dz(k))/(dz(k-1)+dz(k))
      enddo
! 3rd order interpolation to get wi
      fa1 = 9./16.
      fa2 = 1./16.
      wi(kts) = ww(kts)
      wi(kts+1) = 0.5*(ww(kts+1)+ww(kts))
      do k=kts+2,kte-1
        wi(k) = fa1*(ww(k)+ww(k-1))-fa2*(ww(k+1)+ww(k-2))
      enddo
      wi(kte) = 0.5*(ww(kte)+ww(kte-1))
      wi(kte+1) = ww(kte)
!
! terminate of top of raingroup
      do k=kts+1,kte
        if( ww(k).eq.0.0 ) wi(k)=ww(k-1)
      enddo
!
! diffusivity of wi
      con1 = 0.05
      do k=kte,kts,-1
        decfl = (wi(k+1)-wi(k))*dt/dz(k)
        if( decfl .gt. con1 ) then
          wi(k) = wi(k+1) - con1*dz(k)/dt
        endif
      enddo
! compute arrival point
      do k=kts,kte+1
        za(k) = zi(k) - wi(k)*dt
      enddo
!
      do k=kts,kte
        dza(k) = za(k+1)-za(k)
      enddo
      dza(kte+1) = zi(kte+1) - za(kte+1)
!
! computer deformation at arrival point
      do k=kts,kte
        qa(k) = qq(k)*dz(k)/dza(k)
        qa2(k) = qq2(k)*dz(k)/dza(k)
        qr(k) = qa(k)/den(k)
        qr2(k) = qa2(k)/den(k)
      enddo
      qa(kte+1) = 0.0
      qa2(kte+1) = 0.0
!     call maxmin(kte-kts+1,1,qa,' arrival points ')
!
! compute arrival terminal velocity, and estimate mean terminal velocity
! then back to use mean terminal velocity
      if( n.le.iter ) then
        call slope_snow(qr,den,denfac,tk,tmp,tmp1,tmp2,tmp3,wa,kts,kte)
        call slope_graup(qr2,den,denfac,tk,tmp,tmp1,tmp2,tmp3,wa2,kts,kte)
        do k = kts, kte
          tmp(k) = max((qr(k)+qr2(k)), 1.E-15)
          IF ( tmp(k) .gt. 1.e-15 ) THEN
            wa(k) = (wa(k)*qr(k) + wa2(k)*qr2(k))/tmp(k)
          ELSE
            wa(k) = 0.
          ENDIF
        enddo
        if( n.ge.2 ) then
          do k=kts,kte
            wa(k)=0.5*(wa(k)+was(k))
          enddo
        endif
        do k=kts,kte
!#ifdef DEBUG
!        print*,' slope_wsm3 ',qr(k)*1000.,den(k),denfac(k),tk(k),tmp(k),tmp1(k),tmp2(k), &
!           ww(k),wa(k)
!#endif
! mean wind is average of departure and new arrival winds
          ww(k) = 0.5* ( wd(k)+wa(k) )
        enddo
        do k=kts,kte
          was(k) = wa(k)
        enddo
        n=n+1
        go to 100
      endif
      ist_loop : do ist = 1, 2
      if (ist.eq.2) then
        do k=kts,kte
          qa(k) = qa2(k)
        enddo
      endif
!
      precip(i) = 0.
!
! estimate values at arrival cell interface with monotone
      do k=kts+1,kte
        dip=(qa(k+1)-qa(k))/(dza(k+1)+dza(k))
        dim=(qa(k)-qa(k-1))/(dza(k-1)+dza(k))
        if( dip*dim.le.0.0 ) then
          qmi(k)=qa(k)
          qpi(k)=qa(k)
        else
          qpi(k)=qa(k)+0.5*(dip+dim)*dza(k)
          qmi(k)=2.0*qa(k)-qpi(k)
          if( qpi(k).lt.0.0 .or. qmi(k).lt.0.0 ) then
            qpi(k) = qa(k)
            qmi(k) = qa(k)
          endif
        endif
      enddo
      qpi(kts)=qa(kts)
      qmi(kts)=qa(kts)
      qmi(kte+1)=qa(kte+1)
      qpi(kte+1)=qa(kte+1)
!
! interpolation to regular point
      do k = kts,kte
        qn(k) = 0.0
      enddo
      kb=kts
      kt=kts
      intp : do k=kts,kte
             kb=max(kb-1,kts)
             kt=max(kt-1,kts)
! find kb and kt
             if( zi(k).ge.za(kte+1) ) then
               exit intp
             else
               find_kb : do kk=kb,kte
                         if( zi(k).le.za(kk+1) ) then
                           kb = kk
                           exit find_kb
                         else
                           cycle find_kb
                         endif
               enddo find_kb
               find_kt : do kk=kt,kte
                         if( zi(k+1).le.za(kk) ) then
                           kt = kk
                           exit find_kt
                         else
                           cycle find_kt
                         endif
               enddo find_kt
               kt = kt - 1
! compute q with piecewise constant method
               if( kt.eq.kb ) then
                 tl=(zi(k)-za(kb))/dza(kb)
                 th=(zi(k+1)-za(kb))/dza(kb)
                 tl2=tl*tl
                 th2=th*th
                 qqd=0.5*(qpi(kb)-qmi(kb))
                 qqh=qqd*th2+qmi(kb)*th
                 qql=qqd*tl2+qmi(kb)*tl
                 qn(k) = (qqh-qql)/(th-tl)
               else if( kt.gt.kb ) then
                 tl=(zi(k)-za(kb))/dza(kb)
                 tl2=tl*tl
                 qqd=0.5*(qpi(kb)-qmi(kb))
                 qql=qqd*tl2+qmi(kb)*tl
                 dql = qa(kb)-qql
                 zsum  = (1.-tl)*dza(kb)
                 qsum  = dql*dza(kb)
                 if( kt-kb.gt.1 ) then
                 do m=kb+1,kt-1
                   zsum = zsum + dza(m)
                   qsum = qsum + qa(m) * dza(m)
                 enddo
                 endif
                 th=(zi(k+1)-za(kt))/dza(kt)
                 th2=th*th
                 qqd=0.5*(qpi(kt)-qmi(kt))
                 dqh=qqd*th2+qmi(kt)*th
                 zsum  = zsum + th*dza(kt)
                 qsum  = qsum + dqh*dza(kt)
                 qn(k) = qsum/zsum
               endif
               cycle intp
             endif
!
       enddo intp
!
! rain out
      sum_precip: do k=kts,kte
                    if( za(k).lt.0.0 .and. za(k+1).lt.0.0 ) then
                      precip(i) = precip(i) + qa(k)*dza(k)
                      cycle sum_precip
                    else if ( za(k).lt.0.0 .and. za(k+1).ge.0.0 ) then
                      precip(i) = precip(i) + qa(k)*(0.0-za(k))
                      exit sum_precip
                    endif
                    exit sum_precip
      enddo sum_precip
!
! replace the new values
      if(ist.eq.1) then
        do k=kts,kte
          rql(i,k) = qn(k)
        enddo
        precip1(i) = precip(i)
      else
        do k=kts,kte
          rql2(i,k) = qn(k)
        enddo
        precip2(i) = precip(i)
      endif
      enddo ist_loop
!
! ----------------------------------
      enddo i_loop
#ifdef CRAYOPENACC
!$acc end parallel
#endif
!ACC$REGION END
!
  END SUBROUTINE nislfv_rain_plm6
END MODULE module_mp_wsm6
