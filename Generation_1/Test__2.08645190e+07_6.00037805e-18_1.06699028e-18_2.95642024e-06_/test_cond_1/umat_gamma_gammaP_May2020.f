!==================================================================
!
!  This is a version of the Kalidindi/Bronkhorst/Anand UMAT that has
!  been updated to remove COMMON blocks so that it can be run in 
!  parallel. 
!  
!  Temporary note:  This is a cleaned up version of the file located
!  in : ~/projects/curtsKinetics_mcdUmat/debug/run14
!
!  3/1/2016 Jason R. Mayeur
!  Further modifications have been made to replace the phenomenological
!  hardening model with a dislocation density based model.  A Hall-Petch
!  term has also been added to define the initial CRSS.
!
!  Additional code modifications/cleanup efforts have also been carried
!  out in the updating of the code, some of which has been documented
!  elsewher.
!------------------------------------------------------------------
      include "Read_euler.f" 
c +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c +                                                                         +
c +   Uexternaldb...                                                        +
c +       lop = 0 beginning of analysis                                     +
c +             1 start of increment                                        +
c +             2 end of increment                                          +
c +             3 end of analysis                                           +
c +             4 beginning of restart                                      +
c +                                                                         +
c +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine uexternaldb(lop,lrestart,time,dtime,kstep,kinc)
      use mod_read_data
      implicit none
      integer lop,lrestart,kstep,kinc
      real(8) time(2),dtime
      integer ios
      character(len=200) :: CurrentPath=''
      if(lop==0)then   !ini
         print*, kinc,dtime,time(2), 'Starting the analysis...'
         call read_euler()                                      ! subroutine used to read euler angles from a text file as a function of element number
         call getoutdir(CurrentPath, ios)
         print*, 'Running in directory',CurrentPath
      end if
      if(lop==1)then
         print*, kinc,dtime,time(2), 'start of the increment'
      end if

      if(lop==2)then
              print*, kinc,dtime,time(2), 'end of the increment'
      end if
      if(lop==3) then
        call getoutdir(CurrentPath, ios) 
        print*, 'End of the analysis...'       
        open(unit=25, file=trim(CurrentPath)//'/done.txt',
     & status='replace', action='write')
        close(unit=25)
        call XIT
      end if
      return
      end

      subroutine umat (stress,  statev,  ddsdde,  sse,     spd,
     &          scd,     rpl,     ddsddt,  drplde,  drpldt,
     &          stran,  dstran, time,    dtime,   temp,
     &          dtemp,   predef,  dpred,   cmname,  ndi,
     &          nshr,    ntens,   nstatv,  props,   nprops,
     &          coords,  drot,    pnewdt,  celent,  dfgrd0,
     &          dfgrd1,  noel,    npt,     layer,   kspt,
     &          kstep,   kinc)
      use mod_read_data                            !module in Read_euler.f file for reading euler angles as a function of element number using a txt file 
      implicit real*8(a-h,o-z)
      character*80 cmname 

!------------------------------------------------------------------
!  Dimension UMAT argument arrays
!------------------------------------------------------------------

      dimension
     &  coords(3),           ! Coordinates of Gauss point
     &  ddsdde(ntens,ntens), ! Tangent Stiffness Matrix
     &  ddsddt(ntens),       ! Derivative of stress wrt temperature
     &  dfgrd0(3,3),         ! Deformation gradient at beginning of step
     &  dfgrd1(3,3),         ! Deformation gradient at end of step
     &  dpred(1),            ! Predefined state variable array increment
     &  drplde(ntens),       ! Derivative of heat generation wrt strain
     &  drot(3,3),           ! Rotation matrix increment 
     &  dstran(ntens),       ! Strain vector increment
     &  predef(1),           ! Predefined variable array
     &  props(nprops),       ! Material properties array
     &  statev(nstatv),      ! State variables array
     &  stran(ntens),        ! Strain vector at beginning of increment
     &  stress(ntens),       ! Cauchy stress vector
     &  time(2)              ! Step time (1) and total time (2)

!------------------------------------------------------------------           
!  Dimension UMAT local arrays
!------------------------------------------------------------------           

      dimension
     &  array4  (6,6),       ! Full material Jacobian to be passed back to Abaqus 
     &  C_c  (3,3,3,3),      ! Elastic stiffness tensor in crystal basis
     &  C   (3,3,3,3),       ! Elastic stiffness tensor in reference config
     &  del (3,3),           ! Kronecker delta tensor
     &  q_0 (3,3),           ! Rotation tensor from crystal to reference config 
     &  Fe    (3,3),         ! Elastic deformation map
     &  FpInv_n(3,3),        ! Inverse plastic deformation map at beginning of step
     &  FpInv(3,3),          ! Inverse plastic deformation map at end of step
     &  s_n  (props(31)),    ! CRSS at beginning of step
     &  s   (props(31)),     ! CRSS at end of step
     &  bs_n (props(31)),    ! backstress at beginning of step
     &  bs  (props(31)),     ! backstress at end of step
     &  gammadot(props(31)), ! Slip strain rate array
     &  eulangs (3),         ! Euler angles as a function of element number
     &  sig (3,3),           ! Cauchy stress tensor
     &  Spk2    (3,3),       ! 2nd Piola-Kirchhoff stress
     &  tau (props(31)),     ! RSS
     &  xs0 (3,props(31)),   ! Slip vectors in intermediate config
     &  xm0 (3,props(31)),   ! Slip normal vectors in intermediate config
     &  xmC   (3,props(31)), ! Slip normal vectors in crystal basis 
     &  xsC   (3,props(31)), ! Slip vectors in crystal basis
     &  del4th(3,3,3,3),     ! Fourth rank symmetric identity tensor
     &  Fe_tr(3,3),          ! Trial elastic deformation map
     &  Ce_tr(3,3),          ! Trial elastic right stretch tensor
     &  Ee_tr(3,3),          ! Trial elastic Green strain tensor
     &  Spk2_tr(3,3),        ! Trial PK2 stress tensor
     &  Smat(3,3,props(31)), ! Schmid tensor matrix in reference config
     &  Bmat(3,3,props(31)), ! B matrix
     &  Cmat(3,3,props(31)), ! C matrix
     &  rLpBar(3,3),         ! Plastic velocity gradient in intermediate config
     &  smatc(props(31),3,3),! Schmid tensor matrix in crystal basis 
     &  spk2v_n(6),          ! PK2 stress vector at beginning of step
     &  spk2_n(3,3),         ! PK2 stress tensor at beginning of step
     &  spk2v_tr(6),         ! Trial PK2 stress vector
     &  pmatr(props(31),6),  ! Symmetric part of Schmid tensor in reduced form
     &  cmatr(props(31),6),  ! C matrix in reduced form
     &  spk2v(6),            ! PK2 stress vector at end of step
     &  dgamma(props(31)),   ! Plastic slip increment vector
     &  pmat(3,3,props(31)), ! Symmetric part of the Schmid tensor
     &  Re(3,3),             ! Elastic rotation tensor
     &  Ue(3,3),             ! Right elastic stretch tensor
     &  Q_n1(3,3),           ! Total lattice rotation tensor from crystal to current config

! Arrays added for dislocation density hardening law
     &  rhossd_n(props(31)), ! SSD density array at beginning of step
     &  rhossd(props(31)),   ! SSD density array at end of step
     &  amat(props(31),props(31)) ! Slip interaction matrix appearing in definition of CRSS

        pi = 4.d0 * atan(1.d0)
        gamma_acc = 0.d0

!------------------------------------------------------------------
!  Assign props() array to logical variable names
!------------------------------------------------------------------

        THETA0=PROPS(1)              ! Initial temperature [K]
        RHO=PROPS(2)                 ! Initial density [Mg/mm^3]
        CP=PROPS(3)                  ! Specific Heat [mJ/(Mg.K)]
        A11=PROPS(4)                 ! Thermal expansion coefficient [1/K]
        ETA=PROPS(5)                 ! Plastic work heat fraction

        rM11_G=PROPS(6)                ! Slope of C11 vs temp line for gamma [MPa/K]
        c11_0_G=PROPS(7)               ! C11 at 0 K for gamma[MPa]
        rM12_G=PROPS(8)                ! Slope of C12 vs temp line for gamma[MPa/K]
        c12_0_G=PROPS(9)               ! C12 at 0 K for gamma[MPa]
        rM44_G=PROPS(10)               ! Slope of C44 vs temp line for gamma[MPa/K]
        c44_0_G=PROPS(11)              ! C44 at 0 K for gamma[MPa]
        
! Structure independent CRSS [MPa],
        s_inf0=PROPS(12)             ! scaling for gamma' microstructural strength [MPa.mm^1/2], evaluated with the help of J.Lapin et al. Kovove Mater. 2009    
        burg = props(13)             ! Burgers vector [mm]
        rhossd_0 = props(14)         ! Initial SSD density [mm^-2]
        y_c0 = props(15)*burg        ! Capture distance for SSD annihilation [mm]
        A_rec = props(16)            ! Recovery exponent for SSD rate and temp dependence [J]
        a_self = props(17)           ! Dislocation interaction coefficient [-]
        a_copl = props(18)           ! Dislocation interaction coefficient [-]
        a_hirth = props(19)          ! Dislocation interaction coefficient [-]
        a_colli = props(20)          ! Dislocation interaction coefficient [-]
        a_gliss = props(21)          ! Dislocation interaction coefficient [-]
        a_lomer = props(22)          ! Dislocation interaction coefficient [-]
        rk_inter = props(23)         ! Dislocation trapping const for latent interactions [-]
        rk_copl = props(24)          ! Dislocation trapping const for self interactions [-]

        gdot0=PROPS(25)              ! Reference slip rate [1/s]
        f0=PROPS(26)                 ! Activation free energy [J] 
        BK=PROPS(27)                 ! Boltzmann's const. [J/K]
        sl=props(28)                 ! Intrinsic lattice resistance [MPa]
        PF=PROPS(29)                 ! Flow rule barrier exponent [-]
        QF=PROPS(30)                 ! Flow rule barrier exponent [-]
        nslip=int(PROPS(31))         ! Number of slip systems
        
! hall petch parameters 			 
        rk_hp=props(32)              ! Hall-Petch slope [MPa.mm^1/2],taken from Goodfellow et al. Materials & Design(2019)
        d_gr=props(33)               ! Grain diameter [mm]
        
! softening parameter 				       ! softening at very high temperature (close to melting)
        theta_ref=props(34)          ! reference temp for softening [K] (at this temperature softening and hardening balance out like ideal plasticity)
        gamma_ref=props(35)          ! ref strain for softening [-]
        soft_n=props(36)             ! softening exponent
        
! orientation reading flag
!        ebsd_read=int(props(37))    ! flag to decide, orientations are read from exp ebsd file or via props 1:on, 0:off

! APB shearing parameters 	         ! addition in initial slip resistance due to gamma'
        mu_gp100= props(37)          ! shear modulus of gamma prime
        APB_gp111=props(38)          ! in mj/mm2 for 111 plane
        f_gp=props(39)               ! volume fraction of gamma prime
        alpha_oro=props(40)          ! Orowan hardening for narrow gamma channel
        r_gp=props(41)               ! size of gamma prime precipitate
        
! Parameters for back stress        
        r_g=props(42)                ! size of gamma channel
        rb=props(43)                 ! fitting parameter in back stress
        flag_BS=int(props(44))       ! flag for switchon/off back stress 
        
!Bunge Euler angles in degree
        th1 = props(45)              ! Euler angle theta [degree] (Canova convention)
        phi1 = props(46)             ! Euler angle Phi [degree]
        om1 = props(47)              ! Euler angle omega [degree]
        
!Elasticity for gamma'        
        rM11_GP=PROPS(48)            ! Slope of C11 vs temp line for gamma'[MPa/K]
        c11_0_GP=PROPS(49)           ! C11 at 0 K for gamma'[MPa]
        rM12_GP=PROPS(50)            ! Slope of C12 vs temp line for gamma'[MPa/K]
        c12_0_GP=PROPS(51)           ! C12 at 0 K for gamma' [MPa]
        rM44_GP=PROPS(52)            ! Slope of C44 vs temp line for gamma'[MPa/K]
        c44_0_GP=PROPS(53)           ! C44 at 0 K for gamma' [MPa]
        
!Interaction coefficient for cube slip systems         
        cube_self=PROPS(54)          ! self hardening coefficient for cube slip
        cube_latent=PROPS(55)        ! latent hardening coefficient for cube slip
        
        APB_gp100=PROPS(56)          ! APB energy for cube slip
        sl_cube=PROPS(57)            ! Intrinsic lattice resistance for cube slip [MPa]     
        scale_APB= props(58)         ! fitting parameter for APB stress 
        debug = int(props(59))       ! 1: all the print/debug statements active 0: print/debug statements inactive

!------------------------------------------------------------------
!  Reading euler angles from text file: element by element
!------------------------------------------------------------------
!  		if(ebsd_read == 1) then  
            ! th = ang1(noel)*pi/180.d0     ! Euler angle theta [rad] (Canova convention)
            ! phi = ang2(noel)*pi/180.d0    ! Euler angle Phi [rad]
            ! om = ang3(noel)*pi/180.d0     ! Euler angle omega [rad]   
!  		else
        th = (180.d0-om1)*pi/180.d0     ! Euler angle theta [rad] (Canova convention)
        phi = phi1*pi/180.d0    ! Euler angle Phi [rad]
        om = (180.d0-th1)*pi/180.d0     ! Euler angle omega [rad]
!  		end if
        if (noel==20 .AND. npt==1 .AND. kinc==1 .AND. time(2)==0.0) then
           print*,'euler angles for element1 is',th,phi,om
           print*,'value of THETA0 is',THETA0
           print*,'value of RHO',RHO
           print*,'value of CP',CP
           print*,'value of A11',A11
           print*,'value of ETA',ETA
           print*,'value of rM11_G',rM11_G
           print*,'value of C11_0_G',C11_0_G
           print*,'value of rM12_G',rM12_G
           print*,'value of C12_0_G',C12_0_G
           print*,'value of rM44_G',rM44_G
           print*,'value of C44_0_G',C44_0_G
           print*,'value of s_inf0',s_inf0
           print*,'value of Burgers',burg
           print*,'value of rhossd_0',rhossd_0
           print*,'value of y_c0',y_c0
           print*,'value of A_rec',A_rec
           print*,'value of a_self',a_self
           print*,'value of a_copl',a_copl
           print*,'value of a_hirth',a_hirth
           print*,'value of a_colli',a_colli
           print*,'value of a_gliss',a_gliss
           print*,'value of a_lomer',a_lomer
           print*,'value of rk_inter',rk_inter
           print*,'value of rk_copl',rk_copl
           print*,'value of gdot0',gdot0
           print*,'value of f0',f0
           print*,'value of BK',BK
           print*,'value of sl',sl
           print*,'value of PF',PF
           print*,'value of QF',QF
           print*,'value of nslip',nslip
           print*,'value of rk_hp',rk_hp
           print*,'value of d_gr',d_gr
           print*,'value of theta_ref',theta_ref
           print*,'value of gamma_ref',gamma_ref
           print*,'value of soft_n',soft_n
           print*,'value of ebsd_read',ebsd_read
           print*,'value of mu_gp100',mu_gp100
           print*,'value of APB_gp111',APB_gp111
           print*,'value of f_gp',f_gp
           print*,'scaling factor in orowan stress',alpha_oro
           print*,'value of r_gp',r_gp
           print*,'value of r_g',r_g
           print*,'value of rb',rb
           print*,'value of th',th
           print*,'value of phi',phi
           print*,'value of om',om
           print*,'value of rM11_GP',rM11_GP
           print*,'value of c11_0_GP',c11_0_GP
           print*,'value of rM12_GP',rM12_GP
           print*,'value of c12_0_GP',c12_0_GP
           print*,'value of rM44_GP',rM44_GP
           print*,'value of c44_0_GP',c44_0_GP
           print*,'value of cube_self',cube_self
           print*,'value of cube_latent',cube_latent
           print*,'value of APB_gp100',APB_gp100
           print*,'value of sl_cube',sl_cube
           print*,'value of scale_APB',scale_APB
           print*,'value of debug',debug
        end if
                ! if (noel == 80 .AND. npt==1 .AND. kinc==1 .AND. time(2)==0.0) then
                !   print*,'euler angles for element80 is',th,phi,om
                ! print*,'value of the ebsd_read is',ebsd_read
                ! end if
!------------------------------------------------------------------
!  Initialize identity matrices/tensors
!------------------------------------------------------------------
      del = reshape( (/1.,0.,0.,0.,1.,0.,0.,0.,1. /), (/3,3/) )

      do i = 1,3
       do j = 1,3
        do k = 1,3
         do l = 1,3
          del4th(i,j,k,l) = 0.5d0 * 
     &      (del(i,k) * del(j,l) + del(i,l) * del(j,k))
         end do
        end do
       end do
      end do

!------------------------------------------------------------------
!  Assign slip system normals and slip directions for an FCC.
!           plane normal           slip direction
!------------------------------------------------------------------

        xmC(:,1)= (/-1,  1,  1/);  xsC(:,1)=(/0, -1, 1/)   !A2
        xmC(:,2)= (/-1,  1,  1/);  xsC(:,2)=(/1, 0, 1/)    !A3
        xmC(:,3)= (/-1,  1,  1/);  xsC(:,3)=(/1, 1, 0/)    !A6
        xmC(:,4)= (/1,   1,  1/);  xsC(:,4)=(/0, -1,  1/)  !B2
        xmC(:,5)= (/1,   1,  1/);  xsC(:,5)=(/-1,  0, 1/)  !B4   
        xmC(:,6)= (/1 ,  1,  1/);  xsC(:,6)=(/-1,  1, 0/)  !B5
        xmC(:,7)= (/-1,  -1, 1/);  xsC(:,7)=(/0, 1,  1/)   !C1
        xmC(:,8)= (/-1,  -1, 1/);  xsC(:,8)=(/1,  0,  1/)  !C3
        xmC(:,9)= (/-1,  -1, 1/);  xsC(:,9)=(/-1,  1,  0/) !C5
        xmC(:,10)=(/1,  -1,  1/);  xsC(:,10)=(/0,  1,  1/) !D1
        xmC(:,11)=(/1,  -1,  1/);  xsC(:,11)=(/-1,  0, 1/) !D4
        xmC(:,12)=(/1,  -1,  1/);  xsC(:,12)=(/1,  1,  0/) !D6
        xmC(:,13)=(/1,   0,  0/);  xsC(:,13)=(/0, 1, 1/) 
        xmC(:,14)=(/1,   0,  0/);  xsC(:,14)=(/0, 1, -1/) 
        xmC(:,15)=(/0,   1,  0/);  xsC(:,15)=(/1, 0, 1/) 
        xmC(:,16)=(/0,   1,  0/);  xsC(:,16)=(/1, 0, -1/)  
        xmC(:,17)=(/0,   0,  1/);  xsC(:,17)=(/1, 1, 0/)    
        xmC(:,18)=(/0,   0,  1/);  xsC(:,18)=(/1, -1, 0/)  

! extra slip systems with -ve Burgers vector (to implement backstress formulation)         
        ! xmC(:,13)= (/-1,  1,  1/);  xsC(:,13)=(/0, 1, -1/)   !A2
        ! xmC(:,14)= (/-1, 1, 1/);    xsC(:,14)=(/-1, 0, -1/)    !A3
        ! xmC(:,15)= (/-1,  1,  1/);  xsC(:,15)=(/-1, -1, 0/)    !A6
        ! xmC(:,16)= (/1,  1,  1/);   xsC(:,16)=(/0, 1, -1/)  !B2
        ! xmC(:,17)= (/1,  1,  1/);   xsC(:,17)=(/1, 0, -1/)  !B4   
        ! xmC(:,18)= (/1 , 1,  1/);   xsC(:,18)=(/1, -1, 0/)  !B5
        ! xmC(:,19)= (/-1,  -1, 1/);  xsC(:,19)=(/0, -1, -1/)   !C1
        ! xmC(:,20)= (/-1,  -1, 1/);  xsC(:,20)=(/-1, 0, -1/)  !C3
        ! xmC(:,21)= (/-1,  -1, 1/);  xsC(:,21)=(/1, -1,  0/) !C5
        ! xmC(:,22)=(/1, -1,  1/);   xsC(:,22)=(/0, -1, -1/) !D1
        ! xmC(:,23)=(/1, -1,  1/);   xsC(:,23)=(/1, 0, -1/) !D4
        ! xmC(:,24)=(/1, -1,  1/);   xsC(:,24)=(/-1,  -1,  0/) !D6

!  Normalize the vectors
      do i = 1,nslip
        xmCmag = sqrt(dot_product(xmC(:,i),xmC(:,i)))
        xmC(:,i) = xmC(:,i)/xmCmag
        xsCmag = sqrt(dot_product(xsC(:,i),xsC(:,i)))
        xsC(:,i) = xsC(:,i)/xsCmag
      end do

!------------------------------------------------------------------
!  Initialize state at t = t_n
!------------------------------------------------------------------

      if (time(2) .eq. 0.0) then ! t = 0 
        
        theta = props(1)
        eulangs(1) = th 
        eulangs(2) = phi 
        eulangs(3) = om
        rhossd_n(:) = rhossd_0
        FpInv_n = del
        spk2_n = 0.d0
        bs_n = 0.d0
      else  ! time /= 0

        theta = statev(1)

        eulangs(1:3) = statev(2:4)

        n = 4
        do j = 1,3
          do i = 1,3
            n = n + 1
            FpInv_n(i,j) = statev(n)
          end do
        end do

        spk2_n(1,1) = statev(14)
        spk2_n(2,2) = statev(15)
        spk2_n(3,3) = statev(16)
        spk2_n(1,2) = statev(17)
        spk2_n(1,3) = statev(18)
        spk2_n(2,3) = statev(19)
        
        n = 19
        do i = 1,nslip
          n = n + 1
          rhossd_n(i) = statev(n)
        end do
        do i = 1,nslip
          n = n + 1
          bs_n(i) = statev(n)
        end do

        ! n = 31
        ! do i = 1,nslip
        !   n = n + 1
        !   gamma(i) = statev(n)
        ! end do

      end if ! state initialization block

!  Initialize 4th rank elastic stiffness tensor
      call compute_elastConst(props,nrops,theta,c11,c12,c44)
      
      C_c(:,:,:,:) = 0.d0
      do i = 1,3
       do j = 1,3
        do k = 1,3
         do l = 1,3
          C_c(i,j,k,l) = C12 * del(i,j) * del(k,l) +
     &        C44 * (del(i,k)*del(j,l)+del(i,l)*del(k,j))
         end do
        end do
       end do
      end do
      C_c(1,1,1,1) = C11
      C_c(2,2,2,2) = C11
      C_c(3,3,3,3) = C11

!  Rotate stiffness tensor and slip vectors from crystal to reference configuration
      call ROTMAT(TH,PHI,OM,q_0)
      call rotate_4th(q_0,C_c,C)
      xs0 = matmul(q_0,xsC)
      xm0 = matmul(q_0,xmC)
      
!  Compute trial quantities
      Fe_tr = matmul(dfgrd1,FpInv_n)
      Ce_tr = matmul(transpose(Fe_tr),Fe_tr)
      Ee_tr = 0.5d0*(Ce_tr-del) - a11*(theta - theta0)*del
      Spk2_tr = reshape( (/((sum(C(i,j,:,:)*Ee_tr), i = 1,3), j=1,3)/), (/3,3/))

      do i = 1,3
        do j = 1,3
             Smat(i,j,:) = xs0(i,:)*xm0(j,:)
        end do
      end do

      do i = 1,3
        do j = 1,3
             Pmat(i,j,:) = 0.5d0*(Smat(i,j,:)+Smat(j,i,:))
        end do
      end do

      Bmat(:,:,:) = 0.d0
      do ia = 1,nslip
        Bmat(:,:,ia) = matmul(Ce_tr,Smat(:,:,ia)) +
     &                 matmul(transpose(Smat(:,:,ia)),Ce_tr)
      end do

      Cmat(:,:,:) = 0.d0
      do i = 1,3
        do j = 1,3
          do m = 1,3
          do n = 1,3
            Cmat(i,j,:) = Cmat(i,j,:) + 
     &                        0.5d0*C(i,j,m,n)*Bmat(m,n,:)
          end do
          end do
        end do
      end do

!  Package arrays in a manner consistent with what the solver is
!  expecting.  Solver was adapted from Bronkhorst's UMAT.
        spk2_n(2,1) = spk2_n(1,2)
        spk2_n(3,1) = spk2_n(1,3)
        spk2_n(3,2) = spk2_n(2,3) 

        call smatvec(spk2_n,spk2v_n)
        call smatvec(spk2_tr,spk2v_tr)
        do ia = 1,nslip
          call smatvec(cmat(:,:,ia),cmatr(ia,:)) 
          call smatvec(pmat(:,:,ia),pmatr(ia,:))
        end do
        pmatr(:,4) = 2.d0*pmatr(:,4) 
        pmatr(:,5) = 2.d0*pmatr(:,5) 
        pmatr(:,6) = 2.d0*pmatr(:,6)

        C11_0=(1-f_gp)*C11_0_G + f_gp*C11_0_GP
        C12_0=(1-f_gp)*C12_0_G + f_gp*C12_0_GP
        C44_0=(1-f_gp)*C44_0_G + f_gp*C44_0_GP

        xmu0 = sqrt((C11_0-C12_0)/2.d0*C44_0)
        
!  Compute slip resistances from dislocation densities at t_n
      xmu = sqrt((c11-c12)/2.d0*c44)
      call calc_amat(nslip,props,nprops,xsC,xmC,amat)
!  Compute s_inf term used to incorporate softening in Ni-based superalloys at very high temperatures 
      gamma_acc = statev(74)
!      print*, 'accumulated shear',gamma_acc
!      s_inf = -s_inf0*(theta/theta_ref + gamma_acc/gamma_ref)**soft_n
!      print*, 'softening stress',s_inf
!  Compute the flow stress dependence on the size and volume fraction of gamma prime 
!      s_inf = s_inf0*sqrt(f_gp/r_gp)

!      s_n = s_inf + rk_hp/sqrt(d_gr) + xmu*burg*sqrt(matmul(amat,rhossd_n))
      s_n = rk_hp/sqrt(d_gr) + xmu0*burg*sqrt(matmul(amat,rhossd_n))
!      print*, 'total resistance',s_n
      if (noel==20 .AND. npt==1 .AND. kinc==1 .AND. time(2)==0.0 .AND. debug==1) then
       ! if (noel == 1 .AND. npt==1 .AND. kinc==1 .AND. time(2)==0.0) then
         print*,'#######----------#######'
         print*,'#######for element#1'
         print*,'#######----------#######'
         print*,'Interaction matrix is'
         do, i=1,18
            write(*,*) ( amat(i,j), j=1,18 )
         enddo
         print*,'#######----------#######'
         print*,'#######----------#######'
         print*,'#######----------#######'
         print*,'hall-petch hardening',rk_hp/sqrt(d_gr)
         print*,'initial forest hardening',xmu*burg*sqrt(matmul(amat,rhossd_n))
         print*,'current shear modulus',xmu
      end if

!  Call solver      
      dgamma = 0.d0 

      call SOLVCRYS(nslip,props,nprops,noel,npt,theta,spk2v_n,s_n,bs_n,
     &  amat,rhossd_n,spk2v_tr,pmatr,cmatr,DTIME,spk2v,s,bs,rhossd,dgamma,
     &  DGMAX,tau,DOMPSC,ITERK,ITERL,ITERM,ITERERR,xmC)

!  Solver did not converge - cut back and repeat time step
        if (iterErr == 1) then
          pnewdt = 0.5d0
          return
        end if

!  Unpack solver quantities needed for the following state updates
!  and Jacobian calculation
        gammadot = dgamma / dtime

        if (noel==20 .AND. npt==1 .AND. debug==1) then
        ! if (noel == 1 .AND. npt==1 .AND. kinc==1 .AND. time(2)==0.0) then
                print*,'#######----------#######'
                print*,'#######for element#1'
                print*,'#######----------#######'
                ! print*,'back stress at t ',bs_n
                print*,'back stress at t+dt',bs
                ! print*,'slip resistance at t ',s_n
                print*,'slip resistance at t+dt',s
                print*,'slip rate at t+dt',gammadot
        end if
        
        call svec_to_smat(spk2v,spk2)

!  Update Fp
      rLpBar = 0.d0
      do i = 1,3
        do j = 1,3
          do k = 1,nslip
            rLpBar(i,j) = rLpBar(i,j) + 
     &          Smat(i,j,k) * gammadot(k)
          end do
        end do
      end do

      FpInv = matmul(FpInv_n,del-dtime*rLpBar)
      detFp = determinant(FpInv)
      FpInv = FpInv / detFp**(1.d0/3.d0)                                !normalizing Fp tensor for determinent to be 1

!  Compute Fe
      Fe = matmul(dfgrd1,FpInv)

!  Compute Cauchy stress
      detFe = determinant(Fe)
      sig = matmul(matmul(Fe,Spk2),transpose(Fe)) / detFe

!  Compute Jacobian - passed back in array4 in Voigt notation      
      do i = 1,3
      do j = 1,3
        smatc(:,i,j) = xsC(i,:)*xmC(j,:)
      end do
      end do

      call calc_jac(nslip,props,nprops,C,dfgrd0,dfgrd1,FpInv_n,Spk2,
     & gammadot,tau,s,bs,q_0,smatc,theta,del,dtime,array4)

!  Compute the updated temperature
      dWp = dot_product(tau,dgamma)
      dTheta = eta*dWp / (rho*cp)
      theta = theta + dTheta

!  Compute Euler angles at t = t_(n+1) for texture plotting
      call polardecomp(Fe,Re,Ue)
      Q_n1 = matmul(Re,q_0)
      call eulang(Q_n1,th_n1,phi_n1,om_n1,ieulererr)
      if (ieulererr == 1) then
        eulangs(1)  = 0.d0
        eulangs(2) =  0.d0
        eulangs(3)  = 0.d0
      else
        eulangs(1)  = th_n1*180.d0/pi
        eulangs(2) =  phi_n1*180.d0/pi
        eulangs(3)  = om_n1*180.d0/pi
      end if

!  Store the stress tensor in the ABAQUS stress 'vector'
      do i = 1,ndi
         stress(i) = sig(i,i)
      end do
      if (nshr == 1) stress(ndi+1) = sig(1,2)
      if (nshr == 3) then
         stress(4) = sig(1,2)
         stress(5) = sig(1,3)
         stress(6) = sig(2,3)
      end if
      
!  Pass back the Jacobian that Abaqus is looking for
      ddsdde = 0.0

      if (ndi == 1) then          ! 1-D
         ddsdde(1,1) = array4(1,1)
      else if (ndi == 2) then         ! 2-D plane stress & axi
         ddsdde(1:2,1:2) = array4(1:2,1:2)
         ddsdde(1,3) = array4(1,4)
         ddsdde(2,3) = array4(2,4)
         ddsdde(3,1) = array4(4,1)
         ddsdde(3,2) = array4(4,2)
         ddsdde(3,3) = array4(4,4)
      else if (ndi == 3 .and. nshr == 1) then ! plane strain
         ddsdde(1:4,1:4) = array4(1:4,1:4)
      else                  ! Full 3-D
         ddsdde(1:6,1:6) = array4(1:6,1:6)
      end if
!------------------------------------------------------------------
!  Auxillary post-processing calcuations
!------------------------------------------------------------------
! Calculate the misorientation magnitude
      arg = (Re(1,1)+Re(2,2)+Re(3,3)-1.d0)/2.d0
      if (abs(arg-1.d0) < 1e-10) then
       angMis = 0.d0
      else
       angMis = dacos(arg)
      end if
      angMis = angMis*180.d0/pi

!------------------------------------------------------------------
!  Store the state variables
!------------------------------------------------------------------
         n = 0
!  Store the temperature
         n = n + 1
         statev(n) = theta

!  Store the Euler Angles
         do i = 1,3
           n = n + 1
           statev(n) = eulangs(i)
         end do
            
!  Store {F^p}^{-1}
            do j = 1,3
               do i = 1,3
                 n = n + 1
                 statev(n) = FpInv(i,j)
               end do
            end do

!  Store the PK2 stress            
        statev(14)=spk2(1,1)
        statev(15)=spk2(2,2)
        statev(16)=spk2(3,3)
        statev(17)=spk2(1,2)
        statev(18)=spk2(1,3)
        statev(19)=spk2(2,3)

!  Store ssd_density^alpha from 20 to 31 for 12 slip system whereas 20 to 37 for 18 slip systems
           n = 19
           do i = 1,nslip
             n = n + 1
             statev(n) = rhossd(i)
           end do
!  Store backstress^alpha from 32 to 43 for 12 slip system whereas 38 to 55 for 18 slip systems
           do i = 1,nslip
             n = n + 1
             statev(n) = bs(i)
           end do

!  Store the cummulative slips, gamma^alpha from 44 to 55 for 12 slip systems whereas 56 to 73 for 18 slip systems
           cumslip_all = 0.d0
           do i = 1,nslip
             n = n + 1
             statev(n) = statev(n) + abs(gammadot(i)*dtime)
             cumslip_all = cumslip_all + statev(n)
           end do

!  Store the cumulative slip on all systems in 56 for 12 slip whereas in statev 74 for 18 slip systems 
           n = n + 1
           statev(n) = cumslip_all

!  Store the misorientation angle
           n = n + 1
           statev(n) = angMis
           n = n + 1
           statev(n) = s_inf

      return
      end
     
!===================================================================
!===================================================================
!===================== S U B R O U T I N E S =======================
!===================================================================
!===================================================================
!  Store only select outputs from the statev array 
!-------------------------------------------------------------------
      subroutine uvarm(uvar,direct,t,time,dtime,cmname,orname,
     & nuvarm,noel,npt,layer,kspt,kstep,kinc,ndi,nshr,coord,
     & jmac,jmatyp,matlayo,laccfla)
      
      implicit real*8(a-h,o-z)
c
      character*80 cmname,orname
      character*3 flgray(69)
      dimension uvar(nuvarm),direct(3,3),t(3,3),time(2)
      dimension array(69),jarray(69),jmac(*),jmatyp(*),coord(*)

      call getvrm('SDV',array,jarray,flgray,jrcd,jmac,jmatyp,matlayo,
     &            laccfla)

      if (cmname(1:2) == 'CU') then
        uvar(1)    = array(1)
        uvar(2)    = array(44)
        uvar(3)    = array(45)
        uvar(4:6)  = array(2:4) 
        uvar(7:18) = array(20:31)
      end if

      return
      end subroutine
!===================================================================
!  Convert from Voigt notation stress vector to a symmetric 
!-------------------------------------------------------------------
      subroutine svec_to_smat(stressv,stressm)

      implicit none

      real*8 :: stressv(6),stressm(3,3)

      stressm(1,1) = stressv(1)
      stressm(2,2) = stressv(2)
      stressm(3,3) = stressv(3)
      stressm(1,2) = stressv(4)
      stressm(2,1) = stressv(4)
      stressm(1,3) = stressv(5)
      stressm(3,1) = stressv(5)
      stressm(2,3) = stressv(6)
      stressm(3,2) = stressv(6)

      return
      end subroutine
!===================================================================
!===================================================================
!
!  Rotates any 3x3x3x3 tensor by a rotation matrix.
!
!  c(i,j,k,l) = a(i,m) * a(j,n) * a(k,p) * a(l,q) * b(m,n,p,q)
!
!-------------------------------------------------------------------

      subroutine rotate_4th(a,b,c)

      implicit double precision (a-h,o-z)

      dimension a(3,3), b(3,3,3,3), c(3,3,3,3), d(3,3,3,3)

      do m = 1,3
       do n = 1,3
        do k = 1,3
         do l = 1,3
          d(m,n,k,l) = a(k,1) * (a(l,1) * b(m,n,1,1) + 
     &      a(l,2) * b(m,n,1,2) + a(l,3) * b(m,n,1,3)) +
     &      a(k,2) * (a(l,1) * b(m,n,2,1) + 
     &      a(l,2) * b(m,n,2,2) + a(l,3) * b(m,n,2,3)) +
     &      a(k,3) * (a(l,1) * b(m,n,3,1) + 
     &      a(l,2) * b(m,n,3,2) + a(l,3) * b(m,n,3,3))
         end do
        end do
       end do
      end do

      do i = 1,3
       do j = 1,3
        do k = 1,3
         do l = 1,3
          c(i,j,k,l) = a(i,1) * (a(j,1) * d(1,1,k,l) + 
     &      a(j,2) * d(1,2,k,l) + a(j,3) * d(1,3,k,l)) +
     &      a(i,2) * (a(j,1) * d(2,1,k,l) + 
     &      a(j,2) * d(2,2,k,l) + a(j,3) * d(2,3,k,l)) +
     &      a(i,3) * (a(j,1) * d(3,1,k,l) + 
     &      a(j,2) * d(3,2,k,l) + a(j,3) * d(3,3,k,l))
         end do
        end do
       end do
      end do

      return
      end

!===================================================================
!===================================================================
!
!  Calculate the inverse of a 3 x 3 matrix.
!
!-------------------------------------------------------------------

      subroutine inverse_3x3(a,b)

      implicit double precision (a-h,o-z)

      dimension a(3,3), b(3,3)

      b(1,1) = a(2,2) * a(3,3) - a(3,2) * a(2,3)
      b(1,2) = a(3,2) * a(1,3) - a(1,2) * a(3,3)
      b(1,3) = a(1,2) * a(2,3) - a(2,2) * a(1,3)
      b(2,1) = a(3,1) * a(2,3) - a(2,1) * a(3,3)
      b(2,2) = a(1,1) * a(3,3) - a(3,1) * a(1,3)
      b(2,3) = a(2,1) * a(1,3) - a(1,1) * a(2,3)
      b(3,1) = a(2,1) * a(3,2) - a(3,1) * a(2,2)
      b(3,2) = a(3,1) * a(1,2) - a(1,1) * a(3,2)
      b(3,3) = a(1,1) * a(2,2) - a(2,1) * a(1,2)

      det = a(1,1) * b(1,1) + a(1,2) * b(2,1) + a(1,3) * b(3,1)

      do i = 1,3
         do j = 1,3
            b(i,j) = b(i,j) / det
         end do
      end do

      return
      end
!===================================================================
!  Calculate the determinant of a 3 x 3 matrix.
!

      function determinant(a)

      implicit double precision (a-h,o-z)

      dimension a(3,3)

      b1 = a(2,2) * a(3,3) - a(3,2) * a(2,3)
      b2 = a(3,1) * a(2,3) - a(2,1) * a(3,3)
      b3 = a(2,1) * a(3,2) - a(3,1) * a(2,2)

      determinant = a(1,1) * b1 + a(1,2) * b2 + a(1,3) * b3

      return
      end



!===================================================================
!===================================================================
!======================================================================
!***********************************************************************
        SUBROUTINE ROTMAT(TH,PHI,OM,Q)
        IMPLICIT REAL*8(A-H,O-Z)
        DIMENSION Q(3,3)

        STH = DSIN(TH)
        CTH = DCOS(TH)
        SPH = DSIN(PHI)
        CPH = DCOS(PHI)
        SOM = DSIN(OM)
        COM = DCOS(OM)

        Q(1,1) = CPH*COM-SPH*SOM*CTH
        Q(1,2) = SPH*COM+CPH*SOM*CTH
        Q(1,3) = SOM*STH

        Q(2,1) = -CPH*SOM-SPH*COM*CTH
        Q(2,2) = -SPH*SOM+CPH*COM*CTH
        Q(2,3) =  COM*STH

        Q(3,1) = SPH*STH
        Q(3,2) = -CPH*STH
        Q(3,3) = CTH

        RETURN
        END
      
!===================================================================
        SUBROUTINE EULANG(Q,TH,PHI,OM,IEULERERR)
        IMPLICIT REAL*8(A-H,O-Z)
        DIMENSION Q(3,3)
        PI = 4.0*DATAN(1.0D0)
        ICHK=0
        IEULERERR = 0
        IF(DABS(Q(3,3))-1.0 > 1.0D-6)THEN
           IEULERERR = 1
           RETURN
        END IF
        DO I = 1,3
        DO J = 1,3
          IF(DABS(Q(I,J)) < 1.0D-6) Q(I,J) = 0.0
        END DO
        END DO
        IF(DABS(DABS(Q(3,3))-1.0).LT.1.0D-6)THEN
          CALL EULCHECK1(Q,TH,PHI,OM,ICHK)
          IF(ICHK /= 1)GO TO 20
          RETURN
        END IF
        TH = DACOS(Q(3,3))
        STH = DSIN(TH)
        OM = DATAN2(Q(1,3)/STH,Q(2,3)/STH)
        PHI = DATAN2(Q(3,1)/STH,-Q(3,2)/STH)
        CALL EULCHECK(Q,TH,PHI,OM,ICHK)
        IF(ICHK == 1)RETURN
        TH = 2.0*PI-TH
        STH = DSIN(TH)
        OM = DATAN2(Q(1,3)/STH,Q(2,3)/STH)
        PHI = DATAN2(Q(3,1)/STH,-Q(3,2)/STH)
        CALL EULCHECK(Q,TH,PHI,OM,ICHK)
20      IF(ICHK /= 1)THEN
            IEULERERR = 1
            RETURN
        END IF

        RETURN
        END
!===================================================================
        SUBROUTINE EULCHECK(Q,TH,PHI,OM,ICHK)

        IMPLICIT REAL*8(A-H,O-Z)

        DIMENSION Q(3,3)

        TOL=1.0D-3
        ICHK = 0

        A = DCOS(PHI)*DCOS(OM)-DSIN(PHI)*DSIN(OM)*DCOS(TH)
        B = -DSIN(OM)*DCOS(PHI)-DCOS(OM)*DSIN(PHI)*DCOS(TH)
        C = DCOS(OM)*DSIN(PHI)+DSIN(OM)*DCOS(PHI)*DCOS(TH)
        D = -DSIN(PHI)*DSIN(OM)+DCOS(PHI)*DCOS(OM)*DCOS(TH)

        IF(DABS(A-Q(1,1)).LT.TOL.AND.DABS(B-Q(2,1)).LT.TOL.
     &    AND.DABS(C-Q(1,2)).LT.TOL.AND.DABS(D-Q(2,2)).LT.TOL)ICHK=1

        RETURN
        END
C **********************************************************************
        SUBROUTINE EULCHECK1(Q,TH,PHI,OM,ICHK)

        IMPLICIT REAL*8(A-H,O-Z)

        DIMENSION Q(3,3)

        TOL=1.0D-3
        ICHK = 0

        Q(3,3) = 1.0*Q(3,3)/DABS(Q(3,3))
        TH = DACOS(Q(3,3))
        IF(DABS(Q(1,3)).GT.TOL)RETURN
        IF(DABS(Q(2,3)).GT.TOL)RETURN
        IF(DABS(Q(3,1)).GT.TOL)RETURN
        IF(DABS(Q(3,2)).GT.TOL)RETURN
        IF(Q(3,3).EQ.1.0)THEN
           IF(DABS(Q(1,1)-Q(2,2)).GT.TOL)RETURN
           IF(DABS(Q(1,2)+Q(2,1)).GT.TOL)RETURN
        ELSE
           IF(DABS(Q(1,1)+Q(2,2)).GT.TOL)RETURN
           IF(DABS(Q(1,2)-Q(2,1)).GT.TOL)RETURN
        END IF
        PHI = DATAN2(Q(1,2),Q(1,1))
        OM = 0.0
        ICHK = 1

        RETURN
        END
! ********************************************************************
!
!  THIS SUBROUTINE PERFORMS THE RIGHT POLAR DECOMPOSITION
!  [F]=[R][U] OF THE DEFORMATION GRADIENT [F] INTO
!  A ROTATION [R] AND THE RIGHT  STRETCH TENSOR [U].
!  THE EIGENVALUES AND EIGENVECTORS OF [U] AND
!  THE  LOGARITHMIC STRAIN [E] = LN [U]
!  ARE ALSO RETURNED.
!       
!
! *******************************************************************
        SUBROUTINE polardecomp(F,R,U)

        IMPLICIT REAL*8(A-H,O-Z)
!
! ------------------------------------------------------------------
!  VARIABLES
!
        REAL*8 F(3,3),DETF,FTRANS(3,3),
     +   C(3,3), OMEGA(3),
     +   UEIGVAL(3),EIGVEC(3,3), EIGVECT(3,3), 
     +   U(3,3),E(3,3),UINV(3,3),R(3,3),TEMP(3,3)
         
!
!       F(3,3)  -- THE DEFORMATION GRADIENT MATRIX WHOSE
!                          POLAR DECOMPOSITION IS DESIRED.
!       DETF            -- THE DETRMINANT OF [F]; DETF > 0.
!       FTRANS(3,3)     -- THE TRANSPOSE OF [F].
!       R(3,3)  -- THE ROTATION MATRIX; [R]^T [R] = [I];
!                          OUTPUT.
!       U(3,3)  -- THE RIGHT STRETCH TENSOR; SYMMETRIC
!                          AND POSITIVE DEFINITE; OUTPUT.
!       UINV(3,3)       -- THE INVERSE OF [U].
!       C(3,3)  -- THE RIGHT CAUCHY-GREEN TENSOR = [U][U];
!                           SYMMETRIC AND POSITIVE DEFINITE.
!       OMEGA(3)        -- THE SQUARES OF THE PRINCIPAL STRETCHES.
!       UEIGVAL(3)      -- THE PRINCIPAL STRETCHES; OUTPUT.
!       EIGVEC(3,3)     -- MATRIX OF EIGENVECTORS OF [U];OUTPUT.
!       EIGVECT(3,3)-- TRANSPOSE OF THE ABOVE.
!       E(3,3)  -- THE LOGARITHMIC STRAIN TENSOR, [E]=LN[U];
!                          OUTPUT.
! -----------------------------------------------------------------
!       COMPUTATION
       
!  STORE THE IDENTITY MATRIX IN  [R], [U], AND [UINV]

        r = reshape( (/1.,0.,0.,0.,1.,0.,0.,0.,1. /), (/3,3/) )
        u = reshape( (/1.,0.,0.,0.,1.,0.,0.,0.,1. /), (/3,3/) )
        uinv = reshape( (/1.,0.,0.,0.,1.,0.,0.,0.,1. /), (/3,3/) )
        
        

!  STORE THE ZERO MATRIX IN [E]

        E(:,:) = 0.d0

!  CHECK IF THE DETERMINANT OF [F] IS GREATER THAN ZERO.
!  IF NOT, THEN PRINT DIAGONOSTIC AND STOP.

!       detF = determinant(F)
!       IF (DETF .LE. 0.D0) THEN
!          WRITE(91,100)
!          RETURN
!       END IF

!  CALCULATE THE RIGHT CAUCHY GREEN STRAIN TENSOR [C]

       ftrans = transpose(f)
       c = matmul(ftrans,f)
 
!  CALCULATE THE EIGENVALUES AND EIGENVECTORS OF  [C]

         CALL SPECTRAL(C,OMEGA,EIGVEC)

!  CALCULATE THE PRINCIPAL VALUES OF [U] AND [E]

         UEIGVAL(1) = DSQRT(OMEGA(1))
         UEIGVAL(2) = DSQRT(OMEGA(2))
         UEIGVAL(3) = DSQRT(OMEGA(3))

         U(1,1) = UEIGVAL(1)
         U(2,2) = UEIGVAL(2)
         U(3,3) = UEIGVAL(3)

         E(1,1) = DLOG( UEIGVAL(1) )
         E(2,2) = DLOG( UEIGVAL(2) )
         E(3,3) = DLOG( UEIGVAL(3) )

!  CALCULATE THE COMPLETE TENSORS [U] AND [E]
        eigvecT = transpose(eigvec)
        temp = matmul(eigvec,u)
        u = matmul(temp,eigvecT)
        temp = matmul(eigvec,e)
        e = matmul(temp,eigvecT)

!  CALCULATE [UINV]
        call inverse_3x3(U,Uinv)

!  CACULATE [R]
         r = matmul(f,uinv)
! -------------------------------------------------------------------
!       FORMATS

 100  FORMAT(5X,'--ERROR IN KINEMATICS-- THE DETERMINANT OF [F]',
     +       ' IS NOT GREATER THAN 0',/,
     +       10X,'POLAR DECOMPOSITION NOT PERFORMED')
! -------------------------------------------------------------------
        RETURN
        END
! ****************************************************************
!
!  THE FOLLOWING SUBROUTINES CALCULATE THE SPECTRAL
!  DECOMPOSITION OF A SYMMETRIC THREE BY THREE MATRIX
!
! **********************************************************************
        SUBROUTINE SPECTRAL(A,D,V)
!
!  THIS SUBROUTINE CALCULATES THE EIGENVALUES AND EIGENVECTORS OF
!  A SYMMETRIC 3 BY 3 MATRIX [A]. 
!
!  THE OUTPUT CONSISTS OF A VECTOR D CONTAINING THE THREE
!  EIGENVALUES IN ASCENDING ORDER, AND
!  A MATRIX [V] WHOSE COLUMNS CONTAIN THE CORRESPONDING
!  EIGENVECTORS.
! ---------------------------------------------------------------------
        IMPLICIT REAL*8 (A-H,O-Z)
        PARAMETER(NP=3)
        DIMENSION D(NP),V(NP,NP)
        DIMENSION A(3,3),E(NP,NP)
        DO I = 1,3
          DO J= 1,3
            E(I,J) = A(I,J)
          END DO
        END DO
        CALL JACOBI(E,3,NP,D,V,NROT)
        CALL EIGSRT(D,V,3,NP)
        RETURN
        END
! *********************************************************************
        SUBROUTINE JACOBI(A,N,NP,D,V,NROT)

!  COMPUTES ALL EIGENVALUES AND EIGENVECTORS OF A REAL SYMMETRIC
!  MATRIX [A], WHICH IS OF SIZE N BY N, STORED IN A PHYSICAL 
!  NP BY BP ARRAY. ON OUTPUT, ELEMENTS OF [A] ABOVE THE DIAGONAL 
!  ARE DESTROYED, BUT THE DIAGONAL AND SUB-DIAGONAL ARE UNCHANGED
!  AND GIVE FULL INFORMATION ABOUT THE ORIGINAL SYMMETRIC MATRIX.
!  VECTOR D RETURNS THE EIGENVALUES OF [A] IN ITS FIRST N ELEMENTS.
!  [V] IS A MATRIX WITH THE SAME LOGICAL AND PHYSICAL DIMENSIONS AS
!  [A] WHOSE COLUMNS CONTAIN, ON OUTPUT, THE NORMALIZED
!  EIGENVECTORSOF [A]. NROT RETURNS THE NUMBER OF JACOBI ROTATIONS
!  WHICH WERE REQUIRED.
!
!  THIS SUBROUTINE IS TAKEN FROM "NUMERICAL RECIPES", PAGE 346.
 
        IMPLICIT REAL*8 (A-H,O-Z)
        PARAMETER (NMAX =100)
        DIMENSION A(NP,NP),D(NP),V(NP,NP),B(NMAX),Z(NMAX)

!  INITIALIZE [V] TO THE IDENTITY MATRIX

        DO IP = 1,N     
          DO IQ = 1,N
            V(IP,IQ) = 0.D0
          END DO
          V(IP,IP) = 1.D0
        END DO

!  INITIALIZE B AND D TO THE DIAGONAL OF [A], AND Z TO ZERO. THE
!  VECTOR Z WILL ACCUMULATE TERMS OF THE FORM T*A_PQ AS
!  IN EQUATION (11.1.14)

        DO IP = 1,N
          B(IP) = A(IP,IP)
          D(IP) = B(IP)
          Z(IP) = 0.D0
        END DO
        NROT = 0
        DO I = 1,50

!  SUM OFF-DIAGONAL ELEMENTS

        SM = 0.D0
        DO IP = 1, N-1
          DO IQ = IP + 1, N
            SM = SM + DABS ( A(IP,IQ ))
          END DO
        END DO

!  IF SUM = 0., THEN RETURN. THIS IS THE NORMAL RETURN
!  WHICH RELIES ON QUADRATIC CONVERGENCE TO MACHINE 
!  UNDERFLOW.

        IF( SM .EQ. 0.D0) RETURN

!  IF( SM .LT. 1.0D-15) RETURN

!  IN THE FIRST THREE SWEEPS CARRY OUT THE PQ ROTATION ONLY IF
!  |A_PQ| > TRESH, WHERE TRESH IS SOME THRESHOLD VALUE, C
!  SEE EQUATION (11.1.25).  THEREAFTER TRESH = 0.

        IF( I .LT. 4) THEN
           TRESH = 0.2D0*SM/N**2
        ELSE
           TRESH = 0.D0
        END IF
        DO IP = 1, N-1
          DO IQ = IP+1,N
            G = 100.D0*DABS(A(IP,IQ))

!  AFTER FOUR SWEEPS, SKIP THE ROTATION IF THE
!  OFF-DIAGONAL ELEMENT IS SMALL.

        IF( (I .GT. 4) .AND. ( DABS(D(IP))+G  .EQ.  DABS( D(IP)) )
     +               .AND. ( DABS(D(IQ))+G  .EQ.  DABS( D(IQ)) ) )THEN
                       A(IP,IQ) = 0.D0
                   ELSE IF( DABS(A(IP,IQ)) .GT. TRESH) THEN
                       H = D(IQ) - D(IP)
                       IF(DABS(H)+G .EQ. DABS(H)) THEN

!                               T= 1./(2.*THETA), EQUATION(11.1.10)

                                T =A(IP,IQ)/H
                              ELSE
                                  THETA = 0.5D0*H/A(IP,IQ)
                             T =1.D0/(DABS(THETA)+DSQRT(1.D0+THETA**2))
                                  IF(THETA .LT. 0.D0) T = -T
                              END IF
                              C = 1.D0/DSQRT(1.D0 + T**2)
                              S = T*C
                              TAU = S/(1.D0 + C)
                                H = T*A(IP,IQ)
                                Z(IP) = Z(IP) - H
                                Z(IQ) = Z(IQ) + H
                                D(IP) = D(IP) - H
                                D(IQ) = D(IQ) + H
                                A(IP,IQ) = 0.D0

!                               CASE OF ROTATIONS 1<= J < P
                               
                                DO J = 1, IP-1
                                     G = A(J,IP)
                                     H = A(J,IQ)
                                     A(J,IP) = G - S*(H + G*TAU)
                                     A(J,IQ) = H + S*(G - H*TAU)
                                END DO

!                               CASE OF ROTATIONS P < J < Q

                                DO J = IP+1, IQ-1
                                     G = A(IP,J)
                                     H = A(J,IQ)
                                     A(IP,J) = G - S*(H + G*TAU)
                                     A(J,IQ) = H + S*(G - H*TAU)
                                END DO

!                               CASE OF ROTATIONS Q < J <= N

                                DO J = IQ+1, N
                                     G = A(IP,J)
                                     H = A(IQ,J)
                                     A(IP,J) = G - S*(H + G*TAU)
                                     A(IQ,J) = H + S*(G - H*TAU)
                                END DO
                                DO J = 1,N
                                     G = V(J,IP)
                                     H = V(J,IQ)
                                     V(J,IP) = G - S*(H + G*TAU)
                                     V(J,IQ) = H + S*(G - H*TAU)
                                END DO
                            NROT = NROT + 1
                   END IF
                   END DO
            END DO

!  UPDATE D WITH THE SUM OF T*A_PQ, AND REINITIALIZE Z

          DO IP = 1, N
               B(IP) = B(IP) + Z(IP)
               D(IP) = B(IP)
               Z(IP) = 0.D0
          END DO
        END DO

!  IF THE ALGORITHM HAS REACHED THIS STAGE, THEN
!  THERE ARE TOO MANY SWEEPS, PRINT A DIAGNOSTIC
!  AND STOP

!        WRITE(*,'(/1X,A/)') '50 ITERS IN JACOBI SHOULD NEVER HAPPEN'
        RETURN
        END
! **********************************************************************
        SUBROUTINE EIGSRT(D,V,N,NP)
!
!       GIVEN THE EIGENVALUES D AND EIGENVECTORS [V] AS OUTPUT FROM
!       JACOBI, THIS ROUTINE SORTS THE EIGENVALUES INTO ASCENDING ORDER, 
!       AND REARRANGES THE COLUMNS OF [V] ACCORDINGLY.
!
!       THIS SUBROUTINE IS TAKEN FROM "NUMERICAL RECIPES", P. 348.
! ----------------------------------------------------------------------
        IMPLICIT REAL*8 (A-H,O-Z)
        DIMENSION D(NP),V(NP,NP)
        DO I = 1,N-1
                K = I
                P = D(I)
                DO J = I+1,N
                     IF(D(J) .GE. P) THEN
                          K = J
                          P = D(J)
                     END IF
                END DO
                IF(K .NE. I) THEN
                   D(K) = D(I)
                   D(I) = P
                 DO J = 1,N
                        P = V(J,I)
                        V(J,I) = V(J,K)
                        V(J,K) = P
                 END DO
                END IF
        END DO
        RETURN
        END

! **********************************************************************
      SUBROUTINE EMATR(ELMAT,ELMATR)
C      
C  THIS SUBROUTINE REDUCES AN ELASTIC STIFFNESS TYPE MATRIX 
C  3 X 3 X 3 X 3 MATRIX TO A 6 X 6 MATRIX
C          
      IMPLICIT REAL*8(A-H,O-Z) 
      DIMENSION ELMAT(3,3,3,3), ELMATR(6,6)
C      
        DO I = 1,6
        DO J = 1,6
          ELMATR(I,J) = 0.0
        END DO
        END DO
C 
        ELMATR(1,1)  = ELMAT(1,1,1,1)
        ELMATR(1,2)  = ELMAT(1,1,2,2)
        ELMATR(1,3)  = ELMAT(1,1,3,3) 
        ELMATR(1,4)  = 0.5 * ( ELMAT(1,1,1,2) + ELMAT(1,1,2,1) )
        ELMATR(1,5)  = 0.5 * ( ELMAT(1,1,1,3) + ELMAT(1,1,3,1) )
        ELMATR(1,6)  = 0.5 * ( ELMAT(1,1,2,3) + ELMAT(1,1,3,2) )
C
        ELMATR(2,1)  = ELMAT(2,2,1,1)
        ELMATR(2,2)  = ELMAT(2,2,2,2)
        ELMATR(2,3)  = ELMAT(2,2,3,3) 
        ELMATR(2,4)  = 0.5 * ( ELMAT(2,2,1,2) + ELMAT(2,2,2,1) )
        ELMATR(2,5)  = 0.5 * ( ELMAT(2,2,1,3) + ELMAT(2,2,3,1) )
        ELMATR(2,6)  = 0.5 * ( ELMAT(2,2,2,3) + ELMAT(2,2,3,2) )
C       
        ELMATR(3,1)  = ELMAT(3,3,1,1)
        ELMATR(3,2)  = ELMAT(3,3,2,2)
        ELMATR(3,3)  = ELMAT(3,3,3,3) 
        ELMATR(3,4)  = 0.5 * ( ELMAT(3,3,1,2) + ELMAT(3,3,2,1) )
        ELMATR(3,5)  = 0.5 * ( ELMAT(3,3,1,3) + ELMAT(3,3,3,1) )
        ELMATR(3,6)  = 0.5 * ( ELMAT(3,3,2,3) + ELMAT(3,3,3,2) )   
C
        ELMATR(4,1)  = ELMAT(1,2,1,1)
        ELMATR(4,2)  = ELMAT(1,2,2,2)
        ELMATR(4,3)  = ELMAT(1,2,3,3) 
        ELMATR(4,4)  = 0.5 * ( ELMAT(1,2,1,2) + ELMAT(1,2,2,1) )
        ELMATR(4,5)  = 0.5 * ( ELMAT(1,2,1,3) + ELMAT(1,2,3,1) )
        ELMATR(4,6)  = 0.5 * ( ELMAT(1,2,2,3) + ELMAT(1,2,3,2) )        
C
        ELMATR(5,1)  = ELMAT(1,3,1,1)
        ELMATR(5,2)  = ELMAT(1,3,2,2)
        ELMATR(5,3)  = ELMAT(1,3,3,3) 
        ELMATR(5,4)  = 0.5 * ( ELMAT(1,3,1,2) + ELMAT(1,3,2,1) )
        ELMATR(5,5)  = 0.5 * ( ELMAT(1,3,1,3) + ELMAT(1,3,3,1) )
        ELMATR(5,6)  = 0.5 * ( ELMAT(1,3,2,3) + ELMAT(1,3,3,2) ) 
C
        ELMATR(6,1)  = ELMAT(2,3,1,1)
        ELMATR(6,2)  = ELMAT(2,3,2,2)
        ELMATR(6,3)  = ELMAT(2,3,3,3) 
        ELMATR(6,4)  = 0.5 * ( ELMAT(2,3,1,2) + ELMAT(2,3,2,1) )
        ELMATR(6,5)  = 0.5 * ( ELMAT(2,3,1,3) + ELMAT(2,3,3,1) )
        ELMATR(6,6)  = 0.5 * ( ELMAT(2,3,2,3) + ELMAT(2,3,3,2) ) 
        RETURN
        END
!///////////////////////////////////////////////////////////////////////
!  Subroutine for Calculating the Const. Jacobian
!///////////////////////////////////////////////////////////////////////
       SUBROUTINE calc_jac(nslip,props,nprops,elmat,ft,ftau,fpit,tstartau,
     & gammadot,rss,stau,btau,q,smatc,theta,del,DTIME,wjmatr)

        IMPLICIT REAL*8(A-H,O-Z)     
        
        DIMENSION  FTAU(3,3)
        DIMENSION  FPIT(3,3)
        DIMENSION  Q(3,3),QT(3,3),SMATG(nslip,3,3)
        DIMENSION  AMAT(3,3)
        DIMENSION  SALPHA(3,3),PALPHA(3,3),PALPHAV(6),PMAT(nslip,6)
        DIMENSION  CALPHA(3,3),CALPHAV(6),CMAT(nslip,6)
        DIMENSION  TSTARTAU(3,3),STAU(nslip)
        DIMENSION  RSS(nslip),DGAMMA(nslip),BTAU(nslip)
        DIMENSION  FPITAU(3,3),FSTAU(3,3)
        DIMENSION  AUX1(3,3),AUX2(3,3)
C            arrays for analytical jacobian
        DIMENSION  FSTAUINV(3,3)
        DIMENSION  FT(3,3),FTTAU(3,3),RTTAU(3,3),UTTAU(3,3)
        DIMENSION  FSTART(3,3),CJMAT(3,3,3,3),ELMAT(3,3,3,3)
        DIMENSION  DJMAT(3,3,3,3),DJMATR(6,6)
        DIMENSION  GJMAT(nslip,3,3,3,3)
        DIMENSION  AJJMAT(nslip,3,3,3,3)
        DIMENSION  AJJMATR(nslip,6,6)
        DIMENSION  BJMAT(nslip,3,3),BJMATR(nslip,6)
        DIMENSION  AKJMATR(6,6),AKJMATRINV(6,6)
        DIMENSION  AMJMATR(6,6)
        DIMENSION  QJMATR(6,6),QJMAT(3,3,3,3)
        DIMENSION  RJMAT(nslip,3,3)
        DIMENSION  SJMAT(3,3,3,3),WJMAT(3,3,3,3),WJMATR(6,6)
        DIMENSION  INDX(6)
        DIMENSION  AUXMAT1(3,3,3,3),AUXMAT1R(6,6),AUXMAT2(3,3,3,3)
        DIMENSION  AUXV1(6)

        dimension  gammadot(nslip),smatc(nslip,3,3),del(3,3)
        dimension  props(nprops)
C
        f_gp=props(39)
C----------------------------------------------------------------
C             CALCULATE THE JACOBIAN
C----------------------------------------------------------------

!  Intitialization
        dgamma = gammadot*dtime

C  STEP 1. CALCULATE THE MATRIX AMAT
C
        aux1 = matmul(ftau,fpit)
        aux2 = transpose(aux1)
        amat = matmul(aux2,aux1)
        
        DO ISLIP = 1,NSLIP
C
C  GET S MATRIX IN LOCAL COORDINATE SYSTEM
C
        DO I=1,3
        DO J=1,3
          SALPHA(I,J) = SMATC(ISLIP,I,J)
          SMATG(ISLIP,I,J) = 0.0
        END DO
        END DO
C
C  CONVERT S MATRIX TO GLOBAL COORDINATE SYSTEM
        salpha = matmul(q,matmul(salpha,transpose(q)))
        DO I = 1,3
        DO J = 1,3                 
          SMATG(ISLIP,I,J) =SALPHA(I,J)
        END DO
        END DO              
C
        DO I = 1,3
        DO J = 1,3
          PALPHA(I,J) = 0.5*(SALPHA(I,J)+SALPHA(J,I))
        END DO
        END DO
C
        CALL SMATVEC(PALPHA,PALPHAV)
C
C  MULTIPLY THE SHEAR COMPONENTS OF PALPHAV
C  BY A FACTOR OF 2 TO BE CONSISTENT WITH
C  THE DEFINITION OF ENGINEERING SHEAR STRAINS
C
        PALPHAV(4) = 2.0*PALPHAV(4)
        PALPHAV(5) = 2.0*PALPHAV(5)
        PALPHAV(6) = 2.0*PALPHAV(6)
        DO I=1,6
          PMAT(ISLIP,I) = PALPHAV(I)
        END DO
C
        aux1 = matmul(amat,salpha)
        DO I = 1,3
        DO J = 1,3
          AUX2(I,J) = 0.5*(AUX1(I,J) + AUX1(J,I))
        END DO
        END DO

        calpha = reshape( (/((sum(elmat(i,j,:,:)*aux2), i = 1,3), j=1,3)/), (/3,3/))

        CALL SMATVEC(CALPHA,CALPHAV)
        DO I=1,6
          CMAT(ISLIP,I) = CALPHAV(I)
        END DO

        END DO !ISLIP

C  STEP 5. CALCULATE FPI AT TME TAU
C          
        DO I =1,3
        DO J =1,3
          AUX1(I,J)= 0.0
          DO ISLIP=1,NSLIP
            AUX1(I,J)= AUX1(I,J) + 
     &                 DGAMMA(ISLIP)*SMATG(ISLIP,I,J)                   
          END DO
          AUX2(I,J) = del(I,J) - AUX1(I,J)
        END DO
        END DO    

        fpitau = matmul(fpit,aux2)

        det = determinant(fpitau)
        CDET = DET**(1./3.)
        DO I =1,3
        DO J=1,3
          FPITAU(I,J) = FPITAU(I,J)/CDET
        END DO
        END DO
C
C  STEP 6. CALCULATE FSTAU 
        fstau = matmul(ftau,fpitau)

!////////////////////////////////////////////////////////////////
C
C  STEP J1. CALCULATE THE C_IJKL MATRIX
C
        CALL inverse_3x3(FT,AUX1)           

        fttau = matmul(ftau,aux1)

        CALL polardecomp(FTTAU,RTTAU,UTTAU)

        fstart = matmul(ft,fpit)
        DO I = 1,3
        DO J = 1,3
        DO K = 1,3
        DO L = 1,3
          SUM1 = 0.0
          SUM2 = 0.0
          DO M = 1,3
            SUM1 = SUM1 + UTTAU(L,M)  * FSTART(M,J)
            SUM2 = SUM2 + FSTART(M,I) * UTTAU(M,K)
          END DO       
          CJMAT(I,J,K,L) = FSTART(K,I) * SUM1 +
     &                     SUM2 * FSTART(L,J)
        END DO
        END DO
        END DO
        END DO       
   
C  STEP J2. CALCULATE THE L_IJKL MATRIX IN GLOBAL
C           CORDINATES

!  ELMAT is now passed in through the calc_jac subroutine
!  interface
        
             
C  STEP J3. CALCULATE THE D_IJKL MATRIX

        DO I = 1,3
        DO J = 1,3
        DO K = 1,3
        DO L = 1,3
          DJMAT(I,J,K,L) = 0.0
          DO M  = 1,3
          DO N  = 1,3
            DJMAT(I,J,K,L) = DJMAT(I,J,K,L) +
     &                       0.5*ELMAT(I,J,M,N)*CJMAT(M,N,K,L)
          END DO
          END DO
        END DO
        END DO
        END DO
        END DO
C
        CALL EMATR(DJMAT,DJMATR)     
C
C  STEP J4. CALCULATE THE NLIP MATRICES G_MNKL^{\ALPHA}
C
        DO ISLIP = 1,NSLIP
          DO M = 1,3
          DO N = 1,3
          DO K = 1,3
          DO L = 1,3
            GJMAT(ISLIP,M,N,K,L) = 0.0
            DO IP = 1,3
              GJMAT(ISLIP,M,N,K,L) = GJMAT(ISLIP,M,N,K,L) +
     &                     CJMAT(M,IP,K,L)*SMATG(ISLIP,IP,N) +
     &                     SMATG(ISLIP,IP,M)*CJMAT(IP,N,K,L)
            END DO
          END DO
          END DO
          END DO
          END DO
        END DO
C  
C  STEP J5. CALCULATE THE NSLIP MATRICES J_IJKL^{\ALPHA}
C            
        DO ISLIP = 1,NSLIP
          DO I = 1,3
          DO J = 1,3
          DO K = 1,3
          DO L = 1,3
            AJJMAT(ISLIP,I,J,K,L) = 0.0
            DO M  = 1,3
            DO N  = 1,3
              AJJMAT(ISLIP,I,J,K,L) = AJJMAT(ISLIP,I,J,K,L) +
     &                  0.5*ELMAT(I,J,M,N)*GJMAT(ISLIP,M,N,K,L)
            END DO
            END DO
            AUXMAT1(I,J,K,L) = AJJMAT(ISLIP,I,J,K,L)
          END DO
          END DO
          END DO
          END DO
C
          CALL EMATR(AUXMAT1,AUXMAT1R)
          DO I = 1,6
          DO J = 1,6
            AJJMATR(ISLIP,I,J) = AUXMAT1R(I,J)
          END DO
          END DO            
        END DO
C
C  STEP J6. CALCULATE THE NSLIP MATRICES B_IJ^{\ALPHA}

        DO ISLIP = 1,NSLIP
        
          if (ISLIP < 13) then 
            call calc_GAMMADOT_G(props,nprops,noel,npt,theta,rss(islip),
     &      btau(islip),stau(islip),gdotalpha_G,dgdtau_G)
            call calc_GAMMADOT_GP(props,nprops,noel,npt,theta,rss(islip),
     &      btau(islip),stau(islip),gdotalpha_GP,dgdtau_GP)
C           GDOTALPHA_GP = 0.0
C           DGDTAU_GP    = 0.0 
          else
            GDOTALPHA_G = 0.0
            DGDTAU_G    = 0.0 
            call calc_GAMMADOT_GP_cube(props,nprops,noel,npt,theta,rss(islip),
     &      btau(islip),stau(islip),gdotalpha_GP,dgdtau_GP)
          end if
C         if (ISLIP > 12) then
C            GDOTALPHA_G = 0.0
C            DGDTAU_G    = 0.0
C         end if
          DO I     = 1,3
          DO J     = 1,3
            BJMAT(ISLIP,I,J) = 
     &        0.5*(SMATG(ISLIP,I,J)+SMATG(ISLIP,J,I))*
     &        ((1.0-f_gp)*DGDTAU_G + f_gp*DGDTAU_GP)*DTIME
C     &        DGDTAU*DTIME
            AUX1(I,J)=BJMAT(ISLIP,I,J)
          END DO
          END DO
C
          CALL SMATVEC(AUX1,AUXV1)
C
          DO I=1,6
            BJMATR(ISLIP,I) = AUXV1(I)
          END DO
        END DO
C
C  STEP J7. CALCULATE THE  REDUCED MATRIX K_IJ
C
        DO I = 1,6
        DO J = 1,6
          AKJMATR(I,J) = 0.0
          DO ISLIP = 1,NSLIP
            AKJMATR(I,J) = AKJMATR(I,J) +
     &                CMAT(ISLIP,I)*BJMATR(ISLIP,J)
          END DO
          IF(I.EQ.J)THEN
            AKJMATR(I,J) = 1. + AKJMATR(I,J)
          END IF
        END DO
        END DO
C
C  STEP J8. CALCULATE THE INVERSE OF THE MATRIX K_IJ
C
        CALL MATINV(AKJMATR,6,6,INDX,AKJMATRINV)
C
C  STEP J9. CALCULATE THE REDUCED MATRIX M_IJ
C
        DO I = 1,6
        DO J = 1,6
          AMJMATR(I,J) = 0.0
          DO ISLIP= 1, NSLIP
            AMJMATR(I,J) = AMJMATR(I,J) +
     &                     DGAMMA(ISLIP)*AJJMATR(ISLIP,I,J)
          END DO
          AMJMATR(I,J) =  DJMATR(I,J) - AMJMATR(I,J)
        END DO
        END DO
C
C  STEP J10. CALCULATE THE  MATRIX Q_IJKL
C     
        DO I = 1,6
        DO J = 1,6
          QJMATR(I,J) = 0.0
          DO M  = 1,6
            QJMATR(I,J) = QJMATR(I,J)+ AKJMATRINV(I,M)*
     &                    AMJMATR(M,J)
          END DO
        END DO
        END DO
C
          CALL EMATREC(QJMATR,QJMAT)
C
C  STEP J11. CALCULATE THE MATRICES R_IJ^{\ALPHA}
C                 
        DO  ISLIP = 1, NSLIP
        DO  I = 1,3
        DO  J = 1,3
          RJMAT(ISLIP,I,J) = 0.0
          DO K = 1,3
          DO L = 1,3
            RJMAT(ISLIP,I,J) = RJMAT(ISLIP,I,J) +
     &                         BJMAT(ISLIP,K,L)*QJMAT(K,L,I,J)
          END DO
          END DO
        END DO
        END DO
        END DO
C
C  STEP J12. CALCULATE THE MATRIX S_IJKL
C 
        DO IP = 1,3
        DO J = 1,3
        DO K = 1,3
        DO L = 1,3
          AUX1(IP,J)        = 0.0
          AUXMAT1(IP,J,K,L) = 0.0
          DO ISLIP = 1,NSLIP
            AUX1(IP,J) = AUX1(IP,J) 
     &                  + DGAMMA(ISLIP)*SMATG(ISLIP,IP,J)
            AUXMAT1(IP,J,K,L) = AUXMAT1(IP,J,K,L)
     &                  +SMATG(ISLIP,IP,J)*RJMAT(ISLIP,K,L)
          END DO           
        END DO
        END DO
        END DO
        END DO
C                   
        DO L = 1,3
        DO J = 1,3
          AUX2(L,J) = 0.0
          DO IP = 1,3
            AUX2(L,J) = AUX2(L,J) +
     &                 FSTART(L,IP)*AUX1(IP,J)        
          END DO               
        END DO
        END DO               
C
        DO N = 1,3
        DO J = 1,3
        DO K = 1,3
        DO L = 1,3
          AUXMAT2(N,J,K,L) = 0.0
          DO IP = 1,3                
            AUXMAT2(N,J,K,L) =  AUXMAT2(N,J,K,L) +
     &                    FSTART(N,IP)*AUXMAT1(IP,J,K,L)
          END DO
        END DO
        END DO
        END DO
        END DO   
C
        DO I = 1,3
        DO J = 1,3
        DO K = 1,3
        DO L = 1,3
          SJMAT(I,J,K,L) = 0.0
          AUXMAT1(I,J,K,L) = 0.0
          DO N = 1,3                 
            AUXMAT1(I,J,K,L) =  AUXMAT1(I,J,K,L) +
     &                      FTTAU(I,N)*AUXMAT2(N,J,K,L)
          END DO
          SJMAT(I,J,K,L) = SJMAT(I,J,K,L) 
     &           + RTTAU(I,K)*FSTART(L,J)
     &           - RTTAU(I,K)*AUX2(L,J)
     &           - AUXMAT1(I,J,K,L)
        END DO
        END DO
        END DO
        END DO   
C
C  STEP J13. CALCULATE THE JACOBIAN MATRIX FOR EACH CRYSTAL
C
        CALL inverse_3x3(FSTAU,FSTAUINV)
C
        detfstau = determinant(fstau)
C            
        DO I = 1,3
        DO J = 1,3
        DO K = 1,3
        DO L = 1,3
          WJMAT(I,J,K,L) = 0.0
          SUM1 = 0.0
          DO IP = 1,3
          DO IQ = 1,3
            SUM1 = SUM1 + SJMAT(IP,IQ,K,L)*FSTAUINV(IQ,IP)
          END DO
          END DO
          DO M = 1,3
          DO N = 1,3
            WJMAT(I,J,K,L) = WJMAT(I,J,K,L) +
     &                SJMAT(I,M,K,L)*TSTARTAU(M,N)*FSTAU(J,N) +
     &                FSTAU(I,M)*QJMAT(M,N,K,L)*FSTAU(J,N) +
     &                FSTAU(I,M)*TSTARTAU(M,N)*SJMAT(J,N,K,L) -
     &                FSTAU(I,M)*TSTARTAU(M,N)*FSTAU(J,N)*SUM1
          END DO
          END DO
          WJMAT(I,J,K,L) = WJMAT(I,J,K,L)/DETFSTAU
        END DO
        END DO
        END DO
        END DO
C
C  STEP J14. CALCULATE THE REDUCED JACOBIAN MATRIX FOR EACH CRYSTAL
C
        CALL EMATR(WJMAT,WJMATR)
C
        RETURN
        END 
C*********************************************************************
C ****************************************************************
C  THIS SUBROUTINE CONVERTS A SYMMETRIC MATRIX TO A VECTOR
        SUBROUTINE SMATVEC(A,V)
        IMPLICIT REAL*8(A-H,O-Z)       
        DIMENSION A(3,3),V(6)
        DO I = 1,3
          V(I) = A(I,I)
        END DO
        V(4) = A(1,2)
        V(5) = A(1,3)
        V(6) = A(2,3)
        RETURN
        END
C ******************************************************************
C ****************************************************************
      SUBROUTINE MATINV(A,N,NP,INDX,Y)
C
C  Given an NxN matrix [A], with physical dimension NP, this 
C  routine replaces it by the LU decomposition of a rowwise 
C  permutation of itself. [A] and N are input. [A] is output, 
C  arranged in LU form. INDX is an output vector which records
C  the row permutation effected by the partial pivoting; 
C  D is output as +1 or -1 depending on wheter the number of
C  row interchanges was even or odd, respectively.
C
C  Once the LU decomposition is performed, this routine 
C  calculates the inverse of [A] by using subroutine LUBKSB.
C  Note that INDX is input as the permutation vector returned by        
C  LUDCMP. {B} is input as the right-hand side vector {B}, and  
C  returns with the solution vector {X}. [A], N, NP, INDX are not 
C  modified by this routine, and are left in place
C  for succesive calls with different right-hand sides {B}.
C  This routine takes into account that {B} will begin with
C  many zero elements, so it is efficient for use in matrix 
C  inversion.
C
C  The inverse of [A] is calculated using as many unit vectors
C  {B} needed as right hand side vectors. The result is
C  returned as the matrix [Y].
C
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(NP,NP), Y(NP,NP), INDX(NP)
C
C  Set up the identity matrix  
C
        DO I = 1,N
          DO J = 1,N
            Y(I,J) = 0.0
          END DO
          Y(I,I) = 1.
        END DO
C
C  Decompose the matrix just once
C
        CALL LUDCMP(A,N,NP,INDX,D)
C
C  Find the inverse by columns. It is necessary to recognize
C  that FORTRAN stores two dimensional matrices by column, so
C  so that Y(1,J) is the address of the Jth column of Y.
C
        DO J=1,N
          CALL LUBKSB(A,N,NP,INDX,Y(1,J))
        END DO
        RETURN
        END
C**************************************************************
      SUBROUTINE LUDCMP(A,N,NP,INDX,D)
C
C  Given an NxN matrix [A], with physical dimension NP, this 
C  routine replaces it by the LU decomposition of a row-wise 
C  permutation of itself. [A] and N are input. [A] is output, 
C  arranged in LU form. INDX is an output vector which records
C  the row permutation effected by the partial pivoting; 
C  D is output as +1 or -1 depending on wheter the nuber of
C  row interchanges was even or odd, respectively. This routine
C  is used in combination with LUBKSB to solve linear equations 
C  or invert a matrix.
C
        IMPLICIT REAL*8 (A-H,O-Z)
        
        PARAMETER (NMAX=100,TINY=1.0E-20)
        DIMENSION A(NP,NP),INDX(N),VV(NMAX)
        D=1.
        DO I=1,N
        AAMAX=0.
        DO J=1,N
          IF (DABS(A(I,J)).GT.AAMAX) AAMAX=DABS(A(I,J))
        END DO
!        IF (AAMAX.EQ.0.) write(*,*) 'Singular matrix.'
        VV(I)=1./AAMAX
        END DO
        DO J=1,N
        IF (J.GT.1) THEN
          DO I=1,J-1
            rSUM=A(I,J)
            IF (I.GT.1)THEN
              DO K=1,I-1
                rSUM=rSUM-A(I,K)*A(K,J)
              END DO
              A(I,J)=rSUM
            END IF
          END DO
        END IF
        AAMAX=0.
        DO I=J,N
          rSUM=A(I,J)
          IF (J.GT.1)THEN
            DO K=1,J-1
              rSUM=rSUM-A(I,K)*A(K,J)
            END DO
            A(I,J)=rSUM
          END IF
          DUM=VV(I)*DABS(rSUM)
          IF (DUM.GE.AAMAX) THEN
            IMAX=I
            AAMAX=DUM
          END IF
        END DO
        IF (J.NE.IMAX)THEN
          DO K=1,N
            DUM=A(IMAX,K)
            A(IMAX,K)=A(J,K)
            A(J,K)=DUM
          END DO
          D=-D
          VV(IMAX)=VV(J)
        END IF
        INDX(J)=IMAX
        IF(J.NE.N)THEN
          IF(A(J,J).EQ.0.)A(J,J)=TINY
          DUM=1./A(J,J)
          DO I=J+1,N
            A(I,J)=A(I,J)*DUM
          END DO
        END IF
        END DO
        IF(A(N,N).EQ.0.)A(N,N)=TINY
        RETURN
        END
C**********************************************************
       SUBROUTINE LUBKSB(A,N,NP,INDX,B)
C
C  Solves the set of N linear equations [A]{X} = {B}. 
C  Here [A] is input, not as the matrix [A], but as its LU 
C  decomposition, determined by the routine LUDCMP. INDX
C  is input as the permutation vector returned by LUDCMP. {B}
C  is input as the right-hand side vector {B}, and returns
C  with the solution vector {X}. [A], N, NP, INDX are not 
C  modified by this routine, and can be left in place
C  for succesive calls with different right-hand sides {B}.
C  This routine takes into account that {B} will begin with
C  many zero elements, so it is efficient for use in matrix 
C  inversion.
C
        IMPLICIT REAL*8 (A-H,O-Z)
        DIMENSION A(NP,NP),INDX(N),B(N)
       II=0
       DO I=1,N
        LL=INDX(I)
        rSUM=B(LL)
        B(LL)=B(I)
        IF (II.NE.0)THEN
          DO J=II,I-1
            rSUM=rSUM-A(I,J)*B(J)
          END DO
        ELSE IF (rSUM.NE.0.) THEN
          II=I
        END IF
        B(I)=rSUM
       END DO
       DO I=N,1,-1
        rSUM=B(I)
        IF(I.LT.N)THEN
          DO J=I+1,N
            rSUM=rSUM-A(I,J)*B(J)
          END DO
        END IF
        B(I)=rSUM/A(I,I)
       END DO
       RETURN
       END
C ****************************************************************
C******************************************************************************
      SUBROUTINE EMATREC(ELMATR,ELMAT)
C      
C  THIS SUBROUTINE RECONSTRUCTS A STIFFNESS MATRIX ELMATREC
C  AS A 3 X 3 X 3 X 3 MATRIX FROM ITS 
C  REDUCED 6 X 6 FORM ELMATR
C
      IMPLICIT REAL*8(A-H,O-Z) 
      DIMENSION ELMATR(6,6),ELMAT(3,3,3,3)
C
      DO I = 1,3
      DO J = 1,3
      DO K = 1,3
      DO L = 1,3
        ELMAT(I,J,K,L) = 0.0
      END DO
      END DO
      END DO
      END DO
C
C  RECONSTRUCT THE MATRIX ELMAT FROM ITS REDUCED FORM ELMATR
C
        ELMAT(1,1,1,1)  = ELMATR(1,1)
        ELMAT(1,1,2,2)  = ELMATR(1,2)
        ELMAT(1,1,3,3)  = ELMATR(1,3) 
        ELMAT(1,1,1,2)  = ELMATR(1,4)
        ELMAT(1,1,2,1)  = ELMATR(1,4)
        ELMAT(1,1,1,3)  = ELMATR(1,5)
        ELMAT(1,1,3,1)  = ELMATR(1,5)
        ELMAT(1,1,2,3)  = ELMATR(1,6)
        ELMAT(1,1,3,2)  = ELMATR(1,6)
C
        ELMAT(2,2,1,1)  = ELMATR(2,1)
        ELMAT(2,2,2,2)  = ELMATR(2,2)
        ELMAT(2,2,3,3)  = ELMATR(2,3) 
        ELMAT(2,2,1,2)  = ELMATR(2,4)
        ELMAT(2,2,2,1)  = ELMATR(2,4)
        ELMAT(2,2,1,3)  = ELMATR(2,5)
        ELMAT(2,2,3,1)  = ELMATR(2,5)
        ELMAT(2,2,2,3)  = ELMATR(2,6)
        ELMAT(2,2,3,2)  = ELMATR(2,6)
C
        ELMAT(3,3,1,1)  = ELMATR(3,1)
        ELMAT(3,3,2,2)  = ELMATR(3,2)
        ELMAT(3,3,3,3)  = ELMATR(3,3) 
        ELMAT(3,3,1,2)  = ELMATR(3,4)
        ELMAT(3,3,2,1)  = ELMATR(3,4)
        ELMAT(3,3,1,3)  = ELMATR(3,5)
        ELMAT(3,3,3,1)  = ELMATR(3,5)
        ELMAT(3,3,2,3)  = ELMATR(3,6)
        ELMAT(3,3,3,2)  = ELMATR(3,6)
C
        ELMAT(1,2,1,1)  =  ELMATR(4,1)
        ELMAT(2,1,1,1)  =  ELMATR(4,1)
        ELMAT(1,2,2,2)  =  ELMATR(4,2)
        ELMAT(2,1,2,2)  =  ELMATR(4,2)  
        ELMAT(1,2,3,3)  =  ELMATR(4,3)
        ELMAT(2,1,3,3)  =  ELMATR(4,3)
        ELMAT(1,2,1,2)  =  ELMATR(4,4) 
        ELMAT(2,1,1,2)  =  ELMATR(4,4) 
        ELMAT(1,2,2,1)  =  ELMATR(4,4)
        ELMAT(2,1,2,1)  =  ELMATR(4,4)
        ELMAT(1,2,1,3)  =  ELMATR(4,5) 
        ELMAT(2,1,1,3)  =  ELMATR(4,5) 
        ELMAT(1,2,3,1)  =  ELMATR(4,5)
        ELMAT(2,1,3,1)  =  ELMATR(4,5)
        ELMAT(1,2,2,3)  =  ELMATR(4,6) 
        ELMAT(2,1,2,3)  =  ELMATR(4,6) 
        ELMAT(1,2,3,2)  =  ELMATR(4,6)
        ELMAT(2,1,3,2)  =  ELMATR(4,6)
C
        ELMAT(1,3,1,1)  =  ELMATR(5,1)
        ELMAT(3,1,1,1)  =  ELMATR(5,1)  
        ELMAT(1,3,2,2)  =  ELMATR(5,2)
        ELMAT(3,1,2,2)  =  ELMATR(5,2)  
        ELMAT(1,3,3,3)  =  ELMATR(5,3)
        ELMAT(3,1,3,3)  =  ELMATR(5,3)
        ELMAT(1,3,1,2)  =  ELMATR(5,4) 
        ELMAT(3,1,1,2)  =  ELMATR(5,4) 
        ELMAT(1,3,2,1)  =  ELMATR(5,4)
        ELMAT(3,1,2,1)  =  ELMATR(5,4)
        ELMAT(1,3,1,3)  =  ELMATR(5,5) 
        ELMAT(3,1,1,3)  =  ELMATR(5,5) 
        ELMAT(1,3,3,1)  =  ELMATR(5,5)
        ELMAT(3,1,3,1)  =  ELMATR(5,5)
        ELMAT(1,3,2,3)  =  ELMATR(5,6) 
        ELMAT(3,1,2,3)  =  ELMATR(5,6) 
        ELMAT(1,3,3,2)  =  ELMATR(5,6)
        ELMAT(3,1,3,2)  =  ELMATR(5,6)
C
        ELMAT(2,3,1,1)  =  ELMATR(6,1)
        ELMAT(3,2,1,1)  =  ELMATR(6,1)  
        ELMAT(2,3,2,2)  =  ELMATR(6,2)
        ELMAT(3,2,2,2)  =  ELMATR(6,2)  
        ELMAT(2,3,3,3)  =  ELMATR(6,3)
        ELMAT(3,2,3,3)  =  ELMATR(6,3)
        ELMAT(2,3,1,2)  =  ELMATR(6,4) 
        ELMAT(3,2,1,2)  =  ELMATR(6,4) 
        ELMAT(2,3,2,1)  =  ELMATR(6,4)
        ELMAT(3,2,2,1)  =  ELMATR(6,4)
        ELMAT(2,3,1,3)  =  ELMATR(6,5) 
        ELMAT(3,2,1,3)  =  ELMATR(6,5) 
        ELMAT(2,3,3,1)  =  ELMATR(6,5)
        ELMAT(3,2,3,1)  =  ELMATR(6,5)
        ELMAT(2,3,2,3)  =  ELMATR(6,6) 
        ELMAT(3,2,2,3)  =  ELMATR(6,6) 
        ELMAT(2,3,3,2)  =  ELMATR(6,6)
        ELMAT(3,2,3,2)  =  ELMATR(6,6)
        RETURN 
        END  
C **************************************************************
C*************************************************************************
        SUBROUTINE SOLVCRYS(nslip,props,nprops,noel,npt,theta,TSTV,ST,BT,
     &  amat,rhossd_n,TSTRV,PMAT,CMAT,DTIME,TSTAUV,STAU,BTAU,rhossd,DGAMMA,
     &  DGMAX,RSS,DOMPSC,ITERK,ITERL,ITERM,ITERERR,xmC)
C
C       THIS SUBROUTINE SOLVES FOR  (TSTARTAU,STAU)
C       HERE, WE ARE SETTING
C
C       TSTARTV   AS TSTV
C       TSTARTRV  AS TSTRV
C       TSTARTAUV AS TSTAUV
C --------------------------------------------------------------------
        IMPLICIT REAL*8(A-H,O-Z)

        PARAMETER  (MAXIT1=100,MAXIT2=20,NP=6,N=6)
        DIMENSION  ALPHAM(NP,NP),BETA(NP),INDX(NP)
        DIMENSION  TSTV(6),ST(nslip),BT(nslip),TSTRV(6)
        DIMENSION  PMAT(nslip,6)
        DIMENSION  CMAT(nslip,6)
        DIMENSION  TSTAUV(6),STAU(nslip),BTAU(nslip),STRESS(6)
        DIMENSION  RSS(nslip),DGAMMA(nslip),DGAMDTAU(nslip)
        DIMENSION  HAB(nslip,nslip)
        DIMENSION  TAUBTR(nslip)
        DIMENSION  RESOLD(nslip),RESNEW(nslip)
        DIMENSION  TAUBOLD(nslip),TAUBNEW(nslip)
        DIMENSION  TEMPSTRESS(6)
        dimension  props(nprops)
        dimension  rhossd_n(nslip),rhossd(nslip),rhossd_old(nslip)
        dimension  drhossd(nslip)        
        dimension  dtaub(nslip)
        dimension  amat(nslip,nslip)
        dimension  xmC(3,nslip)

!  Initialize iteration counters and other parameters
        ITERK    = 0
        ITERL    = 0
        ITERLNEW = 0
        ITERM    = 0
        ITERMNEW = 0
        TOL      = 1.d-4
        DGMAXTOL = 0.50
        CORRBETA = 4.0
C
        f_gp=props(39)
        debug=int(props(59))
C  INITIAL ESTIMATES FOR THE ROOTS:
C
        DO I=1,6
           STRESS(I) = TSTV(I)
        END DO

        DO ISLIP=1,NSLIP
           RESOLD(ISLIP) = ST(ISLIP)
           rhossd_old(islip) = rhossd_n(islip)
           TAUBOLD(ISLIP) = BT(ISLIP)
        END DO
999     CONTINUE
        ITERERR = 0
C
C  START THE LOOP FOR SOLVING FOR THE COMPONENTS OF TSTAUV
C
        DO K=1,MAXIT1

!        if (K>5) then
!          dgmax = 0.d0
!        end if
C
C  CALCULATE:
C     1. THE RESOLVED SHEAR STRESS ON EACH
C        SLIP SYSTEM CORRESPONDING TO THE ESTIMATE OF TSTAUV.
C     2. CALCULATE THE PLASTIC SHEAR STRAIN INCREMENT FOR
C        EACH SLIP SYSTEM.
C     3. CALCULATE THE MAXIMUM SHEAR STRAIN INCREMENT 
C     4. CALCULATE THE PLASTIC WORK INCREMENT
C
           DGMAX   = 0.0
           DOMPSC = 0.0
           DO ISLIP = 1,NSLIP
             TAUALPHA = 0.0
             DO J=1,6
               TAUALPHA = TAUALPHA + STRESS(J)*PMAT(ISLIP,J)
             END DO
             RSS(ISLIP) = TAUALPHA
             RESALPHA   = RESOLD(ISLIP)
             TAUBALPHA   = TAUBOLD(ISLIP)
           if (noel==20 .AND. npt==1 .AND. debug==1) then
                    print*,'#######----------#######'
                    print*,'#######for element#1'
                    print*,'#######----------#######'
                    print*,'calling flow rule for Gamma'
                    print*,'RSS for slip system',RSS(ISLIP),ISLIP
           end if 
             if (ISLIP < 13) then
               call calc_GAMMADOT_G(props,nprops,noel,npt,theta,TAUALPHA,
     &         TAUBALPHA,RESALPHA,GDOTALPHA_G,DGDTAU_G)
               call calc_GAMMADOT_GP(props,nprops,noel,npt,theta,TAUALPHA,
     &         TAUBALPHA,RESALPHA,GDOTALPHA_GP,DGDTAU_GP)
C              GDOTALPHA_GP = 0.0
C              DGDTAU_GP    = 0.0               
             else
               call calc_GAMMADOT_G(props,nprops,noel,npt,theta,TAUALPHA,
     &         TAUBALPHA,RESALPHA,GDOTALPHA_G,DGDTAU_G)
               GDOTALPHA_G = 0.0
               DGDTAU_G    = 0.0
               call calc_GAMMADOT_GP_cube(props,nprops,noel,npt,theta,TAUALPHA,
     &         TAUBALPHA,RESALPHA,GDOTALPHA_GP,DGDTAU_GP)
             end if
             DGAMMA(ISLIP) = DTIME*((1.0-f_gp)*GDOTALPHA_G + f_gp*GDOTALPHA_GP)

C             DGAMMA(ISLIP) = DTIME*GDOTALPHA
             DGAMDTAU(ISLIP) = DTIME*((1.0-f_gp)*DGDTAU_G + f_gp*DGDTAU_GP)
C             DGAMDTAU(ISLIP) = DTIME*DGDTAU
           if (noel==20 .AND. npt==1 .AND. debug==1) then
            ! if (noel == 1 .AND. npt==1 .AND. kinc==1 .AND. time(2)==0.0) then
                    print*,'#######----------#######'
                    print*,'#######for element#1'
                    print*,'#######----------#######'
                    print*,'Shear rate for gamma',GDOTALPHA_G
                    print*,'Shear rate for gamma prime',GDOTALPHA_GP
           end if

             DGMAX = DMAX1( DGMAX,DABS( DGAMMA(ISLIP) ) )
             TEMP =  DABS ( RSS(ISLIP)*DGAMMA(ISLIP) )
             DOMPSC = DOMPSC + TEMP
           END DO

! No significant plasticity, therefore no change in hardening, dislocation density and back-stress
! Trial stress (tstrv) calculated from elastic strain = current stress
           if (dgmax < 1d-20) then
             tstauv = tstrv
             stau = st
             rhossd = rhossd_n
             btau = bt
             return
           end if
C
C  Calculate BETA_i --- the negative of the values of
C  the functions G_i(T_k)
C
           DO I =1,6
             BETA(I) = 0.0
             DO ISLIP=1,NSLIP
! Beta is the extra stress considerd in the current stress increment due to assumption Fp(t+dt)=Fp(t)             
               BETA(I) = BETA(I) + DGAMMA(ISLIP)*CMAT(ISLIP,I)
             END DO
! Beta should be subtracted from the current stress incerement (TSTRV(I)- STRESS(I))            
             BETA(I) =   TSTRV(I) - BETA(I) - STRESS(I)
           END DO
C
C  Calculate the Jacobian for the Newton corrections
C
           DO I =1,6
           DO J =1,6
             IF (I .EQ. J) THEN      
               ALPHAM(I,J)      = 1.0D0
             ELSE
               ALPHAM(I,J)      = 0.0D0
             END IF
             DO ISLIP=1,NSLIP
               IF(RSS(ISLIP).EQ.0.0)GO TO 60
               ALPHAM(I,J)= ALPHAM(I,J) + 
     &                    DGAMDTAU(ISLIP)*CMAT(ISLIP,I)*PMAT(ISLIP,J)
60           CONTINUE
             END DO
           END DO
           END DO    
C
C  Calculate the Newton corrections
C

          CALL LUDCMP(ALPHAM,N,NP,INDX,D)
          CALL LUBKSB(ALPHAM,N,NP,INDX,BETA)
C
C  BETA  NOW CONTAINS THE NEWTON CORRECTIONS.
C  CHECK TO SEE IF THE NEWTON CORRECTIONS ARE
C  WITHIN SPECIFIED TOLERANCES. IF YES, THEN
C  PROCESS HAS CONVERGED. UPDATE AND GO TO LOOP FOR
C  ITERATING ON THE DEFORMATION RESISTANCES.
C  IF NOT, THEN CONTINUE ITERATING.
C
          IF (   ( DABS( BETA(1) ) .LE. TOL ) .AND.
     +           ( DABS( BETA(2) ) .LE. TOL ) .AND.
     +           ( DABS( BETA(3) ) .LE. TOL ) .AND.
     +           ( DABS( BETA(4) ) .LE. TOL ) .AND.
     +           ( DABS( BETA(5) ) .LE. TOL ) .AND.
     +           ( DABS( BETA(6) ) .LE. TOL)   )  THEN  
C
C  ITERATIVE PROCEDURE HAS CONVERGED. UPDATE 
C
              TSTAUV(1) = STRESS(1)
              TSTAUV(2) = STRESS(2)
              TSTAUV(3) = STRESS(3)
              TSTAUV(4) = STRESS(4)
              TSTAUV(5) = STRESS(5)
              TSTAUV(6) = STRESS(6)           
              ITERK    =  K
              GO TO 550
           END IF
C
C  TEMPORARILY ADD THE CORRECTIONS TO THE VARIABLES
C
           TEMPSTRESS(1) = STRESS(1) + BETA(1)
           TEMPSTRESS(2) = STRESS(2) + BETA(2)
           TEMPSTRESS(3) = STRESS(3) + BETA(3)
           TEMPSTRESS(4) = STRESS(4) + BETA(4)
           TEMPSTRESS(5) = STRESS(5) + BETA(5)
           TEMPSTRESS(6) = STRESS(6) + BETA(6)
C
C  CHECK IF THESE STRESS COMPONENTS GIVE RISE TO LARGE DELTA GAMMAS
C  IF THEY DO THEN REDUCE THE N-R CORRECTIONS BY A FACTOR OF 4
C  AND CONTINUE UNTIL THE DELTA GAMMAS ARE REASONABLE
C
           ICOR = 0
           dgmaxold = 0.d0
85         CONTINUE
           IF(ICOR.GT.10) THEN 
C
C  THE METHOD TO CONSTRAIN THE NEWTON-RAPHSON METHOD 
C  IS NOT WORKING! GET OUT OF THE LOOP! 
C  RETURN TO CUT BACK TIME-STEP AND RECALCULATE
C
!             WRITE(*,*) 'TOO MANY ITERATIONS IN CONSTRAINT....'
             GOTO 500
           END IF
           DGMAX = 0.0
           DO ISLIP = 1,NSLIP
             TAUALPHA = 0.0
             DO J=1,6
               TAUALPHA = TAUALPHA + TEMPSTRESS(J)*PMAT(ISLIP,J)
             END DO
             RSS(ISLIP) = TAUALPHA
             RESALPHA   = RESOLD(ISLIP)
             TAUBALPHA = TAUBOLD(ISLIP)
             if (ISLIP < 13) then
               call calc_GAMMADOT_G(props,nprops,noel,npt,theta,TAUALPHA,
     &         TAUBALPHA,RESALPHA,GDOTALPHA_G,DGDTAU_G)
               call calc_GAMMADOT_GP(props,nprops,noel,npt,theta,TAUALPHA,
     &         TAUBALPHA,RESALPHA,GDOTALPHA_GP,DGDTAU_GP)
C              GDOTALPHA_GP = 0.0
C              DGDTAU_GP    = 0.0               
             else
               GDOTALPHA_G = 0.0
               DGDTAU_G    = 0.0
               call calc_GAMMADOT_GP_cube(props,nprops,noel,npt,theta,TAUALPHA,
     &         TAUBALPHA,RESALPHA,GDOTALPHA_GP,DGDTAU_GP)
             end if
             DGAMMA(ISLIP) = DTIME*((1.0-f_gp)*GDOTALPHA_G + f_gp*GDOTALPHA_GP)

!             DGAMMA(ISLIP) = DTIME*GDOTALPHA
             DGAMDTAU(ISLIP) = DTIME*((1.0-f_gp)*DGDTAU_G + f_gp*DGDTAU_GP)
!             DGAMDTAU(ISLIP) = DTIME*DGDTAU
             DGMAX = DMAX1( DGMAX,DABS( DGAMMA(ISLIP) ) )
           END DO
C
           IF(DGMAX.GT.DGMAXOLD)THEN
             CORRBETA=1./CORRBETA
           END IF
           DGMAXOLD=DGMAX
C
           IF (DGMAX .GT. DGMAXTOL) THEN
              DO I = 1,6
                 BETA(I) = CORRBETA * BETA(I)
                 TEMPSTRESS(I) = STRESS(I) + BETA(I)
              END DO
              ICOR = ICOR +1
              GOTO 85
           ELSE
C
C  ACCEPT THE CORRECTIONS TO THE STRESS AND CONTINUE ITERATING
C
              DO I = 1,6
                 STRESS(I) = TEMPSTRESS(I) 
              END DO
           END IF
C                    
500     CONTINUE
        END DO
C
C  ITERATIVE PROCEDURE FOR THE STRESS 
C  HAS NOT CONVERGED IN MAXIT1 ITERATIONS.
C  RETURN TO CUT BACK TIME-STEP, AND RECALCULATE.
C
        ITERK = MAXIT1
        ITERERR = 1
        RETURN
550   CONTINUE
C
C  ITERATIVE PROCEDURE FOR TSTAUV HAS CONVERGED
C  START ITERATING FOR THE SLIP SYSTEM DEFORMATION 
C  RESISTANCES
C
C  ESTIMATE THE THE SLIP SYSTEM DEFORMATION RESISTANCES
C  AT THE END OF THE STEP
         call calc_HARDENING(props,nprops,nslip,theta,amat,rhossd_old,DGAMMA,
     &        DTIME,HAB,drhossd,xmC)

C  ESTIMATE THE SLIP SYSTEM DEFORMATION RESISTANCES AT THE
C  END OF THE INCREMENT
C       
        DO ISLIP = 1, NSLIP
          RESNEW(ISLIP) = 0.0
          DO JSLIP = 1,NSLIP
            RESNEW(ISLIP) = RESNEW(ISLIP)+
     +                   HAB(ISLIP,JSLIP)*DABS(DGAMMA(JSLIP))
          END DO        
        END DO
        DO ISLIP = 1,NSLIP
          RESNEW(ISLIP) =   RESNEW(ISLIP) + ST(ISLIP)
          rhossd(islip) =   drhossd(islip) + rhossd_n(islip)
        END DO
        ERRRES = 0.0
        DO ISLIP = 1,NSLIP
         ERRRES = DMAX1(DABS(RESNEW(ISLIP)-RESOLD(ISLIP)),ERRRES)
        END DO
C
        TOLS = 1.D-4
        IF (ERRRES .LT. TOLS) THEN
C
C  ITERATIVE PROCEDURE for deformation resistance HAS CONVERGED
C  UPDATE STAU and continue with back stress iteration
C
           DO ISLIP = 1,NSLIP
             STAU(ISLIP) = RESNEW(ISLIP)
           END DO
           GO TO 600
!           RETURN
        ELSE
          DO ISLIP = 1,NSLIP
            RESOLD(ISLIP) = RESNEW(ISLIP)
            rhossd_old(islip) = rhossd(islip)
          END DO
          ITERLNEW = ITERLNEW + 1
          ITERL = MAX0(ITERL ,ITERLNEW)
          IF (ITERL .GT. MAXIT2) THEN
C
C  ITERATIVE PROCEDURE FOR THE DEFORMATION RESISTANCE HAS NOT
C  CONVERGED IN MAXIT2 ITERATIONS. RETURN TO CUT BACK THE TIME
C  STEP, AND RECALCULATE.
C
            ITERL = MAXIT2
            ITERERR = 1
!            WRITE (*,*) 'MAXIT2 EXCEEDED....REPEATING TIME STEP'
            RETURN
          END IF
          GOTO 999         
        END IF
600   	  CONTINUE

C  ITERATIVE PROCEDURE FOR THE DEFORMATION RESISTANCE HAS CONVERGED
C  START ITERATING FOR THE SLIP SYSTEM BACK STRESS
C
C  ESTIMATE THE THE SLIP SYSTEM BACK STRESS
C  AT THE END OF THE STEP
         call calc_BACKSTRESS(props,nprops,nslip,TAUBOLD,DGAMMA,
     &        TAUBTR,DTAUB)
C  ESTIMATE THE SLIP SYSTEM BACK STRESS AT THE
C  END OF THE INCREMENT
C

        DO ISLIP = 1, NSLIP
          TAUBNEW(ISLIP) = 0.0
          TAUBNEW(ISLIP) = TAUBNEW(ISLIP)+ 
     +                     DTAUB(ISLIP)
        END DO
        DO ISLIP = 1,NSLIP
          TAUBNEW(ISLIP) =   TAUBNEW(ISLIP) + BT(ISLIP)
        END DO
        ERRRES = 0.0
        DO ISLIP = 1,NSLIP
         ERRRES = DMAX1(DABS(TAUBNEW(ISLIP)-TAUBOLD(ISLIP)),ERRRES)
        END DO
C
        TOLS = 1.D-4
        IF (ERRRES .LT. TOLS) THEN
C
C  ITERATIVE PROCEDURE HAS CONVERGED
C  UPDATE TAUB AND RETURN
C
           DO ISLIP = 1,NSLIP
             BTAU(ISLIP) = TAUBNEW(ISLIP)
           END DO
           RETURN
        ELSE
          DO ISLIP = 1,NSLIP
            TAUBOLD(ISLIP) = TAUBNEW(ISLIP)
          END DO
          ITERMNEW = ITERMNEW + 1
          ITERM = MAX0(ITERM ,ITERMNEW)
          IF (ITERM .GT. MAXIT2) THEN
C
C  ITERATIVE PROCEDURE FOR THE BACK STRESS HAS NOT
C  CONVERGED IN MAXIT2 ITERATIONS. RETURN TO CUT BACK THE TIME
C  STEP, AND RECALCULATE.
C
            ITERM = MAXIT2
            ITERERR = 1
!            WRITE (*,*) 'MAXIT2 EXCEEDED....REPEATING TIME STEP'
            RETURN
          END IF
          GOTO 999
        END IF
        END
C **********************************************************************
        SUBROUTINE calc_GAMMADOT_G(props,nprops,noel,npt,theta,TAUALPHA,
     &  TAUB,RESALPHA,GDOTALPHA_G,DGDTAU_G)

        IMPLICIT REAL*8(A-H,O-Z)

        dimension props(nprops)

!  Initialize props needed for the calcuation
        C11_0_G= PROPS(7) 
        C12_0_G= PROPS(9)
        C44_0_G= PROPS(11)
        C11_0_GP= PROPS(49) 
        C12_0_GP= PROPS(51)
        C44_0_GP= PROPS(53)

        burg=PROPS(13)

C        call compute_elastConst(props,nprops,theta,c11,c12,c44)
C        xmu = sqrt((c11-c12)/2.d0*c44)
        
        s_inf0= PROPS(12)
        gdot0 = props(25)
        f0    = props(26)
        bk    = props(27)
        sl    = props(28)
        pf    = props(29)
        qf    = props(30)
        alpha_oro=props(40)
        r_g   = props(42)
        flag_BS=props(44)
        f_gp=props(39)
        r_gp=props(41)
        debug=int(props(59))

        C11_0=(1-f_gp)*C11_0_G + f_gp*C11_0_GP
        C12_0=(1-f_gp)*C12_0_G + f_gp*C12_0_GP
        C44_0=(1-f_gp)*C44_0_G + f_gp*C44_0_GP

        xmu0 = sqrt((c11_0-c12_0)/2.d0*c44_0)
C       c12rat = xmu/xmu0               ! works for normal fcc or bcc materials

C       temperature dependent slip resistance is fitted with experimental observation
C       for octahedral slip Y. Ru et al. / Materials and Design 183 (2019)
C       linear fit below 800K and polynomial fit of 6th order above 800K 
        if (THETA < 800.0) then
            c12rat_100 =   0.03055*((THETA - 571.5)/165.5) + 1.045
        else        
            c12rat_100 = - 0.03766*((THETA - 1137.0)/173.1)**6.0
     &                   - 0.05099*((THETA - 1137.0)/173.1)**5.0 
     &                   + 0.2209*((THETA - 1137.0)/173.1)**4.0
     &                   + 0.1929*((THETA - 1137.0)/173.1)**3.0
     &                   - 0.5479*((THETA - 1137.0)/173.1)**2.0
     &                   - 0.3096*((THETA - 1137.0)/173.1) + 1.45
        end if        

        TAUORO = alpha_oro*xmu0*burg/r_g
        s_inf = s_inf0*sqrt(f_gp/r_gp)
C
        ARG1=(DABS(TAUALPHA)-sign(1.d0,TAUALPHA)*TAUB*flag_BS)
     &       /((SL+RESALPHA+s_inf+TAUORO)*c12rat_100)
C        ARG1=(DABS(TAUALPHA)-RESALPHA*C12RAT)/(SL*C12RAT)
C        write(*,*)ARG1,ARG1A
C        write(*,*)TAUALPHA,SL,RESALPHA
C
        IF(noel==20 .AND. npt==1 .AND. debug==1)THEN
        ! if (noel == 1 .AND. npt==1 .AND. kinc==1 .AND. time(2)==0.0) then
                print*,'#######----------#######'
                print*,'#######for element#1'
                print*,'#######----------#######'
                print*,'Driving force TAUALPHA',TAUALPHA
                print*,'Back stress TAUB',TAUB
                print*,'signed Back stress sensedTAUB',sign(1.d0,TAUALPHA)*TAUB
                print*,'orowan stress for gamma',TAUORO
                print*,'gamma prime structure stress ',s_inf
                print*,'slip resistance deacy ratio',c12rat_100
                print*,'current shear modulus',xmu
                ! print*,'shape of TAUB',SHAPE(TAUB)
        END IF

        IF(ARG1 .LE. 0.D0)THEN
          GDOTALPHA_G=0.D0
          DGDTAU_G=0.D0
        ELSE
          ARG2=max(0.d0,1.D0-ARG1**PF)
          
            GDOTALPHA_G=GDOT0*DEXP(-(F0/(BK*THETA))*
     +              (ARG2**QF))
            DGDTAU_G=(F0*QF*PF/(BK*THETA*(SL+RESALPHA+s_inf+TAUORO)*C12RAT_100))*
     +           (ARG2**(QF-1.D0))*(ARG1**(PF-1.D0))*GDOTALPHA_G
            IF(TAUALPHA.LT.0.0)GDOTALPHA_G=-GDOTALPHA_G
        
        END IF

        RETURN
        END

        SUBROUTINE calc_GAMMADOT_GP(props,nprops,noel,npt,theta,TAUALPHA,
     &  TAUB,RESALPHA,GDOTALPHA_GP,DGDTAU_GP)

        IMPLICIT REAL*8(A-H,O-Z)

        dimension props(nprops)
!  Initialize props needed for the calcuation
C        C11_0_G= PROPS(7) 
C        C12_0_G= PROPS(9)
C        C44_0_G= PROPS(11)
C        rM11_GP = PROPS(48)
C        C11_0_GP= PROPS(49) 
C        rM12_GP = PROPS(50)
C        C12_0_GP= PROPS(51)
C        rM44_GP = PROPS(52)
C        C44_0_GP= PROPS(53)

C        C11_0=(1-f_gp)*C11_0_G + f_gp*C11_0_GP
C        C12_0=(1-f_gp)*C12_0_G + f_gp*C12_0_GP
C        C44_0=(1-f_gp)*C44_0_G + f_gp*C44_0_GP

        call compute_elastConst(props,nprops,theta,c11,c12,c44)

        xmu = sqrt((c11-c12)/2.d0*c44)
C        xmu0 = sqrt((C11_0-C12_0)/2.d0*C44_0)

! evaluation of temperature dependent shear modulus of Gamma_prime for APB resistance at {111} plane
C        c11_GP= C11_0_GP +rM11_GP*theta
C        c12_GP= C12_0_GP +rM12_GP*theta
C        c44_GP= C44_0_GP +rM44_GP*theta
C        mu_gp111= sqrt((c11_GP-c12_GP)/2.d0*c44_GP)

        burg = props(13)
        APB_gp111=props(38)          ! in j/mm2 for 111 plane
        f_gp=props(39)
        r_gp=props(41)
        scale_APB= props(58)
        debug= int(props(59))

C        s_APB = scale_APB*mu_gp111*(APB_gp111/(mu_gp111*burg))**1.5 
C     &   *(f_gp*r_gp/burg)**0.5
        s_APB = scale_APB*APB_gp111/(2.d0*burg)

C       c12rat = xmu/xmu0                   ! works for normal fcc or bcc materials

C       temperature dependent slip resistance is fitted with experimental observation
C       for octahedral slip Y. Ru et al. / Materials and Design 183 (2019) 
C       linear fit below 800K and polynomial fit of 6th order above 800K        
        if (THETA < 800.0) then
            c12rat_100 =   0.03055*((THETA - 571.5)/165.5) + 1.045
        else        
            c12rat_100 = - 0.03766*((THETA - 1137.0)/173.1)**6.0
     &                   - 0.05099*((THETA - 1137.0)/173.1)**5.0 
     &                   + 0.2209*((THETA - 1137.0)/173.1)**4.0
     &                   + 0.1929*((THETA - 1137.0)/173.1)**3.0
     &                   - 0.5479*((THETA - 1137.0)/173.1)**2.0
     &                   - 0.3096*((THETA - 1137.0)/173.1) + 1.45
        end if        

        if (noel==20 .AND. npt==1 .AND. debug==1) then
        ! if (noel == 1 .AND. npt==1 .AND. kinc==1 .AND. time(2)==0.0) then
                print*,'#######----------#######'
                print*,'#######for element#1'
                print*,'#######----------#######'
                print*,'Antiphase boundary resistance',s_APB
                print*,'Driving force TAUALPHA',TAUALPHA
                print*,'Back stress TAUB',TAUB
                print*,'signed Back stress sensedTAUB',sign(1.d0,TAUALPHA)*TAUB
C                print*,'shear modulus for 111 plane',mu_gp111
                print*,'slip resistance deacy ratio',c12rat_100
                print*,'current shear modulus',xmu
                ! print*,'shape of TAUB',SHAPE(TAUB)
                ! print*,'Antiphase boundary energy',APB_gp111
                ! print*,'shear modulus of GP',mu_gp111
                ! print*,'Burgers vector mag',burg
                ! print*,'volume fraction GP',f_gp
                ! print*,'approx size of GP precipitate',r_gp
        end if

        gdot0 = props(25)
        f0    = props(26)
        bk    = props(27)
        sl    = props(28)
        pf    = props(29)
        qf    = props(30)
        flag_BS=props(44)


        ARG1=(DABS(TAUALPHA)+sign(1.d0,TAUALPHA)*TAUB*flag_BS-s_APB)
     &       /((SL+RESALPHA)*C12RAT_100)
C        write(*,*)ARG1,ARG1A
C        write(*,*)TAUALPHA,SL,RESALPHA
C
!        IF(ARG1.LE.0.D0 .OR. ARG3.LE.0.D0)THEN
        IF(ARG1.LE.0.D0)THEN
          GDOTALPHA_GP=0.D0
          DGDTAU_GP=0.D0
        ELSE
          ARG2=max(0.d0,1.D0-ARG1**PF)
          
          GDOTALPHA_GP=GDOT0*DEXP(-(F0/(BK*THETA))*
     +              (ARG2**QF))
          DGDTAU_GP=(F0*QF*PF/(BK*THETA*(SL+RESALPHA)*C12RAT_100))*
     +           (ARG2**(QF-1.D0))*(ARG1**(PF-1.D0))*GDOTALPHA_GP
          IF(TAUALPHA.LT.0.0) GDOTALPHA_GP=-GDOTALPHA_GP
        
        END IF

        RETURN
        END

        SUBROUTINE calc_GAMMADOT_GP_cube(props,nprops,noel,npt,theta,TAUALPHA,
     &  TAUB,RESALPHA,GDOTALPHA_GP,DGDTAU_GP)

        IMPLICIT REAL*8(A-H,O-Z)

        dimension props(nprops)
!  Initialize props needed for the calcuation
C        C11_0_G= PROPS(7) 
C        C12_0_G= PROPS(9)
C        C44_0_G= PROPS(11)
C        C11_0_GP= PROPS(49) 
C        C12_0_GP= PROPS(51)
C        C44_0_GP= PROPS(53)

C        C11_0=(1-f_gp)*C11_0_G + f_gp*C11_0_GP
C        C12_0=(1-f_gp)*C12_0_G + f_gp*C12_0_GP
C        C44_0=(1-f_gp)*C44_0_G + f_gp*C44_0_GP

        call compute_elastConst(props,nprops,theta,c11,c12,c44)

        xmu = sqrt((c11-c12)/2.d0*c44)
C        xmu0 = sqrt((C11_0-C12_0)/2.d0*C44_0)

        burg = props(13)
        gdot0  = props(25)
        f0_cube= props(26)*1.25
        bk     = props(27)
        pf     = props(29)
        qf     = props(30)
        mu_gp100=props(37)
        APB_gp100=props(56)          ! in j/mm2 for 100 plane
        sl_cube= props(57)
        scale_APB= props(58)
        f_gp=props(39)
        r_gp=props(41)
        debug=int(props(59))

C        s_APB = scale_APB*mu_gp100*(APB_gp100/(mu_gp100*burg))**1.5 
C     &   *(f_gp*r_gp/burg)**0.5
        s_APB = scale_APB*APB_gp100/(2.d0*burg)

C       c12rat = xmu/xmu0                     ! works for normal fcc or bcc materials

C       temperature dependent slip resistance is fitted with experimental observation
C       for cube slip Y. Ru et al. / Materials and Design 183 (2019)
C       linear fit below 750K and polynomial fit of 6th order above 750K 
        if (THETA < 750.0) then
            c12rat_111 = - 0.01332*((THETA - 511.6)/149.6) + 0.969
        else        
            c12rat_111 = - 0.01784*((THETA - 1120.0)/196.7)**6.0
     &                   - 0.02857*((THETA - 1120.0)/196.7)**5.0 
     &                   + 0.08218*((THETA - 1120.0)/196.7)**4.0
     &                   + 0.09704*((THETA - 1120.0)/196.7)**3.0
     &                   - 0.1759*((THETA - 1120.0)/196.7)**2.0
     &                   - 0.1874*((THETA - 1120.0)/196.7) + 0.9578
        end if
        if (noel==20 .AND. npt==1 .AND. debug==1) then
                print*,'#######----------#######'
                print*,'#######for element#1'
                print*,'#######----------#######'
                print*,'Antiphase boundary resistance',s_APB
                print*,'Driving force TAUALPHA',TAUALPHA
                print*,'Back stress TAUB',TAUB
                print*,'signed Back stress sensedTAUB',sign(1.d0,TAUALPHA)*TAUB
C                print*,'shear modulus for 100 plane',mu_gp100
                print*,'slip resistance deacy ratio',c12rat_111
                print*,'current shear modulus',xmu
                ! print*,'shape of TAUB',SHAPE(TAUB)
                ! print*,'Antiphase boundary energy',APB_gp111
                ! print*,'shear modulus of GP',mu_gp100
                ! print*,'Burgers vector mag',burg
                ! print*,'volume fraction GP',f_gp
                ! print*,'approx size of GP precipitate',r_gp
        end if

        flag_BS= props(44)


        ARG1=(DABS(TAUALPHA)+sign(1.d0,TAUALPHA)*TAUB*flag_BS-s_APB)
     &       /((SL_cube+RESALPHA)*C12RAT_111)
        IF(ARG1.LE.0.D0)THEN
          GDOTALPHA_GP=0.D0
          DGDTAU_GP=0.D0
        ELSE
          ARG2=max(0.d0,1.D0-ARG1**PF)
          
          GDOTALPHA_GP=GDOT0*DEXP(-(F0_cube/(BK*THETA))*
     +              (ARG2**QF))
          DGDTAU_GP=(F0_cube*QF*PF/(BK*THETA*(SL_cube+RESALPHA)*C12RAT_111))*
     +           (ARG2**(QF-1.D0))*(ARG1**(PF-1.D0))*GDOTALPHA_GP
          IF(TAUALPHA.LT.0.0) GDOTALPHA_GP=-GDOTALPHA_GP
        
        END IF

        RETURN
        END

C **********************************************************************
! Subroutine for defining the dislocation interatction matrix that
! appears in the Taylor expression for the 
c **********************************************************************
        subroutine calc_amat(nslip,props,nprops,s0,n0,amat)

        implicit none

        real*8 :: amat(nslip,nslip),props(nprops),s0(3,nslip),
     &           n0(3,nslip),t0(3,nslip)
        integer :: i,j,nslip,nprops,ii,jj,istart,jstart,igliss1,igliss2

        real*8 :: sa_dot_sb,na_dot_nb,na_dot_tb,dp1,dp2,
     &            a_self,a_copl,a_hirth,a_colli,a_gliss,a_lomer,
     &            cube_self,cube_latent

        a_self = props(17)
        a_copl = props(18)
        a_hirth = props(19)
        a_colli = props(20)
        a_gliss = props(21)
        a_lomer = props(22)
        cube_self=PROPS(54)          ! self hardening coefficient for cube slip
        cube_latent=PROPS(55)        ! latent hardening coefficient for cube slip

        amat = 0.d0
        do i = 1,nslip
          do j = 1,nslip
            sa_dot_sb = abs(dot_product(s0(:,i),s0(:,j)))
            na_dot_nb = abs(dot_product(n0(:,i),n0(:,j)))
            ! Self and/or coplanar interaction
            if (na_dot_nb > 0.999d0) then
              if (i==j) then
                amat(i,j) = a_self !0.06d0
              else 
                amat(i,j) = a_copl !0.1d0
              end if
            ! Hirth Lock
            else if ((na_dot_nb < 0.999d0) .and. 
     &               (sa_dot_sb < 0.001d0)) then
              amat(i,j) = a_hirth !0.07d0 
            ! Collinear interaction
            else if ((na_dot_nb < 0.999d0) .and. 
     &               (sa_dot_sb > 0.999d0)) then
              amat(i,j) = a_colli !0.625d0
            else
              if (i < 4) then
               istart = 1
              else if ((i>=4) .and. (i<7)) then
               istart = 4
              else if ((i>=7) .and. (i<10)) then
               istart = 7
              else
               istart = 10
              end if
              if (j < 4) then
               jstart = 1
              else if ((j>=4) .and. (j<7)) then
               jstart = 4
              else if ((j>=7) .and. (j<10)) then
               jstart = 7
              else
               jstart = 10
              end if

              igliss1 = 0
              do jj = jstart,jstart+2
                  dp1 = abs(dot_product(s0(:,i),s0(:,jj)))
                  if (dp1 > 0.999d0) then
                    igliss1 = 1
                    exit
                  end if 
              end do
              igliss2 = 0
              do ii = istart,istart+2
                  dp2 = abs(dot_product(s0(:,ii),s0(:,j)))
                  if (dp2 > 0.999d0) then
                    igliss2 = 1
                    exit
                  end if
              end do 

              ! Glissile Junction
              !if (igliss1 .or. igliss2) then
              if ( (igliss1==1) .or. (igliss2==1) ) then
                  amat(i,j) = a_gliss !0.137d0
              ! Lomer-Cottrell Lock
              else
                  amat(i,j) = a_lomer !0.122d0
              end if
            end if
            if (i > 12 .or. j > 12) then
              if (na_dot_nb > 0.999d0) then
                amat(i,j) = cube_self !1.d0
              else
                amat(i,j) = cube_latent !1.4d0
              end if
            end if
          end do
        end do
                
        return
        end subroutine calc_amat

c **********************************************************************
! Subroutine for defining the dislocation interatction matrix that
! appears in the Taylor expression for the 
        subroutine calc_dmat(nslip,props,nprops,n0,amat,dmat)

        implicit none

        real*8 :: amat(nslip,nslip),dmat(nslip,nslip),props(nprops),
     &            n0(3,nslip),na_dot_nb
        integer :: nslip,nprops,ia,ib

        real*8 :: rk_inter,rk_copl

        rk_inter = props(23)
        rk_copl =  props(24)

        do ia = 1,nslip
        do ib = 1,nslip
            na_dot_nb = abs(dot_product(n0(:,ia),n0(:,ib)))
            if (na_dot_nb > 0.99d0) then
                dmat(ia,ib) = amat(ia,ib) / rk_copl**2.d0
            else
                dmat(ia,ib) = amat(ia,ib) / rk_inter**2.d0
            end if
        end do
        end do

        return
        end subroutine calc_dmat

C **********************************************************************
        SUBROUTINE calc_HARDENING(props,nprops,nslip,theta,amat,rhossd,DGAMMA,
     &  DTIME,HAB,drhossd,xmC)

        IMPLICIT REAL*8(A-H,O-Z)
C
        DIMENSION DGAMMA(nslip),
     +            HAB(nslip,nslip),props(nprops),
     &            rhossd(nslip),rhossd_f(nslip),rl_mfp(nslip),
     &            dmat(nslip,nslip),amat(nslip,nslip),drhossd(nslip),
     &            xmC(3,nslip)
        
        real*8 :: f_gp,c11,c12,c44,burg,theta
        real*8 :: C11_0_G,C12_0_G,C44_0_G,C11_0_GP,C12_0_GP,C44_0_GP
        real*8 :: C11_0,C12_0,C44_0,xmu0

!  Initialize props needed for the calculation
        
        burg = props(13)
        y_c0 = props(15)*burg
        A_rec = props(16)
        gdot0 = props(25)
        bk    = props(27)
        f_gp  = props(39)

        C11_0_G= PROPS(7) 
        C12_0_G= PROPS(9)
        C44_0_G= PROPS(11)
        C11_0_GP= PROPS(49) 
        C12_0_GP= PROPS(51)
        C44_0_GP= PROPS(53)

        C11_0=(1-f_gp)*C11_0_G + f_gp*C11_0_GP
        C12_0=(1-f_gp)*C12_0_G + f_gp*C12_0_GP
        C44_0=(1-f_gp)*C44_0_G + f_gp*C44_0_GP

        xmu0 = sqrt((C11_0-C12_0)/2.d0*C44_0)
        call compute_elastConst(props,nprops,theta,c11,c12,c44)
        xmu = sqrt((c11-c12)/2.d0*c44)

!  Calculate the hardening matrix
        call calc_dmat(nslip,props,nprops,xmC,amat,dmat)

        rhossd_f = sqrt(matmul(amat,rhossd))
        rl_mfp = 1.d0 / sqrt(matmul(dmat,rhossd))

        do ia = 1,nslip
          do ib = 1,nslip
            y_c = y_c0*(abs(dgamma(ib))/gdot0/dtime)**(bk*theta/A_rec)
            hab(ia,ib) = xmu/2.d0/rhossd_f(ia)*amat(ia,ib)*
     &                   (1.d0/rl_mfp(ib) - 2.d0*y_c*rhossd(ib))
          end do
        end do

        drhossd = 0.d0
        do ia = 1,nslip
          y_c = y_c0*(abs(dgamma(ia))/gdot0/dtime)**(bk*theta/A_rec)
          drhossd(ia) = (1.d0/rl_mfp(ia) - 2.d0*y_c*rhossd(ia))/burg*
     &                 abs(dgamma(ia))
        end do

        RETURN 
        END      
C**************************************************************
        subroutine compute_elastConst(props,nprops,theta,c11,c12,c44)

        implicit none

! Global Variables
        integer :: nprops

      	real*8  :: theta

      	real*8  :: props(nprops)

      	real*8  :: c11,c12,c44

! Local Variables
      	real*8  :: rm11_G,rm12_G,rm44_G,c11_0_G,c12_0_G,c44_0_G,f_gp
      	real*8  :: rm11_GP,rm12_GP,rm44_GP,c11_0_GP,c12_0_GP,c44_0_GP
        

        rm11_G  = props(6)
        c11_0_G = props(7)
        rm12_G  = props(8)
        c12_0_G = props(9)
        rm44_G  = props(10)
        c44_0_G = props(11)
        f_gp = props(39)
        rM11_GP=PROPS(48)            
        c11_0_GP=PROPS(49)           
        rM12_GP=PROPS(50)            
        c12_0_GP=PROPS(51)           
        rM44_GP=PROPS(52)            
        c44_0_GP=PROPS(53)

        c11 = (1.0-f_gp)*(c11_0_G + rm11_G*theta) + f_gp*(c11_0_GP + rm11_GP*theta)
        c12 = (1.0-f_gp)*(c12_0_G + rm12_G*theta) + f_gp*(c12_0_GP + rm12_GP*theta)
        c44 = (1.0-f_gp)*(c44_0_G + rm44_G*theta) + f_gp*(c44_0_GP + rm44_GP*theta)

        return
        end subroutine
C************************************************************************
        SUBROUTINE DISP(U,KSTEP,KINC,TIME,NODE,NOEL,JDOF,COORDS)
C
        INCLUDE 'ABA_PARAM.INC'
C
       	REAL*8 L0
        DIMENSION U(3),TIME(2),COORDS(3)
C
        L0=1.D0
        EDOT=1.0D-3
C         EDOT=-2000.
C         EDOT=-2400.
C         EDOT=-2300.
C         EDOT=-2800.
C
        U(1)=L0*(DEXP(EDOT*TIME(2))-1.D0)
C
        RETURN
        END
C************************************************************************
        SUBROUTINE MPC(UE,A,JDOF,MDOF,N,JTYPE,X,U,UINIT,MAXDOF,
     *  LMPC,KSTEP,KINC,TIME,NT,NF,TEMP,FIELD,LTRAN,TRAN)
C
        INCLUDE 'ABA_PARAM.INC'
C
        DIMENSION A(N),JDOF(N),X(6,N),U(MAXDOF,N),UINIT(MAXDOF,N),
     *   TIME(2),TEMP(NT,N),FIELD(NF,NT,N),LTRAN(N),TRAN(3,3,N)
C
        XDIM=0.69993
        YDIM=0.69993
        IF(JTYPE.EQ.1)THEN
         JDOF(1)=1
         JDOF(2)=2
         A(1)=-X(2,2)
         A(2)=XDIM*YDIM/X(2,2)
         UE=((XDIM*YDIM)/(2.*X(2,2)**2))*U(2,2)
        ELSE IF(JTYPE.EQ.2)THEN
         JDOF(1)=1
         JDOF(2)=2
         A(1)=X(2,2)
         A(2)=XDIM*YDIM/X(2,2)
         UE=((-XDIM*YDIM)/(2.*X(2,2)**2))*U(2,2)
        END IF
C
        RETURN
        END
C************************************************************************
        SUBROUTINE calc_BACKSTRESS(props,nprops,nslip,taub_n,DGAMMA,
     &  TAUBTR,DTAUB)

        implicit none

        real*8 :: DTAUB(nslip),taub_n(nslip),TAUBTR(nslip)
        real*8 :: DGAMMA(nslip),props(nprops)
        integer:: nslip,nprops,debug

        real*8 :: f_gp,r_gp,r_g,rb,c11,c12,c44,burg,theta,xmu,pi
        real*8 :: C11_0_G,C12_0_G,C44_0_G,C11_0_GP,C12_0_GP,C44_0_GP
        real*8 :: C11_0,C12_0,C44_0,xmu0

        pi = 4.0*DATAN(1.0D0)
!  Initialize props needed for the calcuation
C       c11_0=PROPS(7)
C       c12_0=PROPS(9)
C       c44_0=PROPS(11)
        C11_0_G= PROPS(7)
        C12_0_G= PROPS(9)
        C44_0_G= PROPS(11)
        C11_0_GP= PROPS(49)
        C12_0_GP= PROPS(51)
        C44_0_GP= PROPS(53)

        burg = props(13)
        f_gp=props(39)
        r_gp=PROPS(41)
        r_g=PROPS(42)
        rb=PROPS(43)
        debug=int(props(59))

        C11_0=(1-f_gp)*C11_0_G + f_gp*C11_0_GP
        C12_0=(1-f_gp)*C12_0_G + f_gp*C12_0_GP
        C44_0=(1-f_gp)*C44_0_G + f_gp*C44_0_GP

        xmu0 = sqrt((C11_0-C12_0)/2.d0*C44_0)
        
        call compute_elastConst(props,nprops,theta,c11,c12,c44)
        xmu = sqrt((c11-c12)/2.d0*c44)

        TAUBTR= taub_n + DGAMMA*xmu*burg/(2.d0*pi*rb)-(r_gp/r_g)*taub_n*DABS(DGAMMA)
        DTAUB = DGAMMA*xmu*burg/(2.d0*pi*rb)-(r_gp/r_g)*TAUBTR*DABS(DGAMMA)
        RETURN
        END
