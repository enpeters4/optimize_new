      module mod_read_data
C        implicit real*8(a-h,o-z)
         implicit none
         integer,parameter:: Numelem = 2500
!         real(8),parameter:: pi      = 4.*atan(1.0)
         real(8) test1, test2, test3
         real(8) ang1(Numelem),ang2(Numelem),ang3(Numelem)
      contains
         subroutine read_euler()
            implicit none
            integer check,i,nemax,noel
            integer, parameter:: pi=3.14159265
            real(8) phi1,phi2,phi3
            character(len=200) :: cwd=''
            call getoutdir(cwd, check)
!            CHARACTER(len=255) :: cwd          
!            call getcwd(cwd)
            open(unit=120,file=trim(cwd)//
     &      '/trial.txt',status='old')
!            open(unit=120,file=
!     &      'C31B_2_XZ_corner.txt',status='old')
            nemax=1
            do i=1,Numelem
               read(120,fmt=*)noel,phi1,phi2,phi3
               ang1(noel)=(180.-phi3)
               ang2(noel)=phi2
               ang3(noel)=(180.-phi1)
               if(noel.gt.nemax)then
                  nemax=noel
               endif
            enddo
!            test1 = ang1(1)*pi/180.d0     ! Euler angle theta [rad] (Canova convention)
!            test2 = ang2(1)*pi/180.d0    ! Euler angle Phi [rad]
!            test3 = ang3(1)*pi/180.d0     ! Euler angle omega [rad]  
!            write(*,*)nemax, test1,test2,test3

            close(120)
            return
         endsubroutine
      endmodule mod_read_data
