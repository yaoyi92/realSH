module run_test_mod
  implicit none
    contains
  subroutine run_test()
    use, intrinsic :: iso_c_binding
    use sqrttable 
    use ISO_FORTRAN_ENV
      implicit none

    interface
      subroutine spherical_harmonics_ref(l_max, theta, phi, Y) bind ( c, name="spherical_harmonics_ref" )
        use iso_c_binding
        integer ( c_int ), value :: l_max
        real ( c_double ), value :: theta,phi
        real ( c_double ), intent(inout):: Y(*)
      end subroutine

      subroutine spherical_harmonics(l_max, theta, phi, Y) bind ( c, name="spherical_harmonics" )
        use iso_c_binding
        integer ( c_int ), value :: l_max
        real ( c_double ), value :: theta,phi
        real ( c_double ), intent(inout):: Y(*)
      end subroutine

    end interface
  
      real*16 :: Y(441)
      real*4 :: Y_single(441)
      real*8 :: Y_double(441)
      real*16 :: Y_ref(441)
      real*8 :: t_init, t_finish
      real*16 :: theta, phi
      real*16 :: sintheta, costheta, sinphi, cosphi
      real*8 :: ylm_times(10,21)
      integer i , j
      integer n
      integer max_n
      integer n_samples

      open(1001,file='ref_data_1000')
      call init_sqrttable(1000)
      open(1002,file='error_sirius_gsl')
      open(1003,file='error_sirius_opt')
      open(1004,file='error_aims_old')
      open(1005,file='error_aims_new')
      open(1006,file='error_SHEval')
      open(1007,file='error_SHEval_quadruple')
      open(1008,file='error_SHEval_single')
      do  n = 1, 1000

        read(1001,*) theta, phi, sintheta,costheta,sinphi,cosphi, Y_ref
        !sintheta = sin(theta)
        !costheta = cos(theta)
        !sinphi = sin(phi)
        !cosphi = cos(phi)
        max_n = 20
        call spherical_harmonics_ref(max_n, dble(theta), dble(phi), Y_double)
        Y = Y_double
        write(1002,'(443ES)') theta,phi,abs(Y) - abs(Y_ref)
        call spherical_harmonics(max_n, dble(theta), dble(phi), Y_double)
        Y = Y_double
        write(1003,'(443ES)') theta,phi,abs(Y) - abs(Y_ref)
        call increment_ylm(dble(sintheta),dble(costheta),dble(sinphi),dble(cosphi),0,max_n,Y_double)
        Y = Y_double
        write(1004,'(443ES)') theta,phi,abs(Y) - abs(Y_ref)
        call increment_ylm_new(dble(sintheta),dble(costheta),dble(sinphi),dble(cosphi),0,max_n,Y_double)
        Y = Y_double
        write(1005,'(443ES)') theta,phi,abs(Y) - abs(Y_ref)
        call SHEval(max_n, dble(sintheta), dble(costheta), dble(sinphi), dble(cosphi), Y_double)
        Y = Y_double
        write(1006,'(443ES)') theta,phi,abs(Y) - abs(Y_ref)
        call SHEvalFqnormalaims(max_n, sintheta, costheta, sinphi, cosphi, Y)
        write(1007,'(443ES)') theta,phi,abs(Y) - abs(Y_ref)
        call SHEvalFsnormalaims(max_n, real(sintheta), real(costheta), real(sinphi), real(cosphi), Y_single)
        Y = Y_single
        write(1008,'(443ES)') theta,phi,abs(Y) - abs(Y_ref)
        !do i = 1, 441
        !  write(*,*) abs(Y(i)) - abs(Y_ref(i))
        !end do

      end do

      

  end subroutine run_test
end module run_test_mod
