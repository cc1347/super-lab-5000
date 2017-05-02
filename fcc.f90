module lattice_stuff
  implicit none

  private
  public::dp,smallest_fcc

  integer,parameter::dp=selected_real_kind(15,300)
  real(kind=dp),dimension(1:16,1:4)::smallest_fcc

  contains
    function define_smallest_fcc()
      real(kind=dp),dimension(1:16,1:4)::define_smallest_fcc
      !    atom       x               y               z
      !     O    0.5000000000    0.7500000000    0.5000000000
      !     O    0.5000000000    0.5000000000    0.0000000000
      !     O    0.0000000000    0.7500000000    0.0000000000
      !     O    0.0000000000    0.5000000000    0.5000000000
      !    Ca    0.0000000000    0.5000000000    0.0000000000
      !    Ca    0.0000000000    0.7500000000    0.5000000000
      !    Ca    0.5000000000    0.5000000000    0.5000000000
      !    Ca    0.5000000000    0.7500000000    0.0000000000
      !     O    0.5000000000    0.2500000000    0.5000000000
      !     O    0.5000000000    0.0000000000    0.0000000000
      !     O    0.0000000000    0.2500000000    0.0000000000
      !     O    0.0000000000    0.0000000000    0.5000000000
      !    Mg    0.0000000000    0.0000000000    0.0000000000
      !    Mg    0.0000000000    0.2500000000    0.5000000000
      !    Mg    0.5000000000    0.0000000000    0.5000000000
      !    Mg    0.5000000000    0.2500000000    0.0000000000
      !     O=0     Ca=1      Mg=2
      define_smallest_fcc(1:4,1)=0.0_dp
      define_smallest_fcc(5:8,1)=1.0_dp
      define_smallest_fcc(9:12,1)=0.0_dp
      define_smallest_fcc(13:16,1)=2.0_dp
      define_smallest_fcc(1:2,2)=0.5_dp
      define_smallest_fcc(3:6,2)=0.0_dp
      define_smallest_fcc(7:10,2)=0.5_dp
      define_smallest_fcc(11:14,2)=0.0_dp
      define_smallest_fcc(15:16,2)=0.5_dp
      define_smallest_fcc(1,3)=0.75_dp
      define_smallest_fcc(2,3)=0.5_dp
      define_smallest_fcc(3,3)=0.75_dp
      define_smallest_fcc(4:5,3)=0.5_dp
      define_smallest_fcc(6,3)=0.75_dp
      define_smallest_fcc(7,3)=0.5_dp
      define_smallest_fcc(8,3)=0.75_dp
      define_smallest_fcc(9,3)=0.25_dp
      define_smallest_fcc(10,3)=0.0_dp
      define_smallest_fcc(11,3)=0.25_dp
      define_smallest_fcc(12:13,3)=0.0_dp
      define_smallest_fcc(14,3)=0.25_dp
      define_smallest_fcc(15,3)=0.0_dp
      define_smallest_fcc(16,3)=0.25_dp
      define_smallest_fcc(1,4)=0.5_dp
      define_smallest_fcc(2:3,4)=0.0_dp
      define_smallest_fcc(4,4)=0.5_dp
      define_smallest_fcc(5,4)=0.0_dp
      define_smallest_fcc(6:7,4)=0.5_dp
      define_smallest_fcc(8,4)=0.0_dp
      define_smallest_fcc(9,4)=0.5_dp
      define_smallest_fcc(10:11,4)=0.0_dp
      define_smallest_fcc(12,4)=0.5_dp
      define_smallest_fcc(13,4)=0.0_dp
      define_smallest_fcc(14:15,4)=0.5_dp
      define_smallest_fcc(16,4)=0.0_dp
    end function define_smallest_fcc

end module lattice_stuff

program fcc
  use lattice_stuff

end program fcc
