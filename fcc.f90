module lattice_stuff
  implicit none

  private
  public::dp,smallest_fcc
  public::define_smallest_fcc
  public::fcc_lattice

  integer,parameter::dp=selected_real_kind(15,300)
  real(kind=dp),dimension(1:16,1:4)::smallest_fcc

  contains
    subroutine define_smallest_fcc()
      !creates the following face centered cubic lattice
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
      smallest_fcc(1:4,1)=0.0_dp
      smallest_fcc(5:8,1)=1.0_dp
      smallest_fcc(9:12,1)=0.0_dp
      smallest_fcc(13:16,1)=2.0_dp
      smallest_fcc(1:2,2)=0.5_dp
      smallest_fcc(3:6,2)=0.0_dp
      smallest_fcc(7:10,2)=0.5_dp
      smallest_fcc(11:14,2)=0.0_dp
      smallest_fcc(15:16,2)=0.5_dp
      smallest_fcc(1,3)=0.75_dp
      smallest_fcc(2,3)=0.5_dp
      smallest_fcc(3,3)=0.75_dp
      smallest_fcc(4:5,3)=0.5_dp
      smallest_fcc(6,3)=0.75_dp
      smallest_fcc(7,3)=0.5_dp
      smallest_fcc(8,3)=0.75_dp
      smallest_fcc(9,3)=0.25_dp
      smallest_fcc(10,3)=0.0_dp
      smallest_fcc(11,3)=0.25_dp
      smallest_fcc(12:13,3)=0.0_dp
      smallest_fcc(14,3)=0.25_dp
      smallest_fcc(15,3)=0.0_dp
      smallest_fcc(16,3)=0.25_dp
      smallest_fcc(1,4)=0.5_dp
      smallest_fcc(2:3,4)=0.0_dp
      smallest_fcc(4,4)=0.5_dp
      smallest_fcc(5,4)=0.0_dp
      smallest_fcc(6:7,4)=0.5_dp
      smallest_fcc(8,4)=0.0_dp
      smallest_fcc(9,4)=0.5_dp
      smallest_fcc(10:11,4)=0.0_dp
      smallest_fcc(12,4)=0.5_dp
      smallest_fcc(13,4)=0.0_dp
      smallest_fcc(14:15,4)=0.5_dp
      smallest_fcc(16,4)=0.0_dp
    end subroutine define_smallest_fcc

    function fcc_lattice(layers)
      integer::layers,n
      real(kind=dp),dimension(1:16*layers,1:4)::fcc_lattice

      call define_smallest_fcc()

      n=1
      do while (n<=layers) !adds the layers of MgO
        fcc_lattice((layers*16)-(8*(n-1)+7),1:4)=smallest_fcc(9,1:4)
        fcc_lattice((layers*16)-(8*(n-1)+6),1:4)=smallest_fcc(10,1:4)
        fcc_lattice((layers*16)-(8*(n-1)+5),1:4)=smallest_fcc(11,1:4)
        fcc_lattice((layers*16)-(8*(n-1)+4),1:4)=smallest_fcc(12,1:4)
        fcc_lattice((layers*16)-(8*(n-1)+3),1:4)=smallest_fcc(13,1:4)
        fcc_lattice((layers*16)-(8*(n-1)+2),1:4)=smallest_fcc(14,1:4)
        fcc_lattice((layers*16)-(8*(n-1)+1),1:4)=smallest_fcc(15,1:4)
        fcc_lattice((layers*16)-(8*(n-1)),1:4)=smallest_fcc(16,1:4)
        fcc_lattice((layers*16)-(8*(n-1)+7):layers*16-((8*n-1)),3)=fcc_lattice(layers*16-((8*n-1)+7):layers*16-((8*n-1)),3)+(0.5_dp*(n-1))
        n=n+1
      end do

      n=1
      do while (n<=layers) !adds the layers of CaO
        fcc_lattice(((layers*16)/2)-(8*(n-1)+7),1:4)=smallest_fcc(1,1:4)
        fcc_lattice(((layers*16)/2)-(8*(n-1)+6),1:4)=smallest_fcc(2,1:4)
        fcc_lattice(((layers*16)/2)-(8*(n-1)+5),1:4)=smallest_fcc(3,1:4)
        fcc_lattice(((layers*16)/2)-(8*(n-1)+4),1:4)=smallest_fcc(4,1:4)
        fcc_lattice(((layers*16)/2)-(8*(n-1)+3),1:4)=smallest_fcc(5,1:4)
        fcc_lattice(((layers*16)/2)-(8*(n-1)+2),1:4)=smallest_fcc(6,1:4)
        fcc_lattice(((layers*16)/2)-(8*(n-1)+1),1:4)=smallest_fcc(7,1:4)
        fcc_lattice(((layers*16)/2)-(8*(n-1)),1:4)=smallest_fcc(8,1:4)
        fcc_lattice(((layers*16)/2)-(8*(n-1)+7):layers*16-((8*n-1)),3)=fcc_lattice(((layers*16)/2)-((8*n-1)+7):((layers*16)/2)-((8*n-1)),3)+((0.5_dp*(n-1))+1.0_dp)
        n=n+1
      end do
    end function fcc_lattice

end module lattice_stuff


module coord_file
  implicit none
  use larrice_stuff

  private
  public::create_file

  contains 

    subroutine create_file(lattice,layers)
      integer :: i,layers,n,istat
      real(kind=dp),dimension(1:16*layers,1:4)::lattice

      n=16*layers !nuber of atoms in the lattice

      open(unit=11,file='data.dat',status='replace', action="write",iostat=istat)
      if (istat .ne. 0) stop "Error opening file"
        do i=1,n
          if(int(lattice(i,1))==0)then
            write(unit=11,fmt=*,iostat=istat) "O",lattice(i,2),lattice(i,3),lattice(i,4)
            if (istat .ne. 0) stop "Error opening file"
          elseif(int(lattice(i,1))==1)then
            write(unit=11,fmt=*,iostat=istat) "Ca",lattice(i,2),lattice(i,3),lattice(i,4)
            if (istat .ne. 0) stop "Error opening file"
          elseif(int(lattice(i,1))==2)then
            write(unit=11,fmt=*,iostat=istat) "Mg",lattice(i,2),lattice(i,3),lattice(i,4)
            if (istat .ne. 0) stop "Error opening file"
          end if
        end do
      close(11)
    end subroutine create_file

end module coord_file

program fcc
  use lattice_stuff
  use coord_file
end program fcc
