program normalmodes
!program that given a equilibrium configuration and the frequencies calculates the evolution of the system

use normal_modes_subroutines

implicit none
!xyzmat: data of the xyz matrix set in a line
!amplitudes: matrix with the values of the amplitudes provided by MOLPRO
!mass: vector with the masses of the atoms of the molecule. The masses are repeated because it's easier to
!work that way

!freq: vector with the frequencies provided by MOLPRO
!ep: phase
!nat: number of atoms of the molecule
!nnm: number of normal modes
!attyp: vector containing the atom type. 
!npoints: number of geometries that the user wants to calculate
!amp: maximum amplitude of a normal mode

real(8), dimension(:), allocatable :: xyzmat, mass, freq
real(8), dimension(:,:), allocatable ::  amplitudes
real(8) :: ep, amp, energy
integer ::  nat, nnm, npoints,nvib, i
character(len=2), dimension(:),allocatable :: attyp

!define the quantum number for vibrations
nvib=4

!the file containing the xyz matrix is open
open(10,file='eq_geometry.xyz', status='old')

!the number of atoms of the molecule is read
read(10,*) nat
read(10,*)

!the total number of points that will be calculated is:
npoints=41

!the phase is set to 0
ep=0

allocate (xyzmat(3*nat))
allocate (mass(3*nat))
allocate (attyp(nat))

!the xyz matrix is read
call xyz_reader(nat,attyp,xyzmat,mass)

!the next read contains the number of normal modes:
read(10,*) nnm
allocate (amplitudes(3*nat,nnm))
allocate (freq(nnm))

!the frequencies and amplitudes are read
call freq_am(nat,nnm,amplitudes,freq)

!some parameters of the potential are calculated
do i=1, nnm
call osc_param(i,nvib,nnm, freq,amp)

amp = 2e0 * amp
call normal_modes(i,nat,nnm,attyp,xyzmat,mass,amplitudes,freq,npoints,ep,amp)
enddo

endprogram
