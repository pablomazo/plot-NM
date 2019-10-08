module normal_modes_subroutines

contains 

subroutine xyz_reader(nat,attyp,xyzmat,mass)
!subroutine that reads the xyz matrix

implicit none

integer :: nat,i, a, b
real(8), dimension(3*nat) :: xyzmat, mass
character(len=2), dimension(nat) :: attyp

!the file containing the xyz matrix is open
open(10,file='eq_geometry.xyz', status='old')

!this is used to be able to write all the xyz matrix in a vector
a=1
b=3
do i=1,nat
 read(10,*) attyp(i), xyzmat(a:b), mass(a)
 mass(a+1)=mass(a)
 mass(b)=mass(a)
 a=a+3
 b=b+3
enddo

endsubroutine
!--------------------------------------------------------------------------------------------

subroutine freq_am(nat,nnm,amplitudes,freq)
!subroutine to read the frequencies and amplitudes from the input file
!norm: vector with the normalization factors for the amplitudes for each NM

implicit none

integer :: nnm, nat, i, k
real(8), dimension(3*nat,nnm) :: amplitudes
real(8), dimension(nnm) :: freq, norm
real(8), parameter :: c=29979245800d0 !cm/s it'll be used to convert the frequencies to s-1


!the file that also contains the amplitudes and frequencies is open
open(10,file='eq_geometry.xyz', status='old')

!the frequencies are read
read(10,*) freq(1:nnm)

!the frecuencies are converted to s-1
do i=1,nnm
freq(i)=freq(i)*c
enddo

!the amplitudes of the movement of each atom are read
do i=1,3*nat
 read(10,*) amplitudes(i,1:nnm)
enddo

!the normalization factors for the NM are calculated

do i=1, nnm
 do k=1,3*nat
  norm(i)=norm(i)+amplitudes(k,i)**2
 enddo
enddo

!we calculate the new normalised amplitudes

do i=1,nnm
 do k=1,3*nat
 amplitudes(k,i)=amplitudes(k,i)/(sqrt(norm(i)))
 enddo
enddo
endsubroutine
!--------------------------------------------------------------------------------------------

subroutine osc_param(i,nvib,nnm, freq,amp)
!subroutine that calculates some parameters related to the potential

implicit none

!ener: energy of a given vibrational state for the ith normal mode
!amp: amplitude of the vibration
!kvib: vibrational constant
!npoints: number of points in which the potential will be calculated

real(8) :: ener, amp, kvib, step
real(8), dimension(nnm) :: freq
integer :: nvib, nnm,i,j, npoints
character(len=18) :: filename
real(8), parameter :: Pi=acos(-1d0), hbar=1.0545718d-34 !J*s
real(8), parameter :: gmoltokg=6.022140857d-26 !converts units of g/mol to units of kg

!the vibrational constant is calculated
kvib=gmoltokg*freq(i)**2  !kg*s-2

!the energy of the oscillator is calculated
ener=hbar*freq(i)*(nvib+0.5)  !J

!the amplitude of the movement is calculated
!in case the freq has imaginary value it will be set as -freq in the input file and the amp would have an imaginary value
! as it is calculated as the root square.

if (ener .gt. 0) then
 amp=sqrt(2*ener/kvib) !m
else 
 amp=sqrt(2*(-ener)/kvib) !m
endif

!the potential will be calculated in the range of values [-2*amp,2*amp]
npoints=301
step=(4*amp)/(npoints-1)

!the files in which the potential info will be saved are created
if (i .lt. 10) then
 write(filename, '("Potential0",I1,".dat")') i
else
 write(filename, '("Potential",I2,".dat")') i
endif
open(11,file=filename,status='replace')

write(11,*) '#         Q/A                           E/J'
do j=1,npoints
 if (ener .gt. 0) then
  write(11,*) (-2*amp+(j-1)*step)*1d10, 0.5*kvib*(-2*amp+(j-1)*step)**2
 else
  write(11,*) (-2*amp+(j-1)*step)*1d10, -0.5*kvib*(-2*amp+(j-1)*step)**2
 endif
enddo
endsubroutine
!--------------------------------------------------------------------------------------------

subroutine normal_modes(i,nat,nnm,attyp,xyzmat,mass,amplitudes,freq,npoints,ep,amp)
!subroutine that calculates the movement of the atoms and saves it in a file

!time: time for the function to cover half of the vibration
!inct: increment of time
!x,y,z. new coordinates of the atom at the time we are making the calculation
!i: the index of the normal mode
implicit none

integer :: nat,nnm, npoints, i,j,k
real(8), dimension(3*nat) :: xyzmat, mass
real(8), dimension(3*nat,nnm) :: amplitudes
real(8), dimension(nnm) :: freq
real(8) :: time, inct, ep, x,y,z, amp
character(len=2), dimension(nat) :: attyp
character(len=18) :: filename
real(8), parameter :: Pi=acos(-1d0)
 

time=Pi/(freq(i))
inct=time/(npoints-1)

!the file in which the info will be save is created
write(filename, '("Normal_modes",I2.2,".xyz")') i
open(12,file=filename,status='replace')


do j=1,npoints
write(12,*) nat
write(12,*) 'Q=', amp*cos(freq(i)*inct*(j-1)+ep)*1e10
 do k=1,3*nat,3
  x=1/(sqrt(mass(k)))*amp*1d10*amplitudes(k,i)*cos(freq(i)*inct*(j-1)+ep)+xyzmat(k)
  y=1/(sqrt(mass(k+1)))*amp*1d10*amplitudes(k+1,i)*cos(freq(i)*inct*(j-1)+ep)+xyzmat(k+1)
  z=1/(sqrt(mass(k+2)))*amp*1d10*amplitudes(k+2,i)*cos(freq(i)*inct*(j-1)+ep)+xyzmat(k+2)

 write(12,'(A3,3(F12.9,X))') attyp((3+k)/3),x,y,z
 enddo
enddo

endsubroutine

endmodule 
