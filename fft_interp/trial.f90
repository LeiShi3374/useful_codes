Program jhp
Implicit None
integer :: i
integer, parameter :: &
     m=7, &    !total number of line
     n=4, &    !line to skip
     p=3      !lines to read
integer,dimension(m)::arr   !file to read

open(12,file='trial.txt',status='old')
do i=1,n
  read(12,*)arr(i)
end do
do i=1,p
  read(12,*)arr(i)
  write(*,*)arr(i)
end do
End Program jhp
