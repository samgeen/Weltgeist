module zone 
!=======================================================================
! (formerly zone.h) global (3D) data arrays
!======================================================================= 
 
 INTEGER :: imax=1000, jmax=1, kmax=1   ! Memory dimensions

 ! DIMENSION (imax,jmax,kmax)
 REAL(kind=8), ALLOCATABLE,DIMENSION(:,:,:) :: zro
 REAL(kind=8), ALLOCATABLE,DIMENSION(:,:,:) :: zpr
 REAL(kind=8), ALLOCATABLE,DIMENSION(:,:,:) :: zux
 REAL(kind=8), ALLOCATABLE,DIMENSION(:,:,:) :: zuy
 REAL(kind=8), ALLOCATABLE,DIMENSION(:,:,:) :: zuz
 REAL(kind=8), ALLOCATABLE,DIMENSION(:,:,:) :: zfl
 REAL(kind=8), ALLOCATABLE,DIMENSION(:,:,:) :: zgr
 
 ! DIMENSION imax, jmax, kmax respectively
 REAL(kind=8), ALLOCATABLE,DIMENSION(:) :: zxa, zdx, zxc
 REAL(kind=8), ALLOCATABLE,DIMENSION(:) :: zya, zdy, zyc
 REAL(kind=8), ALLOCATABLE,DIMENSION(:) :: zza, zdz, zzc

 ! Used only in setup
 REAL(kind=8) :: xmin, xmax, ymin, ymax, zmin, zmax
 
 INTEGER :: ngeomx, ngeomy, ngeomz       ! XYZ Geometry flag
 INTEGER :: nleftx, nlefty, nleftz       ! XYZ Lower Boundary Condition
 INTEGER :: nrightx,nrighty,nrightz      ! XYZ Upper Boundary Condition

end module zone
