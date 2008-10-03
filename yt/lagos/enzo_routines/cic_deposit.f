c=======================================================================
c//////////////////////  SUBROUTINE CIC_DEPOSIT  \\\\\\\\\\\\\\\\\\\\\\\
c
      subroutine cic_deposit(posx, posy, posz, ndim, npositions, 
     &                      mass, field, leftedge, 
     &                      dim1, dim2, dim3, cellsize)
c
c  PERFORMS 1/2/3D CLOUD-IN-CELL INTERPOLATION FROM FIELD TO SUMFIELD
c
c  written by: Greg Bryan
c  date:       January, 1998
c  modified1:
c
c  PURPOSE: This routine performs a three-dimension, second-order
c           interpolation from field to sumfield (without clearing sumfield
c           first) at the positions specified.
c
c  INPUTS:
c     ndim       - dimensionality
c     cellsize   - the cell size of field
c     dim1,2,3   - real dimensions of field
c     leftedge   - the left edge(s) of field
c     npositions - number of particles
c     posx,y,z   - particle positions
c     sumfield   - 1D field (length npositions) of masses
c
c  OUTPUT ARGUMENTS: 
c     field      - field to be deposited to
c
c  EXTERNALS: 
c
c  LOCALS:
c
c-----------------------------------------------------------------------
c
      implicit NONE
c
c-----------------------------------------------------------------------
c
c  argument declarations
c
      integer dim1, dim2, dim3, npositions, ndim
      real*8 posx(npositions), posy(npositions), posz(npositions),
     &        leftedge(3)
      real    mass(npositions), field(dim1, dim2, dim3), cellsize
CF2PY INTENT(INOUT) :: field(dim1, dim2, dim3)
CF2PY INTENT(HIDE),DEPEND(field) :: dim1=shape(field,0)
CF2PY INTENT(HIDE),DEPEND(field) :: dim2=shape(field,1)
CF2PY INTENT(HIDE),DEPEND(field) :: dim3=shape(field,2)
CF@PY INTENT(HIDE),DEPEND(posx) :: nposition=len(posx)
CF2PY INTENT(IN) :: posx(nposition), posy(nposition), posz(npositions)
CF2PY INTENT(IN) :: mass(npositions), cellsize, leftedge(3)
CF2PY INTENT(IN)  :: ndim
c
c  locals
c
      integer iii, jjj, kkk
      integer i1, j1, k1, n
      real    xpos, ypos, zpos, dx, dy, dz, fact
      real*8 edge1, edge2, edge3, half
      parameter (half = 0.5001)
c
c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\///////////////////////////////
c=======================================================================
c

!     write(0,*) npositions, leftedge, dim1, dim2, dim3, cellsize

      fact = 1.0/cellsize
      edge1 = real(dim1) - half
      edge2 = real(dim2) - half
      edge3 = real(dim3) - half
c
c     1D
c
      if (ndim .eq. 1) then
c
         do n=1, npositions
c
c           Compute the position of the central cell
c
            xpos = min(max((posx(n) - leftedge(1))*fact, half), edge1)
c
c           Convert this into an integer index
c
            i1  = int(xpos + 0.5)
c
c           Compute the weights
c
            dx = real(i1) + 0.5 - xpos
c
c           Interpolate from field into sumfield
c
            field(i1  ,1,1) = field(i1  ,1,1) + mass(n)*dx
            field(i1+1,1,1) = field(i1+1,1,1) + mass(n)*(1.0-dx)
c
         enddo
c
      endif
c
c     2D
c
      if (ndim .eq. 2) then
c
         do n=1, npositions
c
c           Compute the position of the central cell
c
            xpos = min(max((posx(n) - leftedge(1))*fact, half), edge1)
            ypos = min(max((posy(n) - leftedge(2))*fact, half), edge2)
c
c           Convert this into an integer index
c
            i1  = int(xpos + 0.5)
            j1  = int(ypos + 0.5)
c
c           Compute the weights
c
            dx = real(i1) + 0.5 - xpos
            dy = real(j1) + 0.5 - ypos
c
c           Interpolate from field into sumfield
c
            field(i1  ,j1  ,1) = field(i1  ,j1  ,1) +
     &                           mass(n)*     dx *     dy
            field(i1+1,j1  ,1) = field(i1+1,j1  ,1) +
     &                           mass(n)*(1.0-dx)*     dy
            field(i1  ,j1+1,1) = field(i1  ,j1+1,1) + 
     &                           mass(n)*     dx *(1.0-dy)
            field(i1+1,j1+1,1) = field(i1+1,j1+1,1) +
     &                           mass(n)*(1.0-dx)*(1.0-dy)
c
         enddo
c
      endif
c
c     3D
c
      if (ndim .eq. 3) then
c
         do n=1, npositions
c
c           Compute the position of the central cell
c
            xpos = min(max((posx(n) - leftedge(1))*fact, half), edge1)
            ypos = min(max((posy(n) - leftedge(2))*fact, half), edge2)
            zpos = min(max((posz(n) - leftedge(3))*fact, half), edge3)
c
c           Convert this into an integer index
c
            i1  = int(xpos + 0.5)
            j1  = int(ypos + 0.5)
            k1  = int(zpos + 0.5)
c
c           Compute the weights
c
            dx = real(i1) + 0.5 - xpos
            dy = real(j1) + 0.5 - ypos
            dz = real(k1) + 0.5 - zpos
c
c           Interpolate from field into sumfield
c     
            field(i1  ,j1  ,k1  ) = field(i1  ,j1  ,k1  ) +
     &                           mass(n)*     dx *     dy *    dz
            field(i1+1,j1  ,k1  ) = field(i1+1,j1  ,k1  ) +
     &                           mass(n)*(1.0-dx)*     dy *    dz
            field(i1  ,j1+1,k1  ) = field(i1  ,j1+1,k1  ) + 
     &                           mass(n)*     dx *(1.0-dy)*    dz
            field(i1+1,j1+1,k1  ) = field(i1+1,j1+1,k1  ) +
     &                           mass(n)*(1.0-dx)*(1.0-dy)*    dz
            field(i1  ,j1  ,k1+1) = field(i1  ,j1  ,k1+1) +
     &                           mass(n)*     dx *     dy *(1.0-dz)
            field(i1+1,j1  ,k1+1) = field(i1+1,j1  ,k1+1) +
     &                           mass(n)*(1.0-dx)*     dy *(1.0-dz)
            field(i1  ,j1+1,k1+1) = field(i1  ,j1+1,k1+1) + 
     &                           mass(n)*     dx *(1.0-dy)*(1.0-dz)
            field(i1+1,j1+1,k1+1) = field(i1+1,j1+1,k1+1) +
     &                           mass(n)*(1.0-dx)*(1.0-dy)*(1.0-dz)
c
         enddo

!        do kkk=1,dim3
!        write(0,'("K = ",i2)') kkk
!        do jjj=1,dim2
!        write(0,'((14(1pe8.1)))')(field(iii,jjj,kkk),iii=1,dim1)
!        end do
!        end do

      endif

      return
      end
