!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!!        FILE: BC-Channel-flow.f90
!!!      AUTHOR: ??
!!!    MODIFIED: Paul Bartholomew
!!! DESCRIPTION: This module describes the channel flow.
!!!   CHANGELOG: [2019-02-19] Making module private by default
!!               [2019-02-19] Turning file into a module
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module channel

  USE decomp_2d
  USE variables
  USE param

  IMPLICIT NONE

  integer :: FS
  character(len=100) :: fileformat
  character(len=1),parameter :: NL=char(10) !new line character

  !probes
  integer, save :: nprobes, ntimes1, ntimes2
  integer, save, allocatable, dimension(:) :: rankprobes, nxprobes, nyprobes, nzprobes

  real(mytype),save,allocatable,dimension(:) :: usum,vsum,wsum,uusum,uvsum,uwsum,vvsum,vwsum,wwsum

  PRIVATE ! All functions/subroutines private by default
  PUBLIC :: init_channel, init_channel_ncdfr, boundary_conditions_channel, postprocess_channel, &
       momentum_forcing_channel, &
       geomcomplex_channel

contains

  subroutine init_channel (ux1,uy1,uz1,ep1,phi1)

    USE decomp_2d
    USE decomp_2d_io
    USE variables
    USE param
    USE MPI

    implicit none

    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux1,uy1,uz1,ep1
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),numscalar) :: phi1

    real(mytype) :: y,r,um,r3,x,z,h,ct
    real(mytype) :: cx0,cy0,cz0,hg,lg
    integer :: k,j,i,fh,ierror,ii,is,it,code
    integer (kind=MPI_OFFSET_KIND) :: disp

    integer, dimension (:), allocatable :: seed
   
    if (iscalar==1) then

       phi1(:,:,:,:) = zero !change as much as you want

    endif
    ux1=zero;uy1=zero;uz1=zero
    if (iin.ne.0) then
       call system_clock(count=code)
       if (iin.eq.2) code=0
       call random_seed(size = ii)
       call random_seed(put = code+63946*nrank*(/ (i - 1, i = 1, ii) /))

       call random_number(ux1)
       call random_number(uy1)
       call random_number(uz1)
    endif

    !modulation of the random noise + initial velocity profile
    do k=1,xsize(3)
       do j=1,xsize(2)
          if (istret.eq.0) y=real(j+xstart(2)-1-1,mytype)*dy-yly/two
          if (istret.ne.0) y=yp(j+xstart(2)-1)-yly/two
          um=exp(-zptwo*y*y)
          do i=1,xsize(1)
             ux1(i,j,k)=init_noise*um*(two*ux1(i,j,k)-one)+one-y*y
             uy1(i,j,k)=init_noise*um*(two*uy1(i,j,k)-one)
             uz1(i,j,k)=init_noise*um*(two*uz1(i,j,k)-one)
          enddo
       enddo
    enddo

    !INIT FOR G AND U=MEAN FLOW + NOISE
    do k=1,xsize(3)
       do j=1,xsize(2)
          do i=1,xsize(1)
             ux1(i,j,k)=ux1(i,j,k)+bxx1(j,k)
             uy1(i,j,k)=uy1(i,j,k)+bxy1(j,k)
             uz1(i,j,k)=uz1(i,j,k)+bxz1(j,k)
          enddo
       enddo
    enddo

#ifdef DEBG
    if (nrank .eq. 0) print *,'# init end ok'
#endif

    return
  end subroutine init_channel
  !############################################################################
  subroutine init_channel_ncdfr (ux1,uy1,uz1,ep1)
    USE decomp_2d
    USE decomp_2d_io
    USE variables
    USE param
    USE MPI

    implicit none
    integer,parameter :: ndim=3

    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux1,uy1,uz1,ep1
   !  character(len=2) :: dimname(ndim)
   !  integer :: dimlen(ndim)

    ux1=zero;uy1=zero;uz1=zero
    call read_restart(f_input,mpid=.TRUE.)

#ifdef DEBG
    if (nrank .eq. 0) print *,'# init end ok'
#endif
   !  dimlen=(/nx,ny,nz/)
   !  dimname=(/'x1','y1','z1'/)
   !  call write_restart('Re550_restarted.nc', dimname, dimlen)

!call write_restart('restart.nc',dimname,dimlen)
    return
  end subroutine init_channel_ncdfr  
  !********************************************************************
  subroutine boundary_conditions_channel (ux,uy,uz,phi)

    USE param
    USE variables
    USE decomp_2d

    implicit none

    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux,uy,uz
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),numscalar) :: phi
!!$  real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ut

    real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: gx
    real(mytype) :: x, y, z
    integer :: i, j, k, is

    call transpose_x_to_y(ux,gx)
    call channel_flrt(gx,two/three)
    call transpose_y_to_x(gx,ux)

    if (iscalar.ne.0) then
       if (nclxS1.eq.2) then
          i = 1
          phi(i,:,:,:) = zero
       endif
       if (nclxSn.eq.2) then
          i = xsize(1)
          phi(i,:,:,:) = phi(i - 1,:,:,:)
       endif

       if ((nclyS1.eq.2).and.(xstart(2).eq.1)) then
          !! Generate a hot patch on bottom boundary
          do k = 1, xsize(3)
             z = real(k + xstart(3) - 2, mytype) * dz - half * zlz
             if (abs(z).lt.zlz/four) then
                j = 1
                do i = 1, xsize(1)
                   x = real(i + xstart(1) - 2, mytype) * dx
                   if ((x.gt.0.1*xlx).and.(x.lt.0.3*xlx)) then
                      do is = 1, numscalar
                         phi(i, j, k, is) = one
                      enddo
                   else
                      do is = 1, numscalar
                         phi(i, j, k, is) = zero
                      enddo
                   endif
                enddo
             endif
          enddo
       endif
    endif

    return
  end subroutine boundary_conditions_channel

  !********************************************************************
  !
  subroutine channel_flrt (ux,constant)
    !
    !********************************************************************

    USE decomp_2d
    USE decomp_2d_poisson
    USE variables
    USE param
    USE var
    USE MPI

    implicit none

    real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: ux
    real(mytype) :: constant

    integer :: j,i,k,code
    real(mytype) :: can,ut3,ut,ut4

    ut3=zero
    do k=1,ysize(3)
       do i=1,ysize(1)
          ut=zero
          do j=1,ny-1
             if (istret.eq.0) then
                ut=ut+dy*(ux(i,j+1,k)-half*(ux(i,j+1,k)-ux(i,j,k)))
             else
                ut=ut+(yp(j+1)-yp(j))*(ux(i,j+1,k)-half*(ux(i,j+1,k)-ux(i,j,k)))
             endif
          enddo
          ut=ut/yly
          ut3=ut3+ut
       enddo
    enddo
    ut3=ut3/(real(nx*nz,mytype))

    call MPI_ALLREDUCE(ut3,ut4,1,real_type,MPI_SUM,MPI_COMM_WORLD,code)

    can=-(constant-ut4)

    if (nrank==0) print *,nrank,'UT',ut4,can

    do k=1,ysize(3)
       do i=1,ysize(1)
          do j=2,ny-1
             ux(i,j,k)=ux(i,j,k)-can
          enddo
       enddo
    enddo

    return
  end subroutine channel_flrt
  !********************************************************************

  !############################################################################
  subroutine init_post(ep1)

    USE MPI

    real(mytype),intent(in),dimension(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3)) :: ep1
    real(mytype) :: x, xprobes, yprobes, zprobes
    integer :: i,j,k,code
    character :: a

#ifdef DEBG
    if (nrank .eq. 0) print *,'# init_post start'
#endif

    allocate(usum(ysize(2)),vsum(ysize(2)),wsum(ysize(2)))
    allocate(uusum(ysize(2)),uvsum(ysize(2)),uwsum(ysize(2)))
    allocate(vvsum(ysize(2)),vwsum(ysize(2)),wwsum(ysize(2)))
    usum=zero;vsum=zero;wsum=zero
    uusum=zero;uvsum=zero;uwsum=zero
    vvsum=zero;vwsum=zero;wwsum=zero
    ntimes1 = 0
    ntimes2 = 0
    nprobes  = 0

    !probes
    !WORK X-PENCILS
    open(10,file='probes.prm',status='unknown',form='formatted')
    read (10,*) nprobes
    read (10,*) a
    if (nprobes .gt. 0) then
       allocate(nxprobes(nprobes), nyprobes(nprobes), nzprobes(nprobes), rankprobes(nprobes))
       rankprobes(:)=0
       do i=1, nprobes
          read (10,*) xprobes, yprobes, zprobes
          !x
          if (nclx) then
             nxprobes(i)=int(xprobes/dx)
          else
             nxprobes(i)=int(xprobes/dx+1)
          end if
          !y
          if (ncly) then
             nyprobes(i)=int(yprobes/dy)
          else
             nyprobes(i)=int(yprobes/dy+1)
          end if
          !z
          if (nclz) then
             nzprobes(i)=int(zprobes/dz)
          else
             nzprobes(i)=int(zprobes/dz+1)
          end if
          if       (xstart(1) .le. nxprobes(i) .and. nxprobes(i) .le. xend(1)) then
             if    (xstart(2) .le. nyprobes(i) .and. nyprobes(i) .le. xend(2)) then
                if (xstart(3) .le. nzprobes(i) .and. nzprobes(i) .le. xend(3)) then
                   rankprobes(i)=1
                endif
             endif
          endif
       enddo
    endif
    close(10)

#ifdef DEBG 
    if (nrank .eq. 0) print *,'# init_post ok'
#endif

  end subroutine init_post
  !############################################################################
  subroutine postprocess_channel(ux1,uy1,uz1,pp3,phi1,ep1) !By Felipe Schuch

    USE MPI
    USE decomp_2d
    USE decomp_2d_io
    USE var, only : umean,vmean,wmean,pmean,uumean,vvmean,wwmean,uvmean,uwmean,vwmean,tmean
    USE var, only : phimean, phiphimean
    USE var, only : ta1, pp1, di1
    USE var, only : ppi3, dip3
    USE var, only : pp2, ppi2, dip2
    
    USE var, ONLY : nxmsize, nymsize, nzmsize
    USE param, ONLY : npress

    real(mytype),intent(in),dimension(xsize(1),xsize(2),xsize(3)) :: ux1, uy1, uz1, ep1
    real(mytype),intent(in),dimension(xsize(1),xsize(2),xsize(3),numscalar) :: phi1
    real(mytype), dimension(ph1%zst(1):ph1%zen(1), ph1%zst(2):ph1%zen(2), nzmsize, npress), intent(in) :: pp3
    character(len=30) :: filename

    integer :: is

    return
  end subroutine postprocess_channel
  !############################################################################
  subroutine write_probes(ux1,uy1,uz1,phi1) !By Felipe Schuch

    real(mytype),intent(in),dimension(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3)) :: ux1, uy1, uz1
    real(mytype),intent(in),dimension(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3),numscalar) :: phi1

    integer :: i
    character(len=30) :: filename
    FS = 1+3+numscalar !Number of columns
    write(fileformat, '( "(",I4,"(E14.6),A)" )' ) FS
    FS = FS*14+1  !Line width

    do i=1, nprobes
       if (rankprobes(i) .eq. 1) then
          write(filename,"('./probe',I4.4)") i
          open(67,file=trim(filename),status='unknown',form='formatted'&
               ,access='direct',recl=FS)
          write(67,fileformat,rec=itime) t,&                         !1
               ux1(nxprobes(i),nyprobes(i),nzprobes(i)),&            !2
               uy1(nxprobes(i),nyprobes(i),nzprobes(i)),&            !3
               uz1(nxprobes(i),nyprobes(i),nzprobes(i)),&            !4
               phi1(nxprobes(i),nyprobes(i),nzprobes(i),:),&         !numscalar
               NL                                                    !+1
          close(67)
       endif
    enddo

  end subroutine write_probes
  !############################################################################

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!
  !!  SUBROUTINE: momentum_forcing
  !!      AUTHOR: Paul Bartholomew
  !! DESCRIPTION: Applies rotation for t < spinup_time.
  !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE momentum_forcing_channel(dux1, duy1, ux1, uy1)

    IMPLICIT NONE

    REAL(mytype), INTENT(IN), DIMENSION(xsize(1), xsize(2), xsize(3)) :: ux1, uy1
    REAL(mytype), DIMENSION(xsize(1), xsize(2), xsize(3), ntime) :: dux1, duy1

    if (itime.lt.spinup_time) then
       if (nrank==0) print *,'Rotating turbulent channel at speed ',wrotation
       dux1(:,:,:,1) = dux1(:,:,:,1) - wrotation*uy1(:,:,:)
       duy1(:,:,:,1) = duy1(:,:,:,1) + wrotation*ux1(:,:,:)
    endif

  ENDSUBROUTINE momentum_forcing_channel

  subroutine geomcomplex_channel(epsi,nxi,nxf,ny,nyi,nyf,nzi,nzf,yp,remp)

    use decomp_2d, only : mytype
    use param, only : zero, one, two
    use ibm

    implicit none

    integer                    :: nxi,nxf,ny,nyi,nyf,nzi,nzf
    real(mytype),dimension(nxi:nxf,nyi:nyf,nzi:nzf) :: epsi
    real(mytype),dimension(ny) :: yp
    real(mytype)               :: remp
    integer                    :: j
    real(mytype)               :: ym
    real(mytype)               :: zeromach
    real(mytype)               :: h

    epsi(:,:,:) = zero
    h = (yly - two) / two

    zeromach=one
    do while ((one + zeromach / two) .gt. one)
       zeromach = zeromach/two
    end do
    zeromach = 1.0e1*zeromach

    do j=nyi,nyf
       ym=yp(j)
       if ((ym.le.h).or.(ym.ge.(h+two))) then
          epsi(:,j,:)=remp
       endif
    enddo

    return
  end subroutine geomcomplex_channel

!call read_restart(filename,mpid=.TRUE.)

!dimlen=(/nx,ny,nz/)
!dimname=(/'x1','y1','z1'/)
!call write_restart('restart.nc',dimname,dimlen)


  subroutine read_restart(file_name, mpid)
   ! -----------------------------------------------------------------------
   ! io : read 3d variable in a netcdf file
   ! -----------------------------------------------------------------------
   ! Hussein RKEIN
   ! 07/2019
   
    use netcdf
    USE var, ONLY : ux1, uy1, uz1
    USE param
    USE variables
    USE decomp_2d
    USE MPI
    implicit none

    character(len=*),intent(in) :: file_name
    integer :: varid(3),i
    integer :: ncid

    integer,parameter :: ndim=3
    integer :: dimid(ndim),dim_len_check
    integer :: dimt(3),coord(3,2)
    integer :: startv(3),countv(3)
    logical,optional :: mpid
    
    !-> open file
    if (mpid) then
            call io_check(nf90_open(path=file_name,&
!!!                       mode=IOR(NF90_WRITE,NF90_MPIPOSIX),ncid=ncid,&
                  mode=IOR(NF90_NOWRITE,NF90_MPIIO),ncid=ncid,&
                  comm=MPI_COMM_WORLD,info=MPI_INFO_NULL))
    else
         call io_check(nf90_open(path=file_name,mode=nf90_nowrite,ncid=ncid))
    endif

!-> get variable id
    call io_check(nf90_inq_varid(ncid,'velocity_x',varid(1)))
    call io_check(nf90_inq_varid(ncid,'velocity_y',varid(2)))
    call io_check(nf90_inq_varid(ncid,'velocity_z',varid(3)))


    !-> recompute dimensions, start and count if mpi
    call mpi_global_coord(dimt,coord,'x')
    do i=1,3
            startv(i)=coord(i,1)
            countv(i)=coord(i,2)
    enddo

    call io_check(nf90_var_par_access(ncid, varid(1),nf90_collective))
    call io_check(nf90_var_par_access(ncid, varid(2),nf90_collective))
    call io_check(nf90_var_par_access(ncid, varid(3),nf90_collective))
   

   !-> read field variable
    call io_check(nf90_get_var(ncid,varid(1),ux1,start=startv,count=countv))
    call io_check(nf90_get_var(ncid,varid(2),uy1,start=startv,count=countv))
    call io_check(nf90_get_var(ncid,varid(3),uz1,start=startv,count=countv))



    !-> close file
    call io_check(nf90_close(ncid))
   
   end subroutine read_restart

   subroutine write_restart(file_name,dim_name,dim_len)
      ! *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
      ! Save restart file in netcdf format
      ! Hussein RKEIN
      ! 11/2019
      !
      use netcdf
      USE var, ONLY : ux1, uy1, uz1
      USE param
      USE variables
      USE decomp_2d
      USE MPI
      
                implicit none
              character(len=*),intent(in) :: file_name !,var_name
              integer,parameter           :: ndim=3
              character(len=*),intent(in) :: dim_name(ndim)
              !type(mpi_data),optional :: mpid
      
              logical :: file_exist
              integer :: varid(3),i
              integer :: ncid
              integer :: dim_len(ndim),dimid(ndim),dim_len_check
              integer :: dimt(3),coord(3,2)
              integer :: start1(3),count1(3)
              integer :: start3(3),count3(3)
      
              call io_check(nf90_create_par(path=file_name,&
      !!!               cmode=IOR(NF90_NETCDF4,NF90_MPIPOSIX),ncid=ncid,&
                      cmode=IOR(NF90_NETCDF4,NF90_MPIIO),ncid=ncid,&
                      comm=MPI_COMM_WORLD,info=MPI_INFO_NULL))
      
              start1 =(/xstart(1),xstart(2),xstart(3)/)
              count1 =(/xend(1),xend(2),xend(3)/)
              count1 = count1-start1+1
              start3 =(/phG%zst(1),phG%zst(2),phG%zst(3)/)
              count3 =(/phG%zen(1),phG%zen(2),phG%zen(3)/)
              count3 = count3-start3+1
      
              !-> create/add dimensions
              do i=1,ndim
                      if (nf90_inq_dimid(ncid,dim_name(i),dimid(i))/=nf90_noerr) then
                              call io_check(nf90_def_dim(ncid,dim_name(i),dim_len(i),dimid(i)))
                      else
                             call io_check(nf90_inquire_dimension(ncid,dimid(i),len=dim_len_check))
                              if (dim_len_check/=dim_len(i)) &
                                    !   call error_stop("NETCDF Error : wrong dimensions")
                                 print *, 'error'
                      endif
              enddo
      
              call io_check(nf90_def_var(ncid,'ux1',nf90_float,dimid,varid(1)))
              call io_check(nf90_def_var(ncid,'uy1',nf90_float,dimid,varid(2)))
              call io_check(nf90_def_var(ncid,'uz1',nf90_float,dimid,varid(3)))
            !   call io_check(nf90_def_var(ncid,'pp1',nf90_float,dimid,varid(4)))
            !   call io_check(nf90_def_var(ncid,'gx1',nf90_float,dimid,varid(5)))
            !   call io_check(nf90_def_var(ncid,'gy1',nf90_float,dimid,varid(6)))
            !   call io_check(nf90_def_var(ncid,'gz1',nf90_float,dimid,varid(7)))
            !   call io_check(nf90_def_var(ncid,'hx1',nf90_float,dimid,varid(8)))
            !   call io_check(nf90_def_var(ncid,'hy1',nf90_float,dimid,varid(9)))
            !   call io_check(nf90_def_var(ncid,'hz1',nf90_float,dimid,varid(10)))
            !   call io_check(nf90_def_var(ncid,'px1',nf90_float,dimid,varid(11)))
            !   call io_check(nf90_def_var(ncid,'py1',nf90_float,dimid,varid(12)))
            !   call io_check(nf90_def_var(ncid,'pz1',nf90_float,dimid,varid(13)))
            !   call io_check(nf90_def_var(ncid,'pp3',nf90_float,dimid,varid(14)))
      
      !	Set collective access on this variable. This will cause all
      !	reads/writes to happen together on every processor.
              call io_check(nf90_var_par_access(ncid, varid(1),nf90_collective))
              call io_check(nf90_var_par_access(ncid, varid(2),nf90_collective))
              call io_check(nf90_var_par_access(ncid, varid(3),nf90_collective))
            !   call io_check(nf90_var_par_access(ncid, varid(4),nf90_collective))
            !   call io_check(nf90_var_par_access(ncid, varid(5),nf90_collective))
            !   call io_check(nf90_var_par_access(ncid, varid(6),nf90_collective))
            !   call io_check(nf90_var_par_access(ncid, varid(7),nf90_collective))
            !   call io_check(nf90_var_par_access(ncid, varid(8),nf90_collective))
            !   call io_check(nf90_var_par_access(ncid, varid(9),nf90_collective))
            !   call io_check(nf90_var_par_access(ncid, varid(10),nf90_collective))
            !   call io_check(nf90_var_par_access(ncid, varid(11),nf90_collective))
            !   call io_check(nf90_var_par_access(ncid, varid(12),nf90_collective))
            !   call io_check(nf90_var_par_access(ncid, varid(13),nf90_collective))
            !   call io_check(nf90_var_par_access(ncid, varid(14),nf90_collective))
      
              !-> end of definition
              call io_check(nf90_enddef(ncid))
      
              !-> write field variable
              call io_check(nf90_put_var(ncid,varid(1),ux1,start=start1,count=count1))
              call io_check(nf90_put_var(ncid,varid(2),uy1,start=start1,count=count1))
              call io_check(nf90_put_var(ncid,varid(3),uz1,start=start1,count=count1))
            !   call io_check(nf90_put_var(ncid,varid(4),tb1,start=start1,count=count1))  !tb1<<pp3
            !   call io_check(nf90_put_var(ncid,varid(5),gx1,start=start1,count=count1))
            !   call io_check(nf90_put_var(ncid,varid(6),gy1,start=start1,count=count1))
            !   call io_check(nf90_put_var(ncid,varid(7),gz1,start=start1,count=count1))
            !   call io_check(nf90_put_var(ncid,varid(8),hx1,start=start1,count=count1))
            !   call io_check(nf90_put_var(ncid,varid(9),hy1,start=start1,count=count1))
            !   call io_check(nf90_put_var(ncid,varid(10),hz1,start=start1,count=count1))
            !   call io_check(nf90_put_var(ncid,varid(11),px1,start=start1,count=count1))
            !   call io_check(nf90_put_var(ncid,varid(12),py1,start=start1,count=count1))
            !   call io_check(nf90_put_var(ncid,varid(13),pz1,start=start1,count=count1))
            !   call io_check(nf90_put_var(ncid,varid(14),pp3,start=start3,count=count3))
      
              !-> close file
              call io_check(nf90_close(ncid))
   end subroutine write_restart

   subroutine io_check(status)
   ! -----------------------------------------------------------------------
   ! io : check netcdf error output and stop code if needed
   ! -----------------------------------------------------------------------
   ! Matthieu Marquillie
   ! 06/2011
   !
      use netcdf
      integer,intent(in) :: status

      if(status /= nf90_noerr) then
               print*,trim(nf90_strerror(status))
               print'(a)',"Netcdf Error : aborting"
               stop
      end if
   end subroutine io_check
   
   subroutine mpi_global_coord(dim,coord,stencil)
   !------------------------------------------------------------------------
   ! md : mpi write coord
   !------------------------------------------------------------------------
   ! Matthieu Marquillie
   ! 09/2012
   ! Modified by Ilkay Solak 10/2015
   !
   USE param
   USE variables
   USE decomp_2d
   USE MPI
   
             implicit none
   !	type(mpi_data) :: mpid
           integer,intent(out) :: dim(3),coord(3,2)
           character(*) :: stencil
           TYPE(DECOMP_INFO) :: decomp
   
           !-> compute total number of points
           dim(1)=nx
           dim(2)=ny
           dim(3)=nz
   
           if (stencil=='x') then
                   !-> compute global indices
                   coord(1,1)=xstart(1)
                   coord(2,1)=xstart(2)
                   coord(3,1)=xstart(3)
   
                   !-> compute dimensions to write
                   coord(1,2)=xend(1)-xstart(1)+1
                   coord(2,2)=xend(2)-xstart(2)+1
                   coord(3,2)=xend(3)-xstart(3)+1
   
           elseif (stencil=='y') then
                   !-> compute global indices
                   coord(1,1)=ystart(1)
                   coord(2,1)=ystart(2)
                   coord(3,1)=ystart(3)
   
                   !-> compute dimensions to write
                   coord(1,2)=yend(1)-ystart(1)+1
                   coord(2,2)=yend(2)-ystart(2)+1
                   coord(3,2)=yend(3)-ystart(3)+1
           elseif (stencil=='z') then
                   !-> compute global indices
                   coord(1,1)=zstart(1)
                   coord(2,1)=zstart(2)
                   coord(3,1)=zstart(3)
   
                   !-> compute dimensions to write
                   coord(1,2)=zend(1)-zstart(1)+1
                   coord(2,2)=zend(2)-zstart(2)+1
                   coord(3,2)=zend(3)-zstart(3)+1
   
                   coord(:,1) = (/phG%zst(1),phG%zst(2),phG%zst(3)/)
                   coord(:,2) = (/phG%zen(1),phG%zen(2),phG%zen(3)/)
                   coord(:,2) = coord(:,2) - coord(:,1) + 1 !count3 = count3-start3+1
   
           endif
   
   end subroutine mpi_global_coord
end module channel



