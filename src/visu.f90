module visu

  implicit none

  private
  public :: postprocessing

contains

  subroutine postprocessing(rho1, ux1, uy1, uz1, pp3, phi1, ep1)

    use decomp_2d, only : mytype, xsize, ph1
    use case, only : postprocess_case

    use stats, only : overall_statistic
    
    use var, only : nzmsize
    use var, only : itime
    use var, only : numscalar, nrhotime, npress
    use var, only : nx, ny, nz
    use param, only : f_output

    integer,parameter :: ndim=3
    real(mytype),dimension(xsize(1),xsize(2),xsize(3)), intent(in) :: ux1, uy1, uz1
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),numscalar), intent(in) :: phi1
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),nrhotime), intent(in) :: rho1
    real(mytype),dimension(xsize(1),xsize(2),xsize(3)), intent(in) :: ep1
    real(mytype), dimension(ph1%zst(1):ph1%zen(1), ph1%zst(2):ph1%zen(2), nzmsize, npress), intent(in) :: pp3
    character(len=2) :: dimname(ndim)
    integer :: dimlen(ndim)

   !  call write_snapshot(rho1, ux1, uy1, uz1, pp3, phi1, ep1, itime)
    dimlen=(/nx,ny,nz/)
    dimname=(/'x1','y1','z1'/)
    call write_snapshot_ncdf(ux1, uy1, uz1, pp3, itime, f_output, dimname, dimlen)
   !  call postprocess_case(rho1, ux1, uy1, uz1, pp3, phi1, ep1)
   !  call overall_statistic(ux1, uy1, uz1, phi1, pp3, ep1)
    
  end subroutine postprocessing

  subroutine write_snapshot(rho1, ux1, uy1, uz1, pp3, phi1, ep1, itime)

    use decomp_2d, only : transpose_x_to_y, transpose_y_to_z, transpose_z_to_y, transpose_y_to_x
    use decomp_2d, only : mytype, xsize, ysize, zsize
    use decomp_2d, only : fine_to_coarsev
    use decomp_2d_io, only : decomp_2d_write_one

    use param, only : ivisu, ioutput, nrhotime, ilmn, iscalar, iibm

    use variables, only : derx, dery, derz 
    use variables, only : ffx, ffxp, fsx, fsxp, fwx, fwxp
    use variables, only : ffy, ffyp, fsy, fsyp, fwy, fwyp, ppy
    use variables, only : ffz, ffzp, fsz, fszp, fwz, fwzp
    use variables, only : sx, cifip6, cisip6, ciwip6, cifx6, cisx6, ciwx6
    use variables, only : sy, cifip6y, cisip6y, ciwip6y, cify6, cisy6, ciwy6
    use variables, only : sz, cifip6z, cisip6z, ciwip6z, cifz6, cisz6, ciwz6
    use variables, only : numscalar

    use var, only : one
    use var, only : uvisu
    use var, only : pp1, ta1, tb1, tc1, td1, te1, tf1, tg1, th1, ti1, di1, nxmsize
    use var, only : pp2, ta2, tb2, tc2, td2, te2, tf2, ppi2, di2, dip2, ph2, nymsize
    use var, only : ppi3, ta3, tb3, tc3, td3, te3, tf3, di3, dip3, ph3, nzmsize
    use var, only : npress

    implicit none

    character(len=30) :: filename

    !! inputs
    real(mytype), dimension(xsize(1), xsize(2), xsize(3)), intent(in) :: ux1, uy1, uz1
    real(mytype), dimension(xsize(1), xsize(2), xsize(3)), intent(in) :: ep1
    real(mytype), dimension(xsize(1), xsize(2), xsize(3), nrhotime), intent(in) :: rho1
    real(mytype), dimension(ph3%zst(1):ph3%zen(1),ph3%zst(2):ph3%zen(2),nzmsize,npress), intent(in) :: pp3
    real(mytype), dimension(xsize(1), xsize(2), xsize(3), numscalar), intent(in) :: phi1
    integer, intent(in) :: itime

    integer :: is

    if ((ivisu.ne.0).and.(mod(itime, ioutput).eq.0)) then
       !! Write velocity
       uvisu=0.
       if (iibm==2) then
          ta1(:,:,:) = (one - ep1(:,:,:)) * ux1(:,:,:)
       else
          ta1(:,:,:) = ux1(:,:,:)
       endif
       call fine_to_coarseV(1,ta1,uvisu)
990    format('ux',I3.3)
       write(filename, 990) itime/ioutput
       call decomp_2d_write_one(1,uvisu,filename,2)

       uvisu=0.
       if (iibm==2) then
          ta1(:,:,:) = (one - ep1(:,:,:)) * uy1(:,:,:)
       else
          ta1(:,:,:) = uy1(:,:,:)
       endif
       call fine_to_coarseV(1,ta1,uvisu)
991    format('uy',I3.3)
       write(filename, 991) itime/ioutput
       call decomp_2d_write_one(1,uvisu,filename,2)

       uvisu=0.
       if (iibm==2) then
          ta1(:,:,:) = (one - ep1(:,:,:)) * uz1(:,:,:)
       else
          ta1(:,:,:) = uz1(:,:,:)
       endif
       call fine_to_coarseV(1,ta1,uvisu)
992    format('uz',I3.3)
       write(filename, 992) itime/ioutput
       call decomp_2d_write_one(1,uvisu,filename,2)

       !! Write pressure
       !WORK Z-PENCILS
       call interzpv(ppi3,pp3(:,:,:,1),dip3,sz,cifip6z,cisip6z,ciwip6z,cifz6,cisz6,ciwz6,&
            (ph3%zen(1)-ph3%zst(1)+1),(ph3%zen(2)-ph3%zst(2)+1),nzmsize,zsize(3),1)
       !WORK Y-PENCILS
       call transpose_z_to_y(ppi3,pp2,ph3) !nxm nym nz
       call interypv(ppi2,pp2,dip2,sy,cifip6y,cisip6y,ciwip6y,cify6,cisy6,ciwy6,&
            (ph3%yen(1)-ph3%yst(1)+1),nymsize,ysize(2),ysize(3),1)
       !WORK X-PENCILS
       call transpose_y_to_x(ppi2,pp1,ph2) !nxm ny nz
       call interxpv(ta1,pp1,di1,sx,cifip6,cisip6,ciwip6,cifx6,cisx6,ciwx6,&
            nxmsize,xsize(1),xsize(2),xsize(3),1)

       uvisu=0._mytype
       if (iibm==2) then
          ta1(:,:,:) = (one - ep1(:,:,:)) * ta1(:,:,:)
       endif
       call fine_to_coarseV(1,ta1,uvisu)
993    format('pp',I3.3)
       write(filename, 993) itime/ioutput
       call decomp_2d_write_one(1,uvisu,filename,2)

       !! LMN - write out density
       if (ilmn) then
          uvisu=0.
          call fine_to_coarsev(1,rho1(:,:,:,1),uvisu)
995       format('rho',i3.3)
          write(filename, 995) itime/ioutput
          call decomp_2d_write_one(1,uvisu,filename,2)
       endif

       !! Scalars
       if (iscalar.ne.0) then
996       format('phi',i1.1,i3.3)
          do is = 1, numscalar
             uvisu=0.
             call fine_to_coarsev(1,phi1(:,:,:,is),uvisu)
             write(filename, 996) is, itime/ioutput
             call decomp_2d_write_one(1,uvisu,filename,2)
          enddo
       endif
    endif
  end subroutine write_snapshot

  subroutine write_snapshot_ncdf(ux1, uy1, uz1, pp3, itime, file_name, dim_name, dim_len)
   use netcdf
   USE MPI
   !use decomp_2d, only : mytype, xsize, ysize, zsize

   !use param, only : ivisu, ioutput, nrhotime, ilmn, iscalar, iibm
   USE param, disabled => itime
   USE variables
   USE decomp_2d

   USE var, only : ta1,tb1,tc1,td1,te1,tf1,tg1,th1,ti1,di1
   USE var, only : ta2,tb2,tc2,td2,te2,tf2,tg2,th2,ti2,di2
   USE var, only : ta3,tb3,tc3,td3,te3,tf3,tg3,th3,ti3,di3

   use var, only : ph3, phG, nzmsize
   use var, only : npress
   use var, only : nx, ny, nz
   use var, only : xstart, xend
   use var, only : pp1, tb1
   use var, only : xnu

   use tools, only : simu_stats
   implicit none

   character(len=*) :: file_name
   character(len=80) :: int2char
   integer,parameter           :: ndim=3
   character(len=*),intent(in) :: dim_name(ndim)
   logical :: file_exist
   integer :: varid(5),i
   integer :: ncid
   integer :: dim_len(ndim),dimid(ndim),dim_len_check
   integer :: dimt(3),coord(3,2)
   integer :: start1(3),count1(3)
   integer :: start3(3),count3(3)
  
   !real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ta1,tb1,tc1,td1,te1,tf1,tg1,th1,ti1,di1
   !real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: ta2,tb2,tc2,td2,te2,tf2,tg2,th2,ti2,tj2
   !real(mytype),dimension(zsize(1),zsize(2),zsize(3)) :: ta3,tb3,tc3,td3,te3,tf3,tg3,th3,ti3,di3
   real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: eps



   !! inputs
   real(mytype), dimension(xsize(1), xsize(2), xsize(3)), intent(in) :: ux1, uy1, uz1
   real(mytype), dimension(ph3%zst(1):ph3%zen(1),ph3%zst(2):ph3%zen(2),nzmsize,npress), intent(in) :: pp3
   integer, intent(in) :: itime

   if ((ivisu.ne.0).and.(mod(itime, ioutput).eq.0)) then
      call simu_stats(5)
      write (int2char, '(I6.6)') itime ! convert timestep to charater for output file name
      call io_check(nf90_create_par(path=trim(file_name)//trim(int2char)//'.nc',&
      !!!               cmode=IOR(NF90_NETCDF4,NF90_MPIPOSIX),ncid=ncid,&
                      cmode=IOR(NF90_NETCDF4,NF90_MPIIO),ncid=ncid,&
                      comm=MPI_COMM_WORLD,info=MPI_INFO_NULL))
      
      start1 =(/xstart(1),xstart(2),xstart(3)/)
      count1 =(/xend(1),xend(2),xend(3)/)
      count1 = count1-start1+1
      start3 =(/phG%zst(1),phG%zst(2),phG%zst(3)/)
      count3 =(/phG%zen(1),phG%zen(2),phG%zen(3)/)
      count3 = count3-start3+1

      !x-derivatives
      call derx (ta1,ux1,di1,sx,ffx,fsx,fwx,xsize(1),xsize(2),xsize(3),0)
      call derx (tb1,uy1,di1,sx,ffxp,fsxp,fwxp,xsize(1),xsize(2),xsize(3),1)
      call derx (tc1,uz1,di1,sx,ffxp,fsxp,fwxp,xsize(1),xsize(2),xsize(3),1)
      !y-derivatives
      call transpose_x_to_y(ux1,td2)
      call transpose_x_to_y(uy1,te2)
      call transpose_x_to_y(uz1,tf2)
      call dery (ta2,td2,di2,sy,ffyp,fsyp,fwyp,ppy,ysize(1),ysize(2),ysize(3),1)
      call dery (tb2,te2,di2,sy,ffy,fsy,fwy,ppy,ysize(1),ysize(2),ysize(3),0)
      call dery (tc2,tf2,di2,sy,ffyp,fsyp,fwyp,ppy,ysize(1),ysize(2),ysize(3),1)
      !!z-derivatives
      call transpose_y_to_z(td2,td3)
      call transpose_y_to_z(te2,te3)
      call transpose_y_to_z(tf2,tf3)
      call derz (ta3,td3,di3,sz,ffzp,fszp,fwzp,zsize(1),zsize(2),zsize(3),1)
      call derz (tb3,te3,di3,sz,ffzp,fszp,fwzp,zsize(1),zsize(2),zsize(3),1)
      call derz (tc3,tf3,di3,sz,ffz,fsz,fwz,zsize(1),zsize(2),zsize(3),0)
      !!all back to x-pencils
      call transpose_z_to_y(ta3,td2)
      call transpose_z_to_y(tb3,te2)
      call transpose_z_to_y(tc3,tf2)
      call transpose_y_to_x(td2,tg1)
      call transpose_y_to_x(te2,th1)
      call transpose_y_to_x(tf2,ti1)
      call transpose_y_to_x(ta2,td1)
      call transpose_y_to_x(tb2,te1)
      call transpose_y_to_x(tc2,tf1)
      !du/dx=ta1 du/dy=td1 and du/dz=tg1
      !dv/dx=tb1 dv/dy=te1 and dv/dz=th1
      !dw/dx=tc1 dw/dy=tf1 and dw/dz=ti1

      ! Compute dissipation
      eps = xnu * (ta1**2+td1**2+tg1**2+tb1**2+te1**2+th1**2+tc1**2+tf1**2+ti1**2)




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
      call io_check(nf90_def_var(ncid,'pp1',nf90_float,dimid,varid(4)))
      call io_check(nf90_def_var(ncid,'eps',nf90_float,dimid,varid(5)))
      ! call io_check(nf90_def_var(ncid,'tb1',nf90_float,dimid,varid(6)))

      call io_check(nf90_var_par_access(ncid, varid(1),nf90_collective))
      call io_check(nf90_var_par_access(ncid, varid(2),nf90_collective))
      call io_check(nf90_var_par_access(ncid, varid(3),nf90_collective))
      call io_check(nf90_var_par_access(ncid, varid(4),nf90_collective))
      call io_check(nf90_var_par_access(ncid, varid(5),nf90_collective))
      ! call io_check(nf90_var_par_access(ncid, varid(6),nf90_collective))

      !-> end of definition
      call io_check(nf90_enddef(ncid))

      !-> write field variable
      call io_check(nf90_put_var(ncid,varid(1),ux1,start=start1,count=count1))
      call io_check(nf90_put_var(ncid,varid(2),uy1,start=start1,count=count1))
      call io_check(nf90_put_var(ncid,varid(3),uz1,start=start1,count=count1))
      call io_check(nf90_put_var(ncid,varid(4),pp1,start=start1,count=count1))
      call io_check(nf90_put_var(ncid,varid(5),eps,start=start1,count=count1))
      ! call io_check(nf90_put_var(ncid,varid(6),tb1,start=start1,count=count1))

      !-> close file
      call io_check(nf90_close(ncid))
      call simu_stats(6)
   endif
 end subroutine write_snapshot_ncdf

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
endmodule visu

!######################################################################################
subroutine mean_plane_x (f1,nx,ny,nz,fm1)

  use param, only : mytype, zero

  implicit none

  integer,intent(in) :: nx, ny, nz
  real(mytype),intent(in),dimension(nx,ny,nz) :: f1
  real(mytype),intent(out),dimension(ny,nz) :: fm1
  integer :: i,j,k

  fm1 = sum(f1,DIM=1)/real(nx,mytype)
  return

end subroutine mean_plane_x
!!######################################################################################
subroutine mean_plane_y (f2,nx,ny,nz,fm2)

  use param, only : mytype, zero

  implicit none

  integer,intent(in) :: nx, ny, nz
  real(mytype),intent(in),dimension(nx,ny,nz) :: f2
  real(mytype),intent(out),dimension(nx,nz) :: fm2
  integer :: i,j,k

  fm2 = sum(f2,DIM=2)/real(ny,mytype)
  return

end subroutine mean_plane_y
!######################################################################################
subroutine mean_plane_z (f3,nx,ny,nz,fm3)

  use param, only : mytype, zero

  implicit none

  integer,intent(in) :: nx, ny, nz
  real(mytype),intent(in),dimension(nx,ny,nz) :: f3
  real(mytype),intent(out),dimension(nx,ny) :: fm3
  integer :: i,j,k

  fm3 = sum(f3,DIM=3)/real(nz,mytype)
  return

end subroutine mean_plane_z
