!23456789012345678901234567890123456789012345678901234567890123456789012
!   Copyright (C) 2021 IAPCM. All rights reserved.
!   
!   文件名称：tstfftw.f90
!   创 建 者：刘占军
!   创建日期：2021年12月29日
!   描    述：对三维数据中的任意两维做FFT，采用FFTW_plan_many_dft
!
!================================================================*/
  use, intrinsic :: iso_c_binding 
  implicit none
  include 'fftw3.f03'
  type(C_PTR) :: planf,planf_many
  integer,parameter::Nx=171,Ny=251,Nz=593
  !integer,parameter::Nx=13,Ny=29,Nz=53
  complex(C_DOUBLE_COMPLEX), dimension(nx,ny,nz) :: in, out,in1,out1
  integer,dimension(2)::N
  integer,dimension(2)::nembed
  integer::i,j,k,howmany
  integer::stride, dist 
   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!1: do 2d fft for the first two dimension 
  planf = fftw_plan_dft_2d(Ny,Nx, in,out, FFTW_FORWARD,FFTW_ESTIMATE)
  call init(nx,ny,nz,in)
  do i=1,nz
     call fftw_execute_dft(planf, in(:,:,i), out(:,:,i))
  enddo

  N=[Ny,Nx]; stride=1; dist=ny*nx; howmany=NZ; nembed=N
  planf_many=fftw_plan_many_dft(2 , N, howmany,  in1,   nembed,     stride,   dist,     out1,  nembed,   stride,  dist, fftw_forward, fftw_estimate);
  call init(nx,ny,nz,in1)
  call fftw_execute_dft(planf_many, in1, out1)
  if(sum(abs(out1-out))<1d-14)then
        print*,"ok1",sum(abs(out1-out))
  else
        print*,'parameters maybe incorrect',sum(abs(out1-out))
  endif
  call fftw_destroy_plan(planf_many)
  call fftw_destroy_plan(planf)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!2:  do 2d fft for the the last two dimension
  planf = fftw_plan_dft_2d(Nz,Ny, in,out, FFTW_FORWARD,FFTW_ESTIMATE)
  call init(nx,ny,nz,in)
  do i=1,nx
     call fftw_execute_dft(planf, in(i,:,:), out(i,:,:))
  enddo

  N= [NZ,Ny]; stride=nx; dist=1; nembed =[Nz, ny]; howmany=NX
  planf_many=fftw_plan_many_dft(2 , N, howmany,  in1,   nembed,     stride,   dist,     out1,  nembed,   stride,  dist, fftw_forward, fftw_estimate);

  call init(nx,ny,nz,in1)
  call fftw_execute_dft(planf_many, in1, out1)
  !if(sum(abs(out1-out))<1d-14)then
  if(maxval(abs(out1-out))<1d-14)then
        print*,"ok2",sum(abs(out1-out))
  else
        print*,'parameters maybe incorrect',sum(abs(out1-out)),maxval(abs(out1-out))
        do k=1,nz
           do j=1,ny
              do i=1,nx
                 !write(7,*)i,j,k,(out1(i,j,k)-out(i,j,k))/(abs(out1(i,j,k))+abs(out(i,j,k)))
              enddo
           enddo
        enddo
  endif
  call fftw_destroy_plan(planf_many)
  call fftw_destroy_plan(planf)



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!3:do fft for the first and last dimension.
  planf = fftw_plan_dft_2d(Nz,Nx, in,out, FFTW_FORWARD,FFTW_ESTIMATE)
  call init(nx,ny,nz,in)
  do i=1,ny
     call fftw_execute_dft(planf, in(:,i,:), out(:,i,:))
  enddo

  N= [NZ,NX]; stride=1; dist=nx;   nembed =[Nz,nx*ny] ;howmany=NY
  planf_many=fftw_plan_many_dft(2 , N, howmany,  in1,   nembed,     stride,   dist,     out1,  nembed,   stride,  dist, fftw_forward, fftw_estimate);

  call init(nx,ny,nz,in1)
  call fftw_execute_dft(planf_many, in1, out1)
  !if(sum(abs(out1-out))<1d-14)then
  if(maxval(abs(out1-out))<1d-14)then
        print*,"ok3",sum(abs(out1-out))
  else
        print*,'parameters maybe incorrect',sum(abs(out1-out)),maxval(abs(out1-out))
  endif

  call fftw_destroy_plan(planf_many)
  call fftw_destroy_plan(planf)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  end

  subroutine init(nx,ny,nz,A)
  implicit none
  integer,intent(in)::nx,ny,nz
  complex(8),dimension(nx,ny,nz)::A
  real(8),dimension(nx,ny,nz)::ReA,Ima
  integer::i,j,k
  integer::seed
  !real(8),parameter::twopi=8d0*datan(1d0)
  !do k=1,nz
  !   do j=1,ny
  !     do i=1,nx
  !       A(i,j,k)=sin(twopi*i/dble(nx))*sin(twopi*j/dble(ny)*2)*sin(twopi*k*3/dble(nz))
  !     enddo
  !   enddo
  !enddo
      call random_init(.true., .true.)
      call random_number(ReA)
      call random_number(Ima)
  a=Rea-0.5d0 +(0d0,1d0)*(Ima-0.5d0)

  end
