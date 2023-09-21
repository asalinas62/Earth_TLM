 program Earth_TLM_MoiT

!! Program for the Earth EM cavity
!! Rev. 2023/06/15
!! English Version
!! gfortran -o temporal -O2 -fopenmp Earth_TLM_MoiT.f95

 use OMP_LIB
 implicit none
 real :: dl
 integer :: ninter,nali,nsal,nlp,ncr,nfro,ncone
 real , allocatable , dimension(:,:) :: alitip
 integer , allocatable , dimension(:,:) :: conec
 integer , allocatable , dimension(:) :: front
 integer , allocatable , dimension(:) :: alinu
 integer , allocatable , dimension(:) :: salnu
 real , allocatable , dimension(:) :: zys
 real , allocatable , dimension(:) :: gaussma
 real :: gauss
 real , allocatable , dimension(:,:,:) :: vsal
 real , allocatable , dimension(:) :: vs, vsd
 real , dimension(7) :: aux
 real , dimension(3) :: campoe
 real :: ys
 real :: rescata
 integer :: i,j,ind,it
 real*8 :: wc1, wc3
 real*8 :: wct
 integer :: nth1, nth2
 integer , parameter :: ncir=6

 open(1, file='geometry1')
 open(2, file='geometry2', form='unformatted')

 read(1,*) dl
 read(1,*) ninter,nali,nsal,nlp,ncr,nfro,ncone
 allocate(alitip(3,nali))
 read(1,*) (alitip(1:3,i),i=1,nali)
 close(1)
 allocate(conec(2,ncone))
 allocate(front(nfro))
 allocate(alinu(nali))
 allocate(salnu(nsal))
 allocate(zys(ncr))
 allocate(gaussma(ninter))
 read(2) conec
 read(2) front
 read(2) alinu
 read(2) salnu
 read(2) zys
 read(2) gaussma
 close(2)

 allocate(vs(ncr*nlp))
 allocate(vsd(nlp))
 allocate(vsal(ncir,nsal,ninter))

 open(unit=3,file="voltages_box",form='unformatted')
 open(unit=4,file="time_box")

 vs=0.0
 wct=0.0
!$ nth1=OMP_get_num_threads()
!$ wc1=omp_get_wtime()


!! TEMPOROAL LOOP
 do it=1,ninter
 gauss=gaussma(it)
 print*,'nt=',it

!! INCIDENT TO REFLECTED 
!$omp parallel private(i,ind,j,vsd,ys,campoe,aux,rescata) & 
!$omp          firstprivate(nlp,ncr,nali,alinu,alitip,gauss,salnu)

!$ nth2=OMP_get_num_threads()

!$omp do
 do i=1,ncr
   ind=(i-1)*nlp
   vsd=vs(ind+1:ind+nlp)
   ys=zys(i)
   do j=1,nali
   if (alinu(j).eq.i) then
     campoe=alitip(1:3,j)*gauss
   else
     campoe=0.
   end if
   end do

!! Lines 1(-), 2(+) V1 I3:
aux(1)=ys*(vsd(1)+vsd(2)+vsd(3)+vsd(4)+campoe(1))
aux(6)=0.5*(-vsd(1)+vsd(2)+vsd(5)-vsd(6))
aux(7)=vsd(1)
vsd(1)=aux(1)+aux(6)-vsd(2)
vsd(2)=aux(1)-aux(6)-aux(7)
!! Lines 3(+), 4(-) V1 I2:
aux(5)=0.5*(vsd(3)-vsd(4)-vsd(7)+vsd(8))
aux(7)=vsd(3)
vsd(3)=aux(1)-aux(5)-vsd(4)
vsd(4)=aux(1)+aux(5)-aux(7)
!! Lines 5(+), 6(-) V2 I3
aux(2)=ys*(vsd(5)+vsd(6)+vsd(9)+vsd(10)+campoe(2))
aux(7)=vsd(5)
vsd(5)=aux(2)-aux(6)-vsd(6)
vsd(6)=aux(2)+aux(6)-aux(7)
!! Lines 7(-), 8(+) V3 I2
aux(3)=ys*(vsd(7)+vsd(8)+vsd(11)+vsd(12)+campoe(3))
aux(7)=vsd(7)
vsd(7)=aux(3)+aux(5)-vsd(8)
vsd(8)=aux(3)-aux(5)-aux(7)
!! Lines 9(-), 10(+) V2 I1
aux(4)=0.5*(-vsd(9)+vsd(10)+vsd(11)-vsd(12))
aux(7)=vsd(9)
vsd(9)=aux(2)+aux(4)-vsd(10)
vsd(10)=aux(2)-aux(4)-aux(7)
!! Lines 11(+), 12(-) V3 I1
aux(7)=vsd(11)
vsd(11)=aux(3)-aux(4)-vsd(12)
vsd(12)=aux(3)+aux(4)-aux(7)

vs(ind+1:ind+nlp)=vsd

 do j=1,nsal
 if (i.eq.salnu(j)) then
    vsal(:,j,it)=aux(1:ncir)
 end if
 end do


 end do
!$omp end do


!! REFLECTED TO INCIDENT 

!$omp do
 do i=1,nfro
 vs(front(i))=vs(front(i))*(-1)
 end do
!$omp end do

!$omp do
 do i=1,ncone
 rescata=vs(conec(1,i))
 vs(conec(1,i))=vs(conec(2,i))
 vs(conec(2,i))=rescata
 end do
!$omp end do

!$omp end parallel

end do

!$ wc3=omp_get_wtime()
 wct=wct+wc3-wc1

 write(3) vsal
 write(4,*) wct
 write(4,*) nth1,nth2
close(4)

end program Earth_TLM_MoiT
