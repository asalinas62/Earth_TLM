 program Earth_TLM_main
 
!! Program for the Earth EM cavity
!! Rev. 2023/07/23
!! English Version
 implicit none 
 integer, parameter :: first=SELECTED_INT_KIND (10)  !! first : aux. par. integer 10^10
 real , parameter :: cluz=2.99793E8  !! clight ; pi : aux. par.
 real , parameter :: pi=3.14159265
 integer , parameter :: ncelmax=178956970    !! Number of nodes max: 2**31/nlp
 integer , parameter :: nlp=12  !! number of link lines 
 integer :: ninter  !! time intervals
 real :: ra1,ra2,dl,dt  !! Earth's radii, cell length ; time increment
 real :: rac
 integer ::ncel !!  number of cells in each dimension
 integer :: sigmatip   !! Atmospheric conductivity 
 integer :: nali,nsal  !! feeding and output points
 real , allocatable , dimension(:,:) :: alicoo, alitip, salcoo
 integer , allocatable , dimension(:) :: alinu,salnu 

 integer (first) :: nct !! nct -first- : total numerical cells
 integer(first) :: ncelb       !! ncelb ; val -first- : aux. var. 
 integer ::  ncr !! ncr : number of cells in the atmosphere
!! GEOBOX
 integer , allocatable , dimension(:,:) :: ngeve !! ngeve(6,ncr) : neighbour indices 
 integer , allocatable , dimension(:) :: ngetot !! ngetot(nct) : real indices BIG matrix
!! Auxiliary variables
 integer , dimension(3) :: celda
 integer :: i,j,k  
 integer , dimension(nlp) :: li
 real , dimension(3) :: coo !! output of the funtion coor
 real , dimension(3) :: cooa , coon
!! Geometry2 matrices
 real :: z0,y0,yt
 real , allocatable , dimension(:) :: zys
 real , allocatable , dimension(:) :: gaussval
 integer , allocatable , dimension(:,:) :: conectados
 integer , allocatable , dimension(:) :: frontera
 integer :: nfro, ncone

 !! Contains
 !!****************************************************************
 !! integer (first) function nudo(celda)   Big integer
 !! function coor(celda) result (coo)
 !! integer function nudoposi(coo)    Needs ngetot
 !! real function sigma(coo)
 !! logical function relleno(coo)
 !! subroutine car2esf(xc, ec)
 !! logical function noche(coo)
 !! logical function terremoto(coo)
 !! subroutine gauss(gaussval)      Gaussian pulse 
!! ********************************************************************

!! READING (see input_Earth_box_MAN.ipynb) 
 open(unit=1, file="input_box.txt")

 read(unit=1, fmt=*) ninter      
 read(unit=1, fmt=*) dl          
 read(unit=1, fmt=*) ra1
 read(unit=1, fmt=*) ra2
 read(unit=1, fmt=*) sigmatip   
 read(unit=1, fmt=*) nali       


 !! Matrix dimensions
 allocate(gaussval(ninter))     
 allocate(alicoo(1:3,nali))     
 allocate(alinu(nali))         
 allocate(alitip(1:3,nali))  
 read(unit=1, fmt=*) (alicoo(1:3,i), i=1,nali) 
 read(unit=1, fmt=*) (alitip(1:3,i), i=1,nali)

 read(unit=1, fmt=*) nsal       
 allocate(salcoo(1:3,nsal))     
 allocate(salnu(nsal))        
 read(unit=1, fmt=*) (salcoo(1:3,i), i=1,nsal)   

 close(unit=1)
!! END READING

 ncel=ceiling(ra2/dl)*2
 rac=ncel/2*dl
 dt=dl/(2.0*cluz)    !! Temporal increment

!! Feeding and output check
 do i=1,nali
 cooa=alicoo(1:3,i)
 call cambia(cooa,coon)
 alicoo(1:3,i)=coon 
 if (relleno(alicoo(1:3,i))) then
         print*,'Feeding ok ',i
 else
         print*,'Feeding problem ',i
         stop
 end if
 end do 
 do i=1,nsal
 cooa=salcoo(1:3,i)
 call cambia(cooa,coon)
 salcoo(1:3,i)=coon 
 if (relleno(salcoo(1:3,i))) then
         print*,'Output ok ',i
 else
         print*,'Output problem ',i
         stop
 end if
 end do
!! ********************************************************************
!! GEOBOX PROCESS
 print*,'GEOBOX'

 ncelb=ncel      !! Big integer
 nct=ncelb**3    !! Big integer
 allocate(ngetot(nct))      !! Big matrix
 print*,ncelmax

 !! Number of cells in the atmosphere; ngeve Big matrix definition
 ncr=0
 do k=1,ncel     !!(3)
 do j=1,ncel     !!(2)
 do i=1,ncel     !!(1)
    celda=(/ i,j,k /)
    coo=coor(celda)
    If (relleno(coo)) then  
        ncr=ncr+1
        ngetot(nudo(celda))=ncr   
    end if
 end do
 end do
 end do
 print*, 'Cells number in the atmosphere=', ncr
 If (ncr.ge.ncelmax) then
     print*,'Overflow in number of nodes ', ncelmax
     stop
 end if
 allocate(ngeve(6,ncr))           !! Neighbour Matrix 

!! ********************************************************************
!! EM PROPERTIES AND NEIGHBOURS 
 print*,'ELECT. DEF. AND NEIGHBOURS'
 z0=120*pi
 y0=1/z0
 allocate(zys(1:ncr))
 alitip=alitip/(2.0*y0)

 ncr=0
 do k=1,ncel    !! (3)
 print*, 'level k=',k
  do j=1,ncel   !! (2)
   do i=1,ncel  !! (1)
     celda=(/ i,j,k /)
     coo=coor(celda)
     If (relleno(coo)) then  
        ncr=ncr+1
        zys(ncr)=2.0/(4.0+sigma(coo)*dl/y0 )

        celda= (/i-1,j,k/) !! x-
        coo=coor(celda)
        If (i==1) then                   
           ngeve(1,ncr)=-1  !! Perfect conductor 
       else if (relleno(coo)) then  
           ngeve(1,ncr)=ngetot(nudo(celda)) !! Big matrix
       else                             
           ngeve(1,ncr)=-1  !! Perfect conductor
       end if

       celda= (/i+1,j,k/) !! x+
       coo=coor(celda)
       If (i==ncel) then
          ngeve(2,ncr)=-1  !! Perfect conductor
       else if (relleno(coo)) then
          ngeve(2,ncr)=ngetot(nudo(celda)) !! Big matrix)
       else 
          ngeve(2,ncr)=-1  !! Perfect conductor
       end if

       celda= (/i,j-1,k/)    !! y-
       coo=coor(celda)
       If (j==1) then
          ngeve(3,ncr)=-1  !! Perfect conductor
       else if (relleno(coo)) then
          ngeve(3,ncr)=ngetot(nudo(celda)) !! Big matrix
       else 
          ngeve(3,ncr)=-1  !! Perfect conductor
       end if

       celda= (/i,j+1,k/)    !! y+
       coo=coor(celda)
       If (j==ncel) then
          ngeve(4,ncr)=-1  !! Perfect conductor
       else if (relleno(coo)) then
          ngeve(4,ncr)=ngetot(nudo(celda)) !! Big matrix
       else 
          ngeve(4,ncr)=-1  !! Perfect conductor
       end if

       celda= (/i,j,k-1/)    !! z-
       coo=coor(celda)
       If (k==1) then
          ngeve(5,ncr)=-1  !! Perfect conductor
       else if (relleno(coo)) then
          ngeve(5,ncr)=ngetot(nudo(celda)) !! Big matrix
       else 
          ngeve(5,ncr)=-1  !! Perfect conductor
       end if

       celda= (/i,j,k+1/)    !! z+
       coo=coor(celda)
       If (k==ncel) then
          ngeve(6,ncr)=-1  !! Perfect conductor
       else if (relleno(coo)) then
          ngeve(6,ncr)=ngetot(nudo(celda)) !! Big matrix
       else 
          ngeve(6,ncr)=-1  !! Perfect conductor
       end if
     end if

    end do     !! (1)
   end do      !! (2)
  end do       !! (3) 

!! ********************************************************************
!! FEEDING AND OUTPUT NUMBERING  
 print*,'FEEDING AND OUPUT'
 do i=1,nali
 alinu(i)=nudoposi(alicoo(1:3,i))     
 end do

 do i=1,nsal
 salnu(i)=nudoposi(salcoo(1:3,i))     
 end do
!! ********************************************************************
!! Memory release:
 deallocate(ngetot)
!! ********************************************************************
!! END GEOMET

!! CONNECTION PROCESS 
 print*,'CONNECTING' 

 nfro=0              !! Boundary line number 
 ncone=0             !! Connected line number

 do i=1,ncr                        
 if(ngeve(1,i).lt.0) then          !!x- 
  nfro=nfro+2
 end if
 if(ngeve(3,i).lt.0) then          !!y-
  nfro=nfro+2
 end if
 if(ngeve(5,i).lt.0) then          !!z-
  nfro=nfro+2
 end if
 if(ngeve(2,i).lt.0) then          !!x+
   nfro=nfro+2
 else
   ncone=ncone+2
 end if
 if(ngeve(4,i).lt.0) then          !!y+
   nfro=nfro+2
 else
   ncone=ncone+2
 end if
 if(ngeve(6,i).lt.0) then          !!z+
   nfro=nfro+2
 else
   ncone=ncone+2
 end if
 end do

 allocate(conectados(1:2,1:ncone))
 allocate(frontera(1:nfro))


!! Borders and conexion 
 nfro=0
 ncone=0

 li= (/ 5,7,11,1,3,9,6,8,12,2,4,10 /)   !! x-,x-,x+,x+,y-,y-,y+,y+,z-,z-,z+,z+
 do i=1,ncr
 j=(i-1)*nlp                   !! Line index 
 if(ngeve(1,i).lt.0) then      !! x-
  frontera(nfro+1)=j+li(1)     
  frontera(nfro+2)=j+li(2)    
  nfro=nfro+2
 end if
 if(ngeve(3,i).lt.0) then     !! y-
  frontera(nfro+1)=j+li(3)
  frontera(nfro+2)=j+li(4)
  nfro=nfro+2
 end if
 if(ngeve(5,i).lt.0) then       !! z-
  frontera(nfro+1)=j+li(5)
  frontera(nfro+2)=j+li(6)
  nfro=nfro+2
 end if
 if(ngeve(2,i).lt.0) then       !! x+
  frontera(nfro+1)=j+li(7)
  frontera(nfro+2)=j+li(8) 
  nfro=nfro+2
 else
  conectados(1,ncone+1)=j+li(7)
  conectados(2,ncone+1)=(ngeve(2,i)-1)*nlp+li(1)
  conectados(1,ncone+2)=j+li(8)
  conectados(2,ncone+2)=(ngeve(2,i)-1)*nlp+li(2)
  ncone=ncone+2
 end if
 if(ngeve(4,i).lt.0) then       !! y+
  frontera(nfro+1)=j+li(9)
  frontera(nfro+2)=j+li(10)
  nfro=nfro+2
 else
  conectados(1,ncone+1)=j+li(9)
  conectados(2,ncone+1)=(ngeve(4,i)-1)*nlp+li(3)
  conectados(1,ncone+2)=j+li(10)
  conectados(2,ncone+2)=(ngeve(4,i)-1)*nlp+li(4)
  ncone=ncone+2
 end if
 if(ngeve(6,i).lt.0) then       !! z+
  frontera(nfro+1)=j+li(11)
  frontera(nfro+2)=j+li(12)
  nfro=nfro+2
 else
  conectados(1,ncone+1)=j+li(11)
  conectados(2,ncone+1)=(ngeve(6,i)-1)*nlp+li(5)
  conectados(1,ncone+2)=j+li(12)
  conectados(2,ncone+2)=(ngeve(6,i)-1)*nlp+li(6)
  ncone=ncone+2
 end if
 end do     !! ncr loop closed

!! ********************************************************************
!! FEEDING
 print*,'Feeding proccess'
 call gauss(gaussval)

!! ********************************************************************
 print*, 'No. nodes = ',ncr
 print*, 'No. line pairs= ', ncone
 print*, 'No. boundary lines= ', nfro
!! ********************************************************************


!! OUTPUTS 
 print*,'OUTPUTS'
 open(unit=1, file="geometry1")
 open(unit=2, file="geometry2",form="unformatted")
 write(1,*) dl
 write(1,*) ninter, nali, nsal, nlp, ncr, nfro, ncone 
 write(1,*) (alitip(1:3,i), i=1,nali)
 close(unit=1)
 write(2) conectados
 write(2) frontera
 write(2) alinu
 write(2) salnu
 write(2) zys
 write(2) gaussval
 close(unit=2)

!! ********************************************************************
!! ********************************************************************


!! CONTAINS SUBROUTINES
!! ********************************************************************    

  contains

!! FUNCTION NUDO
 integer (first) function nudo(celda)
 integer(first) , dimension(3) :: celdab
 integer , intent(in), dimension(3) :: celda
 celdab=celda
 nudo=(celdab(3)-1)*ncelb*ncelb+(celdab(2)-1)*ncelb+celdab(1)
 end function nudo

 subroutine cambia(cooa,coon)   !! Closest node center
 real , dimension(3) , intent(in) :: cooa
 real , dimension(3) :: coon
 coon=(floor(cooa/dl)+0.5)*dl
 end subroutine cambia

 function coor(celda) result (coo)
 real , dimension(3) :: coo 
 integer , intent(in) , dimension(3) :: celda
 coo=(celda-0.5)*dl
 end function coor
 
 !! FUNCTION NUDOPOSI
 integer function nudoposi(coo) 
 real , dimension(1:3) , intent(in) :: coo
 integer , dimension(1:3) :: celdan
 celdan=Floor(coo/dl)+1   
 nudoposi=ngetot(nudo(celdan))
 end function nudoposi

 real function sigma(coo)
 implicit none
 real , dimension(3) , intent(in) :: coo
 real , dimension(3) :: coocentro , cesf
 real :: zkn
 real :: hkn,psia,psib,sigmakn
 real :: psiawhole,psibwhole
 real :: hknn,psian,psibn
 real :: hknd,psiad,psibd

 if (sigmatip.eq.1) then
         sigma=1.0E-10
 else if (sigmatip.eq.10) then
        zkn=sqrt(sum((coo-rac)**2))-ra1  
        hkn=55E3
        psib=2.9E3
        psia=8.3E3
        sigmakn=5.56E-10
        if (zkn.lt.hkn) then
         sigma=sigmakn*exp((zkn-hkn)/psia)
        else if (zkn.ge.hkn) then
         sigma=sigmakn*exp((zkn-hkn)/psib)
        else
         print*,'Error en la conductividad'
         stop
        end if
 else if (sigmatip.eq.11) then
        zkn=sqrt(sum((coo-rac)**2))-ra1  
        hknd=54E3
        hknn=60E3
        psibd=2.7E3
        psibn=3.8E3
        psiad=7.5E3
        psian=9.1E3
        sigmakn=5.56E-10
        if (noche(coo)) then
          if (zkn.lt.hknn) then
            sigma=sigmakn*exp((zkn-hknn)/psian)
          else if (zkn.ge.hknn) then
            sigma=sigmakn*exp((zkn-hknn)/psibn)
          else
            print*,'Error en la conductividad noche'
            stop
          end if
        else  
          if (zkn.lt.hknd) then
            sigma=sigmakn*exp((zkn-hknd)/psiad)
          else if (zkn.ge.hknd) then
            sigma=sigmakn*exp((zkn-hknd)/psibd)
          else
            print*,'Error en la conductividad dia'
            stop
          end if

        end if

 else if (sigmatip.eq.12) then
        zkn=sqrt(sum((coo-rac)**2))-ra1  
        hkn=55E3
        psib=2.9E3
        psia=8.3E3
        psibwhole=2.0E3
        psiawhole=13.0E3
        sigmakn=5.56E-10
        if (terremoto(coo)) then
          if (zkn.lt.hkn) then
            sigma=sigmakn*exp((zkn-hkn)/psiawhole)
          else if (zkn.ge.hkn) then
            sigma=sigmakn*exp((zkn-hkn)/psibwhole)
          else
            print*,'Error in the earthquake conductivity'
            stop
          end if   
        else  
          if (zkn.lt.hkn) then
            sigma=sigmakn*exp((zkn-hkn)/psia)
          else if (zkn.ge.hkn) then
            sigma=sigmakn*exp((zkn-hkn)/psib)
          else
            print*,'Error in the unperturbed conductivity'
            stop
          end if
        end if

 else
         sigma=0.
 end if
 end function sigma

 subroutine car2esf(xc,ec)
 real , dimension(3) , intent(in) :: xc
 real , dimension(3) :: ec
 ec(1)= sqrt(xc(1)**2+xc(2)**2+xc(3)**2)
 ec(2)= atan2(sqrt(xc(1)**2+xc(2)**2),xc(3))
 ec(3)= atan2(xc(2),xc(1))
 end subroutine car2esf

 logical function relleno(coo)
 implicit none
 real , dimension(3) , intent(in) :: coo
 real :: racoo
 racoo=sqrt(sum((coo-rac)**2))  
 relleno=(racoo<=ra2).and.(racoo>=ra1)
 end function relleno

 logical function noche(coo)
 real , dimension (3) , intent(in) :: coo
 real , dimension (3) :: coocentro  
 coocentro = coo-rac
 if (coocentro(1).le.0.) then         !!Sun on the X axis
           noche=.True.
 else if (coocentro(1).gt.0.) then
           noche=.False.
 else
           print*,'Problem in the night'
           stop
 end if
 end function noche

 logical function terremoto(coo)
 real , dimension (3) , intent(in) :: coo
 real , dimension (3) :: coocentro, cooesf
 real :: thetadeg
 coocentro = coo-rac
 call car2esf(coocentro,cooesf)
 thetadeg=cooesf(2)*180.0/pi
 if (thetadeg.le.27.) then            !!Theta <= 27 degrees 
           terremoto=.True.
 else if (thetadeg.gt.27.) then
           terremoto=.False.
 else
           print*,'Problem in the earthquake'
           stop
 end if
 end function terremoto


 subroutine gauss(gaussval)
 real , parameter :: v0=1.0
 real , dimension(ninter) :: gaussval
 real :: gg,t0
 real , parameter :: cluz=2.99793E8
 gg=cluz/(15.0*dl)   !! 2E3
 t0=Sqrt(Log(50.0))/gg   !!  1 ms
 do i=1,ninter
    gaussval(i)=v0*Exp(-(gg**2)*(i*dt-t0)**2)
 end do
 end subroutine gauss

 end program Earth_TLM_main

