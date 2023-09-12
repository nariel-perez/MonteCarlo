program denton
  implicit none
  integer,allocatable:: seed(:)    
  
  integer:: imc,ip,i,j,ig,ntrials,nplot,nha,iq,nq
  integer:: att,acc,counter,nav,nwidom,neg,posit
  integer:: ianew,iatrial,iamin,ipnew, Np
  !integer,allocatable:: ia(:)    
  real(8),allocatable:: g(:),x(:),y(:),z(:),a(:),ha(:),Sq(:,:),alfas(:)
  real(8):: rn,rn4(4)
  real(8):: xnew,ynew,znew,maxd,alfa,alfanew
  real(8):: u,uold,unew,uij,vij,deltau
  real(8):: ratio,wtest,uav,u2av
  real(8):: xtrial,ytrial,ztrial,atrial,utrial,qtrial
  real(8):: v,vb,nid,r,rij,dr,rcut,dq
  real(8):: p,bmu,vol,box
  real(8):: a0,sig,anew,da,ava,a2av,qnew,qaux,maxda
  real(8):: df,bf,bfold,bfnew,duf,bfav,bf2av,bf0
  
  character(8):: chx8
  character(12):: chx12,chx12b
  !!Variables declaradas por mi
  real(8):: x1, x2, y1,y2,z2, z1, a1,a2,dx,dy,dz
  real(8):: b12, u12, v12, u_hz, v_hz, r2
  real(8):: a_max, a_min
  integer:: eq_mc,max_mc

!!!! variables a definir
  
  Np = 10   
  max_mc = 5000
  eq_mc = 500
  alfa = 3.5
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  call random_seed(size=i)
  allocate(seed(i))
  call random_seed(get=seed)

  seed=999
  call random_seed(put=seed)
  
  allocate(x(Np),y(Np),z(Np),a(Np),alfas(Np))


  alfas = alfa
  
  call set_MC_param(Np, box,vol,rcut,dr,nplot,da,nha,maxd,maxda)
  
  call initial_coordinates(x,y,z,Np,box)



  
  do imc=1,max_mc ! loop over the # of MC steps
     
     do ip=1,Np ! MC step --> Np attempts to move randomly chosen particles
          
        call random_number(rn)
          
        i=int(rn*Np)+1
          
        att=att+1
          
        call random_number(rn4)
          
        xnew=x(i)+(rn4(1)-0.5)*maxd ! try new positions
        ynew=y(i)+(rn4(2)-0.5)*maxd
        znew=z(i)+(rn4(3)-0.5)*maxd

        alfanew=alfa + (rn4(4)-0.5)*maxda
        
        if (alfanew == 1d0) alfanew = 1.1
        if (anew < 1d0) alfanew= 1.1
        anew =a(i)*alfanew
         
          
        xnew=xnew-box*anint(xnew/box)
        ynew=ynew-box*anint(ynew/box)
        znew=znew-box*anint(znew/box)


        call df_alfa(alfa,bfold)
        bfold = bfold
        
        call df_alfa(alfanew,bfnew)
        bfnew = bfnew
        
        df=bfnew-bfold
        !write(*,*) 'antes', df
          
        uold=0d0 ! save the energy of particle i in old position
        unew=0d0

          do j=1,Np
             if (i/=j) then
                
                call pair_pot(x(i),y(i),z(i),a(i),&
                     x(j),y(j),z(j),a(j),rij,uij,vij,box,rcut)
                uold=uold+uij
                
                call pair_pot(xnew,ynew,znew,anew,&
                     x(j),y(j),z(j),a(j),rij,uij,vij,box,rcut)
                unew=unew+uij

             endif
          enddo

          deltau=unew-uold
  
          duf=df+deltau          
          !write(*,*) df, deltau, duf
           
          call random_number(rn)

          if (duf<20.) then
             if (duf<=0) then
                u=u+deltau
                bf=bf+df
                acc=acc+1
                neg = neg +1
                x(i)=xnew
                y(i)=ynew
                z(i)=znew
                a(i)=anew
                alfas(i) =  alfanew

             elseif (dexp(-duf)>rn) then
                u=u+deltau
                bf=bf+df
                acc=acc+1
                posit = posit +1
                x(i)=xnew
                y(i)=ynew
                z(i)=znew
                a(i)=anew
                alfas(i) = alfanew
             endif
          endif


       enddo



       write(100,*)imc,(bf+u)/dble(Np),bf/dble(Np),u/dble(Np)

!!!!!!!!!!!!!!!!!!!!!!!!
!!!! NO NECESARIO
!!!!!!!!!!!!!!!!!!!!!

!       if (mod(imc,1000)==0) then          
 
!          write(120,*)Np
!          write(120,*)
!          do ip=1,Np
!             call pbc(x(ip),y(ip),z(ip),box)
             
!             write(120,*)"LJ",x(ip),y(ip),z(ip)
!          enddo
!          
!          
!          write(chx12,"(1pe12.3e2)")real(u/dble(Np))
!          write(chx12b,"(1pe12.3e2)")real(bf/dble(Np))
!
!          write(*,'(4xa,i9,2(1xa,1xa10))')"step:",imc,"/ u:",&
!               trim(adjustl(chx12)),"/ f_MG:",trim(adjustl(chx12b))
!          
!       endif
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!


      !............................................................

       
       ! modificacion de la caja y cambios en alfa
       if (mod(imc,100)==0) then ! correct maxd

          ratio=dble(acc)/dble(att)
          if (ratio>0.45) then
             maxd=maxd*1.05
             maxda=maxda*1.01
          endif
          if (ratio<0.45) then
             maxd=maxd*0.95
             maxda=maxda*0.99
          endif
          if (maxd>box) maxd=box
          !if (maxda>a_max-a_min) maxda=a_max-a_min
       endif
       

       
       !............................................................

!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!
       
       if (imc>=eq_mc) then

!!!!!!!!!!!!
          if (mod(imc,50)==0) then     ! take values every 50 MC steps
             uav=uav+u
             u2av=u2av+u*u
 
             bfav=bfav+bf
             bf2av=bf2av+bf*bf             
 
             ava=ava+sum(a)
             a2av=a2av+sum(a*a)
             
             nav=nav+1             
             
          endif
!!!!!!!!!!!!!!         
         
          
          ! Virial, Pressure & g(r) ...................................
         
          if (mod(imc,50)==0) then
             counter=counter+1
             
             do i=1,Np-1
                do j=i+1,Np
                   
                   call pair_pot(x(i),y(i),z(i),a(i),x(j),y(j),z(j),&
                        a(j),rij,uij,vij,box,rcut)
                   
                   v=v+vij

                   if (rij<box/2.) then              
                      ig=int(rij/dr)           
                      g(ig)=g(ig)+2 ! g(r)

                      Sq(2:nq,2)=Sq(2:nq,2)+sin(Sq(2:nq,1)*rij)/Sq(2:nq,1)/rij
                      Sq(1,2)=Sq(1,2)+1d0
                      
                   endif
                   
                enddo
             enddo


             do i=1,Np
                j=int((a(i)-a_min)/da)+1
                ha(j)=ha(j)+1d0
             enddo
             
          endif



          !............................................................
          
       endif

    enddo

    write(*,*) neg, posit


    close(100)
    close(120)



!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!
!!! save final configuration
    
!    open(unit=121,file="final_config.xyz")    
!    write(121,*)Np
!    write(121,*)
!    do ip=1,Np
!       call pbc(x(ip),y(ip),z(ip),box)
!             
!       write(121,*)"LJ",x(ip),y(ip),z(ip),a(ip)
!    enddo
!    close(121)
!
!
!
!   
!!!................................................
!!!................................................



    

    open(unit=111,file="size.dat")
    write(111,'(a,6xa,15xa)')"#","MG radius","probability"
    
    do i=1,nha
       write(111,*)((i-1)*da+a_min)*sig,ha(i)/dble(Np)/dble(counter)/da/sig
       !write(111,*)((i-1)*da+a_min)*sig,ha(i)/dble(counter)
    enddo
    
    close(111)




    open(unit=110,file="rdf.dat")

    do i=1,nplot
       r=0
       r=dr*(i+0.5)
       vb=((i+1)**3.-i**3.)*dr**3.
       !nid=(4./3.)*pi*vb*dens
       g(i)=g(i)/(dble(counter*Np)*nid)

       !write(110,*)r*sig,g(i)
       write(110,*)r,g(i)
    enddo

    close(110)




   open(unit=112,file="Sq.dat")

    do i=1,nq
       write(112,*)Sq(i,1),1d0+2d0*Sq(i,2)/dble(counter*Np)
    enddo

    close(112)









!!!................................................
!!!................................................
!!!................................................

!!!!!!!!!! check from here !!!!!!!!!!
!!! temp?
!!! cv=.... temp?
!!! bmu=.... temp?

!!! remember that everything (u,bf) is in 4/3*pi*R0^3 kBT units

!!!................................................
!!!................................................
!!!................................................



    wtest=wtest/dble(nwidom)
    !bmu=-temp*dlog(wtest)
    bmu=-dlog(wtest) ! excess 


    v=v/dble(counter)
    !p=dens*temp+v/(3d0*vol)


    uav=uav/dble(nav)
    u2av=u2av/dble(nav)


    !Cv=(u2av-uav*uav)/temp**2/dble(Np) ! residual, in kB units


    uav=uav/dble(Np)
    u2av=u2av/dble(Np*Np)

    bfav=bfav/dble(nav*Np)
    bf2av=bf2av/dble(nav*Np*Np)


    ava=ava/dble(nav*Np)*sig
    a2av=a2av/dble(nav*Np)*sig*sig ! notice Np not Np**2 is correct here!
    

    
    



!!!................................................
!!!................................................
!!!................................................



    write(chx12,"(1pe12.3e2)")maxd
    write(chx12b,"(1pe12.3e2)")box
    write(*,'(/4xa,1xa,1x"("a,")")')"max. displacement* (box):",&
         trim(adjustl(chx12)),trim(adjustl(chx12b))

    write(chx12,"(1pe12.3e1)")maxda*sig
    write(chx12b,"(1pe12.3e1)")(a_max-a_min)*sig
    write(*,'(4xa,1xa," nm (",a," nm)")')"max. swelling:",&
         trim(adjustl(chx12)),trim(adjustl(chx12b))


    write(*,'(/a/)')"--- Average Values:"

    write(chx12,"(1pe12.5e2)")uav
    write(*,'(4xa,1xa,1xa)')"u:",trim(adjustl(chx12)),"iu"

    !write(chx12,"(1pe12.5e2)")u2av
    !write(*,'(4xa,1xa)')"<U^2>:",trim(adjustl(chx12))

    write(chx12,"(1pe12.5e2)")sqrt(u2av-uav*uav)
    write(*,'(4xa,1xa,1xa/)')"Du:",trim(adjustl(chx12)),"iu"



    write(chx12,"(1pe12.5e2)")bfav
    write(*,'(4xa,1xa,1xa)')"f_MG:",trim(adjustl(chx12)),"iu"

    !write(chx12,"(1pe12.5e2)")bf2av
    !write(*,'(4xa,1xa)')"<F_MG^2>:",trim(adjustl(chx12))

    write(chx12,"(1pe12.5e2)")sqrt(bf2av-bfav*bfav)
    write(*,'(4xa,1xa,1xa/)')"Df_MG:",trim(adjustl(chx12)),"iu"


 
    write(chx12,"(1pe12.5e2)")(bfav+uav)
    write(*,'(4xa,1xa,1xa/)')"u+f_MG:",trim(adjustl(chx12)),"iu"




    write(chx12,"(f8.1)")ava
    write(*,'(4xa,1xa,1xa)')"<R>:",trim(adjustl(chx12)),"nm"
   
    write(chx12,"(f8.1)")sqrt(a2av-ava*ava)
    write(chx12b,"(f8.1)")sqrt(a2av-ava*ava)/ava*100
    write(*,'(4xa,2xa,1xa,1xa,1xa)')"DR:",trim(adjustl(chx12)),&
         "nm -",trim(adjustl(chx12b)),"%"
    
    write(chx12,"(f8.1)")2d0*ava
    write(*,'(4xa,1xa,1xa/)')"diameter:",trim(adjustl(chx12)),"nm"
    
  


  
end program denton



!!!..............................................................................
  
  subroutine Herts_Yukawa(x1,y1,z1,a1,x2,y2,z2,a2,r,u12,v12,box)
    
!!!*****************************************************************************
!!!******************************************************************************
!!! This subroutine...
    
    implicit none
    
    real(8),intent(in):: x1,y1,z1,x2,y2,z2,a1,a2,box
    real(8),intent(out):: r,u12,v12
    real(8):: dx,dy,dz,r2,b12
    real(8):: u_hz,v_hz
     

    u12=0d0
    v12=0d0


    dx=x2-x1
    dy=y2-y1
    dz=z2-z1

    dx=dx-box*anint(dx/box)
    dy=dy-box*anint(dy/box)
    dz=dz-box*anint(dz/box)

    r2=dx*dx+dy*dy+dz*dz

    r=sqrt(r2)
  

    


!!! Hertz potential ..........


    u_hz=0d0
    v_hz=0d0
    

    if (r<a1+a2) then

       b12 = 10000
       u_hz=(1d0-r/(a1+a2))**(5d0/2d0)*b12

       v_hz=5d0/2d0*r/(a1+a2)*(1d0-r/(a1+a2))**(3d0/2d0)*b12

    endif
    




    u12=u_hz

    v12=v_hz

  end subroutine Herts_Yukawa



!!!..............................................................................
  
  subroutine pair_pot(x1,y1,z1,a1,x2,y2,z2,a2,r,u12,v12,box,rcut)

!!!*****************************************************************************
!!!******************************************************************************
!!! This subroutine defines the pair potential

    implicit none

    real(8),intent(in):: x1,y1,z1,x2,y2,z2,box,rcut,a1,a2
    real(8),intent(out):: r,u12,v12


    
    call Herts_Yukawa(x1,y1,z1,a1,x2,y2,z2,a2,r,u12,v12,box)

 
  end subroutine pair_pot




  subroutine energy(x,y,z,a,Np,u,box,rcut)
     
!!!*****************************************************************************
!!!******************************************************************************
!!! This subroutine calculates the total potential energy

    implicit none
    
    integer,intent(in):: Np
    real(8),intent(in):: x(Np),y(Np),z(Np),a(Np),box,rcut
     
    real(8),intent(out):: u

    integer:: i,j
    real(8):: rij,uij,vij

    
    
    u=0d0
    
    do i=1,Np-1
       do j=i+1,Np
           
          call pair_pot(x(i),y(i),z(i),a(i),x(j),y(j),z(j),a(j),&
               rij,uij,vij,box,rcut)
           
          u=u+uij

       enddo
    enddo
     
     
     
  end subroutine energy



  subroutine df_alfa( alfa, bf)
    
    implicit none
    real(8)::chi = 0d0
    integer:: Nm, Nch
    real(8)::alfa3,st1, st2, st3
    real(8),intent(in):: alfa
    real(8),intent(out):: bf

    Nm  = 10
    Nch = 10
    alfa3 = alfa**3
    st1 = (alfa3 -1)*LOG(1-1/alfa3)
    st2 = chi*(1- 1/alfa3)
    st3 = 3/2*Nch*(alfa**2 - LOG(alfa) -1)


    bf = Nm*(st1 +st2) + st3 ! ya en unidades de KT

    
   
    
    
  end subroutine df_alfa




  subroutine set_MC_param(Np,box,vol,rcut,dr,nplot,da,nha,maxd,maxda)
    
!!!*****************************************************************************
!!!******************************************************************************
    
    implicit none
    
    integer,intent(in):: Np

    integer,intent(out):: nplot,nha
    
    real(8),intent(out):: box,vol,rcut,dr,maxd,da,maxda
    real(8):: dens


    ! distances in sig units    
    ! energies in 4/3*pi*R0^3 kBT units
    
   
   ! sig=2d0*R0

    dens= 0.5 !now volumen fraction
    
       
    vol=dble(Np)/dens
    !write(*,*) Np*4./3.*pi*(a0/sig)**3/vol

    box=vol**(1d0/3d0)
    
    maxd=0.5 ! in sig units
    !maxd=1.0 ! in sig units

    
    !a_max=a_max/sig
    !a_min=a_min/sig


    !maxda=a0/sig/10d0
    maxda =0.1


    nha=500
    !nha = 100
    !da=(a_max-a_min)/dble(nha-1)
    ! buscar equivalente con alfas
     da = 0.1
    

    nplot=100    
    dr=box/(nplot*2.)


    
    !nq=500
    !dq=2d0*pi/box/10d0
    


    rcut=0d0 !!! obsolete in this code!
    
    


    
    
    
    
  end subroutine set_MC_param




   subroutine initial_coordinates(x,y,z,Np,box)
    
!!!*****************************************************************************
!!!******************************************************************************
    
    implicit none
    
    integer,intent(in):: Np
    
    real(8),intent(in):: box
    real(8),intent(out):: x(Np),y(Np),z(Np)
    
    integer:: i,j,k,l,naux
    real(8):: aux
    
    
    
    naux=int(Np**(1d0/3d0))+1
    aux=box/naux
     
    l=0
    do i=1,naux
       do j=1,naux
          do k=1,naux
             l=l+1

             if (l>Np) cycle

             x(l)=aux/2+aux*(i-1)
             y(l)=aux/2+aux*(j-1)
             z(l)=aux/2+aux*(k-1)

          enddo
       enddo
    enddo


  end subroutine initial_coordinates
