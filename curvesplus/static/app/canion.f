      program Canion
c***********************************************************************
c***** Ion/water analysis           Ver 3.2 RL/KZ/MP/JHM   8/2014 ******
c***********************************************************************
      implicit integer*4 (i-n)
      implicit real*8 (a-h,o-z)
      parameter (cdr=0.017453293d0,crd=57.29577951d0,cba=0.529177249d0)
      parameter (denref=0.0006022d0,pi=3.141592654,vwat=29.9)
      parameter(
     1 n1=200,       !...max bp levels
     1 n2=1200,      !...max bins for histogram
     1 n3=50,        !...max number of ion types
     1 n4=15000,     !...max number of ions per snapshot
     1 n5=45000,     !...max points in axis frames
     1 n6=501,       !...max bins for 3D histogram
     1 n7=500)       !...max no. ions/atoms for rmsf
      character*1 mcha,ct,cs
      character*3 ext
      character*4 inam,imnam,munit,title1*20,title2*20
      character*148 lis,dat,axfrm,solute,prop,type,seq,seqc,
     1 seqin,seqco
      character*150 spr1*(n2),spr2*(n2)
      logical*2 sneg(n1),iacc(n3),rlev(n2),rmsf,circ,series,lprop,
     1 dbrac,arev,kloc,firstr,secstr
      integer*4 intdat(n4*4),idat(n4,2),ken(n7)
      real*8 ionc(n2,61,73),gdat(n4,3),hone(0:n2),his3(n6,n6,n6),
     1 cpts(800,3),cen(n7,3),vac(n7,3),gcor(3),ions
      common/axe/vol(n2,61,73),uvw(n1,4,3),hris(n1),nlev,npt
      common/bin/shell,grida,gridd,gridr,spa,kpd,kpr,nspl,
     1 idx,irx,iax
      common/bisect/bsx(3)
      common/exr/val(n1),voltot,soltot,level(n2),irej,lprop,arev
      common/cha/lis,dat,axfrm,solute,prop,type,seq,seqin
      common/dat/sris,stwi,grid,alow,ahig,dlow,dhig,rlow,rhig,
     1 pmin,pmax,itst,itnd,itdel,istep,rmsf,circ,series,dbrac
      common/ion/ilib(n3),klis(n3),kisa,nion,inam(n3)
      common/scr/uint(n5,4,3),usc(6),var(6),theta,dist
c     NOTE kpd must be an even number
      data shell,spa,kpd,kpr,nspl/30.0,0.2,6,2,200/
c-----------------------------------------------------------------------
      arev=.false.
      call cpu_time(tstart)
c-------------------------------------------------------------setup data
      ext=' '
      lis=' '
      dat=' '
      axfrm=' '
      solute=' '
      prop=' '
      type='*'
      seq='*'
      seqin=' '
      sris=3.38
      stwi=34.5
      grid=1.
      alow=0.
      ahig=360.
      dlow=1.
      dhig=500.
      rlow=0.
      rhig=30.
      pmin=-500.
      pmax=500.
      amaxmol=1.d-10
      aminmol=1.d10
      itst=0
      itnd=0
      itdel=1
      istep=0
      rmsf=.false.
      circ=.false.
      series=.false.
      lprop=.false.
      dbrac=.false.
      gridr=1.d0/kpr
      gridd=1.d0/kpd
      grida=1.d0/spa
c---------------------------------------------------------input question
      call nml
      if(solute.ne.' '.and.axfrm.eq.' ') then
      write(6,*) '---- Need matching axfrm if reading solute pdb ----'
      stop
      endif 
      if(itst.gt.0.and.itnd.eq.0) itnd=itst
      if(alow.gt.ahig) arev=.true.
      kfi=index(dat,' ')-1
         kex=index(dat,'.')+1
         ext=dat(kex:kfi)
         if(kfi-kex.ne.2) then
         write(6,*) '---- Need file type for dat ----'
         stop
         endif
      kfl=index(lis,' ')-1
      open(unit=2,file=dat(:kfi),status='old')
      if(series) open(unit=3,file=lis(:kfl)//'.cser',status='new')
      if(istep.gt.0) open(unit=8,file=lis(:kfl)//'.stp',status='new')
c-------------------------------------------------------initial cub read
      if(ext.eq.'cub'.or.ext.eq.'pts') then
         if(seqin.eq.' ') then
         write(6,*) '---- Need seqin for cub/pts I/P ----'
         stop
         endif
         if(axfrm.eq.' ') then
         write(6,*) '  ---- Need axfrm for cub/pts I/P ----'
         stop
         endif
      nlev=index(seqin,' ')-1
      nst=2
      type='*'
      nion=1
      rmsf=.false.
c-------------------------------------------------------initial cdi read
      else
      read(2,2) nlev,nst,nion
      read(2,3) (seqin(i:i),i=1,nlev)
      read(2,2) (ilib(i),i=1,nion)
      read(2,4) (inam(i),i=1,nion)
      read(2,*)
2     format(40i7)
3     format(100a1)
4     format(40a4)
      endif
c-----------------------------------------------------------------------
      if(dlow.lt.1) dlow=1.
      if(dhig.gt.nlev) then
      if(.not.circ) then
      dhig=nlev
         else
         dhig=nlev+1
         endif
      endif
      tot=0.
      kseen=0
      ksnap=0
      nonzero=0
      maxcount=0
      mincount=1000000
         do i=1,nion
         iacc(i)=.false.
         if(type.eq.'*'.or.type.eq.inam(i)) iacc(i)=.true.
         enddo
      idx=(nlev-1)*kpd
      if(circ) idx=nlev*kpd
      irx=shell*kpr
      if(rhig.lt.shell) irx=rhig*kpr
      iax=360*spa
      irt=2*irx+1
c----------------------------------------------------------sequence scan
      ilen=index(seq,' ')-1
      if(ilen.eq.1.and.seq(:1).eq.'*') then
      do i=1,idx
      level(i)=1
      enddo
      goto 27
      endif
      do i=1,idx
      level(i)=0 
      enddo
      do i=1,nlev
      cs=seqin(i:i)
      if(cs.eq.'G') then
      seqco(i:i)='C'
      else if(cs.eq.'C') then
      seqco(i:i)='G'
      else if(cs.eq.'A') then
      seqco(i:i)='T'
      else if(cs.eq.'T') then
      seqco(i:i)='A'
      else
      stop 'I/P sequnce character not recognized'
      endif
      enddo

      kup=nlev-ilen+1
      if(circ) kup=nlev
      do k=1,kup
      firstr=.false.
      secstr=.false.
      ik=0
      do i=1,ilen
      ig=k+i-1
      if(ig.gt.nlev) ig=ig-nlev
      ct=seq(i:i)
      cs=seqin(ig:ig)
      if(ct.eq.'*') then
      ik=ik+1
      else if(cs.eq.'G') then
      if(ct.eq.'G'.or.ct.eq.'R'.or.ct.eq.'S'.or.ct.eq.'K') ik=ik+1
      else if(cs.eq.'A') then
      if(ct.eq.'A'.or.ct.eq.'R'.or.ct.eq.'W'.or.ct.eq.'M') ik=ik+1
      else if(cs.eq.'C') then
      if(ct.eq.'C'.or.ct.eq.'Y'.or.ct.eq.'S'.or.ct.eq.'M') ik=ik+1
      else if(cs.eq.'T') then
      if(ct.eq.'T'.or.ct.eq.'Y'.or.ct.eq.'W'.or.ct.eq.'K') ik=ik+1
      endif
      enddo
      if(ik.eq.ilen) firstr=.true.
         if(.not.firstr) then
         ik=0
         do i=1,ilen
         ig=k+ilen-i
         if(ig.gt.nlev) ig=ig-nlev
         ct=seq(i:i)
         cs=seqco(ig:ig)
         if(ct.eq.'*') then
         ik=ik+1
         else if(cs.eq.'G') then
         if(ct.eq.'G'.or.ct.eq.'R'.or.ct.eq.'S'.or.ct.eq.'K') ik=ik+1
         else if(cs.eq.'A') then
         if(ct.eq.'A'.or.ct.eq.'R'.or.ct.eq.'W'.or.ct.eq.'M') ik=ik+1
         else if(cs.eq.'C') then
         if(ct.eq.'C'.or.ct.eq.'Y'.or.ct.eq.'S'.or.ct.eq.'M') ik=ik+1
         else if(cs.eq.'T') then
         if(ct.eq.'T'.or.ct.eq.'Y'.or.ct.eq.'W'.or.ct.eq.'K') ik=ik+1
         endif
         enddo
         if(ik.eq.ilen) secstr=.true.
         endif
      if(firstr.or.secstr) then
      imid=ilen/2
      if(imid*2.ne.ilen) then ! odd
      klow=1+(k+imid-1)*kpd-kpd/2
         else ! even
         klow=1+(k+imid-2)*kpd
         endif 
         khig=klow+kpd-1
         print *,'  kbracket= ',klow,' ',khig,' ',khig-klow+1
      do ip=klow,khig
      i=ip
      if(ip.gt.idx) i=ip-idx
      if(firstr) then
      level(i)=1
         else
         level(i)=-1
         endif
      enddo
      endif
      enddo
c--------------------------------------------------------------axis scan
27    if(dbrac) then
      do i=1,idx
      rlev(i)=.true.
      enddo
47    read(5,*,end=48) idlow,idhig
         do i=(idlow-1)*kpd+1,idhig*kpd
         rlev(i)=.false.
         enddo
         goto 47
48       do i=1,idx
         if(rlev(i)) level(i)=0
         enddo
      else
      if(dlow.gt.1) then
      kdlow=max(2,nint((dlow-1)*kpd))
      do i=1,kdlow
      level(i)=0
      enddo
      endif
      if(dhig.lt.nlev) then
      kdhig=1+nint((dhig-1)*kpd)
      do i=kdhig,idx
      level(i)=0
      enddo
      endif
      endif
         k=0
         do i=1,idx
         if(level(i).ne.0) k=k+1
         enddo
         scwt=float(k)/kpd
*--------------------------------------------------------------------O/P
      k=0
      m=0
      do i=1,nlev
      m=m+1
      spr1(m:m)=seqin(i:i)
      spr2(m:m)='|'
      if(i.eq.nlev.and..not.circ) goto 29
      do j=1,kpd
      m=m+1
      spr1(m:m)=' '
      k=k+1
      if(level(k).eq.1) then
      spr2(m:m)='^'
         else if(level(k).eq.-1) then
         spr2(m:m)='v'
             else
             spr2(m:m)='-'
             endif
      enddo
      enddo
29    idprint=m
      nup=(nlev-1)*(kpd+1)
      if(circ) nup=nlev*(kpd+1)
      nline=10*(kpd+1)
      kl=1
      ku=min(kl+nline-1,idprint)
33    if(kl.eq.1) then
      write(6,28) spr1(kl:ku),spr2(kl:ku)
28    format(/2x,'Preselect = ',a,/14x,a)
         else
         write(6,31) spr1(kl:ku),spr2(kl:ku)
31       format(/14x,a,/14x,a)
         endif
      kl=kl+nline
      ku=min(kl+nline-1,idprint)
      if(kl.le.nlev*7) goto 33
      if(scwt.eq.0) goto 500
c----------------------------------------------------------zero matrices
      do i=1,nlev*kpd+1
      do j=1,61
      do k=1,73
      ionc(i,j,k)=0.
      vol(i,j,k)=0.
      enddo
      enddo
      enddo
c----------------------------------------------------start read property
      if(prop.ne.' ') then
      lprop=.true.
      kfi=index(prop,' ')-1
      open(unit=4,file=prop(:kfi)//'.ser',status='old')
      endif
c-------------------------------------------------------read axis frames
      if(axfrm.ne.' ') then
      kfi=index(axfrm,' ')-1
      open(unit=1,file=axfrm(:kfi)//'.afr',status='old')
      read(1,*) nleva
         if(nlev.ne.nleva) then
         write(6,*) '  ---- incompatibility no. BP and axfrm I/P ----'
         stop
         endif
      do i=1,nlev
      read(1,32) ((uvw(i,j,k),k=1,3),j=1,4)
32    format(12f12.7)
      enddo
      close(1)
      call intaxe
      xmin=0.
      ymin=0.
      zmin=0.
      xmax=0.
      ymax=0.
      zmax=0.
      do i=1,npt
      if(uint(i,4,1).lt.xmin) xmin=uint(i,4,1)
      if(uint(i,4,1).gt.xmax) xmax=uint(i,4,1)
      if(uint(i,4,2).lt.ymin) ymin=uint(i,4,2)
      if(uint(i,4,2).gt.ymax) ymax=uint(i,4,2)
      if(uint(i,4,3).lt.zmin) zmin=uint(i,4,3)
      if(uint(i,4,3).gt.zmax) zmax=uint(i,4,3)
         do j=1,3
         x=uint(i,j,1)
         y=uint(i,j,2)
         z=uint(i,j,3)
         r=sqrt(x*x+y*y+z*z)
         uint(i,j,1)=x/r
         uint(i,j,2)=y/r
         uint(i,j,3)=z/r
         enddo
      enddo
      xmin=xmin-shell
      xmax=xmax+shell
      ymin=ymin-shell
      ymax=ymax+shell
      zmin=zmin-shell
      zmax=zmax+shell
c-------------------------------------------standard helical axis frames
      else
      do i=1,nlev
      do j=1,4
      do k=1,3
      uvw(1,j,k)=0.
      enddo
      enddo
      enddo
      uvw(1,1,1)=-1.
      uvw(1,2,2)=-1.
      uvw(1,3,3)= 1.
         do i=2,nlev
         ang=stwi*(i-1)*cdr
         ca=cos(ang)
         sa=sin(ang)
         do j=1,2
         x=uvw(1,j,1)
         y=uvw(1,j,2)
         uvw(i,j,1)=x*ca-y*sa
         uvw(i,j,2)=x*sa+y*ca
         enddo
         uvw(i,3,3)=uvw(1,3,3)
         uvw(i,4,3)=(i-1)*sris
         enddo
      call intaxe
      xmin=-shell
      xmax= shell
      ymin=-shell
      ymax= shell
      zmin=-shell
      zmax=(nlev-1)*sris+shell
      endif
c------------------------------------------------------------read solute
      if(solute.ne.' ') call solvol(solute)
c---------------------------------------------------------initial output
      xsiz=xmax-xmin
      ysiz=ymax-ymin
      zsiz=zmax-zmin
      idimx=int(xsiz/grid)
      idimy=int(ysiz/grid)
      idimz=int(zsiz/grid)
      ishl =int(2*kpr*shell)
         imx=max(idimx,idimy,idimz)
         if(imx.gt.n2) then
         write(6,36) imx
36       format(/2x,'Increase grid, decrease shell, or make n2 >',i4)
         stop
         endif
         do i=1,idimx
         do j=1,idimy
         do k=1,idimz
         his3(i,j,k)=0.
         enddo
         enddo
         enddo

      write(6,38) nlev,nst,nion,seqin(:nlev),
     1 idx,irx,iax,idimx,idimy,idimz,xsiz,ysiz,zsiz
38    format(/2X,'nlev= ',i3,' nst= ',i3,' nion= ',i3,
     1       /2X,'Seq: ',a,
     1      //2x,'Hdim= ',3i7,/2x,'Cdim= ',3i7,/2x,'Csiz= ',3f7.1)
c-------------------------------------------------------------rmsf setup
      if(rmsf) then
      kmin=1.d6
      kmax=0
      do ku=1,n7
      ken(ku)=0.
      do j=1,3
      cen(ku,j)=0.
      vac(ku,j)=0.
      enddo
      enddo
      endif
      flush(6)
c===============================================================cub read
      if(ext.eq.'cub') then
      ksnap=1
      read(2,5) title1
      read(2,5) title2
      read(2,16) idum,xcen,ycen,zcen
      read(2,16) ixd,vx1,vy1,vz1
      read(2,16) iyd,vx2,vy2,vz2
      read(2,16) izd,vx3,vy3,vz3
      do i=1,ixd
      do j=1,iyd
      read(2,*) (his3(i,j,k),k=1,izd)
      enddo
      enddo
      xori=xcen+(vx1+vx2+vx3)/2
      yori=ycen+(vy1+vy2+vy3)/2
      zori=zcen+(vz1+vz2+vz3)/2
      gridin=sqrt(vx1**2+vy1**2+vz1**2)*cba
         if(abs(gridin-grid).gt.1.d-2) then
         grid=gridin
         write(6,41) grid
41       format(/2x,'cub input reset grid to ',f6.2)
         endif
      ifac=4
      grids=grid/ifac
      do isf=0,ifac-1
      xoff=xori+isf*grids
      do jsf=0,ifac-1
      yoff=yori+jsf*grids
      do ksf=0,ifac-1
      zoff=zori+ksf*grids

      do ic=1,ixd
      icm=ic-1
      do jc=1,iyd
      jcm=jc-1
      do kc=1,izd
      kcm=kc-1
      x0=(xoff+icm*vx1+jcm*vx2+kcm*vx3)*cba
      y0=(yoff+icm*vy1+jcm*vy2+kcm*vy3)*cba
      z0=(zoff+icm*vz1+jcm*vz2+kcm*vz3)*cba
      den=his3(ic,jc,kc)*denref*grids**3 ! convert molarity to counts
      do i=1,nlev
      sneg(i)=.false.
      dx=x0-uvw(i,4,1)
      dy=y0-uvw(i,4,2)
      dz=z0-uvw(i,4,3)
      dot=dx*uvw(i,3,1)+dy*uvw(i,3,2)+dz*uvw(i,3,3)
      if(dot.lt.0) sneg(i)=.true.
      enddo
         iminim=0
         rglo=1.d4
         iglo=-1
         do il=1,nlev-1
         if(sneg(il).neqv.sneg(il+1)) then
         imin=-1
         ilow=1+(il-1)*nspl
         ihig=il*nspl+1
         iminim=iminim+1
         bsx(1)=x0
         bsx(2)=y0
         bsx(3)=z0
         call bisection(ilow,ihig,imin)
         dx=bsx(1)-uint(imin,4,1)
         dy=bsx(2)-uint(imin,4,2)
         dz=bsx(3)-uint(imin,4,3)
         rmin=sqrt(dx*dx+dy*dy+dz*dz)
            if(imin.ge.0.and.rmin.lt.rglo) then
            rglo=rmin
            iglo=imin
            endif
         endif
         enddo
         if(rglo.gt.shell) goto 21
         pos=1+float(iglo-1)/nspl
         rad=rglo
         vx=(x0-uint(iglo,4,1))/rglo
         vy=(y0-uint(iglo,4,2))/rglo
         vz=(z0-uint(iglo,4,3))/rglo
         px=uint(iglo,2,1)
         py=uint(iglo,2,2)
         pz=uint(iglo,2,3)
         dot=vx*px+vy*py+vz*pz
         ang=acos(dot)
         dx=py*vz-pz*vy
         dy=pz*vx-px*vz
         dz=px*vy-py*vx
         dot=dx*uint(iglo,3,1)+dy*uint(iglo,3,2)+dz*uint(iglo,3,3)
         if(dot.lt.0) ang=-ang
c-----------------------------------------helical coord. hist.
      id=1+int((pos-1)/gridd)
      if(level(id).eq.-1) ang=pi-ang
      ang3=mod(ang*crd+360.d0,360.d0)
      ir=1+int(rad/gridr)
      ia=1+int(ang3/grida)
      if(vol(id,ir,ia).le.0.) goto 21
      tot=tot+den
      ionc(id,ir,ia)=ionc(id,ir,ia)+den
21    enddo ! kc
      enddo ! jc
      enddo ! ic

      enddo ! ksf
      enddo ! jsf
      enddo ! isf
      goto 10
c===============================================================pts read
      else if(ext.eq.'pts') then
49    read(2,50,end=10) nion
50    format(i10)
      read(2,52) ((cpts(i,j),j=1,3),i=1,nion)
52    format(3f8.3)
      kseen=kseen+1
      if(kseen.lt.itst.and.itst.gt.0) goto 49
      if(kseen.gt.itnd.and.itnd.gt.0) goto 49
      if(mod(kseen-1,itdel).ne.0) goto 49
      ksnap=ksnap+1
      icount=0
      do k=1,nion
      x0=cpts(k,1)
      y0=cpts(k,2)
      z0=cpts(k,3)
      do i=1,nlev
      sneg(i)=.false.
      dx=x0-uvw(i,4,1)
      dy=y0-uvw(i,4,2)
      dz=z0-uvw(i,4,3)
      dot=dx*uvw(i,3,1)+dy*uvw(i,3,2)+dz*uvw(i,3,3)
      if(dot.lt.0) sneg(i)=.true.
      enddo
         iminim=0
         rglo=1.d4
         iglo=-1
         do il=1,nlev-1
         if(sneg(il).neqv.sneg(il+1)) then
         imin=-1
         ilow=1+(il-1)*nspl
         ihig=il*nspl+1
         iminim=iminim+1
         bsx(1)=x0
         bsx(2)=y0
         bsx(3)=z0
         call bisection(ilow,ihig,imin)
         dx=bsx(1)-uint(imin,4,1)
         dy=bsx(2)-uint(imin,4,2)
         dz=bsx(3)-uint(imin,4,3)
         rmin=sqrt(dx*dx+dy*dy+dz*dz)
            if(imin.ge.0.and.rmin.lt.rglo) then
               rglo=rmin
               iglo=imin
            endif
         endif
         enddo
         if(rglo.gt.shell) goto 51
c        if(rglo.gt.shell.or.iglo.eq.1.or.iglo.eq.npt) goto 51
         pos=1+float(iglo-1)/nspl
         rad=rglo
         vx=(x0-uint(iglo,4,1))/rglo
         vy=(y0-uint(iglo,4,2))/rglo
         vz=(z0-uint(iglo,4,3))/rglo
         px=uint(iglo,2,1)
         py=uint(iglo,2,2)
         pz=uint(iglo,2,3)
         dot=vx*px+vy*py+vz*pz
         ang=acos(dot)
         dx=py*vz-pz*vy
         dy=pz*vx-px*vz
         dz=px*vy-py*vx
         dot=dx*uint(iglo,3,1)+dy*uint(iglo,3,2)+dz*uint(iglo,3,3)
         if(dot.lt.0) ang=-ang
c--------------------------------------------------------------store ion
      id=1+int((pos-1)/gridd)
      if(level(id).eq.-1) ang=pi-ang
      ang3=mod(ang*crd+360.d0,360.d0)
      ir=1+int(rad/gridr)
      ia=1+int(ang3/grida)
      if(vol(id,ir,ia).gt.0.) then
      tot=tot+1
      ionc(id,ir,ia)=ionc(id,ir,ia)+1
      icount=icount+1
      endif
c-----------------------------------------------------------------------
      ix=1+int((x0-xmin)/grid)
      if(ix.lt.1.or.ix.gt.idimx) goto 51
      iy=1+int((y0-ymin)/grid)
      if(iy.lt.1.or.iy.gt.idimy) goto 51
      iz=1+int((z0-zmin)/grid)
      if(iz.lt.1.or.iz.gt.idimz) goto 51
      his3(ix,iy,iz)=his3(ix,iy,iz)+1
51    enddo
      if(icount.gt.0) nonzero=nonzero+1
      if(icount.lt.mincount) mincount=icount
      if(icount.gt.maxcount) maxcount=icount
      if(series) write(3,26) kseen,icount
26    format(2i10)
      if(istep.gt.0.and.mod(ksnap,istep).eq.0) write(8,37) ksnap,
     1 tot/ksnap
37    format(i12,2X,f10.3)
      goto 49 ! loop to next snapshot
c===============================================================cdi read
      else
20    read(2,22,end=10) kisa
      read(2,22) (intdat(i),i=1,kisa*4)
22    format(200i7)
      kseen=kseen+1
      if(kseen.lt.itst.and.itst.gt.0) goto 20
      if(kseen.gt.itnd.and.itnd.gt.0) goto 10
      if(mod(kseen-1,itdel).ne.0) goto 20
      ksnap=ksnap+1
            k=0
            do i=1,kisa
            do j=1,3
            k=k+1
            gdat(i,j)=float(intdat(k))/1000.
            enddo
            k=k+1
            idat(i,2)=mod(intdat(k),100)
            idat(i,1)=(intdat(k)-idat(i,2))/100
            if(rmsf) then
            ku=idat(i,1)
            if(ku.lt.kmin) kmin=ku
            if(ku.gt.kmax) kmax=ku
               if(ku.gt.n7) then
               write(6,23)
23             format(/2x,'---- Too many ions/atoms for rmsf ----')
               stop
               endif
            endif
            enddo
c---------------------------------------------------optionally read prop
      if(lprop) read(4,24) (val(i),i=1,nlev)
24    format(8x,450f8.2)
c-------------------------------------------------generate ion positions
      icount=0
      do i=1,kisa
      it=idat(i,2)
      if(.not.iacc(it).or.gdat(i,2).gt.shell) goto 25
      pos=gdat(i,1)
      rad=gdat(i,2)
      ang=gdat(i,3)
c-----------------------------------------helical coord. hist.
      id=1+int((pos-1)/gridd)
      if(level(id).eq.-1) ang=pi-ang
      ang3=mod(ang*crd+360.d0,360.d0)
      ir=1+int(rad/gridr)
      ia=1+int((ang3)/grida)
      if(vol(id,ir,ia).le.0.) goto 25
      ionc(id,ir,ia)=ionc(id,ir,ia)+1
      tot=tot+1
      icount=icount+1
c------------------------------------------------3D with frame
      imin=1+nint((pos-1)*nspl)
      rx=uint(imin,3,1)
      ry=uint(imin,3,2)
      rz=uint(imin,3,3)
      ca=cos(ang)
      sa=sin(ang)
      xx=rad*uint(imin,2,1)
      yy=rad*uint(imin,2,2)
      zz=rad*uint(imin,2,3)
      gx=(rx*rx+(1-rx*rx)*ca)*xx+(rx*ry*(1-ca)-rz*sa)*yy+
     1   (rx*rz*(1-ca)+ry*sa)*zz+uint(imin,4,1)
      gy=(rx*ry*(1-ca)+rz*sa)*xx+(ry*ry+(1-ry*ry)*ca)*yy+
     1   (ry*rz*(1-ca)-rx*sa)*zz+uint(imin,4,2)
      gz=(rx*rz*(1-ca)-ry*sa)*xx+(ry*rz*(1-ca)+rx*sa)*yy+
     1   (rz*rz+(1-rz*rz)*ca)*zz+uint(imin,4,3)
c-------------------------------------------------------------------rmsf
      if(rmsf) then
      ku=idat(i,1)
      ken(ku)=ken(ku)+1
      kion=ken(ku)
      gcor(1)=gx
      gcor(2)=gy
      gcor(3)=gz
      do j=1,3
      del=gcor(j)-cen(ku,j)
      cen(ku,j)=cen(ku,j)+del/kion
      vac(ku,j)=vac(ku,j)+(del**2)*(kion-1)/kion
      enddo
      endif
c-----------------------------------------------------------------------
      ix=1+int((gx-xmin)/grid)
      if(ix.lt.1.or.ix.gt.idimx) goto 25
      iy=1+int((gy-ymin)/grid)
      if(iy.lt.1.or.iy.gt.idimy) goto 25
      iz=1+int((gz-zmin)/grid)
      if(iz.lt.1.or.iz.gt.idimz) goto 25
      his3(ix,iy,iz)=his3(ix,iy,iz)+1
      tmol = ionc(id,ir,ia)/(vol(id,ir,ia)*denref)
      amaxmol = max(tmol,amaxmol)
      aminmol = min(tmol,aminmol)
25    enddo !....i
      if(icount.gt.0) nonzero=nonzero+1
      if(icount.lt.mincount) mincount=icount
      if(icount.gt.maxcount) maxcount=icount
      if(series) write(3,26) kseen,icount
      if(istep.gt.0.and.mod(ksnap,istep).eq.0) write(8,37) ksnap,
     1 tot/ksnap
      goto 20 ! loop to next snapshot
      endif
c=================================================================output
10    close(2)
      if(series) close(3)
      write(6,8) irej
8     format(/2x,'Inaccessible histogram elements = ',i10)
      volsol=voltot-soltot
      aion=tot/ksnap
      amol=aion/(denref*volsol)
      write(6,85) volsol,soltot
      write(6,9) ksnap,aion,amol,amaxmol/ksnap,aminmol/ksnap
85    format(/2x,' Vsolv= ',f8.1,' Vsolu= ',f8.1)
9     format(2x,' Snaps= ',i8,' <Ions>= ',f7.2,
     1 ' <Molarity>= ',f7.2,
     1 ' max(Molarity)= ',f8.1,' min(Molarity)= ',f8.1)
      call cpu_time(tend)
      write(6,30) tend-tstart
30    format(/2x,'Time for run: ',f7.1,' sec')
      if(series) write(6,11) mincount,maxcount,aion,100.*nonzero/ksnap
11    format(/2x,'Series O/P: Min ',i6,' Max ',i6,' Aver ',g12.3,
     1 'Percentage nonzero ',f5.1)
c-------------------------------------------------------------------rmsf
      if(rmsf) then
      open(unit=3,file=lis(:kfl)//'.rmsf',status='new')
      open(unit=7,file=lis(:kfl)//'_cen.pdb',status='new')
      kc=0
      do k=kmin,kmax
      if(ken(k).gt.0) then
      kc=kc+1
      sum=0.
      do j=1,3
      sum=sum+vac(k,j)
      enddo
      write(3,13) k,sqrt(sum/ken(k))
13    format(i4,1x,f8.3)
      write(7,15) kc,type(:4),k,(cen(k,j),j=1,3)
15    format('ATOM  ',i5,1x,a4,1x,'CEN R',i4,4x,3f8.3)
      endif
      enddo
      close(3)
      close(7)
      endif
c----------------------------------------------------------------------D
      do i=1,idx
      ions=0.d0
      vols=0.d0
      do j=1,irx
      do k=1,iax
      ions=ions+ionc(i,j,k)
      vols=vols+vol(i,j,k)
      enddo
      enddo
      if(vols.gt.0.) then
      hone(i)=ions/(vols*denref*ksnap)
         else
         hone(i)=0.
         endif
      enddo
      kfl=index(lis,' ')-1
         open(unit=3,file=lis(:kfl)//'.d',status='new')
         write(3,12) (1+(i-0.5)*gridd,hone(i),i=1,idx)
12       format(2f12.7)
         close(3)
c----------------------------------------------------------------------R
      do j=1,irx
      ions=0.d0
      vols=0.d0
      do i=1,idx
      do k=1,iax
      ions=ions+ionc(i,j,k)
      vols=vols+vol(i,j,k)
      enddo
      enddo
      if(vols.gt.0.) then
      hone(j)=ions/(vols*denref*ksnap)
         else
         hone(j)=0.
         endif
      enddo
         kfl=index(lis,' ')-1
         open(unit=3,file=lis(:kfl)//'.r',status='new')
         write(3,12) ((j-0.5)*gridr,hone(j),j=1,irx)
         close(3)
c--------------------------------------------------------------------RIA
         do j=0,irx
         hone(j)=0.
         enddo
      do j=1,irx
      ions=0.d0
      vols=0.d0
      do i=1,idx
      do k=1,iax
      ions=ions+ionc(i,j,k)
      vols=vols+vol(i,j,k)
      enddo
      enddo
      if(vols.gt.0.) hone(j)=hone(j-1)+(ions/ksnap)
      enddo
         kfl=index(lis,' ')-1
         open(unit=3,file=lis(:kfl)//'.ria',status='new')
         write(3,12) ((j-0.5)*gridr,hone(j),j=1,irx)
         close(3)
c----------------------------------------------------------------------A
      do k=1,iax
      ions=0.d0
      vols=0.d0
      do i=1,idx
      do j=1,irx
      ions=ions+ionc(i,j,k)
      vols=vols+vol(i,j,k)
      enddo
      enddo
      if(vols.gt.0.) then
      hone(k)=ions/(vols*denref*ksnap)
         else
         hone(k)=0.
         endif
      enddo
         kfl=index(lis,' ')-1
         open(unit=3,file=lis(:kfl)//'.a',status='new')
         write(3,12) ((k-0.5)*grida,hone(k),k=1,iax)
         close(3)
c---------------------------------------------------------------------DR
      open(unit=3,file=lis(:kfl)//'.dr',status='new')
      do i=1,idx
      do j=1,irx
      ions=0.d0
      vols=0.d0
      do k=1,iax
      ions=ions+ionc(i,j,k)
      vols=vols+vol(i,j,k)
      enddo
      if(vols.gt.0.) then
      hone(j)=ions/(vols*denref*ksnap)
         else
         hone(j)=0.
         endif
      enddo
      write(3,14) (hone(j),j=1,irx)
14    format(200f10.4)
      enddo
      close(3)
c---------------------------------------------------------------------DA
      open(unit=3,file=lis(:kfl)//'.da',status='new')
      do i=1,idx
      do k=1,iax
      ions=0.d0
      vols=0.d0
      do j=1,irx
      ions=ions+ionc(i,j,k)
      vols=vols+vol(i,j,k)
      enddo
      if(vols.gt.0.) then
      hone(k)=ions/(vols*denref*ksnap)
         else
         hone(k)=0.
         endif
      enddo
      write(3,14) (hone(k),k=1,iax)
      enddo
      close(3)
c---------------------------------------------------------------------RA
      open(unit=3,file=lis(:kfl)//'.ra',status='new')
      do j=1,irx
      do k=1,iax
      ions=0.d0
      vols=0.d0
      do i=1,idx
      ions=ions+ionc(i,j,k)
      vols=vols+vol(i,j,k)
      enddo
      if(vols.gt.0.) then
      hone(k)=ions/(vols*denref*ksnap)
         else
         hone(k)=0.
         endif
      enddo
      ave=(hone(1)+hone(iax))/2
      hone(1)=ave
      hone(iax)=ave
      write(3,14) (hone(k),k=1,iax)
      enddo
      close(3)
c--------------------------------------------------------------------DRA
c     open(unit=3,file=lis(:kfl)//'.dra',status='new')
c     do j=1,irx
c     do k=1,iax
c     do i=1,idx
c     if(vol(i,j,k).gt.0.) then
c     hone(i)=ionc(i,j,k)/(vol(i,j,k)*denref*ksnap)
c        else
c        hone(i)=0.
c        endif
c     enddo
c     write(3,14) (hone(i),i=1,idx)
c     enddo
c     enddo
c     close(3)
c--------------------------------------------------------------cartesian
      if(ext.ne.'cub') then
      facv=denref*ksnap*grid**3
      open(unit=3,file=lis(:kfl)//'.cub',status='new')
      write(3,5) dat(:kfi)
5     format(a)
      write(3,5) 'density inform'
      write(3,16) 0,(xmin+grid/2)/cba,(ymin+grid/2)/cba,
     1 (zmin+grid/2)/cba
      write(3,16) idimx,grid/cba,0.,0.
      write(3,16) idimy,0.,grid/cba,0.
      write(3,16) idimz,0.,0.,grid/cba
16    format(i4,3f10.3)
         do i=1,idimx
         do j=1,idimy
         write(3,18) (his3(i,j,k)/facv,k=1,idimz)
18       format(6e13.5)
         enddo
         enddo
      close(3)
      endif
      close(6)
      if(istep.gt.0) close(8)
500   end
c=======================================================================
      subroutine nml
      implicit integer*4 (i-n)
      implicit real*8 (a-h,o-z)
      parameter (n_real=11,n_int=4,n_log=4,n_cha=8,n_tot=27)
      character*148 lis,dat,axfrm,solute,prop,type,seq,seqin,vc(n_cha)
      character lini*150,input(n_tot)*10,lj*1
      logical*2 rmsf,circ,series,dbrac,vo(n_log),iflag(n_tot),
     1 finml,lanml,start
      integer*4 first,last,skip,vi(n_int),nmls(n_tot)
      dimension vr(n_real)
      common/cha/lis,dat,axfrm,solute,prop,type,seq,seqin
      common/dat/sris,stwi,grid,alow,ahig,dlow,dhig,rlow,rhig,
     1 pmin,pmax,itst,itnd,itdel,istep,rmsf,circ,series,dbrac
      equivalence (vc(1),lis),(vr(1),sris),(vi(1),itst),
     1 (vo(1),rmsf)
c-----------------------------------------------------------------------
      data input/'sris','stwi','grid','alow','ahig','dlow','dhig',
     & 'rlow','rhig','pmin','pmax',
     & 'itst','itnd','itdel','istep',
     & 'rmsf','circ','series','dbrac',
     & 'lis','dat','axfrm','solute','prop','type','seq','seqin'/
      ninr=n_int+n_real
      nlog=n_log+ninr
         do i=1,n_tot
         iflag(i)=.false.
         nmls(i)=index(input(i),' ')-1
         enddo
      finml=.true.
      lanml=.false.
10    read(5,5) lini
5     format(a)
      im=index(lini,'&')
      if(im.gt.0) then
      if(.not.finml) lanml=.true.
      if(finml.and.index(lini(im+1:),'&').ne.0) lanml=.true.
      endif
      do k=1,150
      if(lini(k:k).eq.'=') then
      kl=k
      start=.true.
      do j=k-1,1,-1
      lj=lini(j:j)
      if(start.and.lj.ne.' ') then
      start=.false.
      jh=j
      else if(.not.start.and.(lj.eq.' '.or.lj.eq.',')) then
      jl=j+1
      goto 15
      endif
      enddo
      goto 50
15       do i=1,n_tot
         if(lini(jl:jh).eq.input(i)) then
         iflag(i)=.true.
17          kl=kl+1
            if(lini(kl:kl).eq.' ') goto 17
            do j=kl,100
            lj=lini(j:j)
            if(lj.eq.' '.or.lj.eq.','.or.lj.eq.'&') then
            kh=j-1
            goto 19
            endif
            enddo
19       if(i.le.n_real) then
         read(lini(kl:kh),*,err=50) vr(i)
         goto 25
         else if(i.le.ninr) then
         read(lini(kl:kh),*,err=50) vi(i-n_real)
         goto 25
         else if(i.le.nlog) then
         read(lini(kl:kh),*,err=50) vo(i-ninr)
         goto 25
         else
            if(lini(kl:kl).eq.''''.and.lini(kh:kh).eq.'''') then
            kl=kl+1
            kh=kh-1
            endif
         if(kh.ge.kl) read(lini(kl:kh),5,err=50) vc(i-nlog)(1:kh-kl+1)
         goto 25
         endif
         endif
         enddo
         goto 50
      endif
25    enddo
      finml=.false.
      if(.not.lanml) goto 10
c-----------------------------------------------------------------output
      kfi=index(vc(1),' ')-1
      if(kfi.gt.0) open(unit=6,file=vc(1)(:kfi)//'.lis',status='new')
      write(6,200) 
200   format(
     1/5x,'***************************************',
     1/5x,'****  CANION   Version 3.0 4/2014  ****',
     1/5x,'***************************************'//)
      do i=1,n_tot
      nm=nmls(i)
      if(iflag(i)) then
      do j=1,nm
      ic=ichar(input(i)(j:j))
      if(ic.ge.97.and.ic.le.122) input(i)(j:j)=char(ic-32)
      enddo
      endif
      enddo
      write(6,8) (input(nlog+j),vc(j),j=1,n_cha)
8     format(2x,a5,': ',a32,:,2x,a5,': ',a32)
      write(6,*)
      write(6,12) (input(j),vr(j),j=1,n_real)
12    format(4(:,2x,a6,': ',f8.2))
      write(6,*)
      write(6,14) (input(n_real+j),vi(j),j=1,n_int)
14    format(5(:,2x,a6,': ',i8))
      write(6,*)
      write(6,16) (input(ninr+j),vo(j),j=1,n_log)
16    format(5(:,2x,a6,': ',l8))
c-----------------------------------------------------------------------
      write(6,*)
      return
50    write(6,55) lini(jl:jh)
55    format(/2x,'---- Error in namelist input for ',a,' ----'/)
      stop
      end
c=======================================================================
      function dotdelta(m)
      implicit integer*4 (i-n)
      implicit real*8 (a-h,o-z)
      parameter(n1=200,n5=45000)
      common/scr/uint(n5,4,3),usc(6),var(6),theta,dist
      common/bisect/bsx(3)
      dx=bsx(1)-uint(m,4,1)
      dy=bsx(2)-uint(m,4,2)
      dz=bsx(3)-uint(m,4,3)
      r = sqrt(dx*dx+dy*dy+dz*dz)
      dotdelta = (dx*uint(m,3,1)+dy*uint(m,3,2)+dz*uint(m,3,3))/r
      return
      end
c=======================================================================
      subroutine bisection(a,b,r)
      integer*4 a,b,dtol,m,r,iter,iterMax
      real*8 dotdelta,dx,fa,fm
      parameter(iterMax=1000,dtol=1)
      fa=dotdelta(a)
      fm=dotdelta(b)
c begin bisection
      if(fa*fm.gt.0) return
      if(fa.lt.0) then
         r=a
         dx=b-a
      else
         r=b
         dx=a-b
      endif
c bisections loop
      do iter=1,iterMax
         dx=0.5d0*dx
         m =nint(dx+float(r))
         fm=dotdelta(m)
         if(fm.le.0) r=m
c this check seldom saves 1 iteration
         if(fm.eq.0) return
c decide which edge of step with root
         if(abs(dx).lt.dtol) then
            m =r+sign(1.d0,dx)*dtol
            fm=dotdelta(m)
            fa=dotdelta(r)
            if(abs(fm).lt.abs(fa)) r=m
            return
         endif
      enddo
      write(6,*) 'Too many iterations in routine "bisection"'
      end
c=======================================================================
      subroutine intaxe
      implicit integer*4 (i-n)
      implicit real*8 (a-h,o-z)
      parameter (cdr=0.017453293d0,crd=57.29577951d0,cba=0.529177249d0)
      parameter(n1=200,n2=1200,n5=45000)
      dimension fra(2,4,3)
      logical*2 rmsf,circ,series,lprop,arev
      dimension r1(4,3),r2(4,3),t(3,3),dr(3),v(3)
      common/axe/vol(n2,61,73),uvw(n1,4,3),hris(n1),nlev,npt
      common/bin/shell,grida,gridd,gridr,spa,kpd,kpr,nspl,
     1 idx,irx,iax
      common/exr/val(n1),voltot,soltot,level(n2),irej,lprop,arev
      common/dat/sris,stwi,grid,alow,ahig,dlow,dhig,rlow,rhig,
     1 pmin,pmax,itst,itnd,itdel,istep,rmsf,circ,series,dbrac
      common/scr/uint(n5,4,3),usc(6),var(6),theta,dist
      n=0
      irej=0
      voltot=0.
      iup=nlev
      if(circ) iup=nlev+1
      do i=2,iup
      iu=i
      il=iu-1
      if(i.gt.nlev) iu=1
            do m=1,4
            do j=1,3
            r1(m,j)=uvw(il,m,j)
            r2(m,j)=uvw(iu,m,j)
            enddo
            enddo
            do j=1,3
            v(j)=r2(4,j)-r1(4,j)
            enddo
            call screw(r1,r2,0)
            proj=v(1)*usc(1)+v(2)*usc(2)+v(3)*usc(3)
            hris(il)=var(3) 
      mup=nspl-1
      if(i.eq.iup) mup=mup+1
c-------------------------------------------generate intermediate points
      do k=0,mup
      fth=theta*k/nspl
      fle=proj*k/nspl
      st=dsin(fth)
      ct=dcos(fth)
      cm=ct-1
      t(1,1)= ct       -cm*usc(1)*usc(1)
      t(1,2)=-st*usc(3)-cm*usc(1)*usc(2)
      t(1,3)= st*usc(2)-cm*usc(1)*usc(3)
      t(2,1)= st*usc(3)-cm*usc(2)*usc(1)
      t(2,2)= ct       -cm*usc(2)*usc(2)
      t(2,3)=-st*usc(1)-cm*usc(2)*usc(3)
      t(3,1)=-st*usc(2)-cm*usc(3)*usc(1)
      t(3,2)= st*usc(1)-cm*usc(3)*usc(2)
      t(3,3)= ct       -cm*usc(3)*usc(3)
         do j=1,3
         dr(j)=r1(4,j)-usc(j+3)
         enddo
      n=n+1
      do j=1,3
      uint(n,1,j)=r1(1,1)*t(j,1)+r1(1,2)*t(j,2)+r1(1,3)*t(j,3)
      uint(n,2,j)=r1(2,1)*t(j,1)+r1(2,2)*t(j,2)+r1(2,3)*t(j,3)
      uint(n,3,j)=r1(3,1)*t(j,1)+r1(3,2)*t(j,2)+r1(3,3)*t(j,3)
      uint(n,4,j)=  dr(1)*t(j,1)  +dr(2)*t(j,2)  +dr(3)*t(j,3)
     1            +usc(j+3)+fle*usc(j)
      enddo
      enddo ! k points
c--------------------------------------------------------tangent v and u
      i1=(il-1)*nspl+25
      i2=i1+1
      i3=i1+2
      dx=uint(i3,4,1)-uint(i1,4,1)
      dy=uint(i3,4,2)-uint(i1,4,2)
      dz=uint(i3,4,3)-uint(i1,4,3)
      r=sqrt(dx*dx+dy*dy+dz*dz)
      dx=dx/r
      dy=dy/r
      dz=dz/r
      vz=uint(i2,3,1)*dx+uint(i2,3,2)*dy+uint(i2,3,3)*dz
      dot=dx*usc(1)+dy*usc(2)+dz*usc(3)
      un=(theta/proj)*dot
      unx=usc(1)*un
      uny=usc(2)*un
      unz=usc(3)*un
      ux=  uint(i2,2,1)*unx+uint(i2,2,2)*uny+uint(i2,2,3)*unz
      uy=-(uint(i2,1,1)*unx+uint(i2,1,2)*uny+uint(i2,1,3)*unz)
c-------------------------------------------------generate micro volumes
      dx=r1(4,1)-usc(4)
      dy=r1(4,2)-usc(5)
      dz=r1(4,3)-usc(6)
      dot=dx*usc(1)+dy*usc(2)+dZ*usc(3)
      dx=dx-dot*usc(1)
      dy=dy-dot*usc(2)
      dz=dz-dot*usc(3)
      rh2=dx*dx+dy*dy+dz*dz
      dels=theta*sqrt((proj/theta)**2+rh2)/kpd
      delt=cdr/spa
      delr=1.d0/kpr
      do j=1,shell*kpr
      rad=(j-0.5)/kpr
      if(rad.lt.rlow.or.rad.ge.rhig) goto 21
      rad1=rad-delr/2
      rad2=rad1+delr
      del2r=rad2**2-rad1**2
      del3r=(rad2**3-rad1**3)/3.d0
         do k=1,360*spa
         the1=(k-1)/spa
         the2=k/spa
         the=(the1+the2)/2
         if(.not.arev) then
         if(the.lt.alow.or.the.ge.ahig) goto 31
            else
            if(the.ge.ahig.and.the.lt.alow) goto 31
            endif
         the=the*cdr
         tst=rad*(uy*cos(the)-ux*sin(the))
         if(tst.lt.vz) then
         the1=the1*cdr
         the2=the2*cdr
         uterm=ux*(cos(the2)-cos(the1))+uy*(sin(the2)-sin(the1))
         volume=dels*(delt*del2r*vz/2-del3r*uterm)
            else
            irej=irej+1
            volume=0.
            endif
         ihig=ilow+(kpd-1)
         do m=1,kpd
         id=(il-1)*kpd+m
         if(level(id).ne.0) then
         if(.not.lprop.or.(val(il).ge.pmin.and.val(il).le.pmax)) then
         vol(id,j,k)=volume
         voltot=voltot+volume
         endif
         endif
         enddo
31       enddo
21    enddo
c-----------------------------------------------------------------------
      enddo ! i levels
      npt=n
      return
      end
c=======================================================================
      subroutine screw(r1,r2,key)
      implicit integer*4 (i-n)
      implicit real*8 (a-h,o-z)
      parameter (pi=3.141592654d0,crd=57.29577951d0,n1=450,n5=45000)
      dimension rsc(3),r1(4,3),r2(4,3),urot(4,3),t(3,3),q(3,3),
     1 w(3),v(3),scw(3),a(3,3)
      common/scr/uint(n5,4,3),usc(6),var(6),theta,dist
      do j=1,3
      v(j)=r2(4,j)-r1(4,j)
      enddo
      dist=sqrt(v(1)**2+v(2)**2+v(3)**2)
      do i=1,3
      do j=1,3
      q(i,j)=0.d0
        do k=1,3
        q(i,j)=q(i,j)+r1(k,i)*r2(k,j)
        enddo
      enddo
      enddo
      ct = (q(1,1)+q(2,2)+q(3,3)-1.d0)/2.d0
      if(abs(ct).gt.1.) ct=sign(1.d0,ct)
      theta = dacos(ct)
      if (theta .lt. 1.d-6) then
           if(dist.lt.1.d-3) then
           do j=1,6
           var(j)=0.
           enddo
           return
           endif
 
        usc(1)= v(1)
        usc(2)= v(2)
        usc(3)= v(3)
 
      else if (pi-theta.lt.1.d-6) then
        do i=1,3
           do j=1,3
              a(i,j) = q(i,j)
           enddo
           a(i,i) = a(i,i) + 1.d0
        enddo
        call findaxis(a, 1.d-10, 1.d-6, 40, gamma, usc, it)

      else
        usc(1)=-q(3,2)+q(2,3)
        usc(2)=-q(1,3)+q(3,1)
        usc(3)=-q(2,1)+q(1,2)
      endif
      unorm=dsqrt(usc(1)**2+usc(2)**2+usc(3)**2)
      do j=1,3
      usc(j)=usc(j)/unorm
      scw(j)=usc(j)*theta*crd
      enddo
      hdot=usc(1)*v(1)+usc(2)*v(2)+usc(3)*v(3)
        do k=1,3
        dot=0.d0
        do j=1,3
        dot=dot+scw(j)*r1(k,j)
        enddo
        var(k+3)=dot
        enddo
c------------------------------------------------------------------point
      w(1)=(v(1)-hdot*usc(1))/2.d0
      w(2)=(v(2)-hdot*usc(2))/2.d0
      w(3)=(v(3)-hdot*usc(3))/2.d0
 
      tn=1/tan(theta/2.d0)
      usc(4)=(r2(4,1)+r1(4,1))/2.d0+tn*(-w(2)*usc(3)+w(3)*usc(2))
      usc(5)=(r2(4,2)+r1(4,2))/2.d0+tn*(-w(3)*usc(1)+w(1)*usc(3))
      usc(6)=(r2(4,3)+r1(4,3))/2.d0+tn*(-w(1)*usc(2)+w(2)*usc(1))
 
c--------------------------------------------------------------transform
      if(key.lt.0) return
      st=sin(theta/2)
      ct=cos(theta/2)
      cm=ct-1
      urot(4,1)=(r1(4,1)+r2(4,1))/2
      urot(4,2)=(r1(4,2)+r2(4,2))/2
      urot(4,3)=(r1(4,3)+r2(4,3))/2
 
        t(1,1)= ct       -cm*usc(1)*usc(1)
        t(1,2)=-st*usc(3)-cm*usc(1)*usc(2)
        t(1,3)= st*usc(2)-cm*usc(1)*usc(3)
        t(2,1)= st*usc(3)-cm*usc(2)*usc(1)
        t(2,2)= ct       -cm*usc(2)*usc(2)
        t(2,3)=-st*usc(1)-cm*usc(2)*usc(3)
        t(3,1)=-st*usc(2)-cm*usc(3)*usc(1)
        t(3,2)= st*usc(1)-cm*usc(3)*usc(2)
        t(3,3)= ct       -cm*usc(3)*usc(3)
      do i=1,3
      dot=0.d0
      do j=1,3
        urot(i,j)=0.d0
        do k=1,3
        urot(i,j)=urot(i,j)+r1(i,k)*t(j,k)
        enddo
      dot=dot+v(j)*urot(i,j)
      enddo
      var(i)=dot
      enddo
      return
      end
c=======================================================================
      subroutine findaxis(A, eps, dta, m, gamma, x1, it)
      integer i, it, m, n
      real*8 gamma
      real*8 A(3,3), X(3)
      real*8 eps, dta, X1(3)
      real*8 phi, s, X0(3)
      n = 3
      do i=1,n
      X0(i)=1.d0/dsqrt(dfloat(i))
      enddo
      it=-1
      l=1
      do while (it==-1.and.l<=m)
      gamma=0.d0
         do i=1,n
         X1(i)=0.d0
         do j=1,n
      X1(i)=X1(i)+A(i,j)*X0(j)
         enddo
         if(dabs(X1(i)).gt.dabs(gamma)) gamma=X1(i)
         enddo
      if(dabs(gamma).lt. eps) then
      it=0
      else
      do i=1,n
      X1(i)=X1(i)/gamma
      enddo
      phi=0.d0
         do i=1,n
	 s=dabs(X1(i)-X0(i))
	 if(s.gt.phi) phi=s
        enddo
      if(phi<dta) then
      it=1
      else
         do i=1,n
         X0(i)=X1(i)
         enddo
         l=l+1
         endif
      endif
      enddo
      return
      end
c=======================================================================
      subroutine solvol(solute)
      implicit integer*4 (i-n)
      implicit real*8 (a-h,o-z)
      parameter (cdr=0.017453293d0,crd=57.29577951d0)
      parameter(n1=200,n2=1200,n5=45000,n8=8000)
      logical*2 sneg(n1),lprop,arev
      character*4 mnam(n8),munit(n8),name,mcha(n8)*1
      character*148 solute
      character*80 line
      dimension corm(n8,3),vdw(5),nunit(n8),inx(n8)
      common/axe/vol(n2,61,73),uvw(n1,4,3),hris(n1),nlev,npt
      common/bin/shell,grida,gridd,gridr,spa,kpd,kpr,nspl,
     1 idx,irx,iax
      common/exr/val(n1),voltot,soltot,level(n2),irej,lprop,arev
      common/scr/uint(n5,4,3),usc(6),var(6),theta,dist
      data vdw/1.2,1.6,1.5,1.4,1.9/
      soltot=0.
c------------------------------------------------------------read solute
      kfi=index(solute,' ')-1
      open(unit=1,file=solute(:kfi)//'.pdb',status='old')
      i=0
1     read(1,5,end=114) line
5     format(a)
      if(line(:4).eq.'ATOM') then
      i=i+1
      read(line,10) mnam(i),munit(i),mcha(i),
     1 nunit(i),corm(i,1),corm(i,2),corm(i,3)
10    format(12x,a4,1x,a3,1x,a1,i4,4x,3f8.3)
         inx(i)=2
         name=mnam(i)
         if(index(name,'H').ne.0) then
         inx(i)=1
         else if(index(name,'N').ne.0) then
         inx(i)=3
         else if(index(name,'O').ne.0) then
         inx(i)=4
         else if(index(name,'P').ne.0) then
         inx(i)=5
         endif
      endif
      goto 1
114   kam=i
c---------------------------------------------------define solute volume
      do i=1,idx
      pos=1+(i-0.5)*gridd
      do j=1,irx
      rad=(j-0.5)*gridr
      do k=1,iax
      if(vol(i,j,k).gt.0) then
      ang=(k-0.5)*grida*cdr
      imin=1+nint((pos-1)*nspl)
      rx=uint(imin,3,1)
      ry=uint(imin,3,2)
      rz=uint(imin,3,3)
      ca=cos(ang)
      sa=sin(ang)
      xx=rad*uint(imin,2,1)
      yy=rad*uint(imin,2,2)
      zz=rad*uint(imin,2,3)
      vx=(rx*rx+(1-rx*rx)*ca)*xx+(rx*ry*(1-ca)-rz*sa)*yy+
     & (rx*rz*(1-ca)+ry*sa)*zz+uint(imin,4,1)
      vy=(rx*ry*(1-ca)+rz*sa)*xx+(ry*ry+(1-ry*ry)*ca)*yy+
     & (ry*rz*(1-ca)-rx*sa)*zz+uint(imin,4,2)
      vz=(rx*rz*(1-ca)-ry*sa)*xx+(ry*rz*(1-ca)+rx*sa)*yy+
     & (rz*rz+(1-rz*rz)*ca)*zz+uint(imin,4,3)
         do m=1,kam
         rv=vdw(inx(m))
         dx=abs(vx-corm(m,1))
         if(dx.gt.rv) goto 52
         dy=abs(vy-corm(m,2))
         if(dy.gt.rv) goto 52
         dz=abs(vz-corm(m,3))
         if(dz.gt.rv) goto 52
         rv2=rv*rv
         r2=dx*dx+dy*dy+dz*dz
            if(r2.lt.rv2) then
            soltot=soltot+vol(i,j,k)
            vol(i,j,k)=0.
            goto 51
            endif
52       enddo !...m
      endif
51    enddo !...k
      enddo !...j
      enddo !...i
      return
      end
