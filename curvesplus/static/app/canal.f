      program Canal
c***********************************************************************
c***** Analysis of Curves+ .cda files   Ver 1.3  R.L./K.Z. 9/2014 ******
c                                                                      *
c      Ver 1.2 -> 1.3:                                                 *
c        - implemented circular correlation formulae                   *
c          (see doi: 10.1093/nar/gku855);                              *
c        - added NA to output, with string selection (nastr)           *
c                                                                      *
c***********************************************************************
      implicit integer*4 (i-n)
      implicit real*8 (a-h,o-z)
      parameter(
     & n1=10250, !...maximum number of variables to print (kvar)
     & n2=250,   !...maximum number of base pair levels
     & n3=1800,  !...bins for stocking fine histogram (0.1° / 0.01 ang)
     & n4=41,    !...number of variables per level in Cur+ record
     & n5=100)   !...maximum number of bins in O/P histogram
      parameter (cdr=0.017453293d0,crd=57.29577951d0)
      character*48 lis,seq,sline(100)*482,fils*48,nastr
      character line*482,watson*452,crick*452,olin*48,list(4)*12,
     & key(0:n4)*8,base(14)*1,name*9,namv(0:n1)*9,cstor(0:450)*8,
     & lwat*1,lcrk*1,sugt(10)*4,str1*4,str2*4,slev*3,plev*3,strand*1,
     & selin*48
      logical*2 corr,series,histo,start,seen,even,
     & fop(0:n4),flag(0:n1),range(0:n1,360),ansi,ansj,opened
      integer*4 intdat(n4*n2),lkey(0:n4),kang(0:n4),lenv(0:n1),lev(n2),
     & hstor(0:n1,n5),hist(0:n1,-n3:n3),hsum(0:n1),ks(0:n1),kl(0:n1),
     & kv(0:n1),kstor(0:n1),hfac(0:n4),nbin(0:n4),nlist(4),hmin,hmax
      real*8 cssum(0:n1,2),sum(0:n1),sums(0:n1),
     & sumsq(0:n1),smean(0:n1),vave(0:n1),vdev(0:n1),
     & vmin(0:n1),vmax(0:n1),vstor(0:n1),
     & csval(0:n1,2),smeancs(0:n1,2),
     & cs2sum(0:n1,2),deltacs(2)
      real*4 xcssum(0:n1,0:n1,2,2)
      real*4 bpdat(n2,n4),row(10),delta(0:n1),hstr(0:n4),
     & bin(0:n1,10),val(n2)
      common/cha/lis,seq,nastr
      common/dat/cormin,lev1,lev2,itst,itnd,itdel,corr,series,
     & histo
      data key/'tbend',
     & 'shear','stretch','stagger','buckle','propel','opening',
     & 'xdisp','ydisp','inclin','tip','ax-bend',
     & 'shift','slide','rise','tilt','roll','twist','h-ris','h-twi',
     & 'alpha','beta','gamma','delta','epsil','zeta',
     &  'chi','phase','amp',
     & 'alpha','beta','gamma','delta','epsil','zeta',
     & 'chi','phase','amp',
     & 'minw','mind','majw','majd'/
      data sugt/'C3''n','C4''x','O1''n','C1''x','C2''n',
     &          'C3''x','C4''n','O1''x','C1''n','C2''x'/
      data kang/1,0,0,0,1,1,1,0,0,1,1,1,0,0,0,1,1,1,0,1,18*1,4*0/
      data base/'A','T','G','C','R','Y','*',
     &          'T','A','C','G','Y','R','*'/
      data list/'trajectory','trajectories','level','levels'/
      data nlist/10,12,5,6/,nhist/50/
c-------------------------------------------------------------setup data
      lis=' '
      seq=' '
      nastr='NA'
      cormin=0.6
      lev1=0
      lev2=0
      itst=0
      itnd=0
      itdel=1
      corr=.false.
      series=.false.
      histo=.false.
         do i=0,n4
         j=index(key(i),' ')-1
         lkey(i)=j
         fop(i)=.false.
         enddo
         lenv(0)=lkey(0)
      ksbin=0
c---------------------------------------------------------input question
      call nml
      lseq=index(seq,' ')-1
      even=.false.
      if(lseq.gt.0.and.mod(lseq,2).eq.0) even=.true.
      if(itst.gt.0.and.itnd.eq.0) itnd=itst
      do i=index(nastr,' '),8
       nastr=' '//nastr
      enddo
      
      length=0
      klin=0
1     read(5,5,end=2) line
5     format(a)
      is=index(line,' ')+1
      ie=index(line(is:),' ')-1
      if(length.eq.0) then
      length=ie
         else
         if(ie.ne.length) then
         write(6,60)
60       format(/2X,'---- All oligomers must have same length ----'/)
         stop
         endif
         endif
      klin=klin+1
      sline(klin)=line
      goto 1
c---------------------------------------------------prepare data storage
2     itlev=0
      ittrj=0
      if(lev1.gt.0.and.lev2.eq.0) lev2=lev1
      if(lev1.eq.0) lev1=1
      if(lev2.eq.0) lev2=length
      if(lev2.lt.lev1) stop '  ---- Error: lev2<lev1 ---'
      if(lseq.eq.0) then
      ncol=lev2-lev1+1
         else
         ncol=1
         if(even) ncol=2
         endif
      kvar=ncol*n4
      do k=1,2
       deltacs(k)=0.d0
      enddo
      do i=0,kvar
      ks(i)=0
      flag(i)=.false.
      vmin(i)= 1.d6
      vmax(i)=-1.d6
      sum(i)=0.d0
      sums(i)=0.d0
      smean(i)=0.d0
      sumsq(i)=0.d0
      do k=1,2
       smeancs(i,k)=0.d0
       csval(i,k)=0.d0
       cssum(i,k)=0.d0
       cs2sum(i,k)=0.d0
      enddo
      do j=1,kvar
       do k=1,2
        do l=1,2
         xcssum(i,j,k,l)=0.d0
        enddo
       enddo
      enddo
         do j=1,10
         bin(i,j)=0.d0
         enddo
         do j=1,360
         range(i,j)=.false.
         enddo
         do j=-n3,n3
         hist(i,j)=0
         enddo
      enddo
c------------------------------------------------------choice protection
      if(klin.gt.1) then
      series=.false.
      corr=.false.
48    format(/2x,'Series & corr impossible: ',a64)
      write(6,48) 'multiple trajectories'
      endif
c--------------------------------------------------------find input data
      kln=0
      ksnap=0
      kseen=0
50    kln=kln+1
      if(kln.gt.klin) goto 10
      line=sline(kln)
      ibl=index(line,' ')-1
      olin=line(:ibl)
      watson=line(ibl+2:)
         do i=1,length
         do j=1,4
         if(watson(i:i).eq.base(j)) crick(i:i)=base(j+7)
         enddo
         enddo
      nlev=0
      do i=1,length
      if(lseq.eq.0) then !....................no sequence
      if(i.ge.lev1.and.i.le.lev2) then
      nlev=nlev+1
      lev(nlev)=i
      endif
         else !..............................sequence search
         lact=i-1+nint((lseq+0.1)/2)
         iwat=0
         icrk=0
         if(lact.ge.lev1.and.lact.le.lev2) then
         if(.not.even.or.lact+1.le.lev2) then
         do j=1,lseq
         lwat=watson(i-1+j:i-1+j)
         if((seq(j:j).eq.lwat).or.
     &      (seq(j:j).eq.'R'.and.(lwat.eq.'A'.or.lwat.eq.'G')).or.
     &      (seq(j:j).eq.'Y'.and.(lwat.eq.'T'.or.lwat.eq.'C')).or.
     &      (seq(j:j).eq.'*')) iwat=iwat+1
         lcrk=crick(i+lseq-j:i+lseq-j)
         if((seq(j:j).eq.lcrk).or.
     &      (seq(j:j).eq.'R'.and.(lcrk.eq.'A'.or.lcrk.eq.'G')).or.
     &      (seq(j:j).eq.'Y'.and.(lcrk.eq.'T'.or.lcrk.eq.'C')).or.
     &      (seq(j:j).eq.'*')) icrk=icrk+1
         enddo
         endif
         endif
         if(iwat.eq.lseq) then
         nlev=nlev+1    
         lev(nlev)=lact
            else if(icrk.eq.lseq) then
            nlev=nlev+1
            lev(nlev)=-lact
            if(even) lev(nlev)=-(lact+1) !...RL
            endif
      endif
      enddo
      if(nlev.eq.0) goto 50
      itlev=itlev+nlev
      ittrj=ittrj+1
      write(6,52) ittrj,length,nlev,olin(:10),watson(:min(50,length))
52    format(2x,i2,3x,2i5,3x,a10,2x,a)
         if(length.gt.50) then
         ilow=51
         ihig=min(length,100)
54       write(6,53) watson(ilow:ihig)
53       format(32x,a)
         ilow=ilow+50
         ihig=min(length,ihig+50)
         if(ilow.le.length) goto 54
         endif
c------------------------------------------------------choice protection
      if((series.or.corr).and.lseq.gt.0.and.nlev.gt.1) then
      series=.false.
      corr=.false.
      write(6,48) 'multiple sequence motifs'
      endif
c--------------------------------------------------------------read data
      if(index(olin,'.').eq.0) then
      olin=olin(:ibl)//'.cda'
      ibl=ibl+4
      endif
      inquire(unit=2,opened=opened)
      if(opened) close(2)
      open(unit=2,file=olin(:ibl),status='old')
c     read(2,22) nold,nsto,ntoro,jtot,klig
c     if(nold.ne.length) stop '  ---- I/P seq length incompatible ----'
20    read(2,22,end=50) (intdat(i),i=1,n4*length+1)
22    format(200i6)
      kseen=kseen+1
      if(kseen.lt.itst.and.itst.gt.0) goto 20
      if(kseen.gt.itnd.and.itnd.gt.0) goto 50
      if(mod(kseen-1,itdel).ne.0) goto 20
      ksnap=ksnap+1
            k=0
            do i=1,length
            do j=1,n4
            k=k+1
            bpdat(i,j)=float(intdat(k))/100.
            enddo
            enddo
            k=k+1
            bptot=float(intdat(k))/100.
c------------------------------------------------start loop on variables
      vstor(:)=301.d0
      ivl=1
      do i=0,kvar
      iv=mod(i,n4)
      if(iv.eq.0.and.i.gt.0) iv=n4
      if(iv.le.1.and.i.gt.n4) ivl=ivl+1
      kfn=lkey(iv)
      namv(i)=key(iv)
      if(iv.ge.20.and.iv.le.28) then
      kfn=kfn+1
      namv(i)(kfn:)='W'
         else if(iv.ge.29.and.iv.le.37) then
         kfn=kfn+1
         namv(i)(kfn:)='C'
         endif
      lenv(i)=kfn
      nlow=1
      nhig=nlev
         if(lseq.eq.0.or.i.eq.0) then
         nlow=ivl
         nhig=ivl
         endif
      do nl=nlow,nhig
      isig=1
      idel=0
      ishf=0
      il=abs(lev(nl))
      if(lseq.gt.0) then
      if(lev(nl).gt.0) then
      if(i.gt.n4) il=il+1
         else
         if(i.gt.n4) il=il-1
         if(iv.ge.12.and.iv.le.19) ishf=-1
         if(iv.eq.1.or.iv.eq.4.or.iv.eq.8.or.iv.eq.10
     &   .or.iv.eq.12.or.iv.eq.15) isig=-1
         if(iv.ge.20.and.iv.le.28) idel=9
         if(iv.ge.29.and.iv.le.37) idel=-9
         endif
      endif
      ivd=iv+idel
          if(ks(i).eq.0) then
          kl(i)=il
          kv(i)=ivd
          endif
      if(i.gt.0) then
      v=isig*bpdat(il+ishf,ivd)
         else
         v=bptot
         endif
      if(abs(v).lt.250) then
      ks(i)=ks(i)+1
      sum(i)=sum(i)+v
      if(kang(ivd).eq.0) then !.... distances
      sums(i)=sums(i)+v*v
         else !..... angles
         intv=max(1,nint(v))
         if(v.lt.0) intv=max(1,nint(360.d0+v))
         range(i,intv)=.true.
         vrad=v*cdr
         cssum(i,1)=cssum(i,1)+cos(vrad)
         cssum(i,2)=cssum(i,2)+sin(vrad)
         endif
      if(ivd.lt.20.or.ivd.gt.37)then
         if(v.lt.vmin(i)) vmin(i)=v
         if(v.gt.vmax(i)) vmax(i)=v
      endif
c-------------------------------------------------------------sugar bins
         if(ivd.eq.27.or.ivd.eq.36) then
         vpos=v
         if(v.lt.0.) vpos=360.+v
         jb=1+int(vpos/36)
         bin(i,jb)=bin(i,jb)+1
         endif
c-----------------------------------------------------generate histogram
         if(kang(ivd).eq.0) then
         ipos=nint(v*100)
            else
            ipos=nint(v*10)
            endif
         if(abs(ipos).le.n3) hist(i,ipos)=hist(i,ipos)+1
c------------------------------------------------------epsilon/zeta bins
         if(ivd.eq.24.or.ivd.eq.33) then
         vdif=v-bpdat(il,ivd+1)
         if(abs(vdif).gt.180) vdif=vdif-sign(360.d0,vdif)
         if(vdif.lt.0) then
         bin(i,1)=bin(i,1)+1
            else
            bin(i,2)=bin(i,2)+1
            endif
         endif
c-------------------------------------------------------alpha/gamma bins
         if(ivd.eq.20.or.ivd.eq.29) then
         va=v
         vg=bpdat(il,ivd+2)
         iva=0
         if(abs(va).lt.120) iva=sign(1.d0,va)
         ivg=0
         if(abs(vg).lt.120) ivg=sign(1.d0,vg)
         iva=iva+2
         ivg=ivg+2
         iag=(iva-1)*3+ivg
         bin(i,iag)=bin(i,iag)+1
         endif

         if(series.or.corr) then
         vstor(i)=v
         kstor(i)=ks(i)
         endif
c-----------------------------------------end processing of level/record
            else !....bpdat data doesn't exist
            flag(i)=.true.
            endif
         enddo !....nl
      enddo !....i
c-------------------------------------------------------start series O/P
      if(series) then
      do i=0,n4
         if(.not.fop(i)) then
         kfi=index(lis,' ')-1
         kfn=lenv(i)
         fop(i)=.true.
         open(unit=i+9,file=lis(:kfi)//'_'//namv(i)(:kfn)//'.ser',
     &   status='new')
         endif
      jup=ncol-1
      if(i.eq.0) jup=0
      do j=0,jup
      if(abs(vstor(i+j*n4)).lt.250) then
      write(cstor(j),'(f8.2)') vstor(i+j*n4)
         else
         write(cstor(j),'(a8)') nastr
         endif
      enddo
      write(i+9,32) kseen,(cstor(j),j=0,jup)
32    format(i12,450a8)
      enddo
      endif
c--------------------------------------------------data for correlations
      if(corr) then
       do i=0,kvar
        if(ks(i).ge.1.and..not.flag(i).and.
     &     (i.gt.0.or.(i.eq.0.and.nl.eq.1))) then
         ksi=kstor(i)
         v=vstor(i)
         ivd=kv(i)
         sweep=float(ksi-1)/ksi
         delta(i)=v-smean(i)

         if(kang(ivd).eq.1) then
          tmp = v*cdr
         else
          tmp = atan(v/5.d0)
         endif
         csval(i,1) = cos(tmp)
         csval(i,2) = sin(tmp)
         cs2sum(i,1) = cs2sum(i,1) + cos(2.d0*tmp)
         cs2sum(i,2) = cs2sum(i,2) + sin(2.d0*tmp)
         deltacs(1) = csval(i,1) - smeancs(i,1)
         deltacs(2) = csval(i,2) - smeancs(i,2)
         
         if(ksi.eq.1) then
          smean(i) = v
          smeancs(i,1) = csval(i,1)
          smeancs(i,2) = csval(i,2)
         else
          smean(i) = smean(i)+delta(i)/ksi 
          sumsq(i) = sumsq(i)+delta(i)*delta(i)*sweep
          smeancs(i,1) = smeancs(i,1) + deltacs(1)/ksi
          smeancs(i,2) = smeancs(i,2) + deltacs(2)/ksi
         endif
        endif

        do j=0,i-1
         if(ks(j).eq.ks(i).and..not.flag(j)) then
          jvd=kv(j)
          if(kang(ivd).eq.0.and.kang(jvd).eq.0) then
           xcssum(j,i,1,1)=xcssum(j,i,1,1)+delta(i)*delta(j)*sweep
          else
           do jcs=1,2
            do ics=1,2
             xcssum(j,i,jcs,ics) =
     1        xcssum(j,i,jcs,ics) + csval(j,jcs)*csval(i,ics)
            enddo
           enddo
          endif
         endif
        enddo
        
       enddo
      endif

      goto 20 ! loop to next snapshot
c------------------------------------------------end of processing input
10       if(series) then
         do i=1,n4
         close(unit=i+9)
         enddo
         endif
c-----------------------------------------------------------------------
      itt=1
      if(ittrj.gt.1) itt=2
      itl=3
      if(itlev.gt.1) itl=4
      write(6,9) ittrj,list(itt)(:nlist(itt)),itlev,
     & list(itl)(:nlist(itl))
9     format(/2x,'Data from ',i2,1x,a,' and ',i3,1x,a,' ...')
      if(lseq.gt.0) then
      if(.not.even) then
      write(6,18) 'M',(lev(i),i=1,nlev)
18    format(/2x,a1,') ',15i4,40(:/5x,15i4))
         else
         write(6,18) 'L',(lev(i),i=1,nlev)
         endif
      endif
      write(6,19)
19    format(/5x,'No.',2x,'Lev',6x,'Var',4x,'Aver',4x,'Sdev',4x,
     & 'Range',4x,'Min',5x,'Max',4x,'Ndat')
c------------------------------------------------------------output data
      do i=0,kvar
      iv=kv(i)
      il=kl(i)
      kfn=lenv(i)
         write(slev,'(i3)') il
         if(lseq.gt.0) then
         slev='  M'
         if(even.and.i.le.n4) slev='  L'
         if(even.and.i.gt.n4) slev='  U'
         endif
         if(i.eq.0) slev='  -'
      if(mod(i,n4).eq.1) write(6,24)
24    format(2x,77('='))
      if(ks(i).eq.0) then
      write(6,21) i,slev,namv(i)(:kfn)
      else
      if(kang(iv).eq.0) then
      vave(i)=sum(i)/ks(i)
      vdev(i)=sqrt(sums(i)/ks(i)-vave(i)**2)
      if(klin.gt.1) then
         write(6,21) i,slev,namv(i)(:kfn),vave(i),
     & vdev(i),vmax(i)-vmin(i),vmin(i),vmax(i),ks(i)
      else
         write(6,21) i,slev,namv(i)(:kfn),vave(i),
     & vdev(i),vmax(i)-vmin(i),vmin(i),vmax(i),ks(i)
      endif
21    format(2x,i5,') ',a3,1x,a9,2f8.2,3f8.1,3i12)
         else
         r=sqrt(cssum(i,1)**2+cssum(i,2)**2)
         vave(i)=acos(cssum(i,1)/r)*crd
         if(cssum(i,2).lt.0.) vave(i)=-vave(i)

c        r2=(vcos(i)/ks(i))**2+(vsin(i)/ks(i))**2
c        eps=sqrt(1-r2)
c        vdev2=crd*(1-0.1547*eps**3)*asin(eps)

         vdel=0
         do j=1,360
         if(range(i,j)) vdel=vdel+1
         enddo
            suma=0.
            do j=-n3,n3
            iw=hist(i,j)
            if(iw.gt.0) then
            v=j/10.d0
            del=v-vave(i)
            if(abs(del).gt.180) del=del-sign(360.d0,del)
            suma=suma+iw*del**2
            endif
            enddo
            vdev(i)=sqrt(suma/ks(i))
            if(iv.lt.20) then
      if(klin.gt.1) then
         write(6,26) i,slev,namv(i)(:kfn),vave(i),
     & vdev(i),vdel,vmin(i),vmax(i),ks(i)
      else
         write(6,26) i,slev,namv(i)(:kfn),vave(i),
     & vdev(i),vdel,vmin(i),vmax(i),ks(i)
      endif
26          format(2x,i5,') ',a3,1x,a9,2(f7.1,1x),3f8.1,3i12)
            else
            write(6,27) i,slev,namv(i)(:kfn),vave(i),vdev(i),vdel,
     &      ks(i)
            endif
27       format(2x,i5,') ',a3,1x,a9,2(f7.1,1x),f8.1,16x,i12)
         endif
      if(iv.eq.6.or.iv.eq.11.or.iv.eq.19.or.iv.eq.28.or.iv.eq.37)
     & write(6,23)
23    format(4x,75('-'))
      endif
      enddo

c---------------------------------------------------------O/P sugar bins
         do j=1,10
         row(j)=0.
         enddo
         ksum=0
      ic=0
      do i=1,kvar
      iv=kv(i)
      if(iv.eq.27.or.iv.eq.36) then
      ic=ic+1
      il=kl(i)
      kfn=lenv(i)
         write(slev,'(i3)') il
         if(lseq.gt.0) then
         slev='  M'
         if(even.and.i.le.n4) slev='  L'
         if(even.and.i.gt.n4) slev='  U'
         endif
            ksum=ksum+ks(i)
            do j=1,10
            row(j)=row(j)+bin(i,j)
            bin(i,j)=100*bin(i,j)/ks(i)
            enddo
      if(ic.eq.1) write(6,33) (sugt(j),j=1,10)
33    format(/2x,'Sugar pucker distribution ...',//18x,10a6,
     & /20x,59('-'))
      write(6,345) i,slev,namv(i)(6:kfn),(nint(bin(i,j)),j=1,10)
345   format(2x,i5,') ',a3,1x,a2,1x,10i6)
      endif
      enddo
      if(ic.gt.0) then
      write(6,340) (nint(100*row(j)/ksum),j=1,10)
340   format(20x,59('-'),/11x,'  Tot',10i6)
      endif
c--------------------------------------------------O/P epsilon/zeta bins
         do j=1,2
         row(j)=0.
         enddo
         ksum=0
      ic=0
      do i=1,kvar
      iv=kv(i)
      if(iv.eq.24.or.iv.eq.33) then
      if(ks(i).gt.0.and.ks(i+1).gt.0) then
      ic=ic+1
      il=kl(i)
         write(slev,'(i3)') il
         if(lseq.gt.0) then
         slev='  M'
         if(even.and.i.le.n4) slev='  L'
         if(even.and.i.gt.n4) slev='  U'
         endif
            ksum=ksum+ks(i)
            do j=1,2
            row(j)=row(j)+bin(i,j)
            bin(i,j)=100*bin(i,j)/ks(i)
            enddo
      if(ic.eq.1) write(6,38)
38    format(/2x,'Epsilon/Zeta distribution ...',/
     &       /19x,'  BI',4x,'BII',/19x,12('-'))
      write(6,34) i,slev,namv(i)(6:kfn),(nint(bin(i,j)),j=1,2)
34    format(2x,i5,') ',a3,1x,a2,1x,10i7)
      endif
      endif
      enddo
      if(ic.gt.0) then
      write(6,350) (nint(100*row(j)/ksum),j=1,2)
350   format(19x,12('-'),/11x,'  Tot',2i7)
      endif
c---------------------------------------------------O/P alpha/gamma bins
         do j=1,9
         row(j)=0.
         enddo
         ksum=0
      ic=0
      do i=1,kvar
      iv=kv(i)
      if(iv.eq.20.or.iv.eq.29) then
      if(ks(i).gt.0.and.ks(i+2).gt.0) then
      ic=ic+1
      il=kl(i)
      kfn=lenv(i)
         write(slev,'(i3)') il
         if(lseq.gt.0) then
         slev='  M'
         if(even.and.i.le.n4) slev='  L'
         if(even.and.i.gt.n4) slev='  U'
         endif
            ksum=ksum+ks(i)
            do j=1,9
            row(j)=row(j)+bin(i,j)
            bin(i,j)=100*bin(i,j)/ks(i)
            enddo
      if(ic.eq.1) write(6,36)
36    format(/2x,'Alpha/Gamma distribution ...',/
     & /18x,'  g-/g-','  g-/ t','  g-/g+','  t /g-','  t /t ',
     & '  t /g+','  g+/g-','  g+/t ','  g+/g+',/20x,61('-'))
      write(6,34) i,slev,namv(i)(6:kfn),(nint(bin(i,j)),j=1,9)
      endif
      endif
      enddo
      if(ic.gt.0) then
      write(6,360) (nint(100*row(j)/ksum),j=1,9)
360   format(20x,61('-'),/11x,'  Tot',9i7)
      endif
c----------------------------------------------------------O/P histogram
      if(histo.and.ksnap.gt.50) then
          do i=0,n4
          hmin=0
          hmax=0
          jup=ncol
          if(i.eq.0) jup=1
          do m=-n3,n3 !.............................trim min
          seen=.false.
          do j=1,jup
          iv=i+(j-1)*n4
          if(hist(iv,m).ne.0) seen=.true.
          enddo
          if(seen) then
          hmin=m
          goto 370
          endif
          enddo
370       do m=n3,-n3,-1 !..........................trim max
          seen=.false.
          do j=1,jup
          iv=i+(j-1)*n4
          if(hist(iv,m).ne.0) seen=.true.
          enddo
          if(seen) then
          hmax=m
          goto 375
          endif
          enddo
375       jbin=hmax-hmin+1
          hfac(i)=jbin/nhist
          if(mod(hfac(i),2).eq.0) hfac(i)=hfac(i)+1
          nbin(i)=jbin/hfac(i)
          jdel=(jbin-hfac(i)*nbin(i))/2
          hstr(i)=hmin+jdel-1
          enddo

      do i=0,kvar
      iloc=mod(i,n4)
      if(iloc.eq.0.and.i.gt.0) iloc=n4
      hsum(i)=0
      j=hstr(iloc)
      do k=1,nbin(iloc) !...............................compact 
      hstor(i,k)=0
         do m=1,hfac(iloc)
         j=j+1
         hstor(i,k)=hstor(i,k)+hist(i,j)
         enddo
      hsum(i)=hsum(i)+hstor(i,k)
      enddo !... k compact bins
      if(hsum(i).eq.0) hsum(i)=1
      enddo !... i variables

         kfi=index(lis,' ')-1
         do i=0,n4
         kfn=lenv(i)
         open(unit=9,file=lis(:kfi)//'_'//namv(i)(:kfn)//'.his',
     &   status='new')
         rdel=0.1d0
         if(kang(i).eq.0) rdel=0.01d0
         jd=hstr(i)
         do k=1,nbin(i)
         dist=rdel*(jd+(hfac(i)+1)/2)
         write(9,37) dist,(hstor(i+j*n4,k)*100.0/(hfac(i)*hsum(i+j*n4)),
     &   j=0,ncol-1)
37       format(2x,101f8.2)
         jd=jd+hfac(i)
         enddo
         close(unit=9)
         enddo
      endif
c--------------------------------------------O/P single structure series
      kfi=index(lis,' ')-1
      if(ksnap.eq.1) then
      nout=ncol
      if(lseq.gt.0) nout=nlev
      do i=1,n4
         kfn=lenv(i)
         open(unit=i+9,file=lis(:kfi)//'_'//namv(i)(:kfn)//'.ser',
     &   status='unknown')
      do j=0,nout-1
      if(nout.eq.ncol) then
      il=kl(i+j*n4)
         else
         il=abs(lev(j+1))
         endif
      if(abs(vstor(i+j*n4)).lt.250) then
      write(i+9,'(i8,f8.2)') il,vstor(i+j*n4)
         else
         write(i+9,'(i8,a8)') il, nastr
         endif
      enddo
      close(i+9)
      enddo
      endif
c-------------------------------------------------------O/P correlations
      if(corr) then
       start=.true.
       do i=0,kvar-1
        if(ks(i).gt.0.and..not.flag(i)) then
         iv=kv(i)
         il=kl(i)
         write(slev,'(i3)') il
         if(lseq.gt.0) then
          slev='  M'
          if(even.and.i.le.n4) slev='  L'
          if(even.and.i.gt.n4) slev='  U'
         endif
         if(i.eq.0) slev='  -'
         if(j==5)
     1      write(6,42) i,'','',j,'','',ks(i)
         do j=i+1,kvar
          jv=kv(j)
          jl=kl(j)
          if(ks(i).eq.ks(j).and..not.flag(j)) then
           ksv = ks(i)
           if(kang(iv).eq.0.and.kang(jv).eq.0) then
            cor = xcssum(i,j,1,1)/sqrt(sumsq(i)*sumsq(j))
           else
            A = xcssum(i,j,1,1)
            C = xcssum(i,j,1,2)
            D = xcssum(i,j,2,1)
            B = xcssum(i,j,2,2)
            E = cs2sum(j,1)
            F = cs2sum(j,2)
            G = cs2sum(i,1)
            H = cs2sum(i,2)
            cor = 4.d0*(A*B-C*D)/sqrt(
     1         (ksv**2.d0-E**2.d0-F**2.d0)*
     1         (ksv**2.d0-G**2.d0-H**2.d0))
           endif
           if(abs(cor).gt.cormin) then
            write(plev,'(i3)') jl
            if(lseq.gt.0) then
             plev='  M'
             if(even.and.j.le.n4) plev='  L'
             if(even.and.j.gt.n4) plev='  U'
            endif
            if(iv.eq.18.or.iv.eq.19.or.jv.eq.18.or.jv.eq.19) goto 45
            kin=lenv(i)
            kjn=lenv(j)
            if(start) then
             start=.false.
             write(6,40) cormin
 40          format(/2x,'Correlations (R > ',f4.1,')'/)
            endif
c-----------------------------------------------------Write correlation
            write(6,42) i,slev,namv(i)(:kin),j,plev,namv(j)(:kjn),cor
 42         format(2x,i5,') ',a3,1x,a9,' : ',i5,') ',a3,1x,a9,1x,f5.2)

c--------------------------------------------O/P cross correlation files
          inquire(file=lis(:kfi)//'_'//namv(i)(:kin)//'.ser',exist=ansi)
          inquire(file=lis(:kfi)//'_'//namv(j)(:kjn)//'.ser',exist=ansj)
            if(ansi.and.ansj) then
             open(unit=7,file=lis(:kfi)//'_'//namv(i)(:kin)//'.ser',
     &          status='old')
             isec=7
             if(namv(i).ne.namv(j)) then
              isec=8
              open(unit=8,file=lis(:kfi)//'_'//namv(j)(:kjn)//'.ser',
     &           status='old')
             endif
             write(str1,'(i4)') i
             mi=4
             if(i.gt.0) mi=4-int(log10(float(i))+.0001)
             write(str2,'(i4)') j
             mj=4-int(log10(float(j))+.0001)
             open(unit=9,file=lis(:kfi)//'_'//str1(mi:)//'_'
     &          //str2(mj:)//'.cor',status='new')
             icol=1+(i-1)/n4
 44          read(7,32,end=46) ksn,(vstor(m),m=1,icol)
             vi=vstor(icol)
             jcol=1+(j-1)/n4
             read(isec,32,end=46) ksn,(vstor(m),m=1,jcol)
             vj=vstor(jcol)
             write(9,32) ksn,vi,vj
             goto 44
 46          close(unit=7)
             if(isec.eq.8) close(unit=8)
             close(unit=9)
            endif
c-----------------------------------------------------------------------

           endif
          endif
 45      enddo
        endif
       enddo
      endif
c-----------------------------------------------------------------------
      write(6,*)
      end
c=======================================================================
      subroutine nml
      implicit integer*4 (i-n)
      implicit real*8 (a-h,o-z)
      parameter (n_real=1,n_int=5,n_log=3,n_cha=3,n_tot=12)
      character*48 lis,seq,vc(n_cha),nastr*48
      character lini*100,input(n_tot)*10,lj*1
      logical*2 corr,series,histo,vo(n_log),iflag(n_tot),
     & finml,lanml,start
      integer*4 first,last,skip,vi(n_int),nmls(n_tot)
      dimension vr(n_real)
      common/cha/lis,seq,nastr
      common/dat/cormin,lev1,lev2,itst,itnd,itdel,corr,series,
     & histo
      equivalence (vc(1),lis),(vr(1),cormin),(vi(1),lev1),
     & (vo(1),corr)
c-----------------------------------------------------------------------
      data input /'cormin','lev1','lev2','itst','itnd',
     & 'itdel','corr','series','histo','lis','seq','nastr'/
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
      do k=1,100
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
         if(lini(kl:kl).eq.'''') kl=kl+1
         if(lini(kh:kh).eq.'''') kh=kh-1
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
     1/5x,'****  CANAL   Version 1.3 9/2014  ****',
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
12    format(4(:,2x,a6,': ',f7.2))
      write(6,*)
      write(6,14) (input(n_real+j),vi(j),j=1,n_int)
14    format(5(:,2x,a6,': ',i6))
      write(6,*)
      write(6,16) (input(ninr+j),vo(j),j=1,n_log)
16    format(5(:,2x,a6,': ',l4))
c-----------------------------------------------------------------------
      write(6,*)
      return
50    write(6,55) lini(jl:jh)
55    format(/2x,'---- Error in namelist input for ',a,' ----'/)
      stop
      end
