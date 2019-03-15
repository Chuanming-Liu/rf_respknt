C
C
C     Modified by Chuanming Liu. (CU-Boulder)
C     2019-03-12
C
C
C----------------------------------------------------------------------C
      SUBROUTINE respknt_interface(slowness, dt_in, tduring, nl,
     x rho_in, thick_in, vp_in, vs_in, wavez, waver)
C PARAMETERS
C npts      -    points of waveform
C nl        -    layer of model
C dt_in     -    sampling interval
C tduring   -    signal duration
C waveZ     -    output Z component waveforms
C waveR     -    output R component waveforms
C input parameters
      REAL slowness, dt_in
      INTEGER npts, nl
      REAL rho_in(50), thick_in(50)
      REAL vp_in(50), vs_in(50)
      REAL tduring
C set up limit maxph : dt related
      integer maxph
c     parameter (maxph=4096)
      parameter (maxph=8192)
C output parameters
      REAL wavez(maxph),waver(maxph)
Cf2py intent(out) wavez, waver
      real alfm(50),betm(50),qpm(50),qsm(50),rhom(50),thikm(50),
     *     ta(50),tb(50)
      complex u0(maxph),w0(maxph),u1(maxph),w1(maxph),tn(maxph)
      common /cmparr/u0,w0,u1,w1,tn
      common /innout/ inunit,ounit
      complex dvp,dvs,drp,drs,dts,p,fr
      real*8 wq,t1,t2,qa,qb,qabm,vabm
      character*32 ofil,ofilz,ofilr,ofilt
      character*32 modela,title
C      character*6  comp(3)
      character*1 complt,modcnv
      integer*2 rvb, cnv
      integer inunit,ounit,ipors
C      integer blank
      real dum1(100), dum2(100)
c backup
      integer nft
      real delf
      integer numpts
      integer num1

      real u0_real(maxph),u0_img(maxph),w0_real(maxph),w0_img(maxph)



      include 'kennet.inc'
C      data comp/'_sp.z ','_sp.r ','_sp.t '/
c
C       call inihdr
C       call newhdr
C       inunit=5
C       ounit=6
       twopi = 8.*atan(1.)
c
C----------------------------------------------------------------------C
C                1. Load Model

C      ofil = '                                '
C      ofilr = '                                '
C      ofilz = '                                '
C      ofilt = '                                '

c
C      write(ounit,*) 'Velocity Model Name'
C      read(inunit,'(a)') modela

C      iblank=blank(modela)
C      ofil(1:iblank) = modela(1:iblank)
      nlyrs = nl
      alfm = vp_in
      betm = vs_in
      rhom = rho_in
      thikm = thick_in
C      call rdlyrs(modela,nlyrs,title,alfm,betm,rhom,thikm,
C     *            dum1,dum2,dum1,dum2,-1,ier)

      do 1 i=1,nlyrs
*      qpm(i) = 500.
*      qsm(i) = 225.
       qpm(i) = 125
       qsm(i) = 62.5
       ta(i) = .16
       tb(i) = .26
 1    continue
C----------------------------------------------------------------------C
c
c     	2.  terminal input
c
C      write(ounit,*) 'incident P(1) or S(2) wave'
C      read(inunit,*) ipors
C      write(6,*) 'sampling interval'
C      read(5,*) dt
C      write(6,*) 'signal duration'
C      read(5,*) t
c     write(6,*) 'incident delay'
c     read(5,*) tdelay
c     write(6,*) 'output file base name'
c     read(5,'(a)') ofil
C      write(6,*) ' enter slowness: '
C      read(5,*) pr
C      write(6,*) ' partial(p) or full(f) : '
C      read(5,'(a1)') complt
C      write(6,*) ' mode conversions? (y or n) '
C      read(5,'(a1)') modcnv
      ipors = 1
      dt = dt_in
      t = tduring
      pr = slowness
      complt = 'f'
      modcnv = 'y'
C----------------------------------------------------------------------C
c      3. set up parameters
c      build output filenames
c
C      ofilz(1:iblank+6)=ofil(1:iblank)//comp(1)
C      ofilr(1:iblank+6)=ofil(1:iblank)//comp(2)
C      ofilt(1:iblank+6)=ofil(1:iblank)//comp(3)
c
c     set up the spectral parameters
c
      numpts=ifix(t/dt+1.5)
      nft=npowr2(numpts)
      nfpts=nft/2+1
      fny=1./(2.*dt)
      delf=2.*fny/float(nft)
      t=dt*nft
c
c     set up some computational parameters
c          specifying the type of response
c          requested.
c
      p = cmplx(pr,0.)
      if ( complt(1:1) .eq. 'f' )  then
         rvb = allrvb
       else
         rvb = norvb
       endif
      if ( modcnv(1:1) .eq. 'n' ) then
       cnv = prmphs
       else
       cnv = allphs
      endif
C----------------------------------------------------------------------C
c
c     4.  compute q, alfa, and beta at 1 hz for absorbtion band
c
      t1 = 1.0d04
      wq = twopi
      do 5 i = 1, nlyrs
         qa = qpm(i)
         qb = qsm(i)
         t2 = ta(i)
         alfa(i) = alfm(i) * vabm(wq,t1,t2,qa)
         t2 = tb(i)
         beta(i) = betm(i) * vabm(wq,t1,t2,qb)
         qa = qabm(wq,t1,t2,qa)
         qb = qabm(wq,t1,t2,qb)
         alfa(i) = alfa(i)*( 1. + (0.,0.5)/qa)
         beta(i) = beta(i)*( 1. + (0.,0.5)/qb)
         cnvrsn(i) = cnv
         reverb(i) = rvb
         rho(i) = rhom(i)
 5       thik(i) = thikm(i)
      cnvrsn(0) = cnv
      if ( complt(1:1) .ne. 'f' )  then
         reverb(1) = onervb
      endif
c
C----------------------------------------------------------------------C
c      5. compute kennett's interface matricies for n layer model
      fr = cmplx(1.,0.)
      call ifmat(1,p,fr,nlyrs)
C----------------------------------------------------------------------C
c
c      6. compute receiver function - free surface displacement
      do 10 i = 1, nfpts-1
         fr = cmplx(delf * ( i - 1 ), 0. )
         wq = twopi * fr
         do 6 j = 1, nlyrs
            qa = qpm(j)
            qb = qsm(j)
            t2 = ta(j)
            alfa(j) = alfm(j) * vabm(wq,t1,t2,qa)
            t2 = tb(j)
            beta(j) = betm(j) * vabm(wq,t1,t2,qb)
            qa = qabm(wq,t1,t2,qa)
            qb = qabm(wq,t1,t2,qb)
            alfa(j) = alfa(j)*( 1. + (0.,0.5)/qa)
            beta(j) = beta(j)*( 1. + (0.,0.5)/qb)
 6       continue
         call rcvrfn(p,fr,nlyrs,dvp,dvs,drp,drs,dts)
         u0(i) = dvp * (0.,-1.)*(-1.,0.)
         w0(i) = drp
         u1(i) = dvs
         w1(i) = drs * (0.,1.)
         tn(i) = dts
10    continue
      u0(nfpts) = (0.,0.)
      w0(nfpts) = (0.,0.)
      u1(nfpts) = (0.,0.)
      w0(nfpts) = (0.,0.)
      tn(nfpts) = (0.,0.)
C----------------------------------------------------------------------C
c
c     7. output the responses and ifft (u0, w0)
c
      call dfftr(u0,nft,'inverse',delf)
      call dfftr(w0,nft,'inverse',delf)
c     8. output from complex u0 (z), w0 (r) to real
C	     WRITE(*,*)'u0 =\n',u0(1:numpts)
C						WRITE(*,*)'w0 =\n',w0(1:numpts)
      u0_real = real(u0)
      u0_img = aimag(u0)
      w0_real = real(w0)
      w0_img = aimag(w0)
      num1 = numpts*2-1
      wavez(1:num1:2) = u0_real(1:numpts)
      wavez(2:num1+1:2) = u0_img(1:numpts)
      waver(1:num1:2) = w0_real(1:numpts)
      waver(2:num1+1:2) = w0_img(1:numpts)
C      if(ipors.eq.1)then
C         call dfftr(u0,nft,'inverse',delf)
C         call dfftr(w0,nft,'inverse',delf)
C         call wsac1(ofilz,u0,numpts,0.,dt,nerr)
C 	       call wsac1(ofilr,w0,numpts,0.,dt,nerr)
C      else
C         call dfftr(u1,nft,'inverse',delf)
C         call dfftr(w1,nft,'inverse',delf)
C         call dfftr(tn,nft,'inverse',delf)
C         call wsac1(ofilz,u1,numpts,0.,dt,nerr)
C         call wsac1(ofilr,w1,numpts,0.,dt,nerr)
C         call wsac1(ofilt,tn,numpts,0.,dt,nerr)
C      endif
c
C      stop
      end


C      integer function blank(file)
C      character file*32
C      do 1 i=1,32
C      if(file(i:i).ne.' ') goto 1
C      blank=i-1
C      return
C1     continue
C      write(1,100) file
C100   format(' no blanks found in ',a32)
C      blank = 0
C      return
C      end
