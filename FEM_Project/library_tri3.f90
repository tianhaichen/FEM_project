MODULE library
CONTAINS
! **********************************************
!   该子程序的功能为：输入材料信息、坐标信息和单元信息
SUBROUTINE input(x,y,meo,ae)
	USE global_variables
	IMPLICIT NONE
	REAL(8),INTENT(INOUT)::x(:),y(:),ae(:,:)
	INTEGER,INTENT(INOUT)::meo(:,:)
	INTEGER::i,j,ielem,ip,ie
	!	  输入材料信息，共NM条
	!     每条依次输入：EO-弹性模量(kN/m**2)，   VO-泊松比，
	!                   W-材料重力集度(kN/m**2)，T-单元厚度(m)
	READ(5,*)((ae(i,j),i=1,4),j=1,nm)
	!	  输入坐标信息，共NP条
	!     每条依次输入：结点号，该结点的x坐标和y坐标  
	READ(5,*) (ip,x(ip),y(ip),i=1,np)
	!输入单元信息，共NE条
	DO ielem=1,ne
	!每条依次输入：单元号，该单元的三个结点i、j、m的整体编码及该单元的材料类型号  
		READ(5,*) ie,(meo(j,ie),j=1,4)
	END DO
	WRITE(*,500) ((ae(i,j),i=1,4),j=1,nm)
	WRITE(7,500) ((ae(i,j),i=1,4),j=1,nm)
	WRITE(*,550)
    WRITE(7,550)
	WRITE(*,600) (x(i),y(i),i=1,np)
	WRITE(7,600) (x(i),y(i),i=1,np)
	WRITE(*,*) ' ELEMENT----DATA'
	WRITE(*,650) (ielem,(meo(i,ielem),i=1,4),ielem=1,ne)
	WRITE(7,*)
	WRITE(7,*) ' ELEMENT----DATA'
	WRITE(7,650) (ielem,(meo(i,ielem),i=1,4),ielem=1,ne)
	IF(ni>0) THEN
!	  下4行：如果是平面应变问题，需把E换成E/(1-v**2),v换成v/(1-v)
		DO j=1,nm
			ae(1,j)=ae(1,j)/(1.0-ae(2,j)*ae(2,j))
			ae(2,j)=ae(2,j)/(1.0-ae(2,j))
		ENDDO
	END IF 
500 FORMAT(/1x,'EO=**VO=**W=**T=**'/(1x,4f15.4))
550 FORMAT(/1x,'X--Coordinate  Y--Coordinate')
600 FORMAT(3x,f8.3,7x,f8.3)
650 FORMAT(1x,I5,2x,4i5)
END SUBROUTINE input

! **********************************************
!   该子程序的功能为：形成结点自由度序号矩阵JR(2,NP)
SUBROUTINE MR(jr)
	USE global_variables
	IMPLICIT NONE
	INTEGER,INTENT(OUT)::jr(:,:)
	INTEGER::nn,ir(2)
	INTEGER::i,j,k,l,m
	! 赋初值，初始假设每个结点x向和y向都是自由的 
	jr=1
	!NR>0，输入约束信息，共NR条
	WRITE(7,500)
	DO i=1,nr
	!每条输入：受约束结点的点号NN及该结点x,y方向上的约束信息IR(2)(填1表示自由，填0表示受约束)
		READ(5,*) nn,ir
		jr(:,nn)=ir(:)
		WRITE(7,600) nn,ir
	ENDDO
	n=0
	! N充零，以便累加形成结构的自由度总数
	!下8行：根据每个结点约束信息，形成自由度序号指示矩阵JR(2,NP)
	DO i=1,np
		DO j=1,2
			IF(jr(j,i)>0) THEN
				n=n+1
				jr(j,i)=n
			ENDIF
		ENDDO
	ENDDO
500 FORMAT(/1x,'Constrained message'/6x,'Node No.   State')
600 FORMAT (6X,8(I5,6X,2I1))
	return
END SUBROUTINE MR
!***********************************************
!   该子程序的功能为：形成指示矩阵MA(N)
SUBROUTINE formma(meo,jr,ma)
	USE global_variables
	IMPLICIT NONE
	INTEGER,INTENT(IN)::meo(:,:),jr(:,:)
	INTEGER,INTENT(OUT)::ma(:)
	INTEGER::nn(6),nod(3)
	INTEGER::jb,m,jj,i,ie,l,jp
	!首先主元素指示矩阵MA全部充0
	ma=0
	! 单元循环，循环结束便形成每行的半带宽
element:DO ie=1,ne
	!下8行：取出单元自由度序号数组NN(6)	
		nod(:)=meo(1:3,ie)
		DO i=1,3
			jb=nod(i)
			DO m=1,2
				jj=2*(i-1)+m
				nn(jj)=jr(m,jb)
			ENDDO
		ENDDO
	!下6行：找出该单元非零自由度序号的最小值L
		l=n
		DO i=1,6
			IF(nn(i)>0) THEN
				IF((nn(i)-l)<0) l=nn(i)
			ENDIF
		ENDDO
	!下7行：循环遍历所考察单元的6个自由度，计算相应于该单元的这6个自由度
	!的半带宽JP-L+1,其中JP是自由度序号亦即是方程序号，若JP-L+1大于已存放在
	!MA(JP)中的值，则将它赋值给MA(JP)，对单元的循环结束，则数组MA(JP)存贮了
	!各个方程的半带宽。
		DO m=1,6
			jp=nn(m)
			IF(jp /= 0) THEN
				IF((jp-l+1)>ma(jp)) ma(jp)=jp-l+1
			ENDIF
		ENDDO
	ENDDO element
!     下2行：MX准备存贮最大半带宽，赋初值0；数组MA(N)准备存贮各个方程的主对角元素在
!     一维数组中的序号，全部赋初值1。 
	mx=0
	ma(1)=1
!     下4行：通过对结构的全部自由度循环（第1个自由度除外）找出其中的最大半带宽MX,
!     并将各个方程的半带宽逐个累加获得主对角元素序号指示向量MA(N)
	DO i=2,n
		if(ma(i)>mx) mx=ma(i)
		ma(i)=ma(i)+ma(i-1)
	ENDDO
!   下1行：最后一个方程(即第N个方程)的主对角元素的序号就是按一维存储整体劲度矩阵[K]所需存储的元素总数NH
	nh=ma(n)
	WRITE(*,500) n,nh,mx
	WRITE(7,500) n,nh,mx
	return
500 FORMAT(1x,'Total degrees of freedom N =',i5,/1x,'Total-Storage          NH=',i5,/1x,'Max-Semi-Bandwidth        Mx=',i5)
END SUBROUTINE formma
! *************************************************************
!   该子程序的功能为：形成整体劲度矩阵,并存放在一维数组SK(NH)中
SUBROUTINE mgk(ae,x,y,meo,jr,ma,sk)
	USE global_variables
	IMPLICIT NONE 
	REAL(8),INTENT(IN)::ae(:,:),x(:),y(:)
	REAL(8),INTENT(OUT)::sk(:)
	INTEGER,INTENT(IN)::meo(:,:),ma(:),jr(:,:)
	REAL(8)::ske(6,6)
	INTEGER::ie,nod(3),nn(6),i,j,j2,j3,jj,jk,jl,jm,jn
	REAL(8)::eo,vo,w,t,a,bi(3),ci(3),h11,h12,h21,h22
	!下2行，将劲度矩阵的一维数组赋初值0 
	sk=0.0
element:DO ie=1,ne
	!对单元循环,循环结束便形成SK 
	!下10行：计算单元劲度矩阵SKE(6,6),在对单元循环的循环体内,通过I,J的二重循环
	!来调用子程序KRS,形成单元劲度矩阵式中9个子矩阵,并贮存在数组SKE(6,6) 中   
		CALL div(ie,ae,x,y,meo,nod,eo,vo,w,t,a,bi,ci)
		DO i=1,3
			DO j=1,3
				CALL krs(eo,vo,a,t,bi(i),bi(j),ci(i),ci(j),h11,h12,h21,h22)
				ske(2*i-1,2*j-1)=h11
				ske(2*i-1,2*j)=h12
				ske(2*i,2*j-1)=h21
				ske(2*i,2*j)=h22 
			ENDDO
		ENDDO
		!下7行：将单元自由度序号写入数组NN(6)
		DO i=1,3
			j2=nod(i)
			DO j=1,2
				j3=2*(i-1)+j
				nn(j3)=jr(j,j2)
			ENDDO
		ENDDO
	!下12行：把单元劲度矩阵的元素叠加到整体劲度矩阵SK的相应位置中去
	!由于劲度矩阵的对称性，只存储矩阵的对角元素和下三角部分
	DO i=1,6
		DO j=1,6
			IF((nn(j)/=0) .and. (nn(i)>=nn(j))) THEN
				jj=nn(i)	!JJ整体劲度矩阵的行号
				jk=nn(j)	!JK整体劲度矩阵的列号 
				jl=ma(jj)	!第JJ行对角元素在一维存储数组SK中的位置
				jm=jj-jk 
				jn=jl-jm 	!整体劲度矩阵的第JJ行第JK列的元素K[JJ,JK]在一维存储数组SK中的位置   
				sk(jn)=sk(jn)+ske(i,j)	!把第JJ行第JK列的元素叠加到一维存储数组SK的相应位置中去
			ENDIF
		ENDDO
	ENDDO
	ENDDO element
	return 
END SUBROUTINE mgk 
! **********************************************
!   该子程序的功能为：取出单元单元类型号，并计算单元的BI,CI矩阵等
SUBROUTINE div(jj,ae,x,y,meo,nod,eo,vo,w,t,a,bi,ci)
	IMPLICIT NONE
	INTEGER,INTENT(IN)::jj,meo(:,:)
	INTEGER,INTENT(OUT)::nod(:)
	REAL(8),INTENT(IN)::ae(:,:),x(:),y(:)
	REAL(8),INTENT(OUT)::eo,vo,w,t,a,bi(:),ci(:)
	INTEGER::i,j,m,l 
	!下6行：从单元信息数组MEO中取出JJ号单元的三个结点i,j,m的整体编码，并把它们放在NOD(3)中
	I=meo(1,jj)
	nod(1)=i
	j=meo(2,jj)
	nod(2)=j 
	m=meo(3,jj)
	nod(3)=m 
	l=meo(4,jj) 	!L是该单元的材料类型号
	!下6行：BI(1),BI(2),BI(3)分别为bi,bj,bm,
	!       CI(1),CI(2),CI(3)分别为ci,cj,cm,
	bi(1)=y(j)-y(m)
	ci(1)=x(m)-x(j)
	bi(2)=y(m)-y(i)
	ci(2)=x(i)-x(m)
	bi(3)=y(i)-y(j)
	ci(3)=x(j)-x(i)
	a=(bi(2)*ci(3)-bi(3)*ci(3))/2.0 !A为单元面积
	!下4行：根据材料号L取出该种材料的4个参数
	eo=ae(1,l)
	vo=ae(2,l)
	w=ae(3,l)
	t=ae(4,l)
	return	
END SUBROUTINE div 
! *********************************************
!   该子程序的功能为：计算单元劲度矩阵中的子块Krs
SUBROUTINE krs(eo,vo,a,t,br,bs,cr,cs,h11,h12,h21,h22)
	IMPLICIT NONE
	REAL(8),INTENT(IN)::br,bs,cr,cs,eo,vo,a,t
	REAL(8),INTENT(OUT)::h11,h12,h21,h22 
	REAL(8)::et,v 
	!H11,H12,H21,H22为单元劲度矩阵中子块Krs的4个元素
	et=eo*t/(1.0-vo*vo)/4.0/a 
	v=(1.0-vo)/2.0
	h11=et*(br*bs+v*cr*cs)
	h12=et*(vo*br*cs+v*cr*bs)
	h21=et*(vo*cr*bs+v*br*cs)
	h22=et*(cr*cs+v*br*bs)
	return
END SUBROUTINE krs 
! **********************************************
!   该子程序的功能为：形成整体等效结点荷载列阵{R}
!   程序可以考虑自重和集中力两种荷载情况
SUBROUTINE load(ae,x,y,meo,jr,r)
	USE global_variables
	IMPLICIT NONE 
	REAL(8),INTENT(IN)::ae(:,:),x(:),y(:)
	INTEGER,INTENT(IN)::meo(:,:),jr(:,:)
	REAL(8),INTENT(OUT)::r(:)
	REAL(8)::eo,vo,w,t,a,bi(3),ci(3)
	REAL(8),ALLOCATABLE::fv(:,:)
	INTEGER,ALLOCATABLE::nf(:)
	INTEGER::i,j,j2,j3,ie,jj,m,nod(3)
	r=0.0 ! 将R置零,用来存放等效结点荷载 
	!下12行：计算自重引起的结点等效荷载(对于需要考虑自重荷载的情况)
	!通过对单元的循环将每个单元重量的1/3叠加在R(N) 中与3个结点的y方向自由度相应的位置上   
	IF(ng>0) THEN
		DO ie=1,ne 
			CALL div(ie,ae,x,y,meo,nod,eo,vo,w,t,a,bi,ci)
			DO i=1,3
				j2=nod(i)
				j3=jr(2,j2)
				IF(j3>0) THEN
					r(j3)=r(j3)-t*w*a/3.0
				ENDIF
			ENDDO
		ENDDO
	ENDIF
	!对于有集中载荷作用的情况，输入集中载荷信息
	IF(nl>0) THEN
		!输入荷载信息，共NL条
		!每条依次输入：NF--受集中力的结点号   FV--该结点的x，y方向的荷载分量
		ALLOCATE(fv(2,nl),nf(nl))
		READ(5,*) (nf(j),(fv(i,j),i=1,2),j=1,nl)
		WRITE(*,500) (nf(i),i=1,nl)
		WRITE(7,500) (nf(i),i=1,nl)
		WRITE(*,600) ((fv(i,j),i=1,2),j=1,nl)
		WRITE(7,600) ((fv(i,j),i=1,2),j=1,nl)
		!下8行：通过对NL个受集中荷载作用结点的循环,将集中荷载叠加在R(N)的相应位置上,
		!循环结束形成等效结点荷载列阵R(N)
		DO i=1,nl
			jj=nf(i)
			j=jr(1,jj)
			m=jr(2,jj)
			IF(j>0) r(j)=r(j)+fv(1,i)
			IF(m>0) r(m)=r(m)+fv(2,i)
		ENDDO 
	ENDIF
500 FORMAT(/1x,'Nodes of applied load***nf='/(1x,10i8))
600 FORMAT(/1x,'Constrained-loads***fv='/(1x,6f12.3))
END SUBROUTINE load 
! **********************************************
!   该子程序的功能为：输出结点等效荷载或结点位移
SUBROUTINE output(jr,r)
	use global_variables
	IMPLICIT NONE 
	INTEGER,INTENT(IN)::jr(:,:)
	REAL(8),INTENT(IN)::r(:)
	INTEGER::i,l 
	REAL(8)::s,ss 
	!下16行：循环遍历所有结点，输出等效结点荷载或结点位移
	DO i=1,np 
		l=jr(1,i)
		IF(l>0) THEN
			s=r(l)
		ELSE 
			IF(l==0) s=0.0	!约束自由度方向赋零值 
		ENDIF
		l=jr(2,i)
		IF(l>0) THEN 
			ss=r(l)
		ELSE
			IF(l==0) ss=0.0 !约束自由度方向赋零值 
		ENDIF
		WRITE(*,500) i,s,ss 
		WRITE(7,500) i,s,ss 
	ENDDO
500 FORMAT(5x,i5,2e20.5)
	return 	
END SUBROUTINE output
! ***************************************************************
!     该子程序的主要功能：输入非零已知位移信息，并处理荷载项及对应K中主元素
SUBROUTINE treat(sk,ma,r,jr,ndi,dv)
	USE global_variables
	IMPLICIT NONE
	INTEGER,INTENT(INOUT)::ma(:),jr(:,:),ndi(:)
	REAL(8),INTENT(INOUT)::sk(:),r(:),dv(:,:)
	INTEGER::i,j,l,jn,jj 
	WRITE(7,500)
	!下4行：通过对NDP个已知位移的结点循环，输入非零已知位移的结点号NDI
	!和x,y两个方向的位移分量DV，并输出非零已知位移信息
	DO i=1,ndp
		READ(5,*) ndi(i),(dv(j,i),j=1,2)
		WRITE(7,600) ndi(i),(dv(j,i),j=1,2)
	ENDDO
	!下11行：通过对NDP个已知位移的结点循环，在相应的自由度的主元素赋于
	!一个大数(1E30)，并在荷载的相应行赋于该位移值乘以相同的大数(1E30)
	DO i=1,ndp 
		jj=ndi(i)
		DO j=1,2
			l=jr(j,jj)
			jn=ma(l)
			IF(dv(j,i)/=0.0) THEN
				sk(jn)=1.0e30
				r(l)=dv(j,i)*1.E30
			ENDIF
		ENDDO
	ENDDO
500 FORMAT(/1x,'Known Displacement Nodes And (x,y) values')
600 FORMAT(10x,i8,10x,f10.6,20x,f10.6)
	return
END SUBROUTINE treat

!****************************************************************
!     该子程序的功能为：整体劲度矩阵的分解运算
SUBROUTINE decomp(sk,ma)
	USE global_variables
	IMPLICIT NONE 
	REAL(8),INtent(INOUT)::sk(:)
	INTEGER,INTENT(IN)::ma(:)
	INTEGER::i,j,k,l,m,l1,mp,lp,ip,lpp,ij,jp,ii 
	DO i=2,n
	!下4行：在对行循环的循环体内，计算出第I行第1个非零元素所在的列号L，
	!以及对列循环的下界L1与上界K，若L1 >K，则跳过对列的循环 
		l=i-ma(i)+ma(i-1)+1
		k=i-1
		l1=L+1
		IF(l1<=k) THEN
			DO j=l1,k 
	!	  下5行：在对列循环的循环体内，首先计算出[K]的第I行和第J列的元素在一维数组SK(NH)
	!     中的序号IJ以及第J行第一个非零元素所在列号M，并将L与M比较，较大者存贮在M中，于
	!     是M成为求总和循环时的下界，MP是上界，当M>MP则跳过求总和的循环   
				ij=ma(i)-i+j 
				m=j-ma(j)+ma(j-1)+1
				IF (l .GT. m) m=l 
				mp=j-1 
				IF(m <= mp) THEN
				!下5行：在求总和的循环体内计算出[K]中第I行第P列（程序内用LP）的元素以及第J行第P
				!列的元素在一维数组SK(NH)中的序号IP与JP，总和循环结束Kij成为lij 。
					DO lp=m,mp 
						ip=ma(i)-i+lp 
						jp=ma(j)-j+lp 
						sk(ij)=sk(ij)-sk(ip)*sk(jp)
					ENDDO 
				ENDIF
			ENDDO
		ENDIF
		IF(l<=k)THEN
			DO lp=l,k 
			!下4行：将第I行的各个副元素lip变成uip(程序P用LP表示）,
			!IP与LPP分别是第I行第P列的元素与第P行的主对角线元素在一维数组SK（NH)中的序号
				ip=ma(i)-i+lp 
				lpp=ma(lp) 
				IF(ABS(sk(lpp))<1.0e-12) WRITE(*,*) 'Decom zero error ip=',ip !主对角线元素为零，程序报错
				sk(ip)=sk(ip)/sk(lpp)
				!下2行：将第I行的主对角线元素Kii变成Lii对行循环结束,分解完毕
				ii=ma(i) 
				sk(ii)=sk(ii)-sk(ip)*sk(ip)*sk(lpp)
			ENDDO
		ENDIF
	ENDDO
500 FORMAT(/10X,'SK='/(1X,6F12.4))
    RETURN 
END SUBROUTINE decomp
! *********************************************
!     该子程序的功能为：前代，回代求出未知结点位移并存放在R(N)中
SUBROUTINE foba(sk,ma,r)
	USE global_variables
	IMPLICIT NONE 
	REAL(8),INTENT(IN)::SK(:)
	REAL(8),INTENT(INOUT)::R(:)
	INTEGER,INTENT(IN)::MA(:)
	INTEGER::i,j,k,l,lp,ii,ip,j1,ij 
	!前代过程
	DO i=2,n 
	!L与IP分别为第I行的第一个非零元素所在列号以及
	!第P（程序中用LP）个元素在一维数组SK(NH) 中的序号
	l=i-ma(i)+ma(i-1)+1
	k=i-1 
	IF(l <= k) THEN
		DO lp=l,k 
			ip=ma(i)-i+lp 
			r(i)=r(i)-sk(ip)*r(lp)
		ENDDO
	ENDIF
	ENDDO
	!下5行：实现式 gi=fi/lii (i=1,2,…,n)
	DO i=1,n 
		ii=ma(i)
		IF(ABS(sk(ii))<1.0e-12) WRITE(*,*) 'FOBA ZERO ERROR i=',i 
		r(i)=r(i)/sk(ii)
	ENDDO
	!回代过程
	DO j1=2,n 
	!L与IJ 分别为第I行第一个非零元素所在列号与第J个元素在一维数组SK(NH)中的序号，
	!这里对I从N到2的循环是通过对J1从2到N循环以及I与J1之间的关系I=2+N-J1来实现的
		i=2+n-j1 
		l=i-ma(i)+ma(i-1)+1
		k=i-1
		IF(l<=k) THEN
			do j=l,k 
				ij=ma(i)-i+j 
				r(j)=r(j)-sk(ij)*r(i)
			ENDDO
		ENDIF
	ENDDO
	return
END SUBROUTINE foba
! **********************************************
!   该子程序的功能为：计算并输出单元应力
SUBROUTINE ces(ae,x,y,meo,jr,r) 
	USE global_variables
	IMPLICIT NONE 
	REAL(8),INTENT(INOUT)::ae(:,:),x(:),y(:),r(:)
	REAL(8)::et,b(6),h11,h12,h21,h22,h1,h2,h3,a1,a2,a3,bi(3),ci(3),eo,vo,w,t,a 
	INTEGER,INTENT(IN)::meo(:,:),jr(:,:)
	INTEGER::nod(3),ie,i,j2,i2,i3 
	DO ie=1,ne 
		CALL div(ie,ae,x,y,meo,nod,eo,vo,w,t,a,bi,ci) 
		et=eo/(1.0-vo*vo)/a/2.0
		!在单元循环体中把所考察单元的3个结点的结点位移值送入数组B(6)
		DO i=1,3 
			j2=nod(i)
			i2=jr(1,j2)
			i3=jr(2,j2)
			IF(i2>0)THEN
				B(2*i-1)=r(i2)
			ELSE 
				IF(i2==0) b(2*i-1)=0.0  !约束方向位移为零
			ENDIF 
			IF(i3>0)THEN
				B(2*i)=r(i3)
			ELSE 
				IF(i3==0) b(2*i)=0.0  !约束方向位移为零
			ENDIF 
		ENDDO
		!下11行：计算应力sig-x、sig-y和sig-xy
		h1=0.0
		h2=0.0
		h3=0.0
		DO i=1,3
			h1=h1+bi(i)*b(2*i-1)
			h2=h2+ci(i)*b(2*i)
			h3=h3+bi(i)*b(2*i)+ci(i)*b(2*i-1)
		ENDDO
		a1=et*(h1+vo*h2)
		a2=et*(h2+vo*h1)
		a3=et*(1.0-vo)*h3/2.0
		!下4行：计算主应力并存储在B(4)和B(5)中
		h1=a1+a2 
		h2=sqrt((a1-a2)*(a1-a2)+4.0*a3*a3)
		b(4)=(h1+h2)/2.0 
		b(5)=(h1-h2)/2.0 
		!下9行：计算应力主向并存储在B(6)中
		IF(abs(a3)<=1.0e-4) THEN 
			IF(a1 <= a2) THEN
				b(6)=90.0
			ELSE
				b(6)=0.0
			ENDIF
		ELSE 
			B(6)=atan((B(4)-a1)/a3)/57.29578
		ENDIF
		!下3行：将应力sig-x、sig-y和sig-xy存储在B(1)、B(2)和B(3)中
		b(1)=a1 
		b(2)=a2 
		b(3)=a3 
		!下2行：输出应力分量、主应力和应力主向
		WRITE(*,500) ie,b 
		WRITE(7,500) ie,b 
	ENDDO
500 FORMAT(6x,i4,3x,6f11.3)
return
END SUBROUTINE ces 
! **********************************************
!   该子程序的功能为：计算并输出支座反力
SUBROUTINE erfac(ae,x,y,meo,jr,r)
	USE global_variables
	IMPLICIT NONE
	REAL(8),INTENT(INOUT)::ae(:,:),x(:),y(:),r(:)
	INTEGER,ALLOCATABLE::nci(:),nce(:,:)
	INTEGER,INTENT(INOUT)::meo(:,:),jr(:,:)
	INTEGER::nod(3),i,j,jj,l,ie,m,im,ip,k,ji,jp,flag 
	REAL(8)::fx,fy,s,ss,h11,h12,h21,h22,bi(3),ci(3),eo,vo,w,t,a 
	ALLOCATE(nci(nc),nce(4,nc))
	!下6行：输入并输出数组NCI以及NCE
	READ(5,*) (nci(j),j=1,nc)
	READ(5,*) ((nce(i,j),i=1,4),j=1,nc)
	WRITE(*,500) (nci(j),j=1,nc)
	WRITE(*,600) ((nce(i,j),i=1,4),j=1,nc)
	WRITE(7,500) (nci(j),j=1,nc)
	WRITE(7,600) ((nce(i,j),i=1,4),j=1,nc)
	WRITE(*,700)
	WRITE(7,700)
	DO jj=1,nc 
	!下2行：FX、FY充零，准备存贮所考察结点的两个坐标方向的支座反力
		fx=0.0
		fy=0.0
		l=nci(jj) !所考察的的计算支座反力结点的点号
		DO m=1,4
!在对4个相关单元做循环的循环体内，将所考察的相关单元的单元号存贮在IE中，并调用子程序DIV计算该单元的有关数据
			IF(nce(m,jj)<=0) cycle
			ie=nce(m,jj)
			CALL div(ie,ae,x,y,meo,nod,eo,vo,w,t,a,bi,ci)
!对所考察相关单元的3个结点循环，确定哪个结点就是所考察的计算支座反力结点,将该结点的局部编码序号存贮在K中。
			DO im=1,3
				k=im 
				flag=0
				IF(l==nod(im)) THEN
					flag=1
					EXIT
				ENDIF 
			ENDDO
			IF(flag==0) WRITE(0,750) l !!若找不到这个结点，则相关单元号有错，打印出错信息
!   下18行：对所考察相关单元的3个结点循环，计算每个结点的结点位移贡献给计算支座反力结点的结点力
			DO ip=1,3
				CALL krs(eo,vo,a,t,bi(k),bi(ip),ci(k),ci(ip),h11,h12,h21,h22)
				nl=nod(ip)
				ji=jr(1,nl) ! 所考察结点的x方向的自由度序号
				jp=jr(2,nl) ! 所考察结点的y方向的自由度序号
				IF(ji==0) THEN
					s=0.0  !约束自由度方向位移为零
				ELSE IF(ji>0) THEN
					S=R(JI) !存贮所考察结点的x方向的位移值
				ENDIF 
				IF(jp<=0) THEN 
					ss=0.0 !约束自由度方向位移为零
				ELSE 
					ss=r(jp) !存贮所考察结点的y方向的位移值
				ENDIF
				fx=fx+h11*s+h12*ss 
				fy=fy+h21*s+h22*ss 
			ENDDO
		ENDDO
		!下2行：无约束方向没有支座反力
		IF(jr(1,l)/=0) fx=0.0
		IF(jr(2,l)/=0) fy=0.0
		!下2行：输出支座反力
		WRITE(*,800) l,fx,fy 
		WRITE(7,800) l,fx,fy 
	ENDDO
500 FORMAT(30x,'Node NO.**NCI='/(1x,10i8))
600 FORMAT(30x,'Element-No.**NCE='/(1x,4i6))
700 FORMAT(30x,'Nodal Reactions'/8x,'Node',14x,'X_Comp',14x,'Y_Comp')
750 FORMAT(/10x,'Error of Element Message','****Node number',i5)
800 FORMAT(6x,i5,2f20.3)
return
END SUBROUTINE erfac 
END MODULE