MODULE library
CONTAINS 
! ***************************************************************
!     该子程序的主要功能：输入9个控制变量
SUBROUTINE control 
	USE global_variables 
	IMPLICIT NONE 
	!下1行：输入结点总数 NP ，单元总数 NE ，材料类型数目 NM ，约束结点总数 NR ，水荷载类型数 NW ，
	!受集中荷载结点的数目 NF ，问题性质标识 NI (平面应力填0，平面应变填1)，非零已知位移结点总数 NDP
	!受分布力或者水压力的单元数目 NESW
    READ(15,*)NP,NE,NM,NR,NW,NF,NI,NDP,NESW
    WRITE(16,600)NP,NE,NM,NR,NW,NF,NI,NDP,NESW
	WRITE(*,600)NP,NE,NM,NR,NW,NF,NI,NDP,NESW
600 FORMAT(//1X,'NUMBER OF NODE----------------------------------------------NP=',I5/  &
             1X,'NUMBER OF ELEMENT-------------------------------------------NE=',I5/  &
             1X,'NUMBER OF MATERIAL------------------------------------------NM=',I5/  &
             1X,'NUMBER OF CONSTRAINT----------------------------------------NR=',I5/  &
             1X,'NUMBER OF WATER PRESS KIND----------------------------------NW=',I5/  &
             1X,'NUMBER OF CONCENTRATE LOAD----------------------------------NF=',I5/  &
             1X,'PLANE STRESS OR PLANE STAIN---------------------------------NI=',I5/  &
             1X,'NUMBER OF KNOWN-DISPLACEMENT--------------------------------NDP=',I5/ &
             1X,'NUMBER OF ELEMENT WITH DISTRIBUTING FORCE OR HYDROPRESS----NESW=',I5/)
      RETURN
END SUBROUTINE control 
! ***************************************************************
!     该子程序的主要功能：输入原始数据，形成COOR等数组，处理平面应变问题的E和v
SUBROUTINE input(coor,meo,ae,wg)
	USE global_variables
	IMPLICIT NONE 
	REAL(8),INTENT(INOUT)::coor(:,:),ae(:,:),wg(:,:)
	INTEGER,INTENT(INOUT)::meo(:,:)
	REAL(8)::xy(2),ir(2),e,u
	INTEGER::ie,i,j,nn,l,nee 
	!下9行：对结点循环输入结点号NN及该点的x,y坐标值，并形成结点坐标数组COOR 
	DO i=1,np 
		read(15,*) nn,xy 
		IF(nn/=i) THEN 
			WRITE(16,750) nn,i 
			WRITE(*,750) nn,i 
			stop
		ENDIF
		coor(:,nn)=xy(:)
	ENDDO
	WRITE(16,800) (nn,(coor(j,nn),j=1,2),nn=1,np)
	WRITE(*,800) (nn,(coor(j,nn),j=1,2),nn=1,np)
	!下4行：输入单元信息，共NE条 
	!输入单元序号NEE，该单元的材料类型号NME，是否计算自重信息NET(1-计算，0-不计算)，
	!单元的4个结点号
	DO ie=1,ne 
		READ(15,*) nee,(meo(j,nee),j=1,6)
		WRITE(16,850) nee,(meo(j,nee),j=1,6)
		WRITE(*,850) nee,(meo(j,nee),j=1,6)
	ENDDO
	!下2行：对材料类型数循环，输入材料特性常数：弹性模量、泊松比、材料重力集度和单元厚度
	!其中弹性模量、泊松比要输入平面应力的值
	READ(15,*) ((ae(i,j),i=1,4),j=1,nm)
	WRITE(16,910) (j,(ae(i,j),i=1,4),j=1,nm)
	WRITE(*,910) (j,(ae(i,j),i=1,4),j=1,nm)
	!下8行：若是平面应变问题，E换成E/(1.0-v**2)和v换成v/(1.0-v)
	IF(ni/=0) THEN 
		DO j=1,nm
			u=ae(1,j)
			u=ae(2,j)
			ae(1,j)=e/(1.0-u*u)
			ae(2,j)=u/(1.0-u)
		ENDDO 
	ENDIF 
	!下4行：若有水荷载作用，则输入水荷载特性数组WG(2,NW)，每种水荷载信息依次为水位y坐标值及水的容重
	IF(nw/=0) THEN 
		READ(15,*) ((wg(i,j),i=1,2),j=1,nw)
		WRITE(16,960) (j,(wg(i,j),i=1,2),j=1,nw)
		WRITE(*,960) (j,(wg(i,j),i=1,2),j=1,nw)
	ENDIF 
750 FORMAT(1x,'***Fatal Error***',/,'CARDS ',&
			'Input',i5,'is not equal to',i5)
800 FORMAT(1x,'Node NO.',2x,'X-Coordinat',2x,'Y-Coordinat',&
			/(1x,i5,4x,2f14.4))
850 FORMAT(3x,' NEE= NOD= NME= NET='/&
			3x,i5,4i5,2i5)
910 FORMAT(/20x,'Material Properties',/,&
			2x,'N.M.',5x,'Youngs Moudlus',5x,'Poision ratio',&
			4x,'Unit Weight    Width'/&
			(1x,i5,4e16.4))
960 FORMAT(/5x,'Parameteres of water and ',&
			'silt pressure'/2x,'N.P.',2x,'Zero-Pressure',&
			'Surface',8x,'Unit Weight'/(1x,i5,2f15.5))
END SUBROUTINE input 
END MODULE