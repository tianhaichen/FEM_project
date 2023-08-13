!******************************************************************************
!     1、本程序能计算弹性力学中的平面应力问题和平面应变问题；
!     2、采用三结点三角形单元；
!     3、能考虑自重和结点集中力两种荷载的作用，在计算自重时y轴取垂直向上为正。
!******************************************************************************
PROGRAM FEM_TRI3
! 使用程序模块library和global_variables
USE library
USE global_variables
IMPLICIT NONE

! in_filename -- 输入文件名(最长15个字符)；out_filename --输出文件名(最长12个字符)
CHARACTER(LEN=15)::in_filename,out_filename

! x,y -- 结点坐标  ae -- 材料参数  sk -- 整体刚度矩阵 r -- 等效结点载荷列阵或结点位移列阵
! dv -- 非零已知位移数值
REAL(8),ALLOCATABLE::x(:),y(:),ae(:,:),dv(:,:),sk(:),r(:)

! meo -- 单元结点信息  jr -- 结点自由度序号  ma -- 指示矩阵  ndi -- 非零已知位移结点序号
INTEGER,ALLOCATABLE::meo(:,:),jr(:,:),ndi(:),ma(:)

! istat,ostat -- I/O 状态，0表示打开成功
INTEGER::istat,ostat
! 输入提示信息
in_filename='fem3dat.txt'
WRITE(*,'(A)') 'Please input file name of data=fem3dat.txt'
!READ(*,'(A15)') in_filename

! 输出提示信息
out_filename='output.txt'
WRITE(*,'(A)') 'Please output file name of data=output.txt'
!READ(*,'(A15)') out_filename

OPEN(UNIT=5,FILE=in_filename,STATUS='OLD',ACTION='READ',IOSTAT=istat)
OPEN(UNIT=7,FILE=out_filename,STATUS='REPLACE',ACTION='WRITE',IOSTAT=ostat)
! 输入文件是否打开成功？
IF(istat == 0) THEN
	! 输入9个控制参数
	! np -- 结点总数  
	! ne -- 单元总数  
	! nm -- 材料类型总数 
	! nr -- 约束结点总数
	! ndp -- 非零已知位移结点总数
	! ni -- 问题类型标识，0为平面应力问题，1为平面应变问题
	! nl -- 受集中力的结点数目  
	! ng -- 考虑自重作用为1，不计自重为0
	! nc -- 计算支座反力结点的数目

	READ (5,*) NP,NE,NM,NR,NDP,NI,NL,NG,NC 
	WRITE (*,"(/1X,9(A,I3,2X))") '结点总数=',NP,'单元总数=',NE,'材料类型总数=',NM,'约束结点总数=',NR,'非零已知位移结点总数=',NDP,'0为平面应力问题，1为平面应变问题=',NI,'受集中力的结点数目=',NL,'考虑自重作用为1，不计自重为0=',NG,'计算支座反力结点的数目=',NC                          
	WRITE (7,"(/1X,9(A,I3,2X))") 'NP=',NP,'NE=',NE,'NM=',NM,'NR=',NR,'NDP=',NDP,'NI=',NI,'NL=',NL,'NG=',NG,'NC=',NC 
	! 为数组分配存储空间
	ALLOCATE(x(np),y(np),meo(4,ne),ae(4,nm),jr(2,np))
	! 调用input子程序输入结点坐标，单元信息和材料参数
	CALL input(x,y,meo,ae)
	!调用MR子程序形成结点自由序号矩阵
	CALL mr(jr)
	!为数组MA分配存储空间        
	ALLOCATE(ma(n))
	!调用FORMMA子程序形成指示矩阵MA(N)并调用其他功能子程序
	CALL formma(meo,jr,ma)
	!为劲度矩阵数组SK和等效结点荷载列阵数组R分配存储空间
	ALLOCATE(sk(nh),r(n))
	!调用子程序MGK，形成整体劲度矩阵，并按一维存贮在SK中
	CALL mgk(ae,x,y,meo,jr,ma,sk)
	!调用子程序LOAD，形成整体等效结点荷载列阵 
	CALL load(ae,x,y,meo,jr,r)
	WRITE(*,1020)
	WRITE(7,1020)
	!输出整体等效结点荷载
	CALL output(jr,r)
	!调用TREAT子程序，输入非零已知位移信息，处理荷载项及对应K中主元素
	IF(ndp .GT. 0) THEN
		ALLOCATE(ndi(ndp),dv(2,ndp))
		CALL treat(sk,ma,r,jr,ndi,dv)
	ENDIF
	!整体劲度矩阵的分解运算
	CALL decomp(sk,ma)
	!前代、回代求出未知结点位移
	CALL foba(sk,ma,r)
	WRITE(*,1030)
	WRITE(7,1030)
	!调用子程序OUTPUT输出结点位移
	CALL output(jr,r)
	WRITE(*,1040)
	WRITE(7,1040)
	!调用子程序CES输出单元应力 
	CALL ces(ae,x,y,meo,jr,r)
	!调用子程序ERFAC输出支座反力
	IF(nc>0) CALL erfac(ae,x,y,meo,jr,r)
ELSE
    WRITE(*,1010) in_filename,istat
    WRITE(7,1010) in_filename,istat
ENDIF
1010 FORMAT(' 错误的打开文件',A,': IOSTAT = ',I6)
1020 FORMAT(30x,'Nodal Forces'/8x,'Node',11x,'X-Comp.',14x,'Y-Comp.')
1030 FORMAT(/30x,'Nodal Displacements'/3x,'Node',13x,'X-Comp.',12x,'Y-Comp.')
1040 FORMAT(/30x,'Element Stresses'/5x,'Element',5x,'X_Stress',3x,'Y_Stress',2x,'XY_Stress',1x,'Max_Stress',1x,'Min_Stress',6x,'Angle'/)

END PROGRAM

