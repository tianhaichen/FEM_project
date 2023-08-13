!******************************************************************************
!     1、本程序采用四结点等参单元计算弹性力学平面问题；
!     2、能考虑在自重，结点集中荷载，以及分布面荷载（水压或法向分布压力）的作用；
!     3、在计算自重时y轴取垂直向上为正；
!     4、能处理非零已知位移的边界条件；
!     5、采用压缩稀疏行存储格式存储劲度矩阵；
!     6、采用预条件共轭梯度法求解有限元系统方程组；
!     7、主要输出内容包括：结点位移、单元应力、结点应力、主应力以及主应力与x轴的夹角。
!******************************************************************************
PROGRAM FEM_TRI3
! 使用程序模块library和global_variables
USE library
USE global_variables
!USE csr
!USE pcg
IMPLICIT NONE

! in_filename -- 输入文件名(最长15个字符)；out_filename --输出文件名(最长12个字符)
CHARACTER(LEN=15)::in_filename,out_filename
INTEGER,ALLOCATABLE::jr(:,:),ndi(:),nn(:),meo(:,:),med(:,:),elcount(:),elnd(:,:),ndcon(:),ndptr(:)
REAL(8),ALLOCATABLE::coor(:,:),ae(:,:),wg(:,:),r(:),dv(:,:),sk(:),snod(:,:)
INTEGER::ie,nby

! istat,ostat -- I/O 状态，0表示打开成功
INTEGER::istat,ostat
! 输入提示信息
in_filename='fem4dat.txt'
WRITE(*,'(A)') 'Please input file name of data=fem4dat.txt'
!READ(*,'(A15)') in_filename

! 输出提示信息
out_filename='output.txt'
WRITE(*,'(A)') 'Please output file name of data=output.txt'
!READ(*,'(A15)') out_filename

OPEN(UNIT=15,FILE=in_filename,STATUS='OLD',ACTION='READ',IOSTAT=istat)
OPEN(UNIT=16,FILE=out_filename,STATUS='REPLACE',ACTION='WRITE',IOSTAT=ostat)
! 输入文件是否打开成功？
IF(istat == 0) THEN
!调用CONTROL子程序得到控制变量的具体数据
CALL control
!为数组分配存储空间
ALLOCATE(jr(2,np),coor(2,np),meo(9,ne),med(8,ne),ae(4,nm),wg(2,nw),nn(8),snod(6,np))
!输入结点坐标、单元信息、材料特性常数和水荷载特性常数等原始数据，处理平面应力问题的E和v
CALL input(coor,meo,ae,wg)
ELSE
    WRITE(*,1010) in_filename,istat
    WRITE(16,1010) in_filename,istat
ENDIF
1010 FORMAT(' 错误的打开文件',A,': IOSTAT = ',I6)

END PROGRAM

