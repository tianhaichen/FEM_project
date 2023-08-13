!******************************************************************************
!     1���������ܼ��㵯����ѧ�е�ƽ��Ӧ�������ƽ��Ӧ�����⣻
!     2����������������ε�Ԫ��
!     3���ܿ������غͽ�㼯�������ֺ��ص����ã��ڼ�������ʱy��ȡ��ֱ����Ϊ����
!******************************************************************************
PROGRAM FEM_TRI3
! ʹ�ó���ģ��library��global_variables
USE library
USE global_variables
IMPLICIT NONE

! in_filename -- �����ļ���(�15���ַ�)��out_filename --����ļ���(�12���ַ�)
CHARACTER(LEN=15)::in_filename,out_filename

! x,y -- �������  ae -- ���ϲ���  sk -- ����նȾ��� r -- ��Ч����غ��������λ������
! dv -- ������֪λ����ֵ
REAL(8),ALLOCATABLE::x(:),y(:),ae(:,:),dv(:,:),sk(:),r(:)

! meo -- ��Ԫ�����Ϣ  jr -- ������ɶ����  ma -- ָʾ����  ndi -- ������֪λ�ƽ�����
INTEGER,ALLOCATABLE::meo(:,:),jr(:,:),ndi(:),ma(:)

! istat,ostat -- I/O ״̬��0��ʾ�򿪳ɹ�
INTEGER::istat,ostat
! ������ʾ��Ϣ
in_filename='fem3dat.txt'
WRITE(*,'(A)') 'Please input file name of data=fem3dat.txt'
!READ(*,'(A15)') in_filename

! �����ʾ��Ϣ
out_filename='output.txt'
WRITE(*,'(A)') 'Please output file name of data=output.txt'
!READ(*,'(A15)') out_filename

OPEN(UNIT=5,FILE=in_filename,STATUS='OLD',ACTION='READ',IOSTAT=istat)
OPEN(UNIT=7,FILE=out_filename,STATUS='REPLACE',ACTION='WRITE',IOSTAT=ostat)
! �����ļ��Ƿ�򿪳ɹ���
IF(istat == 0) THEN
	! ����9�����Ʋ���
	! np -- �������  
	! ne -- ��Ԫ����  
	! nm -- ������������ 
	! nr -- Լ���������
	! ndp -- ������֪λ�ƽ������
	! ni -- �������ͱ�ʶ��0Ϊƽ��Ӧ�����⣬1Ϊƽ��Ӧ������
	! nl -- �ܼ������Ľ����Ŀ  
	! ng -- ������������Ϊ1����������Ϊ0
	! nc -- ����֧������������Ŀ

	READ (5,*) NP,NE,NM,NR,NDP,NI,NL,NG,NC 
	WRITE (*,"(/1X,9(A,I3,2X))") '�������=',NP,'��Ԫ����=',NE,'������������=',NM,'Լ���������=',NR,'������֪λ�ƽ������=',NDP,'0Ϊƽ��Ӧ�����⣬1Ϊƽ��Ӧ������=',NI,'�ܼ������Ľ����Ŀ=',NL,'������������Ϊ1����������Ϊ0=',NG,'����֧������������Ŀ=',NC                          
	WRITE (7,"(/1X,9(A,I3,2X))") 'NP=',NP,'NE=',NE,'NM=',NM,'NR=',NR,'NDP=',NDP,'NI=',NI,'NL=',NL,'NG=',NG,'NC=',NC 
	! Ϊ�������洢�ռ�
	ALLOCATE(x(np),y(np),meo(4,ne),ae(4,nm),jr(2,np))
	! ����input�ӳ������������꣬��Ԫ��Ϣ�Ͳ��ϲ���
	CALL input(x,y,meo,ae)
	!����MR�ӳ����γɽ��������ž���
	CALL mr(jr)
	!Ϊ����MA����洢�ռ�        
	ALLOCATE(ma(n))
	!����FORMMA�ӳ����γ�ָʾ����MA(N)���������������ӳ���
	CALL formma(meo,jr,ma)
	!Ϊ���Ⱦ�������SK�͵�Ч��������������R����洢�ռ�
	ALLOCATE(sk(nh),r(n))
	!�����ӳ���MGK���γ����徢�Ⱦ��󣬲���һά������SK��
	CALL mgk(ae,x,y,meo,jr,ma,sk)
	!�����ӳ���LOAD���γ������Ч���������� 
	CALL load(ae,x,y,meo,jr,r)
	WRITE(*,1020)
	WRITE(7,1020)
	!��������Ч������
	CALL output(jr,r)
	!����TREAT�ӳ������������֪λ����Ϣ������������ӦK����Ԫ��
	IF(ndp .GT. 0) THEN
		ALLOCATE(ndi(ndp),dv(2,ndp))
		CALL treat(sk,ma,r,jr,ndi,dv)
	ENDIF
	!���徢�Ⱦ���ķֽ�����
	CALL decomp(sk,ma)
	!ǰ�����ش����δ֪���λ��
	CALL foba(sk,ma,r)
	WRITE(*,1030)
	WRITE(7,1030)
	!�����ӳ���OUTPUT������λ��
	CALL output(jr,r)
	WRITE(*,1040)
	WRITE(7,1040)
	!�����ӳ���CES�����ԪӦ�� 
	CALL ces(ae,x,y,meo,jr,r)
	!�����ӳ���ERFAC���֧������
	IF(nc>0) CALL erfac(ae,x,y,meo,jr,r)
ELSE
    WRITE(*,1010) in_filename,istat
    WRITE(7,1010) in_filename,istat
ENDIF
1010 FORMAT(' ����Ĵ��ļ�',A,': IOSTAT = ',I6)
1020 FORMAT(30x,'Nodal Forces'/8x,'Node',11x,'X-Comp.',14x,'Y-Comp.')
1030 FORMAT(/30x,'Nodal Displacements'/3x,'Node',13x,'X-Comp.',12x,'Y-Comp.')
1040 FORMAT(/30x,'Element Stresses'/5x,'Element',5x,'X_Stress',3x,'Y_Stress',2x,'XY_Stress',1x,'Max_Stress',1x,'Min_Stress',6x,'Angle'/)

END PROGRAM

