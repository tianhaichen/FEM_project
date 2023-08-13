!******************************************************************************
!     1������������Ľ��Ȳε�Ԫ���㵯����ѧƽ�����⣻
!     2���ܿ��������أ���㼯�к��أ��Լ��ֲ�����أ�ˮѹ����ֲ�ѹ���������ã�
!     3���ڼ�������ʱy��ȡ��ֱ����Ϊ����
!     4���ܴ��������֪λ�Ƶı߽�������
!     5������ѹ��ϡ���д洢��ʽ�洢���Ⱦ���
!     6������Ԥ���������ݶȷ��������Ԫϵͳ�����飻
!     7����Ҫ������ݰ��������λ�ơ���ԪӦ�������Ӧ������Ӧ���Լ���Ӧ����x��ļнǡ�
!******************************************************************************
PROGRAM FEM_TRI3
! ʹ�ó���ģ��library��global_variables
USE library
USE global_variables
!USE csr
!USE pcg
IMPLICIT NONE

! in_filename -- �����ļ���(�15���ַ�)��out_filename --����ļ���(�12���ַ�)
CHARACTER(LEN=15)::in_filename,out_filename
INTEGER,ALLOCATABLE::jr(:,:),ndi(:),nn(:),meo(:,:),med(:,:),elcount(:),elnd(:,:),ndcon(:),ndptr(:)
REAL(8),ALLOCATABLE::coor(:,:),ae(:,:),wg(:,:),r(:),dv(:,:),sk(:),snod(:,:)
INTEGER::ie,nby

! istat,ostat -- I/O ״̬��0��ʾ�򿪳ɹ�
INTEGER::istat,ostat
! ������ʾ��Ϣ
in_filename='fem4dat.txt'
WRITE(*,'(A)') 'Please input file name of data=fem4dat.txt'
!READ(*,'(A15)') in_filename

! �����ʾ��Ϣ
out_filename='output.txt'
WRITE(*,'(A)') 'Please output file name of data=output.txt'
!READ(*,'(A15)') out_filename

OPEN(UNIT=15,FILE=in_filename,STATUS='OLD',ACTION='READ',IOSTAT=istat)
OPEN(UNIT=16,FILE=out_filename,STATUS='REPLACE',ACTION='WRITE',IOSTAT=ostat)
! �����ļ��Ƿ�򿪳ɹ���
IF(istat == 0) THEN
!����CONTROL�ӳ���õ����Ʊ����ľ�������
CALL control
!Ϊ�������洢�ռ�
ALLOCATE(jr(2,np),coor(2,np),meo(9,ne),med(8,ne),ae(4,nm),wg(2,nw),nn(8),snod(6,np))
!���������ꡢ��Ԫ��Ϣ���������Գ�����ˮ�������Գ�����ԭʼ���ݣ�����ƽ��Ӧ�������E��v
CALL input(coor,meo,ae,wg)
ELSE
    WRITE(*,1010) in_filename,istat
    WRITE(16,1010) in_filename,istat
ENDIF
1010 FORMAT(' ����Ĵ��ļ�',A,': IOSTAT = ',I6)

END PROGRAM

