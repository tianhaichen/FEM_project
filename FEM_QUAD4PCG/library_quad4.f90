MODULE library
CONTAINS 
! ***************************************************************
!     ���ӳ������Ҫ���ܣ�����9�����Ʊ���
SUBROUTINE control 
	USE global_variables 
	IMPLICIT NONE 
	!��1�У����������� NP ����Ԫ���� NE ������������Ŀ NM ��Լ��������� NR ��ˮ���������� NW ��
	!�ܼ��к��ؽ�����Ŀ NF ���������ʱ�ʶ NI (ƽ��Ӧ����0��ƽ��Ӧ����1)��������֪λ�ƽ������ NDP
	!�ֲܷ�������ˮѹ���ĵ�Ԫ��Ŀ NESW
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
!     ���ӳ������Ҫ���ܣ�����ԭʼ���ݣ��γ�COOR�����飬����ƽ��Ӧ�������E��v
SUBROUTINE input(coor,meo,ae,wg)
	USE global_variables
	IMPLICIT NONE 
	REAL(8),INTENT(INOUT)::coor(:,:),ae(:,:),wg(:,:)
	INTEGER,INTENT(INOUT)::meo(:,:)
	REAL(8)::xy(2),ir(2),e,u
	INTEGER::ie,i,j,nn,l,nee 
	!��9�У��Խ��ѭ���������NN���õ��x,y����ֵ�����γɽ����������COOR 
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
	!��4�У����뵥Ԫ��Ϣ����NE�� 
	!���뵥Ԫ���NEE���õ�Ԫ�Ĳ������ͺ�NME���Ƿ����������ϢNET(1-���㣬0-������)��
	!��Ԫ��4������
	DO ie=1,ne 
		READ(15,*) nee,(meo(j,nee),j=1,6)
		WRITE(16,850) nee,(meo(j,nee),j=1,6)
		WRITE(*,850) nee,(meo(j,nee),j=1,6)
	ENDDO
	!��2�У��Բ���������ѭ��������������Գ���������ģ�������ɱȡ������������Ⱥ͵�Ԫ���
	!���е���ģ�������ɱ�Ҫ����ƽ��Ӧ����ֵ
	READ(15,*) ((ae(i,j),i=1,4),j=1,nm)
	WRITE(16,910) (j,(ae(i,j),i=1,4),j=1,nm)
	WRITE(*,910) (j,(ae(i,j),i=1,4),j=1,nm)
	!��8�У�����ƽ��Ӧ�����⣬E����E/(1.0-v**2)��v����v/(1.0-v)
	IF(ni/=0) THEN 
		DO j=1,nm
			u=ae(1,j)
			u=ae(2,j)
			ae(1,j)=e/(1.0-u*u)
			ae(2,j)=u/(1.0-u)
		ENDDO 
	ENDIF 
	!��4�У�����ˮ�������ã�������ˮ������������WG(2,NW)��ÿ��ˮ������Ϣ����Ϊˮλy����ֵ��ˮ������
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