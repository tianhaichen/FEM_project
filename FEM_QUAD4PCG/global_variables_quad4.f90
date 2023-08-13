MODULE global_variables
! 本程序模块用于声明程序全局变量
IMPLICIT NONE
INTEGER:: NP,NE,NM,NR,NW,NF,NI,NDP,NESW
! N:结构自由度总数  NH:按CSR存贮的整体劲度矩阵的总容量
INTEGER:: N,NH 
! RSTG(2): 高斯积分点坐标  
REAL(8)::RSTG(2)
END MODULE