!変わらない条件
    
module mod_condition
    implicit none
    
    !温度
    real(8),parameter::temp=50.d0+273.15d0 ![℃→K]
    !温度(F)
    real(8),parameter::Temp_F=1.8d0*temp-459.67d0
    !温度(℃)
    real(8),parameter::temp_C=50.0d0

    !初期圧力
    real(8),parameter::iniPressure=10.0d0*10.0d0**6.0d0 ![MPa→Pa]
    
    !円周率
    real(8),parameter::pi=4.0d0*atan(1.0d0)
    
    !気体定数
    real(8),parameter::R = 8.3144598d0
    
    !H2Oの物質量
    real(8),parameter::M1=0.018015d0 !kg/mol!!
    
    end module