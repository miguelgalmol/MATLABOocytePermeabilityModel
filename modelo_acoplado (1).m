%*******************************************************************************
% Programa que permite calcular las dos etapas, primero la de deshidrataci�n 
% y luego la rehidrataci�n con respecto al tiempo al a�adir crioprotector.
%*******************************************************************************
r=63;%Radio del oocito[micras]
inactive=0.19;%Fraccion volumen celula osmoticamente inactivo, 81% agua
Vcelli=4/3*pi*r^3*1e-9;%volumen total inicial de oocito [mm3]
Vi=Vcelli*(1-inactive);%volumen inicial agua osmoticamente activo en oocito [mm3]
A=4*pi*r^2*1e-6;%area membrana celular oocito [mm2]

T=298.15;%temperatura a la que se a�ade el crioprotector [K]

MA=18.01528;%peso molecular del agua pura [g/mol]
VmA=18.61*1e-3;%volumen molar del agua pura a 25�C [litros/mol]
roA=0.9970479; %Densidad del agua pura a 25�C [Kg/litro]=[mg/mm3]


CM2=0.154;%concentracion isotonica inicial de ClNa intracelular[mol/litro] 
v2=1.6836;%coeficiente de disociacion del ClNa a la concentracion isotonica
Vm2=16.68*1e-3;%volumen molar parcial a dilucion infinita de ClNa a 25�C[litros/mol]
n2=CM2*Vi*1e-6;%moles iniciales de ClNa en el medio intracelular
nsol=n2*v2;%osmoles iniciales de de ClNa en el medio intracelular

Vm_sucrose=21.5412*1e-3;%volumen molar sacarosa [litros/osmol]=[cm3/mosmol]---PENDIENTE BUSCAR DATO REFERENCIADO---
Vm_EG=55.37*1e-3;%volumen molar parcial a dilucion infinita de EG a 25�C[litros/osmol]=[cm3/mosmol]
Vm_DMSO=64.323*1e-3; %volumen molar parcial a dilucion infinita de DMSO a 25�C[litros/osmol]=[cm3/mosmol]

%**************************************************************************
%DATOS DE PERMEABILIDAD DE LA MEMBRANA CELULAR A CPAs Y AGUA
%**************************************************************************
%Permeabilidad celular al agua---------------
Lpo=0.69;%Permeabilidad hidr�ulica del agua a 25�C [um/min atm]
Ea_w=14.42;%Energia de activacion de la permeabilidad del agua [kcal/mol]
To_w=298.15;%Temperatura de referencia de la permeabilidad hidra�lica del agua [K]
Lp=Lpo*exp((-Ea_w*1e3/1.987207)*(1/T-1/To_w))*9.869232667;%Permeabilidad de la membrana al agua en presencia de EG [mm4/(min J)]
%Permeabilidad celular al CPA DMSO---------------
Po_DMSO=15;%Coeficiente de permeabilidad del DMSO a 24�C[um/min] 
Ea_DMSO=23.52;%Energia de activacion de la permeabilidad del DMSO [kcal/mol]
To_DMSO=297.15;%Temperatura de referencia de la permeabilidad DMSO [K]
P_DMSO=Po_DMSO*exp((-Ea_DMSO*1e3/1.987207)*(1/T-1/To_DMSO))*1e-3;%Permeabilidad de la membrana al DMSO [mm/min]
%Permeabilidad celular al CPA EG---------------
Po_EG=9.16;%Coeficiente de permeabilidad del EG a 25�C[um/min]
Ea_EG=21.20;%Energia de activacion de la permeabilidad del EG [kcal/mol]
To_EG=298.15;%Temperatura de referencia de la permeabilidad EG [K]
P_EG=Po_EG*exp((-Ea_EG*1e3/1.987207)*(1/T-1/To_EG))*1e-3;%Permeabilidad de la membrana al EG [mm/min]

%##########################################################################
%##########################################################################
%ETAPA 0 ETAPA 0 ETAPA 0 ETAPA 0 ETAPA 0 ETAPA 0 ETAPA 0 ETAPA 0 ETAPA 0 
%##########################################################################
%##########################################################################

%**************************************************************************
%CONCENTRACIONES DE CRIOPROTECTORES Y SACAROSA DEL MEDIO EXTRACELULAR
%**************************************************************************
CM_DMSO=1.056/2;%concentraci�n molar de DMSO [osmol/litro]=[moslmol/cm3]
CM_EG=1.345/2;%concentraci�n molar de Etilenglicol[osmol/litro]=[mosmol/cm3]
CM_sucrose=0;%concentraci�n molar de sacarosa [osmol/litro]=[mosmol/cm3]
%**************************************************************************
%CALCULO MOLALIDADES Y FRACCIONES MOLARES DE LOS SOLUTOS EN EL MEDIO EXTRACELULAR
%**************************************************************************
nA=(1-CM2*Vm2-CM_DMSO*Vm_DMSO-CM_EG*Vm_EG-CM_sucrose*Vm_sucrose)/VmA;%moles de agua en un medio extracelular de un litro [moles/litro]
masaA=nA*MA*1e-3;%[Kg/litro]
%Osmolalidades en el medio extracelular----------------------
mB=v2*CM2/masaA;
m_DMSO=CM_DMSO/masaA;
m_EG=CM_EG/masaA;
m_sucrose=CM_sucrose/masaA;
Smex=mB+m_DMSO+m_EG+m_sucrose;%suma de las osmolalidades de los solutos del medio extracelular[osmol/kg]

%Fracciones molares en el medio extracelular---------------
ntot=CM2*v2+CM_DMSO+CM_EG+CM_sucrose+nA;%moles totales por litro de solucion extracelular[moles/litro]
XB=CM2/ntot;
X_DMSO=CM_DMSO/ntot;
X_EG=CM_EG/ntot;
X_sucrose=CM_sucrose/ntot;
XA=1-XB-X_DMSO-X_EG-X_sucrose;

%**************************************************************************
%DEFINICION DEL SISTEMA DE ECUACIONES DIFERENCIALES DE DESHIDRATACION 
%**************************************************************************
%Las ecuaciones se escriben en un vector columna.
%la primera fila corresponde a la ec. de flujo de agua. El volumen de agua
%intracelular [mm3] es x(1)
%La segunda fila corresponde a la ec. de flujo de moles de CPA DMSO. Los
%mosmoles de CPA DMSO intracelular es x(2)
%La tercera fila corresponde a la ec. de flujo de moles de CPA EG. Los
%mosmoles de CPA intracelular[mol] es x(3)
R=8.314475;%[J/molK]
f=@(t,x) [-Lp*A*R*1e-3*T*(roA*Smex*1e-3-(nsol*1e3+x(2)+x(3))/x(1));P_DMSO*A*(roA*m_DMSO*1e-3-x(2)/x(1));P_EG*A*(roA*m_EG*1e-3-x(3)/x(1))];
%se asume una densidad del agua de 1 litro=1Kg
%el factor 1e-3 que multiplica a R es para pasar los moles a mmoles
%el factor 1e-3 que afecta a Smex es para pasar de mol/litro a mosmol/mm3
%el factor 1e3 que afecta a nsol es para pasar de osmol a mosmol
%el factor 1e-3 que afecta a m_DMSO es para pasar de mol/litro a mosmol/mm3
%el factor 1e-3 que afecta a m_EG es para pasar de mol/litro a mosmol/mm3

%**************************************************************************
%INTEGRACION DEL SISTEMA DE ECUACIONES DIFERENCIALES QUE DEFINEN LOS FLUJOS
%DE AGUA [MM3] Y DE MOLES DE CPA [MOSMOL]
%**************************************************************************
tspan=[0:.01:1];%intervalo de tiempo en minutos
xo=zeros(1,3);
xo(1)=Vi;%El volumen de agua intracelular inicial es Vi[mm3]
xo(2)=0;%Los moles iniciales de CPA DMSO son cero dentro de la celula.
xo(3)=0;%Los moles iniciales de CPA EG son cero dentro de la celula.
[t,x]=ode45(f,tspan,xo);
%Resultados numericos en una grafica---------------
Vn=[];
nfull=[];
Vn=(Vcelli-Vi+x(:,1)+x(:,2)*Vm_DMSO*1e3+x(:,3)*Vm_EG*1e3+nsol*1e3*Vm2*1e3)/Vcelli;%Fraccion del volumen celular total con respecto al inicial
nfull=x(:,2)+x(:,3)+nsol*1e3;%mosmoles intracelulares totales:DMSO+EG+ClNa
flux_w=Lp*A*R*1e-3*T*(Smex*1e-3-(nsol*1e3+x(:,2)+x(:,3))./x(:,1));%caudal de agua instantaneo a traves de la membrana [mm3/min]
flux_EG=P_EG*A*(m_EG*1e-3-x(:,3)./x(:,1));%caudal de mosmoles de Etilenglicol instantaneo a traves de la membrana [mosmol/min]
flux_DMSO=P_DMSO*A*(m_DMSO*1e-3-x(:,2)./x(:,1));%caudal de mosmoles de DMSO instantaneo a traves de la membrana [mosmol/min]

disp('VolumenRel minimo(%)///Tiempo minimo(s)///Caudal Inicial deshidrat(mm3/min)///VolumenRel Final(%)')
[minVol posVol]=min(Vn);
[minFlux posFlux]=min(flux_w);
minTime=t(posVol)*60;
[minVol*100 minTime flux_w(1) Vn(end)*100]

%Representacion grafica de los %resultados---------------------------------
f1=figure;
plot(t,Vn,'b')
hold on;
plot(t(posVol),minVol,'*r')
title('Evolucion deshidratacion celular: Ratio del volumen geometrico celular respecto al volumen total inicial')
ylabel('Fraccion volumen celular geometrico[adim]')
xlabel('Tiempo[min]')
grid on;

s0.time=t
s0.volumen_celular_geometrico = Vn

f2=figure;
plot(t,flux_w,'b')
hold on;
title('Tasa deshidratacion celular: Caudal de agua que atraviesa la membrana celular')
ylabel('Caudal volumentrico [mm3/min]')
xlabel('Tiempo[min]')
grid on;

s0.caudal_volumentrico = flux_w

f3=figure;
plot(t,flux_DMSO,'b',t,flux_EG,'g')
hold on;
title('Tasa de entrada de crioprotector DMSO y EG al interior celular')
ylabel('Caudal molar [mosmol/min]')
xlabel('Tiempo[min]')
grid on;

s0.caudal_molar_DMSO = flux_DMSO
s0.caudal_molar_EG = flux_EG

f4=figure;
hold on;
plot(t,x(:,2)./x(:,1),'b',t,x(:,3)./x(:,1),'g',t,nsol*1e3./x(:,1),'r',t,nfull./x(:,1),'k',t,Smex*1e-3,'*m')
title('Osmolalidad de crioprotector DMSO,EG y sal en el interior celular')
ylabel('Osmolalidad [mosmol/mg]')
xlabel('Tiempo[min]')
grid on;

s0.osmolalidad_DMSO = x(:,2)./x(:,1)
s0.osmolalidad_EG = x(:,3)./x(:,1)
s0.osmolalidad_sal = nsol*1e3./x(:,1)

f5=figure;
plot(t,x(:,1)/Vi*100, 'r')
hold on;
title('Ratio del volumen de agua intracelular con respecto al agua inicialmente presente')
ylabel('Fraccion volumen agua intracelular [%]')
xlabel('Tiempo[min]')
grid on;

s0.fraccion_agua = x(:,1)/Vi*100

f6=figure;
plot(t,nsol*1e3./x(:,1),'r',t,mB*1e-3,'*k')
hold on;
title('Osmolalidad de ClNa en el interior celular')
ylabel('Osmolalidad [mosmol/mg]')
xlabel('Tiempo[min]')
grid on;

s0.osmolalidad_ClNa = nsol*1e3./x(:,1)

f7=figure;
plot(t,x(:,2)./x(:,1),'b',t,m_DMSO*1e-3,'*k')
hold on;
title('Osmolalidad de DMSO en el interior celular')
ylabel('Osmolalidad [mosmol/mg]')
xlabel('Tiempo[min]')
grid on;

s0.osmolalidad_DMSO_interior_celular = x(:,2)./x(:,1)

f8=figure;
plot(t,x(:,3)./x(:,1),'g',t,m_EG*1e-3,'*k')
hold on;
title('Osmolalidad de EG en el interior celular')
ylabel('Osmolalidad [mosmol/mg]')
xlabel('Tiempo[min]')
grid on;

s0.osmolalidad_EG_interior_celular = x(:,3)./x(:,1)

%##########################################################################
%##########################################################################
%ETAPA 1 ETAPA 1 ETAPA 1 ETAPA 1 ETAPA 1 ETAPA 1 ETAPA 1 ETAPA 1 ETAPA 1
%##########################################################################
%##########################################################################

%**************************************************************************
%CONCENTRACIONES DE CRIOPROTECTORES Y SACAROSA DEL MEDIO EXTRACELULAR
%**************************************************************************
CM_DMSO=1.056;%concentraci�n molar de DMSO [osmol/litro]=[moslmol/cm3]
CM_EG=1.345;%concentraci�n molar de Etilenglicol[osmol/litro]=[mosmol/cm3]
CM_sucrose=0;%concentraci�n molar de sacarosa [osmol/litro]=[mosmol/cm3]

%**************************************************************************
%CALCULO MOLALIDADES Y FRACCIONES MOLARES DE LOS SOLUTOS EN EL MEDIO EXTRACELULAR
%**************************************************************************
nA=(1-CM2*Vm2-CM_DMSO*Vm_DMSO-CM_EG*Vm_EG-CM_sucrose*Vm_sucrose)/VmA;%moles de agua en un medio extracelular de un litro de volumen
masaA=nA*MA*1e-3;%Kg de agua(disolvente) en un litro de disolucion
%Osmolalidades en el medio extracelular----------------------
mB=v2*CM2/masaA;
m_DMSO=CM_DMSO/masaA;
m_EG=CM_EG/masaA;
m_sucrose=CM_sucrose/masaA;
Smex=mB+m_DMSO+m_EG+m_sucrose;%suma de las osmolalidades de los solutos del medio extracelular[osmol/kg]

%Fracciones molares en el medio extracelular---------------
ntot=CM2*v2+CM_DMSO+CM_EG+CM_sucrose+nA;%moles totales por litro de solucion extracelular[moles/litro]
XB=CM2/ntot;
X_DMSO=CM_DMSO/ntot;
X_EG=CM_EG/ntot;
X_sucrose=CM_sucrose/ntot;
XA=1-XB-X_DMSO-X_EG-X_sucrose;

%**************************************************************************
%DEFINICION DEL SISTEMA DE ECUACIONES DIFERENCIALES DE DESHIDRATACION 
%**************************************************************************
%Las ecuaciones se escriben en un vector columna.
%la primera fila corresponde a la ec. de flujo de agua. El volumen de agua
%intracelular [mm3] es x(1)
%La segunda fila corresponde a la ec. de flujo de moles de CPA DMSO. Los
%mosmoles de CPA DMSO intracelular es x(2)
%La tercera fila corresponde a la ec. de flujo de moles de CPA EG. Los
%mosmoles de CPA intracelular es x(3)
R=8.314475;%[J/molK]
f=@(t,x) [-Lp*A*R*1e-3*T*(roA*Smex*1e-3-(nsol*1e3+x(2)+x(3))/x(1));P_DMSO*A*(roA*m_DMSO*1e-3-x(2)/x(1));P_EG*A*(roA*m_EG*1e-3-x(3)/x(1))];
% TO DO: FORMULA

%var1 = f(293, 3)
%se asume una densidad del agua 1 es decir, que 1 litro=1Kg de agua.
%el factor 1e-3 que multiplica a R es para pasar los moles a mmoles(estan
%en el denominador)
%el factor 1e-3 que afecta a Smex es para pasar de osmol/Kg a mosmol/mg
%el factor 1e3 que afecta a nsol es para pasar de osmol a mosmol
%el factor 1e-3 que afecta a m_DMSO es para pasar de osmol/Kg a mosmol/mg
%el factor 1e-3 que afecta a m_EG es para pasar de osmol/Kg a mosmol/mg

%**************************************************************************
%INTEGRACION DEL SISTEMA DE ECUACIONES DIFERENCIALES QUE DEFINEN LOS FLUJOS
%DE AGUA [MM3] Y DE MOLES DE CPA [MOSMOL]
%**************************************************************************
tspan=[1:.01:2];%intervalo de tiempo en minutos
xo=zeros(1,3);
xo(1)=x(end,1);%Volumen de agua intracelular inicial es el final del anterior paso[mm3]
xo(2)=x(end,2);%Los osmoles iniciales de CPA DMSO dentro de la celula son los osmoles finales del paso anterior.
xo(3)=x(end,3);%Los osmoles iniciales de CPA EG dentro de la celula son los osmoles finales del paso anterior.
[t,x]=ode45(f,tspan,xo); % integrar a formula f -> permeabilidade; t = tempo, x = osmolaridade
%Resultados numericos-----------------------------------------------------
Vn=[];%inicializamos el vector de resultados con un vector nulo.
nfull=[];
Vn=(Vcelli-Vi+x(:,1)+x(:,2)*Vm_DMSO*1e3+x(:,3)*Vm_EG*1e3+nsol*1e3*Vm2*1e3)/Vcelli;%Fraccion del volumen celular total con respecto al inicial
nfull=x(:,2)+x(:,3)+nsol*1e3;%mosmoles intracelulares totales:DMSO+EG+ClNa
flux_w=Lp*A*R*1e-3*T*(Smex*1e-3-(nsol*1e3+x(:,2)+x(:,3))./x(:,1));%caudal (fluxo) de agua instantaneo a traves de la membrana [mm3/min]
flux_EG=P_EG*A*(m_EG*1e-3-x(:,3)./x(:,1));%caudal de mosmoles de Etilenglicol instantaneo a traves de la membrana [mosmol/min]
flux_DMSO=P_DMSO*A*(m_DMSO*1e-3-x(:,2)./x(:,1));%caudal de mosmoles de DMSO instantaneo a traves de la membrana [mosmol/min]

disp('VolumenRel minimo(%)///Tiempo(segundos)///Caudal inicial agua(mm3/min)///VolumenRel final(%)')
[minVol posVol]=min(Vn);
[minFlux posFlux]=min(flux_w);
minTime=t(posVol)*60;
[minVol*100 minTime flux_w(1) Vn(end)*100]


%Representacion grafica de los %resultados---------------------------------
figure(f1);
plot(t,Vn,'b')
hold on;
plot(t(posVol),minVol,'*r')
title('Evolucion deshidratacion celular: Ratio del volumen geometrico celular respecto al volumen total inicial')
ylabel('Fraccion volumen celular geometrico[adim]')
xlabel('Tiempo[min]')
grid on;

s1.time=t
s1.volumen_celular_geometrico = Vn

figure(f2);
plot(t,flux_w,'b')
hold on;
title('Tasa deshidratacion celular: Caudal de agua que atraviesa la membrana celular')
ylabel('Caudal volumentrico [mm3/min]')
xlabel('Tiempo[min]')
grid on;

s1.caudal_volumentrico = flux_w

figure(f3);
plot(t,flux_DMSO,'b',t,flux_EG,'g')
hold on;
title('Tasa de entrada de crioprotector DMSO y EG al interior celular')
ylabel('Caudal molar [mosmol/min]')
xlabel('Tiempo[min]')
grid on;

s1.caudal_molar_DMSO = flux_DMSO
s1.caudal_molar_EG = flux_EG

figure(f4);
plot(t,x(:,2)./x(:,1),'b',t,x(:,3)./x(:,1),'g',t,nsol*1e3./x(:,1),'r',t,nfull./x(:,1),'k',t,Smex*1e-3,'*m')
hold on;
title('Osmolalidad de crioprotector DMSO,EG y sal en el interior celular')
ylabel('Osmolalidad [mosmol/mg]')
xlabel('Tiempo[min]')
grid on;

s1.osmolalidad_DMSO = x(:,2)./x(:,1)
s1.osmolalidad_EG = x(:,3)./x(:,1)
s1.osmolalidad_sal = nsol*1e3./x(:,1)

figure(f5);
plot(t,x(:,1)/Vi*100, 'r')
hold on;
title('Ratio del volumen de agua intracelular con respecto al agua inicialmente presente')
ylabel('Fraccion volumen agua intracelular [%]')
xlabel('Tiempo[min]')
grid on;

s1.fraccion_agua = x(:,1)/Vi*100

figure(f6);%concentracion de sal intracelular
hold on;
plot(t,nsol*1e3./x(:,1),'r',t,mB*1e-3,'*k')
title('Osmolalidad de ClNa en el interior celular')
ylabel('Osmolalidad [mosmol/mg]')
xlabel('Tiempo[min]')
grid on;

s1.osmolalidad_ClNa = nsol*1e3./x(:,1)

figure(f7);%concentracion de DMSO intracelular
hold on;
plot(t,x(:,2)./x(:,1),'b',t,m_DMSO*1e-3,'*k')
title('Osmolalidad de DMSO en el interior celular')
ylabel('Osmolalidad [mosmol/mg]')
xlabel('Tiempo[min]')
grid on;

s1.osmolalidad_DMSO_interior_celular = x(:,2)./x(:,1)

figure(f8);%concentracion de EG intracelular
hold on;
plot(t,x(:,3)./x(:,1),'g',t,m_EG*1e-3,'*k')
title('Osmolalidad de EG en el interior celular')
ylabel('Osmolalidad [mosmol/mg]')
xlabel('Tiempo[min]')
grid on;

s1.osmolalidad_EG_interior_celular = x(:,3)./x(:,1)

%##########################################################################
%##########################################################################
%ETAPA 2 ETAPA 2 ETAPA 2 ETAPA 2 ETAPA 2 ETAPA 2 ETAPA 2 ETAPA 2 ETAPA 2
%##########################################################################
%##########################################################################

%**************************************************************************
%CONCENTRACIONES DE CRIOPROTECTORES Y SACAROSA DEL MEDIO EXTRACELULAR
%**************************************************************************
CM_DMSO=2.112;%concentraci�n molar de DMSO [osmol/litro]=[moslmol/cm3]
CM_EG=2.69;%concentraci�n molar de Etilenglicol[osmol/litro]=[mosmol/cm3]
CM_sucrose=0.5;%MIGUEL:CAMBIO LA SUCROSA DE 0 A 0,5%concentraci�n molar de sacarosa [osmol/litro]=[mosmol/cm3]

%**************************************************************************
%CALCULO MOLALIDADES Y FRACCIONES MOLARES DE LOS SOLUTOS EN EL MEDIO EXTRACELULAR
%**************************************************************************
nA=(1-CM2*Vm2-CM_DMSO*Vm_DMSO-CM_EG*Vm_EG-CM_sucrose*Vm_sucrose)/VmA;%moles de agua en un medio extracelular de un litro [moles/litro]
masaA=nA*MA*1e-3;%[Kg/litro]
%Osmolalidades en el medio extracelular----------------------
mB=v2*CM2/masaA;
m_DMSO=CM_DMSO/masaA;
m_EG=CM_EG/masaA;
m_sucrose=CM_sucrose/masaA;
Smex=mB+m_DMSO+m_EG+m_sucrose;%suma de las osmolalidades de los solutos del medio extracelular[osmol/kg]

%Fracciones molares en el medio extracelular---------------
ntot=CM2*v2+CM_DMSO+CM_EG+CM_sucrose+nA;%moles totales por litro de solucion extracelular[moles/litro]
XB=CM2/ntot;
X_DMSO=CM_DMSO/ntot;
X_EG=CM_EG/ntot;
X_sucrose=CM_sucrose/ntot;
XA=1-XB-X_DMSO-X_EG-X_sucrose;

%**************************************************************************
%DEFINICION DEL SISTEMA DE ECUACIONES DIFERENCIALES DE DESHIDRATACION 
%**************************************************************************
%Las ecuaciones se escriben en un vector columna.
%la primera fila corresponde a la ec. de flujo de agua. El volumen de agua
%intracelular [mm3] es x(1)
%La segunda fila corresponde a la ec. de flujo de moles de CPA DMSO. Los
%mosmoles de CPA DMSO intracelular es x(2)
%La tercera fila corresponde a la ec. de flujo de moles de CPA EG. Los
%mosmoles de CPA intracelular[mol] es x(3)
R=8.314475;%[J/molK]
f=@(t,x) [-Lp*A*R*1e-3*T*(roA*Smex*1e-3-(nsol*1e3+x(2)+x(3))/x(1));P_DMSO*A*(roA*m_DMSO*1e-3-x(2)/x(1));P_EG*A*(roA*m_EG*1e-3-x(3)/x(1))];
%se asume una densidad del agua de 1 litro=1Kg
%el factor 1e-3 que multiplica a R es para pasar los moles a mmoles
%el factor 1e-3 que afecta a Smex es para pasar de mol/litro a mosmol/mm3
%el factor 1e3 que afecta a nsol es para pasar de osmol a mosmol
%el factor 1e-3 que afecta a m_DMSO es para pasar de mol/litro a mosmol/mm3
%el factor 1e-3 que afecta a m_EG es para pasar de mol/litro a mosmol/mm3

%**************************************************************************
%INTEGRACION DEL SISTEMA DE ECUACIONES DIFERENCIALES QUE DEFINEN LOS FLUJOS
%DE AGUA [MM3] Y DE MOLES DE CPA [MOSMOL]
%**************************************************************************
tspan=[2:.01:3];%intervalo de tiempo en minutos
xo=zeros(1,3);
xo(1)=x(end,1);%Volumen de agua intracelular inicial es el final del anterior paso[mm3]
xo(2)=x(end,2);%Los osmoles iniciales de CPA DMSO dentro de la celula son los osmoles finales del paso anterior.
xo(3)=x(end,3);%Los osmoles iniciales de CPA EG dentro de la celula son los osmoles finales del paso anterior.
[t,x]=ode45(f,tspan,xo);
%Resultados numericos en una grafica---------------
Vn=[];
nfull=[];
Vn=(Vcelli-Vi+x(:,1)+x(:,2)*Vm_DMSO*1e3+x(:,3)*Vm_EG*1e3+nsol*1e3*Vm2*1e3)/Vcelli;%Fraccion del volumen celular total con respecto al inicial
nfull=x(:,2)+x(:,3)+nsol*1e3;%mosmoles intracelulares totales:DMSO+EG+ClNa
flux_w=Lp*A*R*1e-3*T*(Smex*1e-3-(nsol*1e3+x(:,2)+x(:,3))./x(:,1));%caudal de agua instantaneo a traves de la membrana [mm3/min]
flux_EG=P_EG*A*(m_EG*1e-3-x(:,3)./x(:,1));%caudal de mosmoles de Etilenglicol instantaneo a traves de la membrana [mosmol/min]
flux_DMSO=P_DMSO*A*(m_DMSO*1e-3-x(:,2)./x(:,1));%caudal de mosmoles de DMSO instantaneo a traves de la membrana [mosmol/min]

disp('VolumenRel minimo(%)///Tiempo minimo(s)///Caudal Inicial deshidrat(mm3/min)///VolumenRel Final(%)')
[minVol posVol]=min(Vn);
[minFlux posFlux]=min(flux_w);
minTime=t(posVol)*60;
[minVol*100 minTime flux_w(1) Vn(end)*100]

%Representacion grafica de los %resultados---------------------------------
figure(f1);
plot(t,Vn,'b')
plot(t(posVol),minVol,'*r')
title('Evolucion deshidratacion celular: Ratio del volumen geometrico celular respecto al volumen total inicial')
ylabel('Fraccion volumen celular geometrico[adim]')
xlabel('Tiempo[min]')
grid on;
saveas(gcf,'f1.png')

s2.time=t
s2.volumen_celular_geometrico = Vn

figure(f2);
plot(t,flux_w,'b')
title('Tasa deshidratacion celular: Caudal de agua que atraviesa la membrana celular')
ylabel('Caudal volumentrico [mm3/min]')
xlabel('Tiempo[min]')
grid on;
saveas(gcf,'f2.png')

s2.caudal_volumentrico = flux_w

figure(f3);
plot(t,flux_DMSO,'b',t,flux_EG,'g')
title('Tasa de entrada de crioprotector DMSO y EG al interior celular')
ylabel('Caudal molar [mosmol/min]')
xlabel('Tiempo[min]')
lgd = legend('Fluxo de DMSO', 'Fluxo de EG');
grid on;
saveas(gcf,'f3.png')

s2.caudal_molar_DMSO = flux_DMSO
s2.caudal_molar_EG = flux_EG

figure(f4);
plot(t,x(:,2)./x(:,1),'b',t,x(:,3)./x(:,1),'g',t,nsol*1e3./x(:,1),'r',t,nfull./x(:,1),'k',t,Smex*1e-3,'*m')
title('Osmolalidad de crioprotector DMSO,EG y sal en el interior celular')
ylabel('Osmolalidad [mosmol/mg]')
xlabel('Tiempo[min]')
grid on;
saveas(gcf,'f4.png')

s2.osmolalidad_DMSO = x(:,2)./x(:,1)
s2.osmolalidad_EG = x(:,3)./x(:,1)
s2.osmolalidad_sal = nsol*1e3./x(:,1)

figure(f5);
plot(t,x(:,1)/Vi*100, 'r')
title('Ratio del volumen de agua intracelular con respecto al agua inicialmente presente')
ylabel('Fraccion volumen agua intracelular [%]')
xlabel('Tiempo[min]')
grid on;
saveas(gcf,'f5.png')

s2.fraccion_agua = x(:,1)/Vi*100

figure(f6);%concentracion de sal intracelular
plot(t,nsol*1e3./x(:,1),'r',t,mB*1e-3,'*k')
title('Osmolalidad de ClNa en el interior celular')
ylabel('Osmolalidad [mosmol/mg]')
xlabel('Tiempo[min]')
grid on;
saveas(gcf,'f6.png')

s2.osmolalidad_ClNa = nsol*1e3./x(:,1)

figure(f7);%concentracion de DMSO intracelular
plot(t,x(:,2)./x(:,1),'b',t,m_DMSO*1e-3,'*k')
title('Osmolalidad de DMSO en el interior celular')
ylabel('Osmolalidad [mosmol/mg]')
xlabel('Tiempo[min]')
grid on;
saveas(gcf,'f7.png')

s2.osmolalidad_DMSO_interior_celular = x(:,2)./x(:,1)

figure(f8);%concentracion de EG intracelular
plot(t,x(:,3)./x(:,1),'g',t,m_EG*1e-3,'*k')
title('Osmolalidad de EG en el interior celular')
ylabel('Osmolalidad [mosmol/mg]')
xlabel('Tiempo[min]')
grid on;
saveas(gcf,'f8.png')

s2.osmolalidad_EG_interior_celular = x(:,3)./x(:,1)

%%% Putting all structs together

Alldata = [s0, s1, s2];
s_all = struct;  %final structure
for field = fieldnames(Alldata)'
   fname = field{1};
   s_all.(fname) = vertcat(Alldata.(fname));
end

% octave não possui struct2table
% writetable(struct2table(s_all), 'all_data.xlsx')

% Diferença de / e ./
% / - Divisão de matriz
% ./ - Divisão por elemento (element-wise)

% x(:,3) -> coletar todas as linhas da coluna 3
% x(2,:) -> coletar todas as colunas na linha 2

% [y,x] = min(a) -> retorna o valor minimo da matriz a na variável y e o
% índice do valor mínimo de a na variável x

