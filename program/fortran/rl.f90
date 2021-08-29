PROGRAM RegresionLineal

!==========================
!Programita para hacer una regesion lineal y el calculo
!de las incertidumbres
!VCastor
!Programa hecho en 2016 cuando recien iniciaba en la programacion
!esta siendo actualizado en 2021
!==========================

Integer npuntos
real m,b,sx,sx2,sy,sy2,sxy,x,y,r,symxb2,ub,um

WRITE(*,*) "The data should be in a file with the name:"
WRITE(*,*) "data.dat"

OPEN(8,FILE="data.dat")
READ(8,*) npuntos

dimension x(npuntos),y(npuntos)

do i=1,npuntos
 read(8,*)x(i),y(i)
end do

!==========================
!PROCEDIMIENTO, sumas
!==========================
sx=0.0
sy=0.0
sx2=0.0
sy2=0.0
sxy=0.0
do i=1,npuntos
 sx=sx+x(i)
 sy=sy+y(i)
 sx2=sx2+x(i)*x(i)
 sy2=sy2+y(i)*y(i)
 sxy=sxy+x(i)*y(i)
end do

!===========
!PENDIENTE
!============

m=(npuntos*sxy-sx*sy)/(npuntos*sx2-sx*sx)
write(*,*)"m=",m

!===================================
!PROCEDIMIENTO ORDENADA DE ORIGEN
!===================================

b=((sy*sx2)-(sx*sxy))/(npuntos*sx2-(sx*sx))
write(*,*)"b=",b

!===========================
!COEFICIENTE DE CORELACIÓN
!===========================

r=(sxy-(sx*sy/npuntos))/(SQRT((sx2-(sx*sx/npuntos))*(sy2-(sy*sy/npuntos))))
write(*,*)"r=",r
write(*,*)"r^2=",r*r

!========================
!INCERTIDUMBRES
!========================

symxb2=0.0
do i=1,npuntos
 symxb2=symxb2+(y(i)-m*x(i)-b)*(y(i)-m*x(i)-b)
end do

um=SQRT((npuntos*symxb2)/((npuntos-2)*(npuntos*sx2-sx*sx)))
ub=SQRT((sx2*symxb2)/((npuntos-2)*(npuntos*sx2-sx*sx)))

write(*,*)"U(m)=",um
write(*,*)"U(b)=",ub

ENDPROGRAM RegresionLineal
