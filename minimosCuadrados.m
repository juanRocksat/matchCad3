disp("-----------------------------------------------------");
abcisas = input("ingrese abcisas \n");
ordenadas = input ("ingrese ordenadas \n");
%%abcisas = [ 1 2 3 4];
%%ordenadas=[1 2 1.2 2];
cantAbcisas = length(abcisas);
hold on ;
sprintf("cantidad de puntos ingresados es : %f",cantAbcisas)
sx1=0;
sx2=0;
sx3=0;
sx4=0;
sy1=0;
sy2=0;
sy3=0;
sy4=0;
sfxPorx=0;
sfxPorx2=0;
mediaX=0;
mediaY=0;
sLnDeLaImagen=0;
sLnDex=0;
sLnDeLaImagenPorLaImagen=0;
sLnDexAl2=0;
sLnDeLaImagenPorLnDex=0;
sDexiSobreyi = 0;
sDeInversaDeyi=0;
cantDePuntos=length(abcisas);

for i=1:cantDePuntos
  xi=abcisas(i);
  yi=ordenadas(i);
    sx1=sx1+xi;
    sy1=sy1+yi;
    sx2=sx2+(xi^2);
    sx3=sx3+(xi^3);
    sx4=sx4+(xi^4);
    sfxPorx=sfxPorx+xi*yi;
    sfxPorx2=sfxPorx2+xi*xi*yi;
    sLnDeLaImagen=sLnDeLaImagen+log(yi);
    sLnDeLaImagenPorLaImagen=sLnDeLaImagenPorLaImagen+log(yi)*yi;
    sLnDeLaImagenPorLnDex=sLnDeLaImagenPorLnDex+log(yi)*log(xi);
    sLnDex = sLnDex+log(xi);
    sLnDexAl2 = sLnDexAl2+(log(xi)*log(xi));
    sDexiSobreyi = sDexiSobreyi+ (xi/yi);
    sDeInversaDeyi = sDeInversaDeyi +(1/yi);
endfor
matrizNormalParaLaParabola = [sx4 sx3 sx2 ;
                              sx1 sx2 sx3  ;
                              sx2 sx3  cantDePuntos];
matrizNormalParaLaExponencial=[sx2 sx1;
                              sx1 cantDePuntos]; 
matrizNormalParaLaRecta = [sx2 sx1;
                           sx1 cantDePuntos];
matrizNormalParaLaPotencial = [sLnDexAl2 sLnDex;
                                sLnDex cantDePuntos];
matrizNormalParaLaHomografica = [sx2 sx1;
                                  sx1 cantDePuntos];

terminosIndependientesDeLaParabola = [sfxPorx2;sfxPorx;sx1] ; %matriz de un vector columna
terminosIndependientesDeLaExponencial = [sLnDeLaImagenPorLaImagen;sLnDeLaImagen];
terminosIndependientesDeRecta = [sfxPorx;sy1];
terminosIndependientesDeLaHomografica = [sDexiSobreyi;sDeInversaDeyi];
terminosIndependientesDeLaPotencial = [sLnDeLaImagenPorLnDex;sLnDeLaImagen];
coeficientesDeLaParabola  = inv(matrizNormalParaLaParabola)*terminosIndependientesDeLaParabola;%[a;b;c] para ax^2+bx+c
coeficientesDeLaExponencial = inv(matrizNormalParaLaExponencial)*terminosIndependientesDeLaExponencial; % [a;ln(b)]
coeficientesDeLaRecta = inv(matrizNormalParaLaRecta)*terminosIndependientesDeRecta;%[a;b]
coeficientesDeLaPotencial = inv(matrizNormalParaLaPotencial)*terminosIndependientesDeLaPotencial;%[a,ln(b)]
coeficientesDeLaHomografica = inv(matrizNormalParaLaHomografica)*terminosIndependientesDeLaHomografica; %[1/a;b/a]
parabola_a = coeficientesDeLaParabola(1,1);
parabola_b = coeficientesDeLaParabola(2,1);
parabola_c = coeficientesDeLaParabola(3,1);
exponencial_a = coeficientesDeLaExponencial(1,1);
exponencial_b  = e^(coeficientesDeLaExponencial(2,1));                    
recta_a = coeficientesDeLaRecta(1,1);
recta_b = coeficientesDeLaRecta(2,1);
potencial_a = coeficientesDeLaPotencial(1,1);
potencial_b = e.^coeficientesDeLaPotencial(2,1);      
homografica_a = 1/coeficientesDeLaHomografica(1,1);
homografica_b = homografica_a*coeficientesDeLaHomografica(2,1);
x = (abcisas(1,1):0.1:abcisas(1,cantDePuntos));
x_2=x;x_3=x_2;x_4=x_3;x_5=x_4;
parabola_polinomio = parabola_a*(x.^2)+(parabola_b*x)+parabola_c;
exponencial_funcion = exponencial_b*(e.^(exponencial_a*x));
recta_funcion=recta_a*x+recta_b;
potencial_funcion = potencial_b*(x.^(potencial_a));
homografica_funcion = homografica_a./(x.+homografica_b);
%plot(abcisas,ordenadas);
%%minimosCuadrados = figure
%clearplot;
axis;
plot(x,parabola_polinomio,'r','linewidth',1);
plot(x,exponencial_funcion,'g','linewidth',1);% g de green
plot(x,recta_funcion,'k','linewidth',2);%negro
plot(x,potencial_funcion,'m','linewidth',1);%magenta,, c cian
plot(x,homografica_funcion,'b','linewidth',2);%blue
title("Aproximaciones ");
xlabel("<----------Abcisas------------------>");
ylabel("<------------ordenadas------------->");
legend([num2str(parabola_a),'x^2 +',num2str(parabola_b),'x+',num2str(parabola_c)],
      [num2str(exponencial_b),'e**(',num2str(exponencial_a),'x)'],
      [num2str(recta_a),'x+',num2str(recta_b)],
      [num2str(potencial_b),'x**',num2str(potencial_a)],
      [num2str(homografica_a),'/','(x+(',num2str(homografica_b),'))']);
print('minimos cuadrados.png','-dpng');
%---------- errores
disp("---------by Jhon Daniel ,Olmedo Paco ");
errorDeRecta = 0;
errorDePotencial = 0;
errorDeExponencial  = 0;
errorDeParabola=0;
errorDeHomografica=0;
for i=1:cantDePuntos
  xi=abcisas(i);
  yi=ordenadas(i);
  errorDeExponencial=errorDeExponencial+((exponencial_funcion(xi)-yi).^2);
  errorDeParabola=errorDeParabola+((parabola_polinomio(xi)-yi).^2);
  errorDeRecta=errorDeRecta+((recta_funcion(xi)-yi).^2);
  errorDeHomografica=errorDeHomografica+((homografica_funcion(xi)-yi).^2);
  errorDePotencial=errorDePotencial+((potencial_funcion(xi)-yi).^2);
endfor
tabla_de_errores = ['error de recta :        ' num2str(errorDeRecta);
                    'error de parabola:      ' num2str(errorDeParabola);
                    'error de la potencial:  ' num2str(errorDePotencial);
                    'error de la homografica:' num2str(errorDeHomografica);
                    'error de la exponencial:' num2str(errorDeExponencial) ]
%closeplot;
disp("-----------------------------------------------------");
