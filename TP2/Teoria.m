%Si tengo dos procesos independientes, las potencias son aditivas.
%Ejemplo de clase con un proceso continuo y uno complejo.

%Z(t)=x(t)+y(t)
%PSD(Z(t))=PSD(x(t)) + PSD(y(t))
%var(z)=var(x)+var(y)

%RUIDO CIRCULARMENTE SIMETRICOS
%Si se lo multiplica por cualquier exponencial compleja, y se lo 
%desplaza en frecuencia, queda siempre igual


%RUIDO REAL -> PSD(1side)=2*PSD(2side)
%Para ruido complejo no existe la PSD(1side) porque no podemos
%asegurar que sea par

%Hacer el promedio de la fft es un estimador de la FFT asi como la media
%es un estimador de la esperanza. Este proceso da mucho ruido.
%Una forma de mejorarlo es hacerlo en un for e ir promediando las pxx (Unicamente
%para ruidos blacos ).

%Cuando los ruidos son reales se suele decir que la PSD es, por ejemplo
%10MW/Hz

%% Ruidos Filtrados

% Para un ruido blanco la autocorrelacion es practicamente cero, al pasar por un filtro
% se le da a cada muestra una cierta relacion con las anteriores por la
% convolucion. Por lo tanto la salida ya no es blanca.

% La salida depende de varias entradas pasadas.

%%

% Ry(z) = Rx(z)*h(z)*h*(-z) con h(z) rta. al impulso del filtro

% Sy(f) = Sx(f).H(f).H*(f) con H(f) espectro del filtro 

%% Esperanza de la salida del filtro

%E{y(t)}=E{int(-inf;inf)(h(u)x(t-u)du)}
%Metemos la Esperanza dentro de la integral
%E(y(t))= int(-inf;inf)(E(h(u)x(t-u)du))
%Distribuyo la esperanza
%E(y(t))=int(-inf;inf)(E(h(u).E(x(t-u).E(du))))
%La integral de E(hu) es h(u) por ser un proceso deterministico
%es una constante.
%La integral de E(du) es du.
%La integral de E(x(t-u)) es la u(x)=media

%La int(-inf;inf )(h(u)du)=int(-inf;inf)(h(u)e(-jwu)) con w=0
%eso es igual a H(O)= ganancia en continua


%% Filtro anti-aliasing

% Es MUY IMPORTANTE tener un buen filtro anti-aliasing en mi circuito
% cuando hay ruido blanco.