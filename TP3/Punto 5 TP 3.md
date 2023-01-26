# Punto 5 TP 3

## ROC

1. Dividir el disparo en celdas de ancho $\Delta R = \Delta t = \tau = (NOS\ muestras)$

2. Se toma la muestra del medio de la celda y se descarta el resto de las muestras.
   
   Donde la salida del match filter la tomo en intervalos de NOS muestras:
   
   $$
   Y_mf\left(1:NOS:\text{end}\right)
   $$
   
   Y obtengo las muestras del medio con:
   
   $$
   Y_mf\left(1+\frac{NOS}{2}:NOS:\text{end}\right)
   $$
   
   Si no llega a funcionar manualmente podemos hacer:
   
   $$
   Y_mf\left(1+\frac{NOS}{2}+offset:NOS:\text{end}\right)
   $$

3. Con esta nueva señal armamos la ROC. Ponemos un umbral $U_0$ y vemos si cada muestra supera o no el umbral, anotando TP, FN, FP, TN según corresponda.

4. Subimos el umbral y realizamos el punto 3 nuevamente.

5. Se cambia el rango y se realizan los puntos anteriores nuevamente, con el fin de calcular $P_D = \frac{TP}{COI*N_{exp}}$  y $P_{FA} = \frac{FP}{NO\ COI*N_{exp}}$para cada umbral.

6. Con esto se obtiene una curva ROC para cada rango donde se grafica la $PD$ en fuinción de la $PFA$.
























