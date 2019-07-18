# TwoGPLA
Implementación del método Particle-in-Cell en dos dimensiones, con la dinámica de los iones en el sistema.
Los códigos tienen que estar en la misma carpeta.
En el código TwoGPLAmain hay una bandera llamada DINAMICA, con la cual se puede hace 
una grafica dinámica mientras se ejecuta el código
En el código TwoGPLAseed se asignan todos los parámetros necesarios para la semilla.
El orden de ejecución es: 
1-TwoGPLAseed
2-TwoGPLAmain
3-TwoGPLAplots
#Nota: El código TwoGPLAfunctions no se ejecuta directamente, se llama de los otros programas.

Se necesitan acompañar los códigos con las siguientes carpetas:
-Campo
-Densidad
-Energías (con tilde)
-EspacioDeFase
-FuncionDistribucion
-Graficas
-Posiciones
-Potencial
-Semilla
-Velocidades

En ellas se almacena la información.
#Nota: Para el caso de los datos en:
Potencial, Energías, Densidad y Campo se puede dejar 
solamente el último archivo porque en este se contiene 
la información de todos los anteriores.
