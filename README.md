# Actividad-1-6-Multiprocesadores

Se muestran 2 programas que resuelvan 4 ecuaciones diferenciales de forma numérica usando el método de Runge-Kutta de orden 4:

y′=te3t−2y
y′=1+(t−y)2
y′=1+y/t
y′=cos(2ty)+sen(3ty)

En el intervalo [0,π] (con 50, mil puntos en el intervalo de evaluación). Para todas las ecuaciones la condición inicial es y(0)=π/4. 

La ejecución brinda  un archivo por cada una de las soluciones determinadas en el programa.
Estos archivos tambien se adjuntan si no se desea correr el programa

Uno de los programas está hecho de forma secuencial y el otro en forma paralela usando task level parallelism (TLP). 

También se imprime en consola el tiempo de ejecución de ambos programas. Los cuales tambien se muestran a continuación:

Secuencial:
0.871 Segundos

Paralelo:
0.144 Segundos
