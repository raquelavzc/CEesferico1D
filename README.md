# CEesferico1D
En este repositorio se encuentran el código CEesferico1D 


# Requisitos
Para ejecutar los programas de este repositorio se necesita un compilador de Fortran. Se recomienda usar GNU Fortran

Un compilador de Fortran compatible con Fortran 2008, por ejemplo GNU Fortran (gfortran). En Windows, una opción práctica es MSYS2/MinGW-w64.

Python 3 (recomendado 3.10 o superior).

Paquetes de Python:  py -m pip install numpy matplotlib pillow

Pillow es necesario para que ecuacion_de_onda.py guarde la animación como GIF.

# Ejecución del código
También es posible descargarlo con Code → Download ZIP, descomprimirlo y abrir una terminal dentro de la carpeta resultante.
Los archivos que se encuentran en la carpeta CODE son los siguientes:

| Archivo | Función |
|---|---|
| `CEesferico1D.f90` | Este es el programa principal para la evolución esférica en una dimensión. Para correr el programa primero se ejecuta `gfortran .\evolucion.f90 .\CEesferico1D.f90 -o .\CEesferico1D.exe` y el código se compila con la instrucción `.\CEesferico1D.exe 640 64 20 0.05 ` donde `640` es el número de nodos de la malla, `64` es el radio máximo, `20` es el número de pasos temporales y `0.05` es la amplitud inicial del perfil gaussiano. |
| `biseccion.f90` | Este programa calcula mejores aproximaciones de las soluciones tanto críticas como subcrítticas muy cerca del umbral mediante el método de bisección. Se ejecuta primero la intrucción `gfortran .\biseccion_phi0.f90 -o .\biseccion_phi0.exe` y se compila con `.\biseccion_phi0.exe 640 64 20 0.3529 0.3530 1e-13 80` donde `0.35290` es el valor subcrítico calculado manualmente más cercano al umbral, `0.3530` es el valor supercrítico calculado tmabién manualmente más cercnao al umbral, `1e-13` es el orden de nuestra aproximación en el cual queremos que el programa se detenga y `80` es el número de iteraciones que debe hacer el programa antes de detenerse.|
| `evolucion.f90` | Rutinas o programa para evolución temporal. |
| `gamma_rho.f90` | Cálculos necesarios para obtener la ley de potencia poder medio de `gamma` y `rho`. |
| `valores_centrales.dat` | Archivo de datos de entrada con valores centrales para phi0_real y `gamma` y `rho` para obtener la gráfica log-log|
Algunos programas pueden requerir que valores_centrales.dat permanezca en la misma carpeta que el ejecutable.
