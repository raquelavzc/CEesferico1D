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
| `CEesferico1D.f90` | Programa principal para la evolución esférica en una dimensión. |
| `biseccion.f90` | Calcula soluciones mediante el método de bisección para hallar el valor crítico. |
| `evolucion.f90` | Rutinas o programa para evolución temporal. |
| `gamma_rho.f90` | Cálculos necesarios para obtener la ley de potencia poder medio de `gamma` y `rho`. |
| `valores_centrales.dat` | Archivo de datos de entrada con valores centrales para phi0_real y `gamma` y `rho` para obtener la gráfica log-log|
