# CEesferico1D
En este repositorio se encuentran el código CEesferico1D 


# Requisitos
Para ejecutar los programas de este repositorio se necesita un compilador de Fortran. Se recomienda usar GNU Fortran

Un compilador de Fortran compatible con Fortran 2008, por ejemplo GNU Fortran (gfortran). En Windows, una opción práctica es MSYS2/MinGW-w64.

Python 3 (recomendado 3.10 o superior).

Paquetes de Python:  py -m pip install numpy matplotlib pillow

Pillow es necesario para que ecuacion_de_onda.py guarde la animación como GIF.

# Ejecución del código (Carpeta Code)
También es posible descargarlo con Code → Download ZIP, descomprimirlo y abrir una terminal dentro de la carpeta resultante.
Los archivos que se encuentran en la carpeta CODE son los siguientes:

| Archivo | Función |
|---|---|
| `CEesferico1D.f90` | Este es el programa principal para la evolución esférica en una dimensión. Para correr el programa primero se ejecuta `gfortran .\evolucion.f90 .\CEesferico1D.f90 -o .\CEesferico1D.exe` y el código se compila con la instrucción `.\CEesferico1D.exe 640 64 20 0.05 ` donde `640` es el número de nodos de la malla, `64` es el radio máximo, `20` es el número de pasos temporales y `0.05` es la amplitud inicial del perfil gaussiano. |
| `biseccion.f90` | Este programa calcula mejores aproximaciones de las soluciones tanto críticas como subcrítticas muy cerca del umbral mediante el método de bisección. Se ejecuta primero la intrucción `gfortran .\biseccion_phi0.f90 -o .\biseccion_phi0.exe` y se compila con `.\biseccion_phi0.exe 640 64 20 0.3529 0.3530 1e-13 80` donde `0.35290` es el valor subcrítico calculado manualmente más cercano al umbral, `0.3530` es el valor supercrítico calculado tmabién manualmente más cercnao al umbral, `1e-13` es el orden de nuestra aproximación en el cual queremos que el programa se detenga y `80` es el número de iteraciones que debe hacer el programa antes de detenerse.|
| `evolucion.f90` | Este programa contiene las subrutinas del programa principal |
| `gamma_rho.f90` | Cálculos necesarios para obtener la ley de potencia poder medio de `gamma` y la densidad de materia/energía `rho`. |
| `valores_centrales.dat` | Archivo de datos de entrada con valores centrales para phi0_real y `gamma` y `rho` para obtener la gráfica log-log|

Algunos programas pueden requerir que valores_centrales.dat permanezca en la misma carpeta que el ejecutable.

Los programas pueden generar archivos de datos con extensión .dat. Estos contienen los resultados numéricos de la simulación o cálculo y pueden utilizarse posteriormente para hacer gráficas o análisis con los programas de la carpeta Plot.

# Scripts de Python para visualización (Carpeta Plot)
Estos scripts procesan los archivos .dat producidos por las simulaciones Fortran. Generan gráficas, paneles o animaciones en formato GIF o PNG.

| Archivo | Función | Datos requeridos |
|---|---|---|
| `animacion_CEesferico1D.py` | Genera animaciones, series temporales o paneles de las variables de la evolución esférica. | `CEesferico1D_*.dat` |
| `animacion_masa.py` | Genera animaciones o paneles de la función de masa de Misner–Sharp y la relación `2m/r`. | `mass_*.dat` |
| `alpha_critica.py` / `graficar_alpha_critica.py` | Compara la evolución de `alpha` para las soluciones subcrítica y supercrítica cercanas al valor crítico. | `biseccion.csv`, ejecutable de evolución y archivos `CEesferico1D_*.dat` |
| `gamma_rho.py` / `graficar_gamma_rho.py` | Ajusta y grafica el exponente crítico `γ` a partir de datos centrales. | `valores_centrales.dat` |


Para obtener las animaciones de cualquier variable del código principal `animacion_CEesferico1D.py` se utiliza la siguiente intrucción:

`Python .\animacion_CEesferico1D.py --mode profile --variable scalar --output animacion_scalar.gif`

`Python .\animacion_CEesferico1D.py --mode profile --variable alpha --output animacion_alpha.gif`

`Python .\animacion_CEesferico1D.py --mode profile --variable a --output animacion_a.gif`

Además, se puede ajustar la escala de las animaciones con una intrucción como esta:

`Python .\animacion_CEesferico1D.py --mode profile --variable scalar --xmax 5 --ymin -0.1 --ymax 0.5 --interval 80 --output animacion_scalar.gif`

Esta instrucción de ajuste de escala se puede modificar para cualquier variable a observar.

Para la función de masa, se compila el programa  `animacion_masa.py` con la isntrucción: 

`python .\animacion_masa.py`

Para obtener la gráfica de `alpha central` se obtiene con una isntrucción como la siguiente:

`python .\graficar_alpha_critica.py --nr 640 --rmax 64.0 --tfinal 20.0 --output alpha_central_critica.png`


El script `gamma_rho.py` ajusta la relación:

ln(ρ_c^max) = C - 2γ ln(φ* - φ₀)

donde φ* es el valor crítico 

El archivo valores_centrales.dat debe contener los valores de φ₀ y la densidad central máxima ρ_c^max.

Por lo que para ejecutar este programa se debe proporcionar el valor crítico phi_c con la siguiente instrucción

python graficar_gamma_rho.py 0.35294658892608627
