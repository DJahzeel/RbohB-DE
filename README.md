# Efecto del gen RbohB en la expresión diferencial de raíces de frijol

Fecha de inicio: 31/08/2025

**Participantes**:

- Pedro Daniel Pineda Martinez
  <email:[pedropm@lcg.unam.mx]>
  
- Dara Jahzeel Palafox Collado
  <email:[darapc@lcg.unam.mx]>

## Planteamiento del problema

La interacción simbiótica en Phaseolus vulgaris con Rhizobium tropici (nodulación) y Rhizophagus irregularis (micorrización) es crucial para la adquisición de nitrógeno y fósforo, respectivamente. Sin embargo, el gen PvRbohB, que regula la producción de especies reactivas de oxígeno (ROS), ejerce un control diferencial y aparentemente contradictorio sobre estas simbiosis: su silenciamiento (PvRbohB-RNAi) compromete severamente la nodulación, pero promueve la colonización micorrízica.

Ante esta dualidad, surge la siguiente pregunta de investigación: ¿Cómo el silenciamiento del gen PvRbohB en Phaseolus vulgaris modula diferencialmente las redes de interacción genética y los procesos fisiológicos subyacentes, como la homeostasis de ROS, la remodelación de la pared celular y la señalización de fitohormonas, para impactar de manera opuesta la nodulación con Rhizobium tropici y la colonización micorrízica con Rhizophagus irregularis en las etapas tempranas de la simbiosis?

## Calendario de trabajo

[Definir de manera general la actividades que se requerirán para el proyecto. Por ejemplo:]

| Actividad | Fecha   | Responsable  | Entregable |
|----------|----------|----------|----------|
| Descripción de proyecto    | 31-08-25  | Daniel Pineda  | Documento markdown |
| Especificación de requisitos    | Septiembre   | Por definir  | Documento markdown   |
| Análisis y diseño   | Septiembre  | Por definir  | Documento markdown |
| Construcción   | Octubre-Noviembre |  Por definir   | Scripts |
| Pruebas   | Noviembre  | Por definir    | Documento markdown |
| Reporte de resultados  | Noviembre  |  Por definir   | Documentos markdown |
| Presentación del proyecto   | Diciembre  |  Por definir  | repositorio GitHub (release)|



## Metodología
[Descripción general de los pasos a realizar para el proyecto, por ejemplo:]

**Preguntas de investigación**



Paso 1: Localización de fuente de datos  
Paso 2: Descarga de archivos de datos  
Paso 3: Inspección de datos  
Paso 4: Limpieza de datos
Paso 5: Descripción de los datos  
Paso 5: Análsis de datos depurados  
Paso 6: Obtención de resultados  



## Resultados esperados

No se que se espera.







## Especificación de Requisitos

Requisitos funcionales

[Un requisito funcional define una función del sotware. Determina que debe hacer el software: Por ejemplo: Leer números de un archivo dado, Calcular la suma de todos los números leídos del archivo, Producir un mensaje de error si el archivo no existe, etc.]


Requisitos no funcionales

[Estos requisitos se ocupan de aspectos como la seguridad, el rendimiento, la facilidad de uso, la fiabilidad y la escalabilidad. Por ejemplo: El script deberá estar escrito en Python, El tiempo de respuesta debe ser rápido, incluso con archivos de gran tamaño, La entrada del archivo debe ser flexible (i.e. se acepta a través de la línea de comandos), etc.]




## Análisis y Diseño



Para resolver este problema, se utilizarán varias funciones incorporadas en Python, así como el manejo de excepciones para la validación de datos y archivo. A continuación, se muestra un pseudocódigo simple para ilustrar la lógica básica del script:

```
Función principal (Suma_Numero):
    Intentar:
        datos_archivo = Obtener_Datos_Archivo(ruta_archivo)
        numeros = Validar_Datos(datos_archivo)
        resultado = Calcular_Suma(numeros)
        Imprimir_Resultado(numeros, resultado)
    Atrapar cualquier excepción como error:
        Imprimir el error

Función Obtener_Datos_Archivo(ruta_archivo):
    Si la ruta del archivo no existe:
        Levantar un error de "archivo no encontrado"
    Leer y retornar las líneas del archivo

Función Validar_Datos(data):
    Intentar:
        Convertir todos los datos a formato flotante y retornar como una lista
    Atrapar ValorError:
        Levanta un error de "¡Introducir solamente números de base 10!"

Función Calcular_Suma(numeros):
    Retornar la suma de los números

Función Imprimir_Resultado(numeros, resultado):
    Imprimir los números y su suma total en el formato especificado
```

El formato de los datos de entrada será simplemente un archivo, con un número por línea. Los números pueden estar en formato entero o decimal. La salida será una línea de texto que muestra los números sumados y la suma total, en el formato: n1 + n2 + n3 + ... = suma. Los mensajes de error se imprimirán en la consola.


#### Caso de uso: Sumar Números

```
         +---------------+
         |   Usuario     |
         +-------+-------+
                 |
                 | 1. Proporciona archivo de entrada
                 v
         +-------+-------+
         |   Sumador de  |
         |   Números en  |
         |   Archivo     |
         | (Sistema)     |
         +---------------+
```

- **Actor**: Usuario
- **Descripción**: El actor proporciona un archivo de entrada con números a sumar. El sistema valida el archivo y los datos de entrada, calcula la suma de los números y muestra el resultado.
- **Flujo principal**:

	1. El actor inicia el sistema proporcionando el archivo de entrada con los números a sumar.
	2. El sistema valida el archivo y los datos de entrada.
	3. El sistema calcula la suma de los números.
	4. El sistema muestra el resultado.
	
- **Flujos alternativos**:
	- Si el archivo proporcionado no existe
		1. El sistema muestra un mensaje de error diciendo que el archivo no se encuentra.
	- Si los datos de entrada no son números en base 10
		1. El sistema muestra un mensaje de error diciendo que se deben introducir números en base 10.
                
