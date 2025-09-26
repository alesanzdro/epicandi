# **Metodología Sistemática para la Caracterización Exhaustiva de Lecturas Largas de Nanopore No Clasificadas**

## **Parte I: Un Pipeline Estratégico para la Interrogación de Lecturas Largas No Clasificadas**

Esta sección inicial establece el flujo de trabajo fundamental para el procesamiento de datos. La filosofía subyacente es progresar desde evaluaciones amplias y rápidas hacia análisis profundos y computacionalmente intensivos, creando un embudo que enriquece progresivamente el conjunto de datos con lecturas que requieren una inspección detallada. El pipeline se presenta como un script modular en Bash, promoviendo la claridad y facilitando la depuración, de acuerdo con la solicitud.

### **1.1. Introducción: El Contenedor "Unclassified" como Herramienta de Diagnóstico**

El problema central no debe ser visto simplemente como una pérdida de datos, sino como una oportunidad para obtener información diagnóstica valiosa sobre la calidad de la secuenciación y la preparación de la librería. El contenedor de lecturas no clasificadas ("unclassified") es un archivo rico en pistas que pueden desvelar problemas sistémicos en el protocolo experimental.  
La observación de que entre el 10-20% del total de lecturas no se clasifican, pero esta cifra se dispara al 80-90% para lecturas de más de 75,000 nucleótidos, es una pista fundamental. Esta desproporción no es un evento aleatorio; sugiere un problema sistemático directamente relacionado con las propiedades físicas de los fragmentos de ADN de gran longitud o con los procesos que los generan. Las causas potenciales incluyen:

1. **Degradación o baja calidad en los extremos:** Las secuencias de los códigos de barras (barcodes), situadas en los extremos de los fragmentos, pueden ser más susceptibles a la degradación física o a una menor calidad de secuenciación en moléculas extremadamente largas.  
2. **Generación de quimeras pre-barcoding:** Los protocolos que minimizan la fragmentación para obtener lecturas ultra-largas podrían, inadvertidamente, ser más propensos a generar moléculas quiméricas (fusiones de fragmentos de ADN distintos) antes de la ligación de los adaptadores con códigos de barras.  
3. **Precisión del basecaller:** El software de basecalling (p. ej., Guppy o Dorado) podría tener una precisión reducida en la región del código de barras para estas lecturas específicas, impidiendo una asignación correcta.

Por lo tanto, la investigación propuesta no es meramente una tarea de recuperación de datos; es un análisis forense del proceso experimental completo. El pipeline se estructurará siguiendo una estrategia de "embudo": comenzando con un Control de Calidad (QC) global, seguido de una clasificación taxonómica rápida (Triage Taxonómico), y culminando en un análisis estructural profundo (Structural Deep Dive).

### **1.2. Módulo 1: Curación de Datos y Triage de Calidad**

El objetivo de este módulo es crear un subconjunto de lecturas largas de alta confianza, eliminando secuencias de baja calidad que podrían introducir ruido y artefactos en los análisis posteriores, optimizando así el uso de recursos computacionales.

#### **1.2.1. Evaluación Inicial de Calidad con NanoPlot**

Antes de cualquier filtrado, es esencial realizar una evaluación de calidad global sobre el archivo FASTQ de lecturas no clasificadas. Esto proporciona una línea base para comprender las características del conjunto de datos. NanoPlot es una herramienta estándar en la comunidad de Nanopore para este propósito, generando un informe HTML con múltiples gráficos descriptivos.  
Los gráficos clave a examinar en el informe de NanoPlot incluyen la distribución de longitudes de lectura (para confirmar la presencia de la población de interés \>75kb) y el gráfico de dispersión de longitud vs. calidad Phred media. Este último es particularmente útil para visualizar si las lecturas más largas tienden a tener una calidad inferior, lo que podría explicar los fallos en el demultiplexing.  
**Comando de Ejecución:**  
`NanoPlot --fastq unclassified.fastq.gz -o nanoplot_initial --loglength --threads 8`

Este comando procesará el archivo de lecturas no clasificadas, generará los gráficos en el directorio nanoplot\_initial y utilizará una escala logarítmica para el eje de longitud, lo que mejora la visualización de distribuciones amplias.

#### **1.2.2. Filtrado por Longitud y Calidad con Filtlong**

Este es un paso crítico para enfocar los recursos computacionales exclusivamente en las lecturas de interés. Filtlong es una herramienta diseñada específicamente para el filtrado de lecturas largas, ofreciendo métodos más sofisticados que un simple umbral de calidad. En lugar de aplicar un corte de calidad Phred estricto, que puede ser problemático para lecturas largas con calidad variable, se puede optar por conservar un porcentaje de las lecturas de mejor calidad o un número total de bases de las mejores lecturas.  
Para este análisis, se filtrarán las lecturas para retener solo aquellas con una longitud mínima de 75,000 nucleótidos. Simultáneamente, se eliminará el 10% de estas lecturas largas que presenten la peor calidad. Esta estrategia pragmática descarta las secuencias con mayor probabilidad de ser artefactos o de generar resultados de alineamiento espurios, sin ser excesivamente restrictivo.  
**Comando de Ejecución:**  
`filtlong --min_length 75000 --keep_percent 90 unclassified.fastq.gz | gzip > unclassified.long.hq.fastq.gz`

Este comando filtra las lecturas de unclassified.fastq.gz para que cumplan con el requisito de longitud mínima y luego conserva el 90% de esas lecturas que tienen la mejor puntuación de calidad combinada (longitud y calidad Phred). El resultado se comprime y se guarda en unclassified.long.hq.fastq.gz, que será el archivo de entrada para los siguientes módulos.

### **1.3. Módulo 2: Evaluación Taxonómica Rápida**

El objetivo de este módulo es determinar rápidamente el origen biológico de las lecturas. ¿Provienen del organismo de interés (*Candida*), de un contaminante específico, o de una mezcla compleja de especies? La respuesta a esta pregunta es crucial para guiar la estrategia de análisis posterior.

#### **1.3.1. Clasificación Basada en K-meros con Kraken2**

Kraken2 es una herramienta de clasificación taxonómica ultra-rápida que asigna una etiqueta taxonómica a cada secuencia de ADN comparando sus k-meros con una base de datos de genomas de referencia. Su velocidad y eficiencia la hacen ideal para obtener una primera visión de la composición de grandes conjuntos de datos de secuenciación.  
Se recomienda utilizar una base de datos pre-construida y completa, como la base de datos PlusPF (Standard plus protozoa and fungi), que incluye genomas de bacterias, arqueas, virus, protozoos y hongos. Esto asegura una amplia cobertura para detectar posibles contaminantes de diversas fuentes.  
**Comando de Ejecución:**  
`DB_PATH="/ruta/a/la/base_de_datos/PlusPF"`  
`kraken2 --db ${DB_PATH} \`  
        `--threads 16 \`  
        `--report kraken_report.txt \`  
        `unclassified.long.hq.fastq.gz > kraken_output.txt`

Este comando clasificará cada lectura en unclassified.long.hq.fastq.gz utilizando 16 hilos de procesamiento. Generará dos archivos de salida: kraken\_output.txt, con la clasificación por lectura, y kraken\_report.txt, un resumen de la abundancia de cada taxón en el conjunto de datos.

#### **1.3.2. Visualización Interactiva con Krona**

El archivo de reporte de Kraken2 es informativo, pero su formato de texto tabular no es intuitivo para explorar la jerarquía taxonómica. Krona es la herramienta estándar para visualizar estos resultados, creando un gráfico de sol (sunburst plot) interactivo en formato HTML que permite explorar la composición de la muestra de manera jerárquica y visual.  
**Comando de Ejecución:**  
`ktImportTaxonomy -q 2 -t 3 kraken_output.txt -o kraken_interactive.html`

Este comando convierte el archivo de salida por lectura de Kraken2 (kraken\_output.txt) en un archivo HTML interactivo. El resultado kraken\_interactive.html puede abrirse en cualquier navegador web.  
La realización de este perfil taxonómico en una etapa temprana del pipeline es estratégicamente fundamental. Si el gráfico de Krona revela que el 95% de las lecturas pertenecen a *Candida auris*, el problema principal es probablemente un fallo en el reconocimiento del código de barras o la presencia de quimeras estructurales que involucran el genoma diana. En este caso, el análisis debe centrarse intensamente en la caracterización estructural. Por otro lado, si el gráfico muestra una mezcla, por ejemplo, 50% *Escherichia coli* y 40% *Candida*, esto indica un problema de contaminación significativo. Este conocimiento, obtenido antes de ejecutar alineamientos o búsquedas BLAST computacionalmente costosas, permite bifurcar la estrategia: un camino para caracterizar las lecturas de *Candida* y otro para investigar la fuente de la contaminación (p. ej., contaminación de reactivos). Este punto de control temprano ahorra una cantidad considerable de tiempo y recursos, y personaliza el análisis posterior a la naturaleza real del problema.

### **1.4. Módulo 3: Análisis Estructural y de Homología en Profundidad**

El objetivo de este módulo es diseccionar la estructura interna y el origen genómico preciso de cada lectura. Este es el núcleo del análisis, diseñado para responder directamente a las preguntas sobre la presencia de concatémeros, repeticiones internas y quimeras.

#### **1.4.1. Detección de Repeticiones Internas y Concatémeros mediante Auto-alineamiento**

Para identificar estructuras repetitivas dentro de una misma lectura, una técnica poderosa es alinear el conjunto de lecturas contra sí mismo. minimap2 es la herramienta de elección para esta tarea, ya que está altamente optimizada para alinear lecturas largas y ruidosas, como las de Nanopore. Un resultado donde una lectura se alinea múltiples veces consigo misma en diferentes coordenadas es un fuerte indicio de repeticiones en tándem o de un artefacto de concatemerización, donde la misma molécula (p. ej., un plásmido) se ha secuenciado varias veces de forma consecutiva en una sola lectura.  
Se generará la salida en formato PAF (Pairwise mApping Format), un formato tabular simple que es mucho más fácil de procesar mediante scripts que el formato SAM para este tipo de análisis.  
**Comando de Ejecución:**  
`minimap2 -x ava-ont -c --dual=no \`  
         `unclassified.long.hq.fastq.gz unclassified.long.hq.fastq.gz > self_align.paf`

El preset \-x ava-ont está optimizado para la detección de solapamientos entre lecturas de Nanopore. La opción \--dual=no es crucial aquí, ya que evita que se reporte el alineamiento de una lectura consigo misma en ambas direcciones (A vs. B y B vs. A), simplificando el archivo de salida.

#### **1.4.2. Anotación Basada en Homología con BLAST**

Este paso aprovecha el trabajo ya realizado por el usuario, utilizando blastn para buscar homología contra una base de datos exhaustiva como nt (nucleotide collection) del NCBI. El formato de salida tabular (-outfmt 6\) es esencial para el procesamiento programático posterior. El archivo de ejemplo proporcionado ilustra perfectamente el valor de este análisis: una única lectura (b5c2375c-...) muestra alineamientos significativos tanto con el cromosoma 5 como con el cromosoma 6 de *Candidozyma auris*, una evidencia directa de un posible evento quimérico.  
**Comando de Ejecución (Ejemplo):**  
`gunzip -c unclassified.long.hq.fastq.gz | awk 'NR%4==1 |`

`| NR%4==2' | sed 's/^@/>/' > unclassified.long.hq.fasta`

`blastn -query unclassified.long.hq.fasta \`  
       `-db nt \`  
       `-out blast_vs_nt.tsv \`  
       `-outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle" \`  
       `-num_threads 16 \`  
       `-max_target_seqs 5`

Primero, se convierte el archivo FASTQ a FASTA. Luego, blastn se ejecuta con un formato de salida personalizado que incluye el título del sujeto (stitle), lo cual es útil para la interpretación. Limitar el número de alineamientos reportados por lectura (-max\_target\_seqs 5\) mantiene el tamaño del archivo de salida manejable mientras se captura la información de homología más relevante.

#### **1.4.3. Detección Dirigida de Quimeras con VSEARCH**

Para complementar la evidencia de quimerismo obtenida a través de BLAST, se utilizará una herramienta específicamente diseñada para esta tarea. VSEARCH incluye el algoritmo uchime\_ref, que es un estándar para la detección de quimeras basada en referencia. Este método intenta explicar una secuencia de consulta como una combinación lineal de dos o más secuencias de referencia "parentales".  
Se puede utilizar una base de datos de referencia curada, como una colección de genomas fúngicos de alta calidad o el propio genoma de referencia de *C. albicans* del usuario. Esto proporciona una capa adicional y formal de evidencia para clasificar las lecturas como quiméricas o no.  
**Comando de Ejecución:**  
`REF_DB="/ruta/a/genomas_fungicos.fasta"`  
`vsearch --uchime_ref unclassified.long.hq.fasta \`  
        `--db ${REF_DB} \`  
        `--uchimeout vsearch_chimeras.tsv \`  
        `--nonchimeras nonchimeras.fasta \`  
        `--chimeras chimeras.fasta \`  
        `--threads 16`

Este comando evalúa cada lectura en unclassified.long.hq.fasta contra la base de datos de referencia. El archivo vsearch\_chimeras.tsv contendrá la tabla de resultados, mientras que las secuencias se separarán en dos archivos FASTA, chimeras.fasta y nonchimeras.fasta, para una fácil inspección.

## **Parte II: Un Marco Unificado para la Clasificación Algorítmica de Lecturas**

Esta sección aborda el principal cuello de botella identificado: la síntesis de los diversos archivos de salida en una clasificación coherente y procesable. El concepto central que se introduce es la "Tabla Maestra de Anotación", que servirá como la base para todo el análisis y la visualización posteriores.

### **2.1. La Tabla Maestra de Anotación: La Clave para la Síntesis**

El objetivo es consolidar toda la información generada en la Parte I en una única estructura de datos, una tabla o *data frame* donde cada fila representa una única lectura larga y las columnas contienen los resultados de cada análisis. Esta transformación convierte un problema de manejo de múltiples archivos en un problema de análisis de datos estándar, que puede ser abordado eficientemente con herramientas como dplyr en R o pandas en Python.  
El proceso para construir esta tabla implica:

1. **Procesar la salida de Kraken2:** Extraer la asignación taxonómica principal y su puntuación de confianza para cada read\_id.  
2. **Procesar la salida de BLAST:** Agrupar los resultados por qseqid (el read\_id) para derivar métricas de resumen, como el número de cromosomas únicos a los que alinea, la lista de las mejores especies encontradas, y el porcentaje de la lectura cubierto por alineamientos. El uso de dplyr en R o pandas en Python es ideal para estas operaciones de agregación.  
3. **Procesar la salida de auto-alineamiento de minimap2:** Para cada read\_id, calcular la longitud total de sus auto-alineamientos. Una métrica derivada clave es la relación entre esta longitud total y la longitud de la propia lectura, que será un indicador potente de estructuras repetitivas.  
4. **Procesar la salida de VSEARCH:** Extraer la clasificación binaria (quimera/no quimera) para cada read\_id.  
5. **Unir todas las tablas:** Combinar la información de los pasos anteriores en una única tabla, utilizando el read\_id como clave de unión.

Esta tabla es el pilar de todo el análisis posterior. Centraliza la evidencia de métodos ortogonales (basados en k-meros, en alineamiento de homología y en auto-alineamiento), lo que permite conclusiones mucho más robustas. Además, sirve como la única fuente de datos para los resúmenes estadísticos y la aplicación interactiva, garantizando la coherencia en todos los resultados.  
**Tabla 1: Estructura de la Tabla Maestra de Anotación**

| Columna | Descripción | Origen del Dato |
| :---- | :---- | :---- |
| read\_id | Identificador único de la lectura. | FASTQ |
| read\_length | Longitud de la lectura en nucleótidos. | FASTQ |
| mean\_q\_score | Calidad Phred media de la lectura. | FASTQ |
| kraken\_taxon | Taxón principal asignado por Kraken2. | Kraken2 |
| kraken\_confidence | Puntuación de confianza de la asignación de Kraken2. | Kraken2 |
| blast\_hit\_count | Número total de alineamientos BLAST significativos. | BLAST |
| blast\_unique\_subjects | Número de especies/secuencias de referencia únicas en los hits. | BLAST |
| blast\_unique\_chromosomes | Número de cromosomas únicos en los hits (para una especie dada). | BLAST |
| is\_chimeric\_vsearch | Clasificación de VSEARCH ('Y' o 'N'). | VSEARCH |
| self\_align\_total\_len | Suma de las longitudes de todos los auto-alineamientos. | minimap2 |
| self\_align\_ratio | self\_align\_total\_len / read\_length. | Derivado |
| structural\_category | Categoría final asignada (ver sección 2.2). | Derivado |

### **2.2. Un Motor Basado en Reglas para la Categorización de Lecturas**

El siguiente paso es automatizar la clasificación de cada lectura en una categoría estructural y biológicamente significativa, basándose en los datos sintetizados en la Tabla Maestra de Anotación.  
Los artefactos y las secuencias biológicas genuinas a menudo presentan firmas distintivas que son multivariadas. Una sola pieza de evidencia puede ser ambigua: una lectura con dos hits de BLAST podría ser una quimera o simplemente abarcar una región repetitiva que existe en dos cromosomas. El poder analítico surge de la combinación de señales. Por ejemplo:

* Una **Quimera Intra-Genómica** se caracterizará por tener blast\_unique\_chromosomes \> 1 y is\_chimeric\_vsearch \== 'Y', pero kraken\_taxon apuntará a una única especie.  
* Una **Quimera Inter-Especies** tendrá blast\_unique\_subjects \> 1 con nombres de especies diferentes.  
* Un **Concatémero** mostrará un self\_align\_ratio muy elevado (p. ej., \> 2.0), indicando que la longitud total de sus repeticiones internas supera con creces la longitud de la propia lectura.  
* Una lectura **Canónica** ideal tendrá blast\_unique\_chromosomes \<= 1, un self\_align\_ratio cercano a 1.0, y is\_chimeric\_vsearch \== 'N'.

Esta lógica multifactorial es inherentemente más robusta que confiar en la salida de una sola herramienta. Se implementará mediante una función en R o Python que aplique una serie de reglas priorizadas (una cascada de sentencias if/elif/else) a cada fila de la Tabla Maestra de Anotación para rellenar la columna structural\_category.  
**Categorías Estructurales Propuestas:**

1. **Canónica:** Lectura que alinea de forma continua y única a una sola localización genómica. Representa el ideal de una secuencia biológica intacta.  
2. **Repetitiva/Concatémero:** Lectura con una alta puntuación de auto-alineamiento (self\_align\_ratio \> 1.8). Probablemente un artefacto de la secuenciación de moléculas circulares o regiones con repeticiones en tándem de alto número de copias.  
3. **Quimera Intra-Genómica:** Partes de la lectura mapean a diferentes localizaciones (cromosomas distintos o regiones distantes del mismo cromosoma) dentro de la *misma* especie.  
4. **Quimera Inter-Especies:** Partes de la lectura mapean a genomas de *diferentes* especies. Un claro artefacto de librería o contaminación.  
5. **Contaminante:** Lectura que alinea de forma canónica pero a un organismo no esperado (p. ej., *E. coli*).  
6. **Fragmentada/Alineamiento Pobre:** La lectura no alinea de forma coherente en una porción significativa de su longitud. Puede deberse a muy baja calidad o a una estructura extremadamente compleja.  
7. **Novedosa/No Anotada:** Lectura de alta calidad que no obtiene hits significativos en BLAST contra la base de datos nt. Podría representar regiones genómicas nuevas o muy divergentes.

### **2.3. Generación de Salidas de Resumen**

Una vez que cada lectura ha sido clasificada, se pueden generar resúmenes de alto nivel para el informe final. Estos resúmenes cuantifican el alcance y la naturaleza del problema de las lecturas no clasificadas.  
**Salidas Clave:**

1. **Tabla de Resumen:** Una tabla que cuenta el número de lecturas y el total de pares de bases para cada structural\_category. Esto proporciona una cuantificación directa del impacto de cada tipo de artefacto o secuencia.  
2. **Gráfico de Barras:** Una visualización de la tabla de resumen, mostrando la proporción de lecturas en cada categoría. Esto ofrece una visión general inmediata de la composición del conjunto de datos.  
3. **Archivo CSV Final:** La Tabla Maestra de Anotación completa, con la columna de structural\_category añadida, se guardará como un archivo CSV. Este archivo es un resultado final valioso que puede ser utilizado para análisis adicionales o como material suplementario.

## **Parte III: Una Suite de Visualización Interactiva y Multicapa**

Esta sección proporciona el diseño y el código para los productos analíticos finales, centrándose en el panel de control interactivo y las visualizaciones detalladas por lectura, que son el objetivo final de la consulta del usuario para lograr una interpretación clara y profunda.

### **3.1. El Panel de Control Global: Una Aplicación R Shiny**

Para facilitar la exploración interactiva de todo el conjunto de lecturas largas no clasificadas, se propone el desarrollo de una aplicación web utilizando el framework R Shiny. Shiny permite crear interfaces de usuario web interactivas directamente desde R, utilizando data frames como fuente de datos, lo que lo hace perfecto para este caso de uso.  
**Diseño de la Interfaz de Usuario (UI) y Lógica del Servidor:**  
La aplicación se estructurará en varias pestañas para una navegación lógica:

* **Pestaña 1: Resumen General (Overview).**  
  * Mostrará el gráfico de barras de resumen de las categorías estructurales (generado en la sección 2.3).  
  * Incrustará el gráfico interactivo de Krona (kraken\_interactive.html) dentro de un iframe para una exploración taxonómica directa.  
* **Pestaña 2: Explorador de Datos (Data Explorer).**  
  * Presentará la Tabla Maestra de Anotación completa en una tabla interactiva utilizando el paquete DT::datatable.  
  * Los usuarios podrán ordenar por cualquier columna, buscar por read\_id, y filtrar los datos para seleccionar subconjuntos de interés (p. ej., "mostrar todas las lecturas clasificadas como 'Quimera Intra-Genómica' con una longitud superior a 90,000 nt").  
* **Pestaña 3: Inspector de Lectura Individual (Single Read Inspector).**  
  * Contendrá un menú desplegable poblado con los read\_ids. La selección de un read\_id (que puede ser alimentada desde la tabla en la Pestaña 2\) activará la generación y visualización del gráfico detallado multi-panel descrito en la siguiente sección.

### **3.2. El Inspector de Lectura Individual: Una Visualización Profunda**

El objetivo de esta visualización es crear un único gráfico, denso en información, para cualquier lectura seleccionada. Este gráfico sintetizará todos los datos disponibles, permitiendo una interpretación visual inmediata de su estructura y origen.  
La clave de esta visualización es la síntesis. En lugar de examinar gráficos separados para la calidad, los alineamientos y la composición, se integran todos en una única vista de "navegador genómico", donde la propia lectura actúa como secuencia de referencia. Este enfoque permite al analista reconocer patrones complejos de forma intuitiva. Por ejemplo, una caída en el perfil de calidad puede correlacionarse visualmente con un punto de unión quimérico en la pista de alineamientos BLAST. Un patrón regular y repetitivo en los arcos de auto-alineamiento es una firma inequívoca de un concatémero. Esta vista integrada elimina la carga cognitiva de tener que conectar mentalmente múltiples gráficos y facilita la obtención de conclusiones directas.  
Para la implementación en R, se utilizará el poder de ggplot2 en combinación con paquetes de Bioconductor especializados en visualización genómica, como ggbio o trackViewer , que están diseñados para crear este tipo de gráficos con pistas apiladas.  
**Diseño del Gráfico Multi-Panel:**  
El gráfico para una lectura individual se compondrá de las siguientes pistas apiladas, todas compartiendo el mismo eje x que representa las coordenadas a lo largo de la lectura (desde 0 hasta su longitud total).

* **Panel 1: Perfil de Composición de Secuencia.**  
  * **Contenido:** Un gráfico de líneas que muestra el contenido de GC calculado en una ventana deslizante (p. ej., ventana de 500 pb con un paso de 100 pb).  
  * **Implementación:** Se puede calcular fácilmente en R utilizando funciones como zoo::rollapply sobre la secuencia de la lectura.  
  * **Propósito:** Revelar cambios abruptos en la composición de bases, que pueden ser indicativos de uniones quiméricas entre fragmentos de ADN de diferentes orígenes genómicos (p. ej., de regiones ricas en GC a regiones ricas en AT).  
* **Panel 2: Mapa de Anotación Estructural (Panel Central).**  
  * **Contenido:** Este es el panel principal que visualiza la estructura de la lectura. La lectura se representa como una barra gris a lo largo del eje x.  
  * **Alineamientos BLAST:** Se dibujan como rectángulos de colores (geom\_rect) sobre la barra de la lectura. La posición y longitud del rectángulo corresponden a las coordenadas del alineamiento en la lectura (qstart, qend). El color puede mapearse al cromosoma o a la especie del sujeto del alineamiento, haciendo que las quimeras sean visualmente obvias. Una funcionalidad interactiva (usando plotly::ggplotly) puede mostrar detalles como el porcentaje de identidad y el E-value al pasar el ratón sobre cada rectángulo. Esto visualiza directamente los datos tabulares de BLAST, como los del archivo de ejemplo.  
  * **Auto-Alineamientos:** Se representan como arcos (geom\_curve) que conectan las regiones de la lectura que se alinean consigo mismas. Un patrón denso y regular de arcos idénticos es la firma visual clásica de un concatémero. Alineamientos esporádicos indican la presencia de familias de genes repetidos o elementos transponibles.  
* **Panel 3: Perfil de Puntuación de Calidad.**  
  * **Contenido:** Un gráfico de líneas o de área que muestra la puntuación de calidad Phred a lo largo de la lectura.  
  * **Implementación:** Requiere extraer los valores de calidad del archivo FASTQ original para la lectura seleccionada.  
  * **Propósito:** Las caídas bruscas en la calidad a menudo se correlacionan con regiones problemáticas de la secuencia, como homopolímeros largos, o pueden coincidir con los puntos de ruptura de alineamientos quiméricos, proporcionando una evidencia convergente del artefacto.

## **Parte IV: Guía de Implementación y Repositorio de Código**

Esta sección final proporciona el código completo y comentado para ejecutar todo el flujo de trabajo, garantizando que el análisis sea reproducible desde el principio hasta el final.

### **4.1. El Script Completo del Pipeline en Bash (characterize\_long\_reads.sh)**

Se proporciona un único script ejecutable en Bash que encadena todos los comandos de la Parte I. El script está parametrizado para facilitar su uso y reutilización, e incluye comprobaciones básicas de errores y de la existencia de archivos de entrada.  
`#!/bin/bash`

`# ==============================================================================`  
`# SCRIPT: characterize_long_reads.sh`  
`# DESCRIPCIÓN: Pipeline para la caracterización de lecturas largas de Nanopore`  
`#              no clasificadas.`  
`# AUTOR: Experto en Bioinformática Computacional`  
`# FECHA: 2024-09-05`  
`# ==============================================================================`

`set -e # El script se detendrá si algún comando falla`  
`set -u # Trata las variables no definidas como un error`  
`set -o pipefail # El código de salida de un pipeline es el del último comando que falló`

`# --- PARÁMETROS ---`  
`INPUT_FQ="unclassified.fastq.gz"`  
`OUTPUT_DIR="analysis_output"`  
`MIN_LEN=75000`  
`THREADS=16`  
`KRAKEN_DB="/path/to/your/kraken2/PlusPF_db"`  
`BLAST_DB="nt"`  
`REF_DB_VSEARCH="/path/to/your/fungal_ref.fasta"`

`# --- INICIO DEL PIPELINE ---`  
`echo "--- Iniciando el pipeline de caracterización de lecturas largas ---"`  
`mkdir -p ${OUTPUT_DIR}`  
`cd ${OUTPUT_DIR}`

`# --- MÓDULO 1: QC Y FILTRADO ---`  
`echo "[Módulo 1] Ejecutando QC inicial y filtrado..."`

`# 1.1 NanoPlot para QC inicial`  
`if [! -d "1_nanoplot_initial" ]; then`  
    `echo "  -> Ejecutando NanoPlot..."`  
    `NanoPlot --fastq../${INPUT_FQ} -o 1_nanoplot_initial --loglength --threads ${THREADS}`  
`else`  
    `echo "  -> Directorio de NanoPlot ya existe. Saltando."`  
`fi`

`# 1.2 Filtlong para filtrar por longitud y calidad`  
`LONG_HQ_FQ="unclassified.long.hq.fastq.gz"`  
`if [! -f "${LONG_HQ_FQ}" ]; then`  
    `echo "  -> Ejecutando Filtlong..."`  
    `filtlong --min_length ${MIN_LEN} --keep_percent 90../${INPUT_FQ} | gzip > ${LONG_HQ_FQ}`  
`else`  
    `echo "  -> Archivo filtrado ${LONG_HQ_FQ} ya existe. Saltando."`  
`fi`

`# Convertir a FASTA para herramientas posteriores`  
`LONG_HQ_FA="unclassified.long.hq.fasta"`  
`if [! -f "${LONG_HQ_FA}" ]; then`  
    `echo "  -> Convirtiendo FASTQ a FASTA..."`  
    `gunzip -c ${LONG_HQ_FQ} | awk 'NR%4==1 |`

`| NR%4==2' | sed 's/^@/>/' > ${LONG_HQ_FA}`  
`else`  
    `echo "  -> Archivo FASTA ${LONG_HQ_FA} ya existe. Saltando."`  
`fi`

`echo "[Módulo 1] Completado."`

`# --- MÓDULO 2: ANÁLISIS TAXONÓMICO ---`  
`echo "[Módulo 2] Ejecutando análisis taxonómico..."`

`# 2.1 Kraken2 para clasificación`  
`KRAKEN_OUT="kraken_output.txt"`  
`KRAKEN_REPORT="kraken_report.txt"`  
`if; then`  
    `echo "  -> Ejecutando Kraken2..."`  
    `kraken2 --db ${KRAKEN_DB} \`  
            `--threads ${THREADS} \`  
            `--report ${KRAKEN_REPORT} \`  
            `${LONG_HQ_FQ} > ${KRAKEN_OUT}`  
`else`  
    `echo "  -> Archivo de salida de Kraken2 ya existe. Saltando."`  
`fi`

`# 2.2 Krona para visualización`  
`KRONA_HTML="kraken_interactive.html"`  
`if; then`  
    `echo "  -> Generando informe de Krona..."`  
    `ktImportTaxonomy -q 2 -t 3 ${KRAKEN_OUT} -o ${KRONA_HTML}`  
`else`  
    `echo "  -> Archivo de Krona ya existe. Saltando."`  
`fi`

`echo "[Módulo 2] Completado."`

`# --- MÓDULO 3: ANÁLISIS ESTRUCTURAL ---`  
`echo "[Módulo 3] Ejecutando análisis estructural en profundidad..."`

`# 3.1 Auto-alineamiento con minimap2`  
`SELF_ALIGN_PAF="self_align.paf"`  
`if; then`  
    `echo "  -> Ejecutando auto-alineamiento con minimap2..."`  
    `minimap2 -x ava-ont -c --dual=no \`  
             `${LONG_HQ_FQ} ${LONG_HQ_FQ} > ${SELF_ALIGN_PAF}`  
`else`  
    `echo "  -> Archivo PAF de auto-alineamiento ya existe. Saltando."`  
`fi`

`# 3.2 BLASTn contra la base de datos nt`  
`BLAST_TSV="blast_vs_nt.tsv"`  
`if; then`  
    `echo "  -> Ejecutando BLASTn..."`  
    `blastn -query ${LONG_HQ_FA} \`  
           `-db ${BLAST_DB} \`  
           `-out ${BLAST_TSV} \`  
           `-outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle" \`  
           `-num_threads ${THREADS} \`  
           `-max_target_seqs 5`  
`else`  
    `echo "  -> Archivo de resultados de BLAST ya existe. Saltando."`  
`fi`

`# 3.3 Detección de quimeras con VSEARCH`  
`VSEARCH_TSV="vsearch_chimeras.tsv"`  
`if; then`  
    `echo "  -> Ejecutando VSEARCH para detección de quimeras..."`  
    `vsearch --uchime_ref ${LONG_HQ_FA} \`  
            `--db ${REF_DB_VSEARCH} \`  
            `--uchimeout ${VSEARCH_TSV} \`  
            `--nonchimeras nonchimeras.fasta \`  
            `--chimeras chimeras.fasta \`  
            `--threads ${THREADS}`  
`else`  
    `echo "  -> Archivo de resultados de VSEARCH ya existe. Saltando."`  
`fi`

`echo "[Módulo 3] Completado."`  
`echo "--- Pipeline finalizado con éxito. Los resultados se encuentran en el directorio ${OUTPUT_DIR} ---"`

### **4.2. El Código de Análisis y Visualización en R (analysis\_and\_visualization.R)**

Se proporciona un script de R bien documentado que realiza todos los pasos de las Partes II y III. Este script está diseñado para ser ejecutado después del pipeline de Bash y asume que se encuentra en el directorio de trabajo que contiene la carpeta analysis\_output.  
`# ==============================================================================`  
`# SCRIPT: analysis_and_visualization.R`  
`# DESCRIPCIÓN: Script en R para procesar, clasificar y visualizar los resultados`  
`#              del pipeline de caracterización de lecturas largas.`  
`# AUTOR: Experto en Bioinformática Computacional`  
`# FECHA: 2024-09-05`  
`# ==============================================================================`

`# --- 0. Carga de Librerías ---`  
`# Asegúrate de tener estas librerías instaladas:`  
`# install.packages(c("tidyverse", "data.table", "Biostrings", "zoo", "DT", "shiny", "plotly", "ggbio"))`  
`library(tidyverse)`  
`library(data.table)`  
`library(Biostrings)`  
`library(zoo)`  
`library(DT)`  
`library(shiny)`  
`library(plotly)`  
`library(ggbio)`

`# --- A. Carga y Procesamiento de Datos ---`  
`message("--- Sección A: Cargando y procesando datos del pipeline ---")`  
`analysis_dir <- "analysis_output"`

`# 1. Cargar secuencias FASTA`  
`sequences <- readDNAStringSet(file.path(analysis_dir, "unclassified.long.hq.fasta"))`  
`seq_info <- tibble(`  
  `read_id = names(sequences),`  
  `read_length = width(sequences),`  
  `sequence = as.character(sequences)`  
`)`

`# 2. Procesar salida de BLAST`  
`message("  -> Procesando resultados de BLAST...")`  
`blast_cols <- c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",`   
                `"qstart", "qend", "sstart", "send", "evalue", "bitscore", "stitle")`  
`blast_results <- fread(file.path(analysis_dir, "blast_vs_nt.tsv"), col.names = blast_cols) %>%`  
  `as_tibble() %>%`  
  `mutate(subject_species = str_extract(stitle, "^[^ ]+ [^ ]+")) %>% # Extraer especie (ej. "Candida auris")`  
  `mutate(chromosome = str_extract(stitle, "chromosome [^,]+|plasmid [^,]+")) # Extraer cromosoma o plásmido`

`blast_summary <- blast_results %>%`  
  `group_by(qseqid) %>%`  
  `summarise(`  
    `blast_hit_count = n(),`  
    `blast_unique_subjects = n_distinct(subject_species),`  
    `blast_unique_chromosomes = n_distinct(chromosome[!is.na(chromosome)])`  
  `) %>%`  
  `rename(read_id = qseqid)`

`# 3. Procesar salida de auto-alineamiento (minimap2 PAF)`  
`message("  -> Procesando resultados de auto-alineamiento...")`  
`paf_cols <- c("q_name", "q_len", "q_start", "q_end", "strand", "t_name", "t_len",`   
              `"t_start", "t_end", "matches", "aln_len", "mapq")`  
`self_align_paf <- fread(file.path(analysis_dir, "self_align.paf"), select = 1:12, col.names = paf_cols) %>%`  
  `as_tibble() %>%`  
  `filter(q_name == t_name) %>% # Asegurarse de que son auto-alineamientos`  
  `filter(q_start!= 0 | q_end!= q_len) # Excluir el alineamiento perfecto de una lectura consigo misma`

`self_align_summary <- self_align_paf %>%`  
  `group_by(q_name) %>%`  
  `summarise(self_align_total_len = sum(aln_len)) %>%`  
  `rename(read_id = q_name)`

`# 4. Procesar salida de VSEARCH`  
`message("  -> Procesando resultados de VSEARCH...")`  
`vsearch_results <- fread(file.path(analysis_dir, "vsearch_chimeras.tsv"), header = FALSE, select = 2, col.names = "read_id") %>%`  
  `as_tibble() %>%`  
  `mutate(is_chimeric_vsearch = "Y")`

`# --- B. Creación de la Tabla Maestra de Anotación ---`  
`message("--- Sección B: Creando la Tabla Maestra de Anotación ---")`  
`master_table <- seq_info %>%`  
  `left_join(blast_summary, by = "read_id") %>%`  
  `left_join(self_align_summary, by = "read_id") %>%`  
  `left_join(vsearch_results, by = "read_id") %>%`  
  `mutate(`  
    `across(where(is.numeric), ~replace_na(., 0)), # Reemplazar NAs numéricos con 0`  
    `is_chimeric_vsearch = replace_na(is_chimeric_vsearch, "N"),`  
    `self_align_ratio = ifelse(read_length > 0, self_align_total_len / read_length, 0)`  
  `)`

`# --- C. Motor de Clasificación de Lecturas ---`  
`message("--- Sección C: Clasificando lecturas con motor de reglas ---")`  
`classify_read <- function(df_row) {`  
  `if (df_row$self_align_ratio > 1.8) {`  
    `return("Repetitiva/Concatémero")`  
  `} else if (df_row$is_chimeric_vsearch == "Y" |`

`| df_row$blast_unique_chromosomes > 1) {`  
    `if (df_row$blast_unique_subjects > 1) {`  
      `return("Quimera Inter-Especies")`  
    `} else {`  
      `return("Quimera Intra-Genómica")`  
    `}`  
  `} else if (df_row$blast_hit_count > 0) {`  
    `# Aquí se podría añadir lógica para detectar contaminantes vs. organismo de interés`  
    `return("Canónica/Fragmento")`  
  `} else {`  
    `return("Novedosa/No Anotada")`  
  `}`  
`}`

`master_table <- master_table %>%`  
  `mutate(structural_category = apply(., 1, classify_read))`

`# Guardar la tabla maestra`  
`fwrite(master_table %>% select(-sequence), "master_annotation_table.csv")`

`# --- D. Funciones de Visualización ---`  
`message("--- Sección D: Definiendo funciones de visualización ---")`

`# 1. Gráfico de resumen`  
`summary_plot <- ggplot(master_table, aes(x = fct_infreq(structural_category), fill = structural_category)) +`  
  `geom_bar() +`  
  `labs(title = "Clasificación Estructural de Lecturas Largas No Clasificadas",`  
       `x = "Categoría Estructural", y = "Número de Lecturas") +`  
  `theme_minimal() +`  
  `theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")`  
`ggsave("summary_classification_plot.png", summary_plot, width = 10, height = 7)`

`# 2. Función para el Inspector de Lectura Individual`  
`plot_single_read <- function(read_id_to_plot, master_df, blast_df, self_align_df) {`  
    
  `# Filtrar datos para la lectura seleccionada`  
  `read_data <- filter(master_df, read_id == read_id_to_plot)`  
  `blast_hits <- filter(blast_df, qseqid == read_id_to_plot)`  
  `self_hits <- filter(self_align_df, q_name == read_id_to_plot)`  
    
  `if (nrow(read_data) == 0) {`  
    `return(ggplot() + annotate("text", x=0, y=0, label="Read ID no encontrado.") + theme_void())`  
  `}`  
    
  `# Panel 1: Contenido GC`  
  `seq <- DNAString(read_data$sequence)`  
  `gc_freq <- letterFrequencyInSlidingView(seq, view.width = 500, letters = "GC", as.prob = TRUE)`  
  `gc_df <- tibble(position = 1:length(gc_freq) * 500, gc_content = gc_freq * 100)`  
    
  `p_gc <- ggplot(gc_df, aes(x = position, y = gc_content)) +`  
    `geom_line(color = "darkgreen") +`  
    `labs(y = "GC (%)") +`  
    `theme_minimal() +`  
    `xlim(0, read_data$read_length)`  
    
  `# Panel 2: Anotación Estructural`  
  `p_structure <- ggplot() +`  
    `# Barra de la lectura`  
    `geom_segment(aes(x = 0, xend = read_data$read_length, y = 0, yend = 0), size = 2, color = "grey") +`  
    `# Alineamientos BLAST`  
    `geom_rect(data = blast_hits, aes(xmin = qstart, xmax = qend, ymin = -0.1, ymax = 0.1, fill = chromosome), alpha = 0.8) +`  
    `# Auto-alineamientos`  
    `geom_curve(data = self_hits, aes(x = q_start, xend = q_end, y = 0.15, yend = 0.15), curvature = -0.5, color = "purple", alpha = 0.5) +`  
    `labs(y = "Estructura", x = "Posición en la Lectura (pb)") +`  
    `theme_minimal() +`  
    `theme(axis.text.y = element_blank()) +`  
    `ylim(-0.5, 0.5)`  
      
  `# Combinar plots usando ggbio`  
  `tracks(p_gc, p_structure, heights = c(1, 2), title = paste("Inspector para:", read_id_to_plot))`  
`}`

`# Ejemplo de uso:`  
`# plot_single_read("b5c2375c-85d3-43ad-bce7-05ba76834e14", master_table, blast_results, self_align_paf)`

`# --- E. Aplicación Shiny ---`  
`message("--- Sección E: Código de la aplicación Shiny ---")`  
`# Este código puede guardarse como app.R y ejecutarse con shiny::runApp()`

`ui <- fluidPage(`  
  `titlePanel("Explorador de Lecturas Largas de Nanopore No Clasificadas"),`  
    
  `tabsetPanel(`  
    `tabPanel("Resumen General",`  
             `fluidRow(`  
               `column(8, plotOutput("summaryPlot")),`  
               `column(4, h4("Explorador Taxonómico (Krona)"),`  
                      `htmlOutput("krona_iframe"))`  
             `)`  
    `),`  
    `tabPanel("Explorador de Datos",`  
             `DT::dataTableOutput("masterTable")`  
    `),`  
    `tabPanel("Inspector de Lectura Individual",`  
             `sidebarLayout(`  
               `sidebarPanel(`  
                 `selectInput("selected_read", "Seleccionar Read ID:", choices = master_table$read_id)`  
               `),`  
               `mainPanel(`  
                 `plotOutput("singleReadPlot", height = "600px")`  
               `)`  
             `)`  
    `)`  
  `)`  
`)`

`server <- function(input, output, session) {`  
    
  `output$summaryPlot <- renderPlot({`  
    `summary_plot`  
  `})`  
    
  `output$krona_iframe <- renderUI({`  
    `tags$iframe(src = "kraken_interactive.html", height = "600px", width = "100%")`  
  `})`  
    
  `output$masterTable <- DT::renderDataTable({`  
    `DT::datatable(master_table %>% select(-sequence), options = list(pageLength = 15), filter = 'top')`  
  `})`  
    
  `output$singleReadPlot <- renderPlot({`  
    `req(input$selected_read)`  
    `plot_single_read(input$selected_read, master_table, blast_results, self_align_paf)`  
  `})`  
    
`}`

`# Para ejecutar la app:`  
`# shinyApp(ui, server)`

### **4.3. Interpretación de los Resultados y Próximos Pasos**

La ejecución de este pipeline y la exploración de sus resultados proporcionarán una imagen clara y multifacética del contenido del contenedor "unclassified". La interpretación final guiará las acciones a tomar:

* **Si la mayoría de las lecturas son "Quimera Intra-Genómica" o "Canónica/Fragmento" del organismo de interés:** Esto es una excelente noticia. Significa que los datos son biológicamente relevantes y valiosos. El problema principal radica en el proceso de demultiplexing. Las lecturas pueden ser recuperadas y añadidas al conjunto de datos de la muestra correspondiente (si se puede inferir su origen) o analizadas como un conjunto agregado. Se podría investigar el uso de herramientas de demultiplexing más sensibles o que operen en el espacio de la señal eléctrica cruda, como Deepbinner.  
* **Si una fracción significativa es "Repetitiva/Concatémero":** Esto apunta a un artefacto específico de la preparación de la librería, posiblemente relacionado con la ligación de moléculas circulares o la amplificación por círculo rodante. Estas lecturas deben ser tratadas con precaución. Aunque pueden ser útiles para ensamblar regiones repetitivas, su estructura no representa la del genoma nativo y deben ser excluidas de análisis de variantes estructurales.  
* **Si hay una alta proporción de "Contaminante" o "Quimera Inter-Especies":** Esto indica un problema de contaminación en el laboratorio. El análisis taxonómico de estas lecturas puede ayudar a identificar la fuente (p. ej., reactivos, ambiente, cross-contaminación entre muestras). Estas lecturas deben ser descartadas, y se deben revisar los protocolos de laboratorio para minimizar la contaminación en futuras secuenciaciones.  
* **Si predominan las lecturas "Novedosa/No Anotada":** Este es un escenario emocionante que podría indicar la presencia de grandes regiones genómicas no ensambladas previamente, elementos móviles masivos, o incluso el genoma de un organismo simbionte o un virus de gran tamaño no presente en las bases de datos. Estas lecturas merecen un análisis de ensamblaje *de novo* por separado.

En conclusión, este pipeline transforma el problema de las lecturas no clasificadas de una simple "pérdida de datos" a una poderosa herramienta de diagnóstico y descubrimiento. Proporciona un camino claro para clasificar cada lectura, tomar decisiones informadas sobre si descartarla o recuperarla, y obtener información valiosa para mejorar los protocolos experimentales y analíticos futuros.

#### **Obras citadas**

1\. Assessing Read Quality, Trimming and Filtering – QC & Assembly \- GitHub Pages, https://cloud-span.github.io/nerc-metagenomics03-qc-assembly/01-QC-quality-raw-reads/index.html 2\. NanoPlot | PSC \- Pittsburgh Supercomputing Center, https://www.psc.edu/resources/software/nanoplot/ 3\. wdecoster/NanoPlot: Plotting scripts for long read sequencing data \- GitHub, https://github.com/wdecoster/NanoPlot 4\. NanoPlot Report, https://physiology.med.cornell.edu/faculty/skrabanek/lab/angsd/lecture\_notes/NanoPlotReport.pdf 5\. NanoPack: visualizing and processing long-read sequencing data \- Oxford Academic, https://academic.oup.com/bioinformatics/article/34/15/2666/4934939 6\. rrwick/Filtlong: quality filtering tool for long reads \- GitHub, https://github.com/rrwick/Filtlong 7\. modules/filtlong \- nf-core, https://nf-co.re/modules/filtlong 8\. Manual · DerrickWood/kraken2 Wiki \- GitHub, https://github.com/DerrickWood/kraken2/wiki/Manual 9\. How do I measure the proportion of different bacteria in a sample from a high-throughput sequencer? \- Bioinformatics Stack Exchange, https://bioinformatics.stackexchange.com/questions/20542/how-do-i-measure-the-proportion-of-different-bacteria-in-a-sample-from-a-high-th 10\. wf-metagenomics \- Oxford Nanopore Technologies, https://nanoporetech.com/document/epi2me-workflows/wf-metagenomics 11\. Microbiome / Pathogen detection from (direct Nanopore) sequencing data using Galaxy \- Foodborne Edition / Hands-on, https://training.galaxyproject.org/training-material/topics/microbiome/tutorials/pathogen-detection-from-nanopore-foodborne-data/tutorial.html 12\. Generating Krona plots from Kraken data | Microbiome binfies \- Andrea Telatin, https://telatin.github.io/microbiome-bioinformatics/Kraken-to-Krona/ 13\. 6\. Taxonomic investigation \- Computational Genomics Tutorial \- Sebastian Schmeier, https://genomics.sschmeier.com/ngs-taxonomic-investigation/index.html 14\. 3\. Metagenomic screening of shotgun data — Physalia Paleogenomics 0.1.0 documentation, https://physalia-paleogenomics-2020.readthedocs.io/en/latest/3\_Metagenomics\_v2.html 15\. Manual Reference Pages \- minimap2 (1), https://lh3.github.io/minimap2/minimap2.html 16\. New strategies to improve minimap2 alignment accuracy | Bioinformatics \- Oxford Academic, https://academic.oup.com/bioinformatics/article/37/23/4572/6384570 17\. minimap2 — Debian testing, https://manpages.debian.org/testing/minimap2/minimap2.1.en.html 18\. VSearch chimera detection \- Galaxy | Tool Shed, https://toolshed.g2.bx.psu.edu/repository/display\_tool?repository\_id=8812a0146423e217\&tool\_config=%2Fsrv%2Ftoolshed-repos%2Fmain%2F001%2Frepo\_1629%2Fchimera.xml\&changeset\_revision=9495df9dd6ef\&render\_repository\_actions\_for=tool\_shed 19\. Reference-based chimera detection | EUKARYOME, https://eukaryome.org/reference-based-chimera-detection/ 20\. Learn dplyr: Learn R: Joining Tables Cheatsheet | Codecademy, https://www.codecademy.com/learn/learn-dplyr/modules/learn-r-multiple-tables/cheatsheet 21\. 19 Joins \- R for Data Science (2e), https://r4ds.hadley.nz/joins.html 22\. Mutating joins \- dplyr, https://dplyr.tidyverse.org/reference/mutate-joins.html 23\. BLAST+6 format (skbio.io.format.blast6), https://scikit.bio/docs/latest/generated/skbio.io.format.blast6.html 24\. Merge Multiple Dataframes \- Pandas \- GeeksforGeeks, https://www.geeksforgeeks.org/pandas/how-to-merge-multiple-dataframes-in-pandas/ 25\. Merge, join, concatenate and compare — pandas 2.3.2 documentation \- PyData |, https://pandas.pydata.org/docs/user\_guide/merging.html 26\. GraphBio: A shiny web app to easily perform popular visualization analysis for omics data, https://www.frontiersin.org/journals/genetics/articles/10.3389/fgene.2022.957317/full 27\. shinyCircos: an R/Shiny application for interactive creation of Circos plot \- Oxford Academic, https://academic.oup.com/bioinformatics/article/34/7/1229/4657077 28\. StructuRly: A novel shiny app to produce comprehensive, detailed and interactive plots for population genetic analysis | PLOS One, https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0229330 29\. ggbio \- Visualize genomic data \- Easy Guides \- Wiki \- STHDA, https://www.sthda.com/english/wiki/ggbio-visualize-genomic-data 30\. ggbio: Visualization tools for genomic data \- bioc \- R-universe, https://bioc.r-universe.dev/ggbio 31\. ggbio: visualization toolkits for genomic data \- Bioconductor, https://www.bioconductor.org/packages/devel/bioc/vignettes/ggbio/inst/doc/ggbio.pdf 32\. trackViewer Vignette: overview \- Bioconductor, https://bioconductor.org/packages/devel/bioc/vignettes/trackViewer/inst/doc/trackViewer.html 33\. jianhong/trackViewer: A bioconductor package with minimalist design for drawing elegant tracks or lollipop plot \- GitHub, https://github.com/jianhong/trackViewer 34\. Genome coverage as sliding window \- Stack Overflow, https://stackoverflow.com/questions/48721332/genome-coverage-as-sliding-window 35\. Demultiplexing barcoded Oxford Nanopore reads with deep convolutional neural networks, https://www.biorxiv.org/content/10.1101/366526v2.full-text