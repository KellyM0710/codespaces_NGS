# NGS Hackathon — Informe

**Universidad Publica de Navarra**  
**Master en Metodos Computacionales en Ciencias**  
**Modulo:** Analisis e interpretacion de datos de alto rendimiento  
**Nombre:** Kelly Mendoza   

---

## Reportes.

Todos los reportes de calidad generados durante este analisis estan disponibles de forma interactiva en los siguientes links:

| Reporte | Link |
|---|---|
| MultiQC — comparacion completa de todos los alineamientos | [Ver reporte](https://kellym0710.github.io/codespaces_NGS/multiqc_report.html) |
| FastQC — Negative.fq (antes del trimming) | [Ver reporte](https://kellym0710.github.io/codespaces_NGS/figures/trimmed_Negative_fastqc.html) |
| FastQC — Negative.fq (despues del trimming) | [Ver reporte](https://kellym0710.github.io/codespaces_NGS/figures/Negative_trimN_fastqc.html) |
| FastQC — BQ.fq (antes del trimming) | [Ver reporte](https://kellym0710.github.io/codespaces_NGS/figures/trimmed_BQ_fastqc.html) |
| FastQC — BQ.fq (despues del trimming) | [Ver reporte](https://kellym0710.github.io/codespaces_NGS/figures/BQ_trimN_fastqc.html) |
| FastQC — BQ.fq optimizado (resultado final Q3) | [Ver reporte](https://kellym0710.github.io/codespaces_NGS/figures/BQ_opt1_fastqc.html) |

---

# Pregunta 1 — Alineacion de Negative.fq

## Problematica encontrada.

Al intentar alinear el archivo `trimmed_Negative.fq` contra el genoma de referencia con Bowtie2, se observo que ninguna de las 1.076.320 lecturas consiguio mapearse, dando una tasa de alineamiento del 0%.

```bash
bowtie2 -x genomes/AFPN02.1/AFPN02.1_merge \
  -U fastq/trimmed_Negative.fq \
  -S results_NGS1/Negative.sam \
  2> results_NGS1/Negative_bowtie_stats.txt
```

```
1076320 reads; of these:
  1076320 (100.00%) aligned 0 times
  0 (0.00%) aligned exactly 1 time
0.00% overall alignment rate
```

---

## Inspeccionando las lecturas

Primero, se observo como eran las lecturas por dentro, al ejecutar `head` sobre el archivo, el patron se hizo visible inmediatamente:

```bash
head -20 fastq/trimmed_Negative.fq
```

```
@AFPN02.1_merge-1076320
NNNNNAAGGTCGCTAATCTCTTTACGCAGATTTTTTATTCCTTCAACTAACAGCGNNNN
+
#####BACCCGFGGGGGGGGGG1GGEGGG0GBGG1GGGGGGGGBCG####
```

Todas las lecturas; tenian, 5 N's al inicio y 4 N's al final. Las N's son bases ambiguas, posiciones donde el secuenciador no fue capaz de identificar que nucleotido habia. En este caso, no eran errores de secuenciacion sino algo introducido durante el procesamiento de los datos.

> 📊 [Ver reporte FastQC de Negative.fq antes del trimming](https://kellym0710.github.io/codespaces_NGS/figures/trimmed_Negative_fastqc.html)

---

## Hipotesis propuestas como solucion a la alineacion. 

**Hipotesis 1: Las N's son el problema**
Se infiere que durante el demultiplexado, las secuencias de barcode (usadas para identificar cada muestra) fueron enmascaradas con N's en lugar de eliminarse. Bowtie2 no sabe que hacer con bases ambiguas y simplemente rechaza esas lecturas.
**Solucion propuesta:** eliminar las N's con `cutadapt --trim-n`.

**Hipotesis 2: Los parametros de Bowtie2 son demasiado estrictos**
Tambien se puede dar el hecho de que quizas el alineador en su modo por defecto es demasiado exigente, y usando `--very-sensitive` podria conseguir alinear las lecturas igualmente.
**Solucion propuesta:** repetir el alineamiento con `--very-sensitive` sin tocar las lecturas.

---

## Prueba de las dos hipotesis.

### Hipotesis 1: Eliminar las N's con cutadapt.

```bash
cutadapt --trim-n -m 30 \
  -o results_Q1/Negative_trimN.fq \
  fastq/trimmed_Negative.fq \
  > results_Q1/Negative_trimN_report.txt 2>&1
```

Resultado: Se elimino aproximadamente el 15% de todas las bases, que eran N's. Al volver a alinear:

```bash
bowtie2 -x genomes/AFPN02.1/AFPN02.1_merge \
  -U results_Q1/Negative_trimN.fq \
  -S results_Q1/Negative_trimN.sam \
  2> results_Q1/Negative_trimN_bowtie_stats.txt
```

```
1076320 reads; of these:
  308 (0.03%) aligned 0 times
  1020259 (94.79%) aligned exactly 1 time
  55753 (5.18%) aligned >1 times
99.97% overall alignment rate
```

Al volver alinear, se paso de un 0% a un **99.97%** con un solo cambio. Los flagstat confirmaron que 1.076.012 lecturas quedaron mapeadas correctamente (los valores de "paired" aparecen a 0 porque son lecturas single-end, lo cual es completamente normal).

> 📊 [Ver reporte FastQC de Negative.fq despues del trimming](https://kellym0710.github.io/codespaces_NGS/figures/Negative_trimN_fastqc.html)

---

### Hipotesis 2: Modo very-sensitive sin tocar las lecturas.

```bash
bowtie2 --very-sensitive \
  -x genomes/AFPN02.1/AFPN02.1_merge \
  -U fastq/trimmed_Negative.fq \
  -S results_Q1/Negative_verysensitive.sam \
  2> results_Q1/Negative_verysensitive_bowtie_stats.txt
```

```
1076320 reads; of these:
  1076320 (100.00%) aligned 0 times
0.00% overall alignment rate
```

En el caso de esta segunda hipotesis, cambiar la sensibilidad del alineador del tuvo ningun ejecto, el problema de alineacion es exactamente igual que antes. 

---

## Comparacion de los resultados de las hipotesis.

| Condicion | Alineamiento | unicas | Multimapping | No alineadas |
|---|---|---|---|---|
| Original (sin cambios) | 0,00% | — | — | 100% |
| **H1: eliminando N's** | **99,97%** | 94,79% | 5,18% | 0,03% |
| H2: modo very-sensitive | 0,00% | — | — | 100% |

> 📊 [Ver comparacion completa en MultiQC](https://kellym0710.github.io/codespaces_NGS/multiqc_report.html)

---

## Conclusion.

Al evaluar el resultado en ambas hipotesis, se determino que la Hipotesis 1 es la correcta, ya que, eliminar las N's transformo un alineamiento fallido en uno casi perfecto, mientras que cambiar los parametros de Bowtie2 no se mostro ningun cambio.

Por consiguiente, se entiende que, lo que ocurrio es que durante el demultiplexado de la libreria multiplexada, las secuencias de barcode (situadas en posiciones fijas al inicio y al final de cada lectura) fueron reemplazadas por N's en lugar de recortarse. Es un artefacto tecnico del pipeline de procesamiento, no un problema con el ADN en si.

Desde el punto de vista biologico, esto es importante, los datos de secuenciacion subyacentes eran de buena calidad. El 94,79% de lecturas unicas y el bajo porcentaje de multimapping tras el trimming confirman que las lecturas provienen genuinamente del organismo de referencia (`AFPN02.1`). El problema era puramente computacional, un paso de procesamiento mal ejecutado que se habria pasado por alto sin inspeccionar las lecturas directamente.

---

# Pregunta 2: Mejorando el alineamiento de BQ.fq.

## Problematica encontrada.

Al igual que con Negative.fq, el archivo `trimmed_BQ.fq` tampoco conseguia alinearse, arrojando una tasa de alineamiento del 0% en 1.086.246 lecturas. Pero al inspeccionar las lecturas, se observo una problematica diferente, encontrandose dos problemas a la vez:

```bash
head -20 fastq/trimmed_BQ.fq
```

```
@chr1-2
NNNNNGTGAAAGAAAAGAAGGAAGAAATATCTGAATTAAGTGTCATCAGGTACAGNNNN
+
#####77777777777777777777777777777777771777777777777777####
```

1. Las N's en los extremos que vimos en Negative.fq.
2. Las puntuaciones de calidad son demasiado uniformes, es decir, se observo que todos los valores son 7, lo cual es raro; ya que, en una secuenciacion normal, la calidad varia a lo largo de la lectura y suele bajar hacia el final, un patron tan plano y homogeneo sugiere que toda la secuencia tuvo problemas, dando a entender porque el bad quality y la falta de alineacion.

> 📊 [Ver reporte FastQC de BQ.fq — calidad baja uniforme](https://kellym0710.github.io/codespaces_NGS/figures/trimmed_BQ_fastqc.html)

---

## Hipotesis propuestas. 

**Hipotesis 1:** Solo eliminar las N's, igual que en el ejercicio anterio, asumiendo que que en este caso tambien podria funcionar. 

**Hipotesis 2:** Eliminar N's y ademas filtrar por calidad (umbral Phred 20). Entendiciendose de que si la baja calidad de las bases tambien esta afectando al alineamiento, hay que limpiar mas a fondo.

---

## Resultados

### Hipotesis 1 — Solo N's

```bash
cutadapt --trim-n -m 30 \
  -o results_Q2/BQ_trimN.fq \
  fastq/trimmed_BQ.fq \
  2> results_Q2/BQ_trimN_report.txt
```

```bash
bowtie2 -x genomes/AFPN02.1/AFPN02.1_merge \
  -U results_Q2/BQ_trimN.fq \
  -S results_Q2/BQ_trimN.sam \
  2> results_Q2/BQ_trimN_bowtie.txt
```

```
82.82% overall alignment rate — 899.607 lecturas mapeadas
```

> 📊 [Ver FastQC de BQ.fq despues del trimming de N's](https://kellym0710.github.io/codespaces_NGS/figures/BQ_trimN_fastqc.html)

---

### Hipotesis 2 — N's + filtro de calidad Phred 20

```bash
cutadapt --trim-n -q 20 -m 30 \
  -o results_Q2/BQ_trimNQ.fq \
  fastq/trimmed_BQ.fq \
  2> results_Q2/BQ_trimNQ_report.txt
```

En este caso, se observo que el filtro de calidad fue demasiado agresivo, como la baja calidad no estaba solo al final de las lecturas sino en toda su longitud, al recortar las bases de baja calidad las lecturas quedaban tan cortas que eran descartadas por el filtro de longitud minima.

```
Reads that were too short: 763.613 (70.3%)  ← se perdio el 70% de las lecturas
Reads written (passing):   322.633 (29.7%)
```

```bash
bowtie2 -x genomes/AFPN02.1/AFPN02.1_merge \
  -U results_Q2/BQ_trimNQ.fq \
  -S results_Q2/BQ_trimNQ.sam \
  2> results_Q2/BQ_trimNQ_bowtie.txt
```

```
87.93% overall alignment rate — pero solo 283.685 lecturas mapeadas
```

---

## Tabla comparativa de las hipotesis. 

| Condicion | Lecturas conservadas | Alineamiento | Total mapeadas |
|---|---|---|---|
| Original | 1.086.246 (100%) | 0,00% | 0 |
| **H1: solo N's** | **1.086.246 (100%)** | **82,82%** | **899.607** |
| H2: N's + Phred 20 | 322.633 (29,7%) | 87,93% | 283.685 |

> 📊 [Ver comparacion en MultiQC](https://kellym0710.github.io/codespaces_NGS/multiqc_report.html)

---

## Conclusion. 

En este caso se presenta un problema entre cantidad y calidad, es decir, H2 tiene un porcentaje de alineamiento un poco mayor, sin embargo esta descarta el 70% de los datos, obteniendo muchas menos lecturas mapeadas que con H1. Entonces, la mejor opcion en este caso, seria la H1, donde la calidad baja es uniforme en todas las posiciones, lo que significa que no es un artefacto localizado que se deba eliminar de manera agresiva es una caracristica de toda la corrida de secuenciacion, descartar el 70% de los datos por un problema sistematico empeoraria la cobertura del genoma y dificultaria analisis posteriores como la llamada de variantes. 

Se puede inferir que el bad quality de este conjunto de datos, posiblemente se deba por degradacion de reactivos, problemas en la celda de flujo o un error de calibracion del secuenciador. El patron de calidad uniforme y bajo, sin el reconocible declive hacia el extremo 3', apunta a un problema global y no a errores puntuales. A pesar de esto, el 82,82% de las lecturas si contienen suficiente informacion como para alinearse correctamente. Esto demuestra que datos de baja calidad no son necesariamente datos inutilisables, con el procesamiento adecuado todavia pueden proporcionar informacion biologica valiosa sobre el organismo secuenciado.

---

# Pregunta 3 — Hackathon: consiguiendo el mejor resultado posible.

## Estrategia.
De acuerdo a lo realizado anteriormente en los ejercicios anteriores, el objetivo en este caso es tratar de alinear lo mas posible el BQ.fq. En el ejercicio dos, se descubrio que filtrando la calidad se pueden destruir datos y a su vez combinado con el modo de alineamiento de Bowtie2; sin embargo, tambien se comprendio que al ser demasiado agresivos, podemos eliminar la mayoria de los datos. Entonces se planteo tres combinaciones de prueba:

| Intento | cutadapt | Bowtie2 | Lecturas conservadas | Alineamiento | Total mapeadas |
|---|---|---|---|---|---|
| Linea base (Q2-H1) | `--trim-n` | defecto | 100% | 82,82% | 899.607 |
| Q2-H2 | `--trim-n -q 20` | defecto | 29,7% | 87,93% | 283.685 |
| **Opt1** | **`--trim-n -q 10`** | **`--very-sensitive`** | **100%** | **89,43%** | **971.428** |
| Opt2 | `--trim-n -q 15` | `--very-sensitive` | 100% | 89,41% | 971.259 |
| Opt3 | `--trim-n` | `--very-sensitive` | 100% | 89,43% | 971.437 |

---

## BOX Final

```bash
# Paso 1 — Trimming
cutadapt --trim-n -q 10 -m 30 \
  -o results_Q3/BQ_opt1.fq \
  fastq/trimmed_BQ.fq \
  2> results_Q3/BQ_opt1_cutadapt.txt

# Paso 2 — Alineamiento
bowtie2 --very-sensitive \
  -x genomes/AFPN02.1/AFPN02.1_merge \
  -U results_Q3/BQ_opt1.fq \
  -S results_Q3/BQ_opt1.sam \
  2> results_Q3/BQ_opt1_bowtie.txt

# Paso 3 — Estadisticas
samtools flagstat results_Q3/BQ_opt1.sam > results_Q3/BQ_opt1_flagstat.txt
```

**Resultado final:**
```
Total lecturas:          1.086.246
Lecturas conservadas:    1.086.246 (100%)
Tasa de alineamiento:    89,43%
Total lecturas mapeadas: 971.428
```

> 📊 [Ver FastQC del resultado final optimizado](https://kellym0710.github.io/codespaces_NGS/figures/BQ_opt1_fastqc.html)
> 📊 [Ver comparacion completa en MultiQC](https://kellym0710.github.io/codespaces_NGS/multiqc_report.html)

---

## Conclusion. 

Al comparar los resultados obtenidos, se puede destacar que, Opt1, Opt2 y Opt3 dan practicamente el mismo resultado, esto indica que el factor determinante fue el modo `--very-sensitive` de Bowtie2, no el umbral de calidad del trimming. El alineador en modo estandar simplemente no es capaz de colocar estas lecturas de baja calidad, pero al darle mas intentos y permitirle ser mas flexible con los mismatches, consigue mapear casi un 7% mas de lecturas. Por consiguiente, se elige Opt1 como solucion final porque, aunque el umbral `-q 10` tiene un efecto minimo en este dataset, añadir un filtro de calidad ligero siempre es buena practica; ya que, elimina las bases genuinamente malas sin sacrificar datos.

El hecho de que `--very-sensitive` mejore tanto el resultado confirma que las lecturas BQ, aunque sea de baja calidad, contienen secuencia genomica real de `AFPN02.1`. No es solo ruido, sino una señal degradada. El alineador necesita mas flexibilidad para reconocerlas, pero una vez que se le da esa flexibilidad, las coloca correctamente en el genoma.

# Pregunta 4 — Explorando el archivo de alineamiento.

## Preparacion.

Se convierte el SAM final a BAM ordenado e indexado para poder trabajar con el de forma eficiente:

```bash
samtools view -bS results_Q3/BQ_opt1.sam | \
  samtools sort -o results_Q3/BQ_aligned.bam

samtools index results_Q3/BQ_aligned.bam
```

```
Total: 1.086.246 lecturas
Mapeadas: 971.428 (89,43%)
No mapeadas: 114.818 (10,57%)
```
## a) Lecturas unicas vs lecturas multimapping

En Bowtie2, cuando una lectura puede alinearse en mas de un sitio del genoma, se le asigna una calidad de mapeo (MAPQ) de 0, indicando que no hay certeza sobre su posicion correcta. Las lecturas que mapean en un unico lugar tienen MAPQ mayor o igual a 1.

```bash
# Lecturas unicas (MAPQ >= 1)
samtools view -b -F 4 -q 1 results_Q3/BQ_aligned.bam \
  > results_Q4/BQ_unique.bam

# Lecturas multimapping (mapeadas pero con MAPQ = 0)
samtools view -h -F 4 results_Q3/BQ_aligned.bam | \
  awk '$1~/^@/ || $5==0' | \
  samtools view -b > results_Q4/BQ_multimapping.bam
```

| Categoria | Lecturas | % del total |
|---|---|---|
| unicas (MAPQ >= 1) | 940.924 | 86,62% |
| Multimapping (MAPQ = 0) | 30.504 | 2,81% |
| No mapeadas | 114.818 | 10,57% |

## b) Lecturas mapeadas vs no mapeadas.

El flag SAM numero 4 marca las lecturas no mapeadas. Con `-F 4`, estas se excluyen (nos quedamos con las mapeadas) y con `-f 4` las seleccionamos directamente.

```bash
# Mapeadas
samtools view -b -F 4 results_Q3/BQ_aligned.bam \
  > results_Q4/BQ_mapped.bam

# No mapeadas
samtools view -b -f 4 results_Q3/BQ_aligned.bam \
  > results_Q4/BQ_unmapped.bam
```

```
Mapeadas:     971.428 lecturas (89,43%)
No mapeadas:  114.818 lecturas (10,57%)
```

## c) Directamente desde Bowtie2.

En lugar de filtrar el BAM despues del alineamiento, Bowtie2 permite separar las lecturas durante el propio alineamiento usando el flag `--un`:

```bash
bowtie2 --very-sensitive \
  -x genomes/AFPN02.1/AFPN02.1_merge \
  -U results_Q3/BQ_opt1.fq \
  -S results_Q3/BQ_aligned.sam \
  --un results_Q3/BQ_unmapped_direct.fastq
```

Los flags mas utiles para esto son:
- `--un` → guarda las lecturas no mapeadas en un FASTQ
- `--al` → guarda las lecturas mapeadas en un FASTQ
- `--un-gz` / `--al-gz` → igual pero comprimido

Este enfoque es mas eficiente porque evita tener que leer el archivo alineado una segunda vez. Sin embargo, el filtrado con samtools es mas flexible cuando ya tienes el archivo alineado y quieres aplicar criterios adicionales como el MAPQ o el cromosoma.

## d) Procedencia de las lecturas que no se mapearon. 

Se extrajeron las lecturas no mapeadas y miraron sus nombres:

```bash
samtools fastq results_Q4/BQ_unmapped.bam \
  > results_Q4/BQ_unmapped.fastq

samtools view results_Q4/BQ_unmapped.bam | \
  awk '{print $1}' | sed 's/-[0-9]*$//' | sort -u | head -10
```

```
AFPN02.1_merge
chr1, chr2, chr3 ... chr22, chrX, chrY, chrM
GL000191.1, JH159131.1, KB021644.1 ...
```
La respuesta estaba ya en los propios nombres de las lecturas: `chr1`, `chr2`, etc. Es decir, son los nombres de los cromosomas humanos.Para confirmarlo, se envio las secuencias a BLAST (NCBI nucleotide nr/nt):

```
>chr1-2
GTGAAAGAAAAGAAGGAAGAAATATCTGAATTAAGTGTCATCAGGTACAG
>chr1-1
GTTTGATGAGAATGATGACTTCCAGCTTCATCCATGTCCCTGCAAAGGAC
```

**Resultado BLAST:**
> **AL592436.7** — *Human DNA sequence from clone RP11-710N8 on chromosome 1*
> Organismo: *Homo sapiens* | Identidad: **100%**

El origen de las lecturas no mapeadas quedo completamente confirmado.

| Origen | Lecturas | Descripcion |
|---|---|---|
| AFPN02.1_merge | 104.892 | Lecturas del organismo de referencia que fallaron por baja calidad |
| Cromosomas humanos (chr1-22, X, Y, M) | ~1.400 | Contaminacion con ADN humano |
| Contigs GRCh38 (GL/JH/KB) | ~8.500 | Fragmentos no colocados del genoma humano |

```
Error rate:     4,83%   (vs ~0,1% tipico en NGS)
Longitud media: 49 pb
```

---

## Conclusion. 

Las lecturas no mapeadas pueden representar dos cosas distintas.

La mayoria (104.892 lecturas, ~91%) provienen del propio `AFPN02.1` pero no se consigio alinearse porque su calidad era demasiado baja incluso despues del trimming. Tienen informacion real sobre el organismo, pero la señal estaba demasiado degradada para ser recuperable.

El resto (9.926 lecturas, 9%), son representado por ADN humano, todos los cromosomas estan representados, incluyendo el mitocondrial, con una uniformidad de 56 lecturas por cromosoma. Esta distribucion tan regular no es tipica de una contaminacion accidental de laboratorio, que normalmente afecta a unas regiones mas que a otras. Encontrar contaminacion humana en datos de secuenciacion de otro organismo seria una señal de alarma importante. Identificarla y eliminarla es esencial para obtener resultados fiables, especialmente si se trabaja con muestras clinicas donde la privacidad del paciente tambien entra en juego.

> 📊 [Ver todos los resultados de alineamiento en MultiQC](https://kellym0710.github.io/codespaces_NGS/multiqc_report.html)
