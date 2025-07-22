nextflow.enable.dsl = 2

/* ─── Importa el proceso desde el módulo ───────────────────────── */
include { DGE_MATRIX } from './modules/DGE_matrix/dropseq.nf'

/* ─── Flujo principal ──────────────────────────────────────────── */
workflow {

    /*
     * Ejecutar solo si se invoca con --cmd dge_matrix
     *   Ejemplo:  nextflow run main.nf --cmd dge_matrix
     */
    if( params.cmd == 'dge_matrix' ) {

        /*
         * Canal que emite:  (read1, read2, gtf, fasta, dict, sample_id)
         */
        Channel
            /* 1️⃣  toma todos los FASTQ de lectura 1 */
            .fromPath('datasets/*_1.fastq.gz')
            /* 2️⃣  empareja con su correspondiente lectura 2 */
            .map { Path r1 ->

                /* sample_id = nombre sin “_1.fastq.gz” */
                def sample_id = r1.name.replaceAll(/_1\.fastq\.gz$/, '')

                /* Path del FASTQ de lectura 2 */
                def r2 = file("datasets/${sample_id}_2.fastq.gz")

                /* chequeo: si falta el R2, aborta con mensaje claro */
                if( !r2.exists() )
                    error "No se encontró el R2 para ${sample_id}: ${r2}"

                /* devuelve el tuple que consumirá el proceso */
                tuple(
                    r1,                              // read1 Path
                    r2,                              // read2 Path
                    file('genomas/genes.gtf'),       // GTF
                    file('genomas/genome.fa'),       // FASTA
                    file('genomas/genome.dict'),     // diccionario
                    sample_id                        // id
                )
            }
            /* 3️⃣  DEBUG opcional: muestra qué se está enviando */
            .view { t -> println "[DEBUG] ${t[5]} ⇒ ${t[0].name}, ${t[1].name}" }
            /* 4️⃣  asigna el canal a la variable samples */
            .set { samples }

        /* 5️⃣  lanza el proceso DGE_MATRIX */
        samples | DGE_MATRIX
    }
}



















