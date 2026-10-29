version 1.0
# samples_tsv columns (no header): sample_id  cram  crai  sex(M|F)
import "str_genotyper.wdl" as STR

workflow MultiSampleSTRGenotyper {
    input {
        File samples_tsv
        File reference_fasta
        File reference_fasta_index
        File variant_catalog_female
        File variant_catalog_male
        Int threads = 8
        Int? memory_gb_override
        Boolean output_detailed_csv = true
        Boolean output_summary_csv = true
    }
    Array[Array[String]] sample_rows = read_tsv(samples_tsv)

    scatter (row in sample_rows) {
        String sample_id = row[0]
        File catalog = if row[3] == "M" then variant_catalog_male else variant_catalog_female
        call STR.STRGenotyper as Genotype {
            input:
                base_file_name = sample_id,
                cram_file = row[1],
                cram_index = row[2],
                reference_fasta = reference_fasta,
                reference_fasta_index = reference_fasta_index,
                variant_catalog = catalog,
                threads = threads,
                memory_gb_override = memory_gb_override,
                output_detailed_csv = output_detailed_csv,
                output_summary_csv = output_summary_csv
        }
    }
    output {
        Array[File] genotypes_vcf = Genotype.genotypes_vcf
        Array[File] genotypes_vcf_index = Genotype.genotypes_vcf_index
        Array[File] genotypes_bed = Genotype.genotypes_bed
        Array[Array[File]] summary_csv_files = Genotype.summary_csv_files
        Array[Array[File]] detailed_csv_files = Genotype.detailed_csv_files
    }
}