#!/opt/common/CentOS_6-dev/R/R-3.4.1/bin/Rscript

suppressPackageStartupMessages({
          library(data.table)
          library(dplyr)
          library(stringr)
          library(purrr)
          library(tidyr)
          library(readr)
          library(jsonlite)
          library(Hmisc)
})
"%nin%" <- Negate("%in%")
#suppressMessages(devtools::source_url('https://raw.githubusercontent.com/kpjonsson/kpjmisc/master/R/hotspot_annotate_maf.R'))
source("./hotspot_annotate_maf.R")
#source('ngs_recurrent_fp_annotate_maf.R')

# Write output ----------------------------------------------------------------------------------------------------
   write_out = function (x, file, col.names = T) {
            write.table(x, file = file, sep = "\t", quote = F, row.names = F, col.names = col.names)
       }

# Harcoded paths
# hotspot_files = c('/ifs/work/taylorlab/jonssonp/hotspots/24k/hotspots_24k_FULL.txt',
#                  '/ifs/work/taylorlab/jonssonp/hotspots/nbt/hotspots_NBT_FULL.txt',
#                  '/ifs/work/taylorlab/jonssonp/hotspots/3d/hotspots.txt')
oncokb_file = '/ifs/work/taylorlab/jonssonp/oncokb/oncokb-all-variants-2018-03-06.txt'
repeatmasker_file = '/home/jonssonp/git/ngs-filters/data/rmsk_mod.bed'
blacklist_file = '/home/jonssonp/git/ngs-filters/data/wgEncodeDacMapabilityConsensusExcludable.bed'
broad_pon_file = '/ifs/res/taylorlab/chavans/WES_QC_filters/annot/Broad_PoN_distinct_by_allele_5cols_with_tag.txt'


coding_mutations = c('Frame_Shift_Del',
                     'Frame_Shift_Ins',
                     'In_Frame_Del',
                     'In_Frame_Ins',
                     'Missense_Mutation',
                     'Nonsense_Mutation',
                     'Nonstop_Mutation',
                     'Splice_Site',
                     'Targeted_Region',
                     'Translation_Start_Site')

truncating_mutations = c('Nonsense_Mutation',
                         'Frame_Shift_Ins',
                         'Frame_Shift_Del',
                         'Splice_Site',
                         'Nonstop_Mutation')

args = commandArgs(TRUE)

if (is.null(args) | length(args)<1) {
    cat("Usage: maf-postprocess.R input.maf\n")
    quit()
}

maf_file = args[1]


outname_flagged = gsub('.maf$', '.postprocessed.maf', maf_file)
outname_filter = gsub('.maf$', '.postprocessed.filter.maf', maf_file)
outname_recurr = gsub('.maf$', '.postprocessed.recurr.maf', maf_file)
maf = fread(maf_file) %>%
    replace_na(list(MUTECT = 0,
                    PINDEL = 0,
                    VARDICT = 0))


if ('t_var_freq' %nin% names(maf)) maf = mutate(maf, t_var_freq = t_alt_count/t_depth)

# Annotate hotspots -----------------------------------------------------------------------------------------------
maf = hotspot_annotate_maf(maf)

# Annotate with OncoKB  -------------------------------------------------------------------------------------------
oncokb = fread(oncokb_file) %>%
    select(Gene, Alteration, Oncogenicity, Mutation_Effect = `Mutation Effect`) %>%
    filter(Alteration %nin% c('Amplification', 'Deletion'),
           !Alteration %like% 'Fusion')

maf = mutate(maf, oncokb_tag = str_replace(HGVSp_Short, 'p.', '')) %>%
    left_join(., oncokb, by = c('Hugo_Symbol' = 'Gene', 'oncokb_tag' = 'Alteration')) %>%
    select(-oncokb_tag)

# Blacklist variants ----------------------------------------------------------------------------------------------
repeatmasker = fread(repeatmasker_file, header = F,
                     col.names = c('Chromosome', 'Start_Position', 'End_Position', 'type')) %>%
    mutate(Chromosome = str_replace(Chromosome, 'chr', ''))
blacklist = fread(blacklist_file, header = F,
                  col.names = c('Chromosome', 'Start_Position', 'End_Position', 'type', 'na1', 'na2')) %>%
    mutate(Chromosome = str_replace(Chromosome, 'chr', ''))

data.table::setkey(data.table::setDT(repeatmasker),
                   Chromosome, Start_Position, End_Position)
data.table::setkey(data.table::setDT(blacklist),
                   Chromosome, Start_Position, End_Position)
data.table::setkey(data.table::setDT(maf),
                   Chromosome, Start_Position, End_Position)

repeatmasker_overlap = data.table::foverlaps(maf[, .(Chromosome, Start_Position, End_Position)],
                                             repeatmasker, type = 'any', mult = 'first')

blacklist_overlap = data.table::foverlaps(maf[, .(Chromosome, Start_Position, End_Position)],
                                          blacklist, type = 'any', mult = 'first')

maf = maf[, repeat_masker := repeatmasker_overlap$type]
maf = maf[, blacklist_region := blacklist_overlap$type]

maf = replace_na(maf, list(repeat_masker = '', blacklist_region = ''))

# Annotate with Broad PoN -----------------------------------------------------------------------------------------
broad_pon = fread(broad_pon_file)
maf = mutate(maf, pon_tag = str_c(Chromosome,
                                  Start_Position,
                                  End_Position,
                                  Reference_Allele,
                                  Tumor_Seq_Allele2,
                                  sep = ':'),
             Broad_PoN = pon_tag %in% broad_pon$TAG)

# Annotate with recurrent filter ----------------------------------------------------------------------------------

#maf = ngs_recurrent_fp(maf)
#write_out(maf, outname_recurr)

# Apply filters ---------------------------------------------------------------------------------------------------
maf = mutate(maf, remove = (t_var_freq < .05 |
                            t_alt_count < 3 |
                            t_depth < 20 |
                            n_depth < 10 |
                            n_alt_count > 3 |
                            FILTER != 'PASS' |
                            blacklist_region != '' |
                            repeat_masker != '' |
                            #is_recurrent == T |
                            (Broad_PoN == T & t_var_freq < .15)),
             whitelist = ((Hotspot == T & sum(snv_hotspot, indel_hotspot) > 0) | # don't whitelist 3D hotspots, too many FPs
                         (Oncogenicity %like% 'Oncogenic' & Variant_Classification %nin% truncating_mutations)))

maf_filter = filter(maf, !(remove == T & whitelist == F) & # removes everything tagged by "remove" unless whitelisted
                         !(Hotspot == T & t_var_freq < .02)) # only retain hotspots if VAF >= .03

write_out(maf, outname_flagged)
write_out(maf_filter, outname_filter)
