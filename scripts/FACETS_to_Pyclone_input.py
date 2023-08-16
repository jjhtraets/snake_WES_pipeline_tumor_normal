# https://github.com/vanallenlab/PyCloneTSVGeneration_FACETS_or_TITAN/blob/master/generate_wgs_pyclone_input.py
# lazy, using only Facets, rest not working

from argparse import ArgumentParser
import os
import pandas as pd
import sys
from collections import defaultdict
import numpy as np

DEFAULT_MAJOR_COPY_NUMBER = 1
DEFAULT_MINOR_COPY_NUMBER = 1
DEFAULT_CELL_COPY_NUMBER = 2


##

#Mut_list = ["Missense_Mutation", "Frame_Shift_Del", "In_Frame_Del", "Nonsense_Mutation", "Splice_Site", "In_Frame_Ins"]
Mut_list = ["Frame_Shift_Del", "Frame_Shift_Ins", "Splice_Site", "Translation_Start_Site","Nonsense_Mutation", "Nonstop_Mutation", "In_Frame_Del","In_Frame_Ins", "Missense_Mutation"]

class FacetsCols:
    CHR = "chrom"
    START = "start"
    END = "end"
    TOTAL_CN = "tcn.em"
    MINOR_CN = "lcn.em"
    PURITY = "Purity"


class MafCols:
    START = "Start_Position"
    END = "End_Position"
    CHR = "Chromosome"
    GENE = "Hugo_Symbol"
    MUT_TYPE = "Variant_Classification"
    REF = "Reference_Allele"
    ALT = "Tumor_Seq_Allele2"
    REF_COUNT = "t_ref_count"
    ALT_COUNT = "t_alt_count"

    @classmethod
    def all(cls):
        return [cls.GENE, cls.CHR, cls.START, cls.END, cls.MUT_TYPE, cls.REF, cls.ALT,
                cls.REF_COUNT, cls.ALT_COUNT]

def mutation_id(row):
    """Generate a unique mutation ID"""
    return 'chr{}:{}({}):{}'.format(row[MafCols.CHR], row[MafCols.START], row[MafCols.GENE], row[MafCols.MUT_TYPE])


def center_position(row):
    """Find the point halfway between the start and end of the mutation -- relevant for SNVs and indels"""
    start_int = int(row[MafCols.START])
    end_int = int(row[MafCols.END])
    length = end_int - start_int
    return end_int - (length / 2)

## unhash them when adding strelka_indel

# def combine_snv_and_indel_mafs(snv_maf_df, indel_maf_df):
#     """Combine data from snv and indel maf dataframes into one dataframe"""
#     snv_df = snv_maf_df[MafCols.all()]
#     subset_indels = indel_maf_df[MafCols.all()]
#     combined_snvs_and_indels_df = pd.concat([snv_df, subset_indels])

#     sys.stdout.write("Adding mutation ID...\n")
#     combined_snvs_and_indels_df['mutation_id'] = combined_snvs_and_indels_df.apply(mutation_id, axis=1)
#     sys.stdout.write("Added mutation ID...\n")

#     sys.stdout.write("Adding section center...\n")
#     combined_snvs_and_indels_df['Center_Position'] = combined_snvs_and_indels_df.apply(center_position, axis=1)
#     sys.stdout.write("Added section center...\n")
#     return combined_snvs_and_indels_df

## Remove below when adding strelka_indel
def prepare_snv_maf(snv_maf_df):
    """Loading data from snv dataframe"""
    snv_df = snv_maf_df[MafCols.all()]

    sys.stdout.write("Adding mutation ID...\n")
    snv_df['mutation_id'] = snv_df.apply(mutation_id, axis=1)
    sys.stdout.write("Added mutation ID...\n")

    sys.stdout.write("Adding section center...\n")
    snv_df['Center_Position'] = snv_df.apply(center_position, axis=1)
    sys.stdout.write("Added section center...\n")

    sys.stdout.write("Filtering oncogenes only...\n")
    snv_df = snv_df[snv_df["Variant_Classification"].isin(Mut_list)]
    sys.stdout.write("Filtered oncogenes...\n")

    return snv_df

# def total_copy_number_based_on_sex_and_chrom(sex, chrom):
#     """Based on whether looking at a sex chromosome or not, determine the normal cell copy number"""
#     # Figure out normal copy number based on sex if this a sex chromosome
#     if chrom in ['X', 'Y']:
#         if sex == 'female':
#             assert(chrom == 'X')
#             return int(2)
#         elif sex == 'male':
#             return int(1)
#     else:
#         return int(DEFAULT_CELL_COPY_NUMBER)


def build_chrom_map(facet_df, input_type='facets'):
    """Given a dataframe of FACET data, build a dictionary where the chromosomes are the keys and the values
    are the segments"""
    sys.stdout.write("Building chromosome map for CNA data...\n")
    chr_map = defaultdict(lambda: defaultdict(dict))
    for index, row in facet_df.iterrows():
        if input_type == 'facets':
            chrom = str(int(row[FacetsCols.CHR]))
            try:
                minor_cn = int(row[FacetsCols.MINOR_CN])
            except:
                minor_cn = int(1)
            total_cn = int(row[FacetsCols.TOTAL_CN])
            major_cn = total_cn - minor_cn
            purity = row[FacetsCols.PURITY]
            start = row[FacetsCols.START]   
            end = row[FacetsCols.END]

            chr_map[chrom][start] = {"end": end,
                                        "major_cn": major_cn,
                                        "minor_cn": minor_cn,
                                        "normal_cn": total_cn,
                                        "purity": purity
                                        }

            # if not (chrom == 'Y' and sex == 'female'):
            #     chr_map[chrom][start] = {"end": end,
            #                             "major_cn": major_cn,
            #                             "minor_cn": minor_cn,
            #                             "normal_cn": total_copy_number_based_on_sex_and_chrom(sex, chrom),
            #                             "purity": purity
            #                             }
    return chr_map



def add_cn_info_to_indel_snv_maf(indel_snv_maf, chrom_map, sample_id):
    """Add the facet allelic copy number information to the combined snv indel maf"""
    all_chrs = chrom_map.keys()
    # To extract purity
    chr_list = list(all_chrs)
    start_list = list(chrom_map[chr_list[0]].keys())
    purity = chrom_map[chr_list[0]][start_list[0]]["purity"]

    def major_and_minor_cns(combined_maf_row):
        # For all rows where the chromosome isn't even in the search tree, set the minor and major alleles to the
        # default, non-aneuploidy values
        chrom = combined_maf_row[MafCols.CHR]
        # print(type(chrom))
        position = combined_maf_row["Center_Position"]
        if chrom in all_chrs:
            starts = list(chrom_map[chrom].keys())
            start_index = np.searchsorted(starts, position)
            relevant_start = starts[start_index - 1]
            seg = chrom_map[chrom][relevant_start]
            if seg.get("end") >= position:
                return sample_id, seg.get("major_cn"), seg.get("minor_cn"), seg.get("normal_cn"), purity
            else:
                return sample_id, int(DEFAULT_MAJOR_COPY_NUMBER), int(DEFAULT_MINOR_COPY_NUMBER), int(DEFAULT_CELL_COPY_NUMBER), purity        
        else:
            return sample_id, int(DEFAULT_MAJOR_COPY_NUMBER), int(DEFAULT_MINOR_COPY_NUMBER), int(DEFAULT_CELL_COPY_NUMBER), purity
            


    cn_info = indel_snv_maf.apply(major_and_minor_cns, axis=1)

    # Get a series with a major_cn and minor_cn column of out this Series of tuples
    cn_info_cols = cn_info.apply(pd.Series, index=["sample_id","major_cn", "minor_cn", "normal_cn", "tumour_content"])

    cn_info_added_to_snv_and_maf_df = pd.concat([indel_snv_maf, cn_info_cols], axis=1)
    return cn_info_added_to_snv_and_maf_df


def get_sex(snv_maf):
    """Based on whether both a Y chromosome is present, make a sex determination"""
    chromosomes = snv_maf['Chromosome'].unique()
    y_in_sample = 'Y' in chromosomes
    if y_in_sample:
        sys.stdout.write('Y chromosome found in sample SNV file. Determining male sex. \n')
        return 'male'
    else:
        sys.stdout.write('Y chromosome not found in sample SNV. Determining female sex. Will ignore any Y chromosome '
                         'information from other files. \n')
        return 'female'


def main():
    parser = ArgumentParser(description='Generate .tsv input for PyClone from SNV maf,'
                                        ' GATK indel file, and TITAN output')
    parser.add_argument('snv_maf', metavar='snv_maf', type=str)
#    parser.add_argument('indel_maf', metavar='indel_maf', type=str)
#    parser.add_argument('--titan_output', metavar='titan_output', type=str)
    parser.add_argument('--facets_output', metavar='facets_output', type=str)
    parser.add_argument('output_dir', metavar='output_dir', type=str)
    parser.add_argument('--sample_id', metavar='handle', type=str)
    args = parser.parse_args()

    snv_maf = args.snv_maf

## unhash below when adding strelka_indel
#    indel_maf = args.indel_maf

    if args.facets_output:
        cn_output = args.facets_output
        input_type = 'facets'
    else:
        sys.exit('A facets output file must be provided with allelic copy number information')

    handle = args.sample_id
    if not handle:
        # If not analysis handle is provided, simply name the file after the titan output file prefix
        handle = os.path.split(cn_output)[-1].split('.')[0]

    output_dir = args.output_dir

    snv_maf_df = pd.read_csv(snv_maf, delimiter='\t', comment='#', header='infer')

    # sex = get_sex(snv_maf_df)
    
## unhash below when adding strelka_indel
#    indel_maf_df = pd.read_csv(indel_maf, delimiter='\t', comment='#',header='infer')

    # if sex == 'female':
    #     # if there are no SNVs on the Y chromosome then let's ignore the data on Y chromosome for indels as well
    #     indel_maf_df = indel_maf_df[indel_maf_df.Chromosome != 'Y']

    ## unhash below when adding strelka_indel
    # sys.stdout.write("Combining SNV and indel maf file data...\n")
    snv_maf = prepare_snv_maf(snv_maf_df)

    sys.stdout.write("Loading CNA data...\n")
    cn_df = pd.read_csv(cn_output, delimiter=',', comment='#', header='infer')

    # Build a data structure that allows for easier searching across the CNA data
    chrom_map = build_chrom_map(cn_df, input_type)

    sys.stdout.write("Merging copy number, SNV, and indel information...\n")

    ## Change input snv_maf_df into indel_snv_maf when adding strelka_indel input
    final_df = add_cn_info_to_indel_snv_maf(snv_maf, chrom_map, handle)

    final_output_filepath = os.path.join(output_dir, '{}_pyclone_maf_filtered_ready.tsv'.format(handle))
    sys.stdout.write("Writing output to {}...\n".format(final_output_filepath))

    #sys.stdout.write("Filtering homozygous deletion sites out of output (major cn == 0)...\n")
    #sys.stdout.write("Length before filtering: {}\n".format(len(final_df)))
    final_df = final_df[final_df['major_cn'] > 0]
    sys.stdout.write("Length: {}\n".format(len(final_df)))

    final_df[MafCols.REF_COUNT] = final_df[MafCols.REF_COUNT].fillna(0)
    final_df[MafCols.REF_COUNT] = final_df[MafCols.REF_COUNT].astype(int)
    final_df[MafCols.ALT_COUNT] = final_df[MafCols.ALT_COUNT].fillna(0)
    final_df[MafCols.ALT_COUNT] = final_df[MafCols.ALT_COUNT].astype(int)

    final_df.to_csv(final_output_filepath, sep='\t',
                    columns=['mutation_id', 'sample_id', MafCols.REF_COUNT, MafCols.ALT_COUNT, 'normal_cn', 'major_cn', 'minor_cn', 'tumour_content'],
                    header=['mutation_id', 'sample_id', 'ref_counts', 'var_counts', 'normal_cn', 'major_cn', 'minor_cn', 'tumour_content'],
                    index=False)

    sys.stdout.write("Done!\n")

if __name__ == '__main__':
    main()

