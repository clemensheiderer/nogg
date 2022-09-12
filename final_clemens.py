import pandas as pd
import argparse
import numpy as np
import regex as re
from pandas import *
import os.path


pd.set_option('mode.chained_assignment', None)




def Main():
    parser = argparse.ArgumentParser(description='details',
                                     usage='use "%(prog)s --help" for more information',
                                     formatter_class=argparse.RawTextHelpFormatter)
    # parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--intersection", help="""
    load it with: python3 filename.py eggnog4.species_list.tsv meNOG.members.tsv meNOG.annotations.tsv -i

    Enter 2 mammal species your choice.                                  
    All genes/ProteinIds of first species are calculated, which have at least one homolog in the second species.            
                      """, action="store_true")

    parser.add_argument("-t", "--three_species", help="""
    load it with: python3 filename.py eggnog4.species_list.tsv meNOG.members.tsv meNOG.annotations.tsv -t

    The relative compliment of Homo sapiens in Mu musculus is made, intersection with Pan troglodytes generates 
    all needed indices.
                      """, action="store_true")
    parser.add_argument("-f", "--family_specific_genes", help="""
    load it with: python3 filename.py eggnog4.species_list.tsv meNOG.members.tsv meNOG.annotations.tsv -i -f

    Mus Musculus and Rattus norvegicus share the same Muridae family. 

                      """, action="store_true")
    parser.add_argument('speciesfile', nargs='+')
    parser.add_argument('membersfile', nargs='+')
    parser.add_argument('annotations', nargs='+')

    args = parser.parse_args()

    print("Loaded files in order: \n")
    for file_name1 in args.speciesfile:
        df_sp = pd.read_csv(file_name1, sep='\t',
                            names=['species', 'taxid', 'core/periphery/adherent', 'source', 'sourceversion'])
    print(df_sp.head())

    for file_name2 in args.membersfile:
        df_me = pd.read_csv(file_name2, sep='\t',
                            names=['TaxonomicLevel', 'GroupName', 'ProteinCount', 'SpeciesCount', 'Functional Category',
                                   'TaxonIDProteinID'])
    print(df_me.head())

    for file_name3 in args.annotations:
        df_anno = pd.read_csv(file_name3, sep='\t',
                              names=['GroupName', 'ProteinCount', 'SpeciesCount', 'Functional Category', 'Description'])
    print(df_anno.head())
    # pattern_sp1 = r'[Mm]us [Mm]usculus'
    # df_row_sp1 = df_sp[df_sp['species'].str.contains(pattern_sp1)]
    # print(df_row_sp1)
    if args.intersection:

        # print("intersection")
        ############################
        # species 1
        ############################

        # pattern_sp1 = r'[Mm]us [Mm]usculus'

        print("Species 1: \n")
        pattern_sp1 = input()

        df_row_sp1 = df_sp[df_sp['species'].str.contains(pattern_sp1)]
        print(df_row_sp1)

        with open('defaultdata.txt', 'w') as f:
            print(
                f" {df_row_sp1} \n",
                file=f)

        pattern_row_sp1 = df_row_sp1['taxid']
        # print(pattern_row_sp1)

        vals_row_sp1 = pattern_row_sp1.to_csv(header=None, index=False)  # .strip('\n').split('\n').replace('\r', '')
        # als_row_sp.vals_row_sp.replace('\n, '')
        vals_row_sp1 = vals_row_sp1.strip()
        # print(vals_row_sp1)

        ##############################################################################################################
        ############################
        # species 2
        ############################
        print("\n")
        print("Species 2: \n")
        pattern_sp2 = input()
        # pattern_sp2 = r'Homo sapiens'

        df_row_sp2 = df_sp[df_sp['species'].str.contains(pattern_sp2)]
        print(df_row_sp2)
        with open('defaultdata.txt', 'a+') as f:
            print(
                f" {df_row_sp2} \n",
                file=f)
        pattern_row_sp2 = df_row_sp2['taxid']
        vals_row_sp2 = pattern_row_sp2.to_csv(header=None, index=False)  # .strip('\n').split('\n').replace('\r', '')
        # als_row_sp.vals_row_sp.replace('\n, '')
        vals_row_sp2 = vals_row_sp2.strip()
        # print(vals_row_sp2)

        # vals_row_sp2 += '.[A-Za-z0-9]{7,}'
        # attern_me2 = vals_row_sp2
        pattern_me2 = vals_row_sp2 + '.[.A-Za-z0-9_-]{3,}'
        my_indices_species_2 = df_me.index[df_me['TaxonIDProteinID'].str.contains(pattern_me2)].tolist()

        # vals_row_sp1 += '.[A-Za-z0-9]{7,}'
        # pattern_me1 = vals_row_sp1
        pattern_me1 = vals_row_sp1 + '.[.A-Za-z0-9_-]{3,}'
        my_indices_species_1 = df_me.index[df_me['TaxonIDProteinID'].str.contains(pattern_me1)].tolist()

        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        df_me['TaxonIDProteinID_s1'] = None
        df_me['ProteinID1'] = None
        index_proteinid1 = df_me.columns.get_loc('ProteinID1')
        # print(index_proteinid1)
        # print(df_me.head())
        index_description_me_1 = df_me.columns.get_loc('TaxonIDProteinID')
        index_date_me_1 = df_me.columns.get_loc('TaxonIDProteinID_s1')
        # print(index_description_me_1, index_date_me_1)

        # ha = r'([A-Za-z0-9]{7,})'
        # id_p = r'([.A-Za-z0-9_-]{3,})'
        id_p = r'\b(?<=.)([A-Za-z0-9-_].{4,}).\b'
        length_1 = np.array([])

        for row in range(0, len(df_me)):
            try:
                number = re.findall(pattern_me1, df_me.iat[row, index_description_me_1])
                length_1 = np.append(length_1, len(number))
                # print(number)
                up = []
                if number != []:
                    for nu in number:
                        # print(nu)
                        m = re.search(id_p, nu)
                        # print(m.group())
                        up.append(m.group())
                # print(up)

                # print(up)
            except AttributeError:
                number = re.findall(pattern_me1, df_me.iat[row, index_description_me_1])
                # up = re.findall(pattern, df_me.iat[row, index_description_me_1])
                # df.iat[row, index_proteinid1] = up

                # print(up)

            df_me.iat[row, index_date_me_1] = number
            df_me.iat[row, index_proteinid1] = up
        print(f"genearray of {pattern_sp1}:  {length_1}")
        print(f"sum: {sum(length_1)}")

        # -------------------------------------------------------------------------

        df_me['TaxonIDProteinID_s2'] = None
        df_me['ProteinID2'] = None
        index_proteinid2 = df_me.columns.get_loc('ProteinID2')

        index_description_me_2 = df_me.columns.get_loc('TaxonIDProteinID')
        index_date_me_2 = df_me.columns.get_loc('TaxonIDProteinID_s2')

        # print(index_description_me_2, index_date_me_2)

        length_2 = np.array([])

        for row in range(0, len(df_me)):
            try:
                number_2 = re.findall(pattern_me2, df_me.iat[row, index_description_me_2])  # .group()
                length_2 = np.append(length_2, len(number_2))
                # print(len(df['number10090']))
                up_2 = []
                if number_2 != []:
                    for nu in number_2:
                        # print(nu)
                        m = re.search(id_p, nu)
                        # print(m.group())
                        up_2.append(m.group())
                # print(up_2)

            except AttributeError:
                number2 = re.findall(pattern_me2, df_me.iat[row, index_description_me_2])  # .group()
                # df['number10090'] = None

            df_me.iat[row, index_date_me_2] = number_2
            df_me.iat[row, index_proteinid2] = up_2
        # print(length_2)
        # print(sum(length_2))
        print(f"genearray of {pattern_sp2}:  {length_2}")
        print(f"sum: {sum(length_2)}")

        def intersection(my_numbers, my_numbers_1):
            inter = [value for value in my_numbers if value in my_numbers_1]
            return inter

        inter_both_species = intersection(my_indices_species_1, my_indices_species_2)

        inter_spezies = np.array([])
        for i in inter_both_species:
            # print(i)
            inter_spezies = np.append(inter_spezies, length_1[i])
        # print(inter_spezies)
        # print(sum(inter_spezies))

        print(
            f"Of species1 {pattern_sp1} {sum(inter_spezies)} genes have at least one homolog in species2 {pattern_sp2}.")

        with open('defaultdata.txt', 'a+') as f:
            print(
                f"Of species1 {pattern_sp1} {sum(inter_spezies)} genes have at least one homolog in species2 {pattern_sp2}.",
                file=f)
        #####################################################################

        list_of_indexes = inter_both_species
        rows_me = df_me.iloc[list_of_indexes, :]

        # rows_me['indices'] = rows_me.index
        new_rows_me = rows_me[['TaxonIDProteinID', 'ProteinCount', 'ProteinID1', 'ProteinID2', 'GroupName']]



        with open('parameter.tsv', 'w') as f:
            with pd.option_context('display.max_rows', None, 'display.max_columns', None):
                print(f"{new_rows_me} ", file=f)

                if args.family_specific_genes:

                    df_family = new_rows_me['TaxonIDProteinID']
                    valst = df_family.to_csv(header=None, index=False).strip('\n').split('\n')

                    new_indices = []
                    for i, val in zip(inter_both_species, valst):

                        fx = val.split(",")
                        if len(fx) < 3:
                            # print(fx)
                            new_indices.append(i)

                    rows_fam = df_me.iloc[new_indices, :]
                    rows_fam['indices'] = rows_fam.index
                    new_rows_me = rows_fam[['indices', 'ProteinCount', 'ProteinID1', 'ProteinID2', 'GroupName']]
                    # new_rows_me
                    df_anno_fam = df_anno.iloc[new_indices, :]
                    df_anno_fam = df_anno_fam[['Description']]

                    df_f = pd.concat([new_rows_me.reset_index(drop=True), df_anno_fam.reset_index(drop=True)], axis=1)
                    print("\n")
                    print(f"{pattern_sp1} share with {pattern_sp1} {len(df_anno_fam)} homolog groups.")
                    print("\n")

                    print(df_f)
                    with open('family.tsv', 'w') as f:
                        with pd.option_context('display.max_rows', None, 'display.max_columns', None):
                            print(f"{df_f} ", file=f)

    if args.three_species:
        ###########################################################################
        ###########################################################################
        ##            File:       eggnog4.species_list.tsv                        #
        ###########################################################################
        ###########################################################################

        print("Choosen species")
        print("\n")
        ##########################################################
        #                 Mus Musculus                           #
        ##########################################################

        pattern_sp1 = r'[Mm]us [Mm]usculus'
        df_row_sp1 = df_sp[df_sp['species'].str.contains(pattern_sp1)]
        print(df_row_sp1)

        with open('default_three_data_t.txt', 'w') as f:
            print(
                f" {df_row_sp1} \n",
                file=f)
        pattern_row_sp1 = df_row_sp1['taxid']
        vals_row_sp1 = pattern_row_sp1.to_csv(header=None, index=False)  # .strip('\n').split('\n').replace('\r', '')
        # als_row_sp.vals_row_sp.replace('\n, '')
        vals_row_sp1 = vals_row_sp1.strip()
        # print(vals_row_sp1)

        ##########################################################
        #                 Homo sapiens                           #
        ##########################################################

        pattern_sp2 = r'Homo sapiens'
        df_row_sp2 = df_sp[df_sp['species'].str.contains(pattern_sp2)]
        print(df_row_sp2)
        with open('default_three_data_t.txt', 'a+') as f:
            print(
                f" {df_row_sp2} \n",
                file=f)
        pattern_row_sp2 = df_row_sp2['taxid']
        pattern_row_sp2
        vals_row_sp2 = pattern_row_sp2.to_csv(header=None, index=False)  # .strip('\n').split('\n').replace('\r', '')
        # als_row_sp.vals_row_sp.replace('\n, '')
        vals_row_sp2 = vals_row_sp2.strip()
        # print(vals_row_sp2)

        ##########################################################
        #                 Pan troglodytes                        #
        ##########################################################

        pattern_sp3 = r'Pan troglodytes'
        df_row_sp3 = df_sp[df_sp['species'].str.contains(pattern_sp3)]
        print(df_row_sp3)

        with open('default_three_data_t.txt', 'a+') as f:
            print(
                f" {df_row_sp3} \n",
                file=f)
        pattern_row_sp3 = df_row_sp3['taxid']
        # pattern_row_sp3
        vals_row_sp3 = pattern_row_sp3.to_csv(header=None, index=False)
        # als_row_sp.vals_row_sp.replace('\n, '')
        vals_row_sp3 = vals_row_sp3.strip()
        # print(vals_row_sp3)

        ###########################################################################
        ###########################################################################
        ##            File:       meNOG.members.tsv                               #
        ###########################################################################
        ###########################################################################

        df_me['TaxonIDProteinID_s1'] = None
        df_me['ProteinID1'] = None

        index_description_me_1 = df_me.columns.get_loc('TaxonIDProteinID')
        index_date_me_1 = df_me.columns.get_loc('TaxonIDProteinID_s1')
        index_proteinid1 = df_me.columns.get_loc('ProteinID1')
        # vals_row_sp2 = '(' + vals_row_sp2
        # vals_row_sp2 += '.[A-Za-z0-9]{7,}'
        pattern_me1 = vals_row_sp2 + '.[.A-Za-z0-9_-]{3,}'
        # pattern_me1 = vals_row_sp2
        # print(pattern_me1)

        # id_p = r'([.A-Za-z0-9_-]{3,})'  # pattern for protein ID
        id_p = r'\b(?<=.)([A-Za-z0-9-_]{5,})\b'
        length_1 = np.array([])

        length_1 = np.array([])

        for row in range(0, len(df_me)):
            try:
                number = re.findall(pattern_me1, df_me.iat[row, index_description_me_1])
                length_1 = np.append(length_1, len(number))
                # print(number)
                up = []
                if number != []:
                    for nu in number:
                        # print(nu)
                        m = re.search(id_p, nu)
                        # print(m.group())
                        up.append(m.group())
                # print(up)

                # print(up)
            except AttributeError:
                number = re.findall(pattern_me1, df_me.iat[row, index_description_me_1])
                # up = re.findall(pattern, df_me.iat[row, index_description_me_1])
                # df.iat[row, index_proteinid1] = up

                # print(up)

            df_me.iat[row, index_date_me_1] = number
            df_me.iat[row, index_proteinid1] = up

        print(f"genearray1 sum: {sum(length_1)}")

        with open('default_three_data_t.txt', 'a+') as f:
            print(
                f"genearray1 sum: {sum(length_1)} \n",
                file=f)
        print("\n")
        # print(df_me)

        ##########################################################
        #                 Mus Musculus                           #
        ##########################################################

        # vals_row_sp1 += '.[A-Za-z0-9]{7,}'
        # vals_row_sp1 + '.[.A-Za-z0-9_-]{3,}'
        # pattern_me2 = vals_row_sp1
        pattern_me2 = vals_row_sp1 + '.[.A-Za-z0-9_-]{3,}'
        df_me['TaxonIDProteinID_s2'] = None
        df_me['ProteinID2'] = None
        index_description_me_2 = df_me.columns.get_loc('TaxonIDProteinID')
        index_date_me_2 = df_me.columns.get_loc('TaxonIDProteinID_s2')
        # print(index_description_me_2, index_date_me_2)
        length_2 = np.array([])

        for row in range(0, len(df_me)):
            try:
                number_2 = re.findall(pattern_me2, df_me.iat[row, index_description_me_2])  # .group()
                length_2 = np.append(length_2, len(number_2))
                # print(len(df['number10090']))
            except AttributeError:
                number2 = re.findall(pattern_me2, df_me.iat[row, index_description_me_2])  # .group()
                # df['number10090'] = None

            df_me.iat[row, index_date_me_2] = number_2
        # print(length_2)

        print(f"genearray2 sum: {sum(length_2)}")
        with open('default_three_data_t.txt', 'a+') as f:
            print(
                f"genearray2 sum: {sum(length_2)} \n",
                file=f)
        # print(df_me)

        #####################################################
        # Indices who contain pattern
        ######################################################
        my_indices_homo_sapiens = df_me.index[df_me['TaxonIDProteinID'].str.contains(pattern_me1)].tolist()
        # print(f"indices of homo sapiens are {my_indices_homo_sapiens}")

        my_indices_mus_musculus = df_me.index[df_me['TaxonIDProteinID'].str.contains(pattern_me2)].tolist()
        # print(f"indices of mus musculus are {my_indices_mus_musculus}")

        ############################################################
        #             complement of 2 species - Homo sapiens
        ########################################################
        my_indices_homo_sapiens
        my_indices_mus_musculus

        def listComplementElements(list1, list2):
            storeResults = []

            for num in list1:
                if num not in list2:  # this will essentially iterate your list behind the scenes
                    storeResults.append(num)

            return storeResults

        compl_spez = listComplementElements(my_indices_homo_sapiens, my_indices_mus_musculus)
        # print(compl_spez) #expected result [1,3,4,7,8]

        compl_homo_sapiens = np.array([])
        for i in compl_spez:
            compl_homo_sapiens = np.append(compl_homo_sapiens, length_1[i])
        # print(compl_homo_sapiens)

        # vals_row_sp3 += '.[A-Za-z0-9]{7,}'
        # pattern_me3 = vals_row_sp3
        pattern_me3 = vals_row_sp3 + '.[.A-Za-z0-9_-]{3,}'
        # pattern_me3
        # pattern_me3 = vals_row_sp3
        my_indices_pan_troglodytes = df_me.index[df_me['TaxonIDProteinID'].str.contains(pattern_me3)].tolist()
        # print(f"indices of pan troglodytes are {my_indices_pan_troglodytes}")

        df_me['TaxonIDProteinID_s3'] = None
        df_me['ProteinID3'] = None

        index_description_me_3 = df_me.columns.get_loc('TaxonIDProteinID')
        index_date_me_3 = df_me.columns.get_loc('TaxonIDProteinID_s3')
        # print(index_description_me_3, index_date_me_3)

        length_3 = np.array([])

        for row in range(0, len(df_me)):
            try:
                number_3 = re.findall(pattern_me3, df_me.iat[row, index_description_me_3])  # .group()
                length_3 = np.append(length_3, len(number_3))
                # print(len(df['number10090']))
            except AttributeError:
                number3 = re.findall(pattern_me2, df_me.iat[row, index_description_me_3])  # .group()
                # df['number10090'] = None

            df_me.iat[row, index_date_me_3] = number_3
        # print(length_3)

        print(f"genearray3 sum: {sum(length_3)}")
        with open('default_three_data_t.txt', 'a+') as f:
            print(
                f"genearray3 sum: {sum(length_3)} \n",
                file=f)

        def intersection(my_numbers, my_numbers_1):
            inter = [value for value in my_numbers if value in my_numbers_1]
            return inter

        inter_homo_pan_species = intersection(compl_spez, my_indices_pan_troglodytes)

        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        #          Final                                  #
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        inter_homo_sapiens = np.array([])
        for i in inter_homo_pan_species:
            # print(i)
            inter_homo_sapiens = np.append(inter_homo_sapiens, length_1[i])
        # print(inter_homo_sapiens)

        # print(sum(inter_homo_sapiens))

        print(
            f"{sum(compl_homo_sapiens)} Homo sapiens genes have not a homolog in Mus musculus, {sum(inter_homo_sapiens)} of these have at least one homolog in Pan troglodytes")

        with open('default_three_data_t.txt', 'a+') as f:
            print(
                f"{sum(compl_homo_sapiens)} Homo sapiens genes have not a homolog in Mus musculus, {sum(inter_homo_sapiens)} of these have at least one homolog in Pan troglodytes",
                file=f)

        rows_interspec = df_me.iloc[inter_homo_pan_species, :]
        rows_interspec
        rows_interspec['indices'] = rows_interspec.index

        rows_pid = rows_interspec[['ProteinID1']]

        print(rows_pid)
        # with np.printoptions(threshold=np.inf):
        with open('three_species_proteinid.tsv', 'w') as f:
            with pd.option_context('display.max_rows', 999999, 'display.max_columns', 999999):
                print(f"{rows_pid} ", file=f)
        # pandas.set_option("display.max_rows", max_rows, "display.max_columns", max_cols)
        rows_interspec_anno = df_anno.iloc[inter_homo_pan_species, :]
        rows_interspec_anno = rows_interspec_anno[['GroupName', 'ProteinCount', 'Description', 'Functional Category']]
        rows_interspec_anno
        with open('three_species_categories.tsv', 'w') as f:
            with pd.option_context('display.max_rows', 999999, 'display.max_columns', 999999):
                print(f"{rows_interspec_anno} ", file=f)


if __name__ == '__main__':
    Main()