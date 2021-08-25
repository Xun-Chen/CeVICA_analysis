# fourth stage CeVICA_analysis, input VHH protein segments, output clusters based on CDR similarity

# required packages
import sqlite3
import xlrd
import multiprocessing
import time
from functools import partial
import pathlib
import glob


def match_test(item, match_to): # input two sets of CDRs and calculate the match score, return the two sets plus Y/N match
    global match_cutoff
    output = [match_to, item]
    score_cdr1 = -5
    score_cdr2 = -5
    score_cdr3 = -5

    if len(match_to[0]) != len(item[0]):
        score_cdr1 = -5
    else:
        for p in range(len(match_to[0])):
            score_cdr1 += matrix_aas[match_to[0][p] + item[0][p]]
    
    if len(match_to[1]) != len(item[1]):
        score_cdr2 = -5
    else:
        for p in range(len(match_to[1])):
            score_cdr2 += matrix_aas[match_to[1][p] + item[1][p]]
    
    if len(match_to[2]) != len(item[2]):
        score_cdr3 = -5
    else:
        for p in range(len(match_to[2])):
            score_cdr3 += matrix_aas[match_to[2][p] + item[2][p]]

    if score_cdr1 + score_cdr2 > match_cutoff or score_cdr1 + score_cdr3 > match_cutoff or score_cdr2 + score_cdr3 > match_cutoff:
        output += ['Y']
    else:
        output += ['N']
    return output


def pick_clusterer(item, list):
    output = [item]
    match_count = 0
    status = ''
    for i in list:
        if match_test(item, i)[2] == 'Y':
            match_count += 1
        if match_count >= cluster_size_min:
            output += ['Y']
            status = 'complete'
            break
    if status != 'complete':
        output += ['N']
    return output


def cluster_size():
    path = f'{name_change}cluster*.txt'
    files = glob.glob(path)
    list_output = []
    for file_path in files:
        count = 0
        for line in open(file_path).readlines(): 
            count += 1
        list_output += [count]
    
    # include the following block if including loners in graph
    path_loner = f'{name_change}Loners.txt' 
    count_loner = 0
    for line in open(path_loner).readlines(): 
        count_loner += 1
    for n in range(count_loner):
        list_output += [1]

    return list_output


matrix_aas = {}
workbook_aas = xlrd.open_workbook("Amino_Acids_Similarity_Matrix.xlsx")
worksheet_aas = workbook_aas.sheet_by_index(0)
for ncols in range(1,23):
    for nrows in range(1,23):
        ax = worksheet_aas.cell_value(0,ncols)
        ay = worksheet_aas.cell_value(nrows,0)
        score_aas = worksheet_aas.cell_value(nrows,ncols)
        codepair_aas = str(ax) + str(ay)
        matrix_aas[codepair_aas] = float(score_aas)


# settings
name_change = 'RC' # scoring strategy note
match_cutoff = 35 # match score cut off value
cluster_size_min = 5


if __name__ == '__main__':
    # Input definition, change accordingly #
    list_sample = ['SR3', 'GP3', 'STR'] # sample identifier for generating filepath for analysis
    for s in list_sample:
        start = time.time()
        filepath = f'{s}_aa_profile.db'
        pathlib.Path(f'{s}_clusters').mkdir(parents=True, exist_ok=True)
        list_remove = ['','RTFSSYA']
        print(filepath)
        print('match_cutoff: ', match_cutoff)
        print('cluster size minimum: ', cluster_size_min)
        print('CDR1 removal:', list_remove)
        db_input = sqlite3.connect(filepath)
        cursor_input = db_input.cursor()
        cursor_input.execute(f'SELECT * FROM parts')

        # extract CDRs from input
        count_input = 0
        list_init = []
        list_loners = []
        list_init_size = 1500000
        while count_input < list_init_size:
            data = cursor_input.fetchone()
            pick_cdrs = ['','','','']
            if data == None:
                break
            if '*' not in data[1] and data[2] != 'bad sequence' and data[3] not in list_remove and data[7] != 'GRGGGS':
                count_input += 1
                pick_cdrs[0] = data[3]
                pick_cdrs[1] = data[5]
                pick_cdrs[2] = data[7]
                pick_cdrs[3] = '#########' + data[9] + '##########' + data[10] + '\n'
                list_init += [pick_cdrs]
        
        # pick out items that form clusters, including this stage reduced time dramatically (10min vs 1min for 10000 seqs)
        pool = multiprocessing.Pool()
        result_select = pool.map(partial(pick_clusterer, list = list_init), list_init)
        list_init = []
        for i in result_select:
            if i[1] == 'Y':
                list_init += [i[0]]
            else:
                list_loners += [i[0]]
        with open(f'{s}_clusters\\{name_change}Loners.txt', 'w') as l:
            for i in list_loners:
                l.write(i[0] + '#' + i[1] + '#' + i[2] + '#' + i[3])

        mid = time.time()
        print(f'{len(list_init)} Clusterer found. Time spent: {int((mid - start)//60)}mins {round((mid - start)%60,2)}s')

        # construct clusters
        num_analyzed = 0
        cluster_id = 0
        while len(list_init) > 1:  
            standard = list_init[0]
            list_init.remove(list_init[0])
            temp_save = [standard[0] + '#' + standard[1] + '#' + standard[2] + '#' + standard[3]]
            result = pool.map(partial(match_test, match_to = standard), list_init) 

            list_init = []
            for i in result:
                if i[2] == 'Y':
                    temp_save += [i[1][0] + '#' + i[1][1] + '#' + i[1][2] + '#' + i[1][3]]
                else:
                    list_init += [i[1]]

            if len(temp_save) >= cluster_size_min:    
                with open(f'{s}_clusters\\{name_change}cluster{cluster_id}.txt', 'w') as save:
                    for i in temp_save:    
                        save.write(i)
                cluster_id += 1
            num_analyzed += 1
            if num_analyzed % 500 == 0:
                mid = time.time()
                print(f'{len(list_init)} sequences remaining... time spent: {int((mid - start)//60)}mins {round((mid - start)%60,2)}s')
        end = time.time()
        print(f'total processing time: {int((end - start)//60)}mins {round((end - start)%60,2)}s')
