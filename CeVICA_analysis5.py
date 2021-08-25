# fifth stage CeVICA_analysis, input cluster files, output excel spreadsheet containing summaries of characteristics of each cluster

# required packages
import glob
import multiprocessing
import xlrd
import xlwt
from statistics import stdev
from collections import Counter


def aa_consensus(seq_list): # input a list of aa sequences, output up to two aa at each position along with percentage (at least 34%, 100 represented by 99), different positions separated by ., C11D22.M11N22...
    # check length distribution, if length varies too much, a single 'B00' is recorded.
    output = ['', 0]
    output_seq = ''
    num_seq = len(seq_list)
    output[1] = num_seq
    list_len = []
    for i in seq_list:
        list_len += [len(i)]
    std_dev = stdev(list_len)
    list_positions = []
    if std_dev < 1:
        best_len = Counter(list_len).most_common(1)[0][0]
        for n in range(best_len):
            str_positions_temp = ''
            for i in seq_list:
                try:
                    str_positions_temp += i[n]
                except IndexError:
                    str_positions_temp += 'B'  # fill in shorter sequences with B
            list_positions += [str_positions_temp]
        for p in list_positions:
            output_temp = ''
            c = Counter(p).most_common(2)
            # print(c)
            if len(c) == 1:
                output_temp = c[0][0] + str(99) +c[0][0] + str(99) + '.' 
            else:
                if c[0][1]/num_seq > 0.34:
                    output_temp += c[0][0] + str(int(c[0][1]/num_seq*100))
                else: 
                    output_temp += 'B00'
                if c[1][1]/num_seq > 0.34:
                    output_temp += c[1][0] + str(int(c[1][1]/num_seq*100)) + '.'
                else: 
                    output_temp += 'B00' + '.'
            output_seq += output_temp
    else:
        output_seq = 'B00'
    output[0] = output_seq
    return output
            

def pick_CDR_most_freq(list_CDR): # from list of CDRs (up to 100?) and pick the most frequent CDR as the representative real sequence of the CDR, when some CDRs has the same frequency, a random one of these CDRs are returned
    input = list_CDR
    CDR_rep = max(set(input), key = input.count)
    return CDR_rep


def CDR_score(input): # use the output of aa_consensus() as input, processes the consensus sequence
    aa_cns = input[0]
    score = 0
    list_positions = aa_cns.rstrip().split('.')
    for n in range(len(list_positions)-1):
        if list_positions[n][0] != 'B':
            if int(list_positions[n][1:3]) > 80:
                score += 3
            elif int(list_positions[n][1:3]) > 50:
                score += 2
            else:
                score += 1
        else:
            score += 1
    return score


def cluster_stats(name): # name is the generated file path
    save_to_excel = [name.split('\\')[1][:2], int(name.split('\\')[1].split('cluster')[1].split('.')[0])]
    with open(name, 'r') as f:
        CDR1s = []
        CDR2s = []
        CDR3s = []
        while True:
            seq = f.readline()
            data_split = seq.rstrip().split('#')
            try:
                CDR3s += [data_split[2]]
                CDR2s += [data_split[1]]
                CDR1s += [data_split[0]]
            except IndexError:
                break
        CDR1_and_cluster_size = aa_consensus(CDR1s)
        cluster_size = CDR1_and_cluster_size[1]
        CDR1_cns = CDR1_and_cluster_size[0]
        CDR2_cns = aa_consensus(CDR2s)[0]
        CDR3_cns = aa_consensus(CDR3s)[0]
        CDR1_score = CDR_score(aa_consensus(CDR1s))
        CDR2_score = CDR_score(aa_consensus(CDR2s))
        CDR3_score = CDR_score(aa_consensus(CDR3s))
        total_score = CDR1_score + CDR2_score + CDR3_score
        CDR1_topseq = CDR1_cns[::7]
        CDR2_topseq = CDR2_cns[::7]
        CDR3_topseq = CDR3_cns[::7]
        CDR1_rep = pick_CDR_most_freq(CDR1s)
        CDR2_rep = pick_CDR_most_freq(CDR2s)
        CDR3_rep = pick_CDR_most_freq(CDR3s)
        CDR1_unique = 'Yes'
        CDR2_unique = 'Yes'
        CDR3_unique = 'Yes'
        if CDR1_topseq in CDRs_bg[0]:
            CDR1_unique = 'No'
        if CDR2_topseq in CDRs_bg[1]:
            CDR2_unique = 'No'
        if CDR3_topseq in CDRs_bg[2]:
            CDR3_unique = 'No' 
        save_to_excel += [cluster_size, CDR1_rep, CDR2_rep, CDR3_rep, CDR1_cns, CDR2_cns, CDR3_cns, CDR1_score, CDR2_score, CDR3_score, total_score, CDR1_unique, CDR2_unique, CDR3_unique]
    return save_to_excel


# optional: read in background CDR sequences
CDRs_bg = [[],[],[]]
try:
    workbook_bg_path = 'background_cluster.xls'
    workbook_bg = xlrd.open_workbook(workbook_bg_path)
    worksheet_bg = workbook_bg.sheet_by_index(0)
    for n in range(3):
        nrows = 1
        while True:
            try: 
                content = worksheet_bg.cell_value(nrows,n+3)
                if 'B' not in content and content not in CDRs_bg[n]:
                    CDRs_bg[n] += [content]
                nrows += 1
            except IndexError:
                break
except:
    pass


if __name__ == '__main__':
    # Input definition, change accordingly #
    list_sample = ['SR3', 'GP3', 'STR'] # sample identifier for generating filepath for analysis
    for s in list_sample:
        inputpath = f'{s}_clusters'
        path = inputpath + '\\*cluster*.txt'
        files = glob.glob(path)
        cluster_score_output_head = ['cluster_type', 'cluster_ID', 'cluster_size', 'CDR1_rep', 'CDR2_rep', 'CDR3_rep', 'CDR1_consensus', 'CDR2_consensus', 'CDR3_consensus', 'CDR1_score', 'CDR2_score', 'CDR3_score', 'total_score', 'CDR1_unique', 'CDR2_unique', 'CDR3_unique']
        workbook_output_cluster_score = xlwt.Workbook(encoding='utf-8')
        worksheet_output_cluster_score = workbook_output_cluster_score.add_sheet('cluster_score', cell_overwrite_ok=True)
        col_num = 0
        row_num = 0
        for i in cluster_score_output_head:
            worksheet_output_cluster_score.write(row_num, col_num, i)
            col_num += 1

        pool_cluster_stats = multiprocessing.Pool()
        output = pool_cluster_stats.map(cluster_stats, files)
        for i in output:
            col_num = 0
            row_num += 1
            for s in i:
                worksheet_output_cluster_score.write(row_num, col_num, s)
                col_num += 1
        save_path = path.split('\\')[0] + f'\\{inputpath}_cluster_score.xls'
        workbook_output_cluster_score.save(save_path)
