# second stage CeVICA_analysis, input index separated merged reads fastq file, output trimmed reads based on Q_score, then split vhhlib DNA into segments

# required packages
import xlrd
import sqlite3
import multiprocessing


def read_fastq(r1): # r1 is the file handle of the fastq file
    lines_4 = []
    lines_4 += [r1.readline().strip()] + [r1.readline().strip()] + [r1.readline().strip()] + [r1.readline().strip()]
    if len(lines_4[0]) > 2:
        return lines_4
    else:
        return 'end of file'


def trim_seq(lines_4):
    list_trim = []
    str_trim = ''
    for n in range(0,len(lines_4[1])-1):
        if qs_encoding[lines_4[3][n]] > qs_cutoff:
            str_trim = str_trim + lines_4[1][n]
        else:
            list_trim = list_trim + [str_trim]
            str_trim = ''
    list_trim = list_trim + [str_trim]
    if len(max(list_trim, key=len)) > len_cutoff:
        lines_4[1] = max(list_trim, key=len)
        # check the timmed sequence still contain full CDR
        seq_scan = 0
        match_score_cutoff = 12
        status_5p = ''
        status_3p = ''
        check_5p = 'ACCATGCAGGTTCAA'
        check_3p = 'GAACAGAAACTGATT'
        for l in range(0,50):
            score = 0
            for n in range(0,15):
                if lines_4[1][n+seq_scan] == check_5p[n]:
                    score += 1
            if score >= match_score_cutoff:
                status_5p = 'complete'
                break
            seq_scan += 1
        seq_scan = 0
        for l in range(0,30):
            score = 0
            for n in range(0,15):
                if lines_4[1][len(lines_4[1])+n-15+seq_scan] == check_3p[n]:
                    score += 1
            if score >= match_score_cutoff:
                status_3p = 'complete'
                break
            seq_scan += -1
        if status_5p == 'complete' and status_3p == 'complete':
            return lines_4
        else:
            return 'CDR incomplete'
    else:
        return 'too short'


def split_vhhDNA(seq): # split vhhDNA into segments, segments shorter than defined length are filled in by blanks
    # whole CDS deletion and insertion check
    match_score_cutoff = 12
    list_CDRref = ['ATGCAGGTTCAACTGCAAGAATCTGGTGGCGGTTTAGTTCAAGCCGGTGG', 'GTCACGGTCTCAAGCGGCGG']
    dic_refseq = {'frame1':'ATGCAGGTTCAACTGCAAGAATCTGGTGGCGGTTTAGTTCAAGCCGGTGGTAGCCTTCGTCTGAGCTGCGCCGCAAGCGGA', 'frame2':'ATGGGTTGGTTTCGCCAGGCACCTGGTAAAGAACGTGAATTTGTGGCAGCCATTTCA', 'frame3':'ACCTACTATGCCGATTCTGTTAAAGGTCGCTTCACTATTTCACGCGATAATGCGAAGAATACCGTGTACTTACAGATGAACTCACTGAAACCTGAGGATACCGCCGTTTATTATTGTGCGGCA', 'frame4':'GACTATTGGGGTCAGGGCACACAAGTCACGGTCTCAAGC'}
    list_marks = [0,0,0,0,0,0]
    seq_scan = 0
    for i in dic_refseq.values():
        score = 0
        for l in range(seq_scan,seq_scan+1):
            no_match = 0
            for n in range(0,len(i)):
                if i[n] == seq[n]:
                    score += 1
                    no_match = 0
                else:
                    no_match += 1
                if no_match > 12:
                    seq_scan += 1
                    break
        if score > len(i)*0.7:
            pass


def find_cdr(lines_4): # finds CDRs and their neigboring sequence
    seq = lines_4[1]
    match_score_cutoff = 15
    list_CDRref = ['TGAGCTGCGCCGCAAGCGGA','ATGGGTTGGTTTCGCCAGGC','AATTTGTGGCAGCCATTTCA','ACCTACTATGCCGATTCTGT','CCGTTTATTATTGTGCGGCA','GACTATTGGGGTCAGGGCAC']
    list_marks = [0,0,0,0,0,0]
    seq_scan = 0
    for n in range(0,6):
        for l in range(seq_scan,len(seq)-50):
            score = 0
            no_match = 0
            for p in range(0,len(list_CDRref[n])):
                if list_CDRref[n][p] == seq[p + seq_scan]:
                    score += 1
                    no_match = 0
                else:
                    no_match += 1
                if no_match > 5:
                    seq_scan += 1
                    break
            if score < match_score_cutoff and no_match <=5:
                seq_scan += 1
            elif score >= match_score_cutoff:
                list_marks[n] = seq_scan
                break
    # retrieve CDR region sequence
    if 0 in list_marks:
        pass
    else:
        CDR1_seq = seq[list_marks[0]:list_marks[1]+20]
        CDR2_seq = seq[list_marks[2]:list_marks[3]+20]
        CDR3_seq = seq[list_marks[4]:list_marks[5]+20]
        output = [CDR1_seq, CDR2_seq, CDR3_seq, lines_4[0]]
        return output


workbook_qs = xlrd.open_workbook("quality_score.xlsx")
worksheet_qs = workbook_qs.sheet_by_index(0)
comp_encoding = {'A':'T', 'T':'A','G':'C','C':'G','N':'N'}
qs_encoding = {}
for nrows in range(1,42):
    try:
        qs_code = int(worksheet_qs.cell_value(nrows,0))
    except:
        qs_code = worksheet_qs.cell_value(nrows,0)
    qs_encoding[str(qs_code)] = worksheet_qs.cell_value(nrows,2)
qs_cutoff = 10 # quality score cutoff value
len_cutoff = 300


if __name__ == '__main__':
# Input definition, change accordingly: no change required if using default output files from CeVICA_analysis1.py #
    list_sample = ['SR3', 'GP3', 'STR'] # sample identifier for generating filepath for analysis
    print(f'qs_cutoff = {qs_cutoff}')
    print(f'reads trimming...')
    for s in list_sample:
        status = ''
        filepath = f'{s}.fastq'
        with open(filepath, 'r') as r1:
            with open(f'{s}_trimmed.fastq', 'w') as w1:
                while True:
                        if status == 'complete':
                            break
                        list_of_blocks = []
                        for n in range(0,100000): 
                            lines_4 = read_fastq(r1)
                            if lines_4 != 'end of file':
                                list_of_blocks += [lines_4]
                            else:    
                                status = 'complete'
                                break
                        pool = multiprocessing.Pool()
                        result = pool.map(trim_seq, list_of_blocks)
                        for block in result:
                            if len(block) == 4:
                                for line in block:
                                    w1.write(line + '\n')

# generate CDRs
    print(f'generating CDRs...')
    for s in list_sample:
        db_CDR = sqlite3.connect(f'{s}_CDRs.db')
        cursor_CDR = db_CDR.cursor()
        try:
            cursor_CDR.execute('CREATE TABLE CDR_sequences(id INTEGER PRIMARY KEY, CDR1 TEXT, CDR1_len INTEGER, CDR2 TEXT, CDR2_len INTEGER, CDR3 TEXT, CDR3_len INTEGER, read_id TEXT)')
        except:
            pass

        with open(f'{s}_trimmed.fastq', 'r') as r1:
            status = ''
            while True:
                if status == 'complete':
                    break
                list_of_blocks = []
                for n in range(0,100000): 
                    block = read_fastq(r1)
                    if len(block) == 4:
                        list_of_blocks += [block]
                    else:    
                        status = 'complete'
                        break
                pool = multiprocessing.Pool()
                result = pool.map(find_cdr, list_of_blocks)

                for block in result:
                    try:
                        cursor_CDR.execute('INSERT INTO CDR_sequences(CDR1, CDR1_len, CDR2, CDR2_len, CDR3, CDR3_len, read_id) VALUES(?,?,?,?,?,?,?)', (block[0], len(block[0]), block[1], len(block[1]), block[2], len(block[2]), block[3]))
                    except:
                        pass
                db_CDR.commit()
