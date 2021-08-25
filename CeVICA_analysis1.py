# first stage CeVICA_analysis: input fastq for paired end sequencing reads, and output merged reads separated by index

# required packages
import xlrd 
import sqlite3
import multiprocessing
import time


def read_paired_end(r1,r2):
    lines_4x2 = []
    lines_4x2 += [r1.readline().strip()] + [r1.readline().strip()] + [r1.readline().strip()] + [r1.readline().strip()]
    lines_4x2 += [r2.readline().strip()] + [r2.readline().strip()] + [r2.readline().strip()] + [r2.readline().strip()]
    if len(lines_4x2[0]) > 2:
        return lines_4x2
    else:
        return 'end of file'


def merge_2X250(lines_4X2): # merges long reads and combine quality scores for overlap region, first 4 lines read1, last 4 lines read2.
    str_r2c = '' # reverse complement read2 sequence
    for n in range(0,len(lines_4X2[5])):
        str_r2c = str_r2c + comp_encoding[lines_4X2[5][n]]
    str_r2rc = str_r2c[::-1]
    lines_4X2[5] = str_r2rc
    lines_4X2[7] = lines_4X2[7][::-1] # reverse read2 quality score

    # sequence alignment for merging
    status = ''
    score = 0
    slice_location = 0
    for n in range(50,100): # adjust this search window depending on how much the two reads overlap, changed from 50,150 8-9-2019
        for seq_location in range(0,16):
            if lines_4X2[1][len(lines_4X2[1])-16+seq_location] == lines_4X2[5][n + seq_location]:
                score +=1
        if score > match_score_cutoff:
            status = 'complete'
            slice_location = n
            break
        score = 0
    if status == 'complete':
        r1_only = lines_4X2[1][:len(lines_4X2[1])-slice_location-16]
        r2_only = lines_4X2[5][slice_location+16:]
        r1_overlap = lines_4X2[1][len(lines_4X2[1])-slice_location-16:]
        r2_overlap = lines_4X2[5][:slice_location+16]
        r1_onlyQS = lines_4X2[3][:len(lines_4X2[1])-slice_location-16]
        r2_onlyQS = lines_4X2[7][slice_location+16:]
        r1_overlapQS = lines_4X2[3][len(lines_4X2[1])-slice_location-16:]
        r2_overlapQS = lines_4X2[7][:slice_location+16]
        overlap_comb = ''
        overlap_combQS = ''
        for n in range(0,slice_location+16):
            if r1_overlap[n] == r2_overlap[n]:
                overlap_comb += r1_overlap[n]
                try:
                    overlap_combQS += list(qs_encoding.keys())[list(qs_encoding.values()).index(qs_encoding[r1_overlapQS[n]] + qs_encoding[r2_overlapQS[n]])]
                except:
                    overlap_combQS += 'I'
            if r1_overlap[n] != r2_overlap[n]:
                if qs_encoding[r1_overlapQS[n]] >= qs_encoding[r2_overlapQS[n]]:
                    overlap_comb += r1_overlap[n]
                    overlap_combQS += list(qs_encoding.keys())[list(qs_encoding.values()).index(qs_encoding[r1_overlapQS[n]])]
                if qs_encoding[r1_overlapQS[n]] < qs_encoding[r2_overlapQS[n]]:
                    overlap_comb += r2_overlap[n]
                    overlap_combQS += list(qs_encoding.keys())[list(qs_encoding.values()).index(qs_encoding[r2_overlapQS[n]])]
        output = [lines_4X2[0]]
        output += [r1_only + overlap_comb + r2_only]
        output += [r1_onlyQS + overlap_combQS + r2_onlyQS]
        output += ['None']
        output += ['None']
        # identify index
        index_neighboring_seq = 'TCAGAAG' # for RDS5 index
        if index_neighboring_seq in output[1][-25:]:
            output[3] = output[1][-25:][output[1][-25:].index(index_neighboring_seq)+7:]
        
        # identify UMI
        if 'CTTTAAG' in output[1][:25]:
            output[4] = output[1][:output[1][:25].index('CTTTAAG')]
        return output
                      

def DNA_match_score(seq1,seq2):
    score = 0
    if len(seq1) == len(seq2):
        for n in range(len(seq1)):
            if seq1[n] == seq2[n]:
                score += 1
    else:
        score = 0
    return score


def split_by_index8(db_path):
    # RDS5 index set
    sampA1 = ['C']
    sampA2 = ['CCAAG']
    sampB1 = ['TT']
    sampB2 = ['CGTTGGT']
    sampC1 = ['CGA']
    sampC2 = ['CTAGATCC']
    sampD1 = ['TACG']
    sampD2 = ['GAGCAGCTA']

    db = sqlite3.connect(db_path)
    cursor = db.cursor()
    with open('sample_A1.fastq', 'w') as a:
        with open('sample_A2.fastq', 'w') as b:
            with open('sample_B1.fastq', 'w') as c:
                with open('sample_B2.fastq', 'w') as d:
                    with open('sample_C1.fastq', 'w') as e:
                        with open('sample_C2.fastq', 'w') as f:
                            with open('sample_D1.fastq', 'w') as g:
                                with open('sample_D2.fastq', 'w') as h:
                                    cursor.execute('SELECT * FROM fastq_processed')
                                    for row in cursor:
                                        if row[5] in sampA1:
                                            a.write(row[3] + '#' + row[4] + '\n')
                                            a.write(row[1] + '\n')
                                            a.write('+' + '\n')
                                            a.write(row[2] + '\n')
                                        if row[5] in sampA2:
                                            b.write(row[3] + '#' + row[4] + '\n')
                                            b.write(row[1] + '\n')
                                            b.write('+' + '\n')
                                            b.write(row[2] + '\n')
                                        if row[5] in sampB1:
                                            c.write(row[3] + '#' + row[4] + '\n')
                                            c.write(row[1] + '\n')
                                            c.write('+' + '\n')
                                            c.write(row[2] + '\n')
                                        if row[5] in sampB2:
                                            d.write(row[3] + '#' + row[4] + '\n')
                                            d.write(row[1] + '\n')
                                            d.write('+' + '\n')
                                            d.write(row[2] + '\n')
                                        if row[5] in sampC1:
                                            e.write(row[3] + '#' + row[4] + '\n')
                                            e.write(row[1] + '\n')
                                            e.write('+' + '\n')
                                            e.write(row[2] + '\n')
                                        if row[5] in sampC2:
                                            f.write(row[3] + '#' + row[4] + '\n')
                                            f.write(row[1] + '\n')
                                            f.write('+' + '\n')
                                            f.write(row[2] + '\n')
                                        if row[5] in sampD1:
                                            g.write(row[3] + '#' + row[4] + '\n')
                                            g.write(row[1] + '\n')
                                            g.write('+' + '\n')
                                            g.write(row[2] + '\n')
                                        if row[5] in sampD2:
                                            h.write(row[3] + '#' + row[4] + '\n')
                                            h.write(row[1] + '\n')
                                            h.write('+' + '\n')
                                            h.write(row[2] + '\n')


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
match_score_cutoff = 15 
len_cutoff = 40 


if __name__ == '__main__':
    start = time.time()
# Input definition, change accordingly on the line below #
    list_sample = ['SR3', 'GP3', 'STR'] # sample identifier for generating filepath for analysis
    for s in list_sample:
        status = ''
        print('merging reads...')
        total_count = 0
        merged_count = 0
        # create a database file to store merged reads before demultiplexing, not required if demultiplexing is already done
        db_fastq = sqlite3.connect(f'{s}_fastq.db')
        cursor_fastq = db_fastq.cursor()
        try:
            cursor_fastq.execute(
                'CREATE TABLE fastq_processed(id INTEGER PRIMARY KEY, sequences TEXT, Q_score TEXT, read_id TEXT, UMI TEXT, barcode TEXT)')
        except:
            pass
        with open(f'{s}_read1.fastq', 'r') as r1: # read1 fastq file 
            with open(f'{s}_read2.fastq', 'r') as r2: # read2 fastq file
                while True:
                    if status == 'complete':
                        break
                    list_of_blocks = []
                    for n in range(0,100000):
                        lines_4x2 = read_paired_end(r1,r2)
                        if len(lines_4x2[1])>225 and len(lines_4x2[5])>225:
                            list_of_blocks += [lines_4x2]
                            total_count += 1
                        elif lines_4x2 == 'end of file':    
                            print('end of file')
                            status = 'complete'
                            break
                    pool = multiprocessing.Pool()
                    result = pool.map(merge_2X250, list_of_blocks)

                    # save merged reads to database
                    for block in result:
                        try:
                            cursor_fastq.execute('INSERT INTO fastq_processed(sequences, Q_score, read_id, UMI, barcode) VALUES(?,?,?,?,?)', (block[1], block[2], block[0], block[4], block[3]))
                            merged_count += 1
                        except:
                            pass
                    db_fastq.commit()

                    # save merged reads to fastq
                    with open(f'{s}.fastq', 'a') as save_fastq:
                        for block in result:
                            if block != None:
                                save_fastq.write(block[0] + '\n')
                                save_fastq.write(block[1] + '\n')
                                save_fastq.write('+' + '\n')
                                save_fastq.write(block[2] + '\n')

        print('merging reads completed')
        print(f'{total_count} sequences total for {s}')
        print(f'{merged_count} sequences merged for {s}')
        end = time.time()
        print(f'total processing time: {int((end - start)//60)}mins {round((end - start)%60,2)}s')

    # optional demultiplexing
    # print('sorting reads by index...')
    # split_by_index8('_fastq.db')
