# third stage CeVICA_analysis, input trimmed fastq files, output VHH DNA translated into VHH protein sequence, then split VHH protein into segments, generate aa_profile

# required packages
import sqlite3
import multiprocessing
import xlrd
import xlwt
import datetime


def DNA_orf_split(DNA_seq): # convert one DNA sequence into vhh peptide filled-in to full length and then split into segments
    output = ['','','','','','','','','','']

    workbook_bdr = xlrd.open_workbook("CDR_Boundries.xlsx")
    worksheet_bdr = workbook_bdr.sheet_by_index(0)
    list_frame_len = [24,17,39,11]

    orf_check = 'SLRLSCAASG'
    seq = DNA_seq[19:]
    total_bases = len(seq)
    status_orf = ''
    for orf_shift in range(0,5):
        if status_orf == 'complete':
            break
        orf_step = 0
        data_aa = matrix_cdo.get(seq[orf_step+orf_shift:orf_step+3+orf_shift])
        while orf_step+3 < 135:
            orf_step = orf_step +3
            data_aa = data_aa + matrix_cdo.get(seq[orf_step+orf_shift:orf_step+3+orf_shift])
        if '*' in data_aa:
            continue
        for orf_position in range(0,26):
            match_score = 0
            for orf_check_position in range(0,9):
                match_score = match_score + matrix_aas.get(data_aa[orf_position]+orf_check[orf_check_position])
                orf_position +=1
            if match_score < 15:
                continue
            else:
                orf_step = 0 
                data_aa = matrix_cdo.get(seq[orf_step+orf_shift:orf_step+3+orf_shift])
                while orf_step+orf_shift < total_bases-25: 
                    orf_step = orf_step +3
                    data_aa = data_aa + matrix_cdo.get(seq[orf_step+orf_shift:orf_step+3+orf_shift])
                if orf_position < 26:
                    data_aa ='          ' + data_aa
                output[0] = data_aa
                output[9] = seq  
                status_orf = 'complete'
                break
    if orf_shift==4:
        for n in range(0,8):
            output[n] = 'bad sequence'
    # split the orf into vhh frames and cdrs
    else:
        data_aa_len = len(data_aa)
        list_bdr = [int(0),int(0),int(0),int(0)]
        for bdr_number in range(0,4):
            scan_seq = list_bdr[bdr_number-1]
            frame_std = worksheet_bdr.cell_value(bdr_number,1)
            status_parts = ''
            match_score = 0
            for scan_seq in range(list_bdr[bdr_number-1],data_aa_len-list_frame_len[bdr_number]):
                for frame_location in range(0,list_frame_len[bdr_number]+1):
                    match_score = match_score + matrix_aas.get(data_aa[scan_seq]+frame_std[frame_location])
                    scan_seq +=1
                if match_score > 20:
                    list_bdr[bdr_number] = scan_seq
                    status_parts = 'complete'
                    break
                else:
                    match_score = 0
            if status_parts != 'complete':
                for n in range(1,8):  
                    output[n] = 'bad sequence'
                break
        if status_parts == 'complete':
            output[1] = data_aa[list_bdr[0]-list_frame_len[0]-1:list_bdr[0]+1]
            output[2] = data_aa[list_bdr[0]+1:list_bdr[1]-list_frame_len[1]-1]
            output[3] = data_aa[list_bdr[1]-list_frame_len[1]-1:list_bdr[1]+1]
            output[4] = data_aa[list_bdr[1]+1:list_bdr[2]-list_frame_len[2]-1]
            output[5] = data_aa[list_bdr[2]-list_frame_len[2]-1:list_bdr[2]+1]
            output[6] = data_aa[list_bdr[2]+1:list_bdr[3]-list_frame_len[3]-1]
            output[7] = data_aa[list_bdr[3]-list_frame_len[3]-1:list_bdr[3]+1]
            output[8] = output[1].lower()+output[2]+output[3].lower()+output[4]+output[5].lower()+output[6]+output[7].lower() 
    return output


def feature_count(sample):
    filepath = f'{sample}_vhhlib_orf_parts.db'
    db_input = sqlite3.connect(filepath)
    cursor_input = db_input.cursor()
    cursor_input.execute(f'SELECT Full, CDR1, CDR2, CDR3 FROM parts')
    total = 0
    unknown = 0 # bad_sequence in Full
    frame_shift = 0 # okay sequence in Full but bad_sequence in all others
    stop_early = 0
    stop_free = 0
    stop_free_good = 0 # CDR1 = 7, CDR2 >= 4, CDR3 >= 4
    stop_free_full_len = 0
    for d in cursor_input.fetchall():
        total += 1
        if d[0] == 'bad sequence':
            unknown += 1
        elif d[1] == 'bad sequence':
            frame_shift += 1
        elif '*' in d[0]:
            stop_early += 1
        else:
            stop_free += 1 
            if len(d[1]) == 7 and len(d[2]) == 5:
                if len(d[3]) == 6 or len(d[3]) == 9 or len(d[3]) == 10 or len(d[3]) == 13:
                    stop_free_full_len += 1
            if len(d[1]) == 7 and len(d[2]) >= 4 and len(d[3]) >= 4:
                stop_free_good += 1 
    print(f'{total} sequences total')
    print(f'{unknown} unknown sequences')
    print(f'{frame_shift} sequences have frame_shift')
    print(f'{stop_free} sequences stop free')
    print(f'{stop_free_good} sequences stop free good')
    print(f'{stop_free_full_len} sequences stop free and full_len')
    output = [sample, total, unknown, frame_shift, stop_early, stop_free, stop_free_good, stop_free_full_len, stop_free_good/total, stop_free_full_len/total]
    return output


def aa_profile(sample): # analyze amino acid profile of segments stored in database, results are stored in a list
    filepath = f'{sample}_vhhlib_orf_parts.db'
    db_input = sqlite3.connect(filepath)
    cursor_input = db_input.cursor()
    list_aa_table = ['P','A','G','M','I','L','V','F','W','C','S','T','N','Q','Y','H','R','K','D','E',' ']
    list_parts_names = ['Frame1','CDR1','Frame2','CDR2','Frame3','CDR3','Frame4']
    list_parts_len = [26,10,19,10,41,30,13]
    set_num_seq_analyzed_at_a_time = 500
    num_seq_analyzed = 0

    # # build a database to store sequences selected for aa_profile
    db_select = sqlite3.connect(f'{sample}_aa_profile.db')
    cursor_select = db_select.cursor()
    try:
        cursor_select.execute('CREATE TABLE parts(id INTEGER PRIMARY KEY, Full TEXT, Frame1 TEXT, CDR1 TEXT, Frame2 TEXT, CDR2 TEXT, Frame3 TEXT, CDR3 TEXT, Frame4 TEXT, CDS TEXT, DNA TEXT)')
        db_select.commit()
    except:
        pass
    cursor_input.execute(f'SELECT * FROM parts')
    for d in cursor_input.fetchall(): 
        if '*' not in d[1] and d[2] != 'bad sequence':  
            cursor_select.execute('INSERT INTO parts(Full, Frame1, CDR1, Frame2, CDR2, Frame3, CDR3, Frame4, CDS, DNA) VALUES(?,?,?,?,?,?,?,?,?,?)', (d[1],d[2],d[3],d[4],d[5],d[6],d[7],d[8],d[9],d[10]))
    db_select.commit()


    # build a list of lists to serve as a table to store intermediate data, do not build by multiplying one list
    list_output = []
    for list_columns in range(0,163):
        list_one_column = []
        for list_rows in range(0,21):
            list_one_column = list_one_column + [0]
        list_output = list_output + [list_one_column]
    
    cell_position = 0
    part_number = 0
    for part_number in range(0,7):
        total_num_analyzed = 0
        nrows_parts = 0
        status = ''
        while True:
            list_aa_profile = []
            while True:
                nrows_parts +=1
                cursor_select.execute('SELECT {colv} FROM parts WHERE id= {nrows}'.\
                        format(colv=list_parts_names[part_number], nrows=nrows_parts))
                content_db = cursor_select.fetchone()
                if content_db == None:
                    status = 'complete'
                    num_seq_analyzed = 0
                    break
                list_aa_add = content_db[0]
                if list_aa_add != 'bad sequence':
                    list_aa_profile = list_aa_profile + [list_aa_add]
                num_seq_analyzed +=1
                if num_seq_analyzed == set_num_seq_analyzed_at_a_time:
                    num_seq_analyzed = 0
                    break
            total_num_analyzed = total_num_analyzed + len(list_aa_profile)
            # fill frames from different antibodies to the length of at least what is defined in list_parts_len
            for nrows_aa_profile in range(0,len(list_aa_profile)):
                if len(list_aa_profile[nrows_aa_profile]) < list_parts_len[part_number]:
                    list_aa_profile[nrows_aa_profile] =list_aa_profile[nrows_aa_profile] + ' '*(list_parts_len[part_number]-len(list_aa_profile[nrows_aa_profile]))
            list_aa_positions = []
            str_aa_positions = ''
            for n_sub_str in range(0,list_parts_len[part_number]):
                for nrows_aa_profile in range(0,len(list_aa_profile)):
                    str_aa_positions = str_aa_positions + list_aa_profile[nrows_aa_profile][n_sub_str]
                list_aa_positions = list_aa_positions + [str_aa_positions]
                str_aa_positions = ''

            for aa_positions_number in range(0,list_parts_len[part_number]):
                cell_position = cell_position + 1
                for aa_number in range(0,21):
                    count = 0
                    for n_sub_str_aa_positions in range(0,len(list_aa_profile)):
                        if list_aa_table[aa_number] == list_aa_positions[aa_positions_number][n_sub_str_aa_positions]:
                            count +=1
                    list_output[cell_position-1][aa_number] = list_output[cell_position-1][aa_number] + count
            cell_position = cell_position - list_parts_len[part_number]
            if status == 'complete':
                break
        cell_position = cell_position + 2 + list_parts_len[part_number]
    
    # write list_output into excel
    workbook_output_aa_profile = xlwt.Workbook(encoding='utf-8')
    worksheet_output_aa_profile = workbook_output_aa_profile.add_sheet('aa_profile', cell_overwrite_ok=True)
    cell_position = 1
    for part_number in range(0,7):
        worksheet_output_aa_profile.write_merge(0,0,cell_position,cell_position+list_parts_len[part_number]-1,list_parts_names[part_number])
        for aa_number in range(0,21):
            worksheet_output_aa_profile.write(aa_number+2,cell_position-1,list_aa_table[aa_number])
        for aa_location in range(1,list_parts_len[part_number]+1):
            worksheet_output_aa_profile.write(1,cell_position,aa_location)
            cell_position +=1
        cell_position = cell_position+2
    cell_position = 0
    for part_number in range(0,7):
        for aa_positions_number in range(0,list_parts_len[part_number]):
            cell_position = cell_position + 1
            for aa_number in range(0,21):
                worksheet_output_aa_profile.write(aa_number+2,cell_position,list_output[cell_position-1][aa_number]/total_num_analyzed*100)
        cell_position = cell_position + 2
    workbook_output_aa_profile.save(f'{sample}_aa_profile.xls')

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
matrix_cdo = {}
workbook_cdo = xlrd.open_workbook("codon_table.xlsx")
worksheet_cdo = workbook_cdo.sheet_by_index(0)
for ncols_cdo in range(1,5):
    for nrows_cdo in range(1,17):
        cx = worksheet_cdo.cell_value(0,ncols_cdo)
        cy = worksheet_cdo.cell_value(nrows_cdo,0)
        cz = worksheet_cdo.cell_value(nrows_cdo,5)
        aa_cdo = worksheet_cdo.cell_value(nrows_cdo,ncols_cdo)
        codon_cdo = str(cy) + str(cx) + str(cz)
        matrix_cdo[codon_cdo] = str(aa_cdo)



if __name__ == '__main__':
# Input definition, change accordingly: no change required if using default output files from CeVICA_analysis2.py #
    list_sample = ['SR3', 'GP3', 'STR'] # sample identifier for generating filepath for analysis
    feature_count_output_head = ['smaple', 'total', 'unknown', 'frame_shift', 'stop_early', 'stop_free', 'stop_free_good', 'stop_free_full_len', 'good percent', 'full_len percent']
    workbook_output_feature_count = xlwt.Workbook(encoding='utf-8')
    worksheet_output_feature_count = workbook_output_feature_count.add_sheet('freature_count', cell_overwrite_ok=True)
    col_num = 0
    row_num = 0
    for i in feature_count_output_head:
        worksheet_output_feature_count.write(row_num, col_num, i)
        col_num += 1
    for s in list_sample:
        print(f'translating and analyzing {s}')
        with open(f'{s}_trimmed.fastq','r') as f:
            db_orf_parts = sqlite3.connect(f'{s}_vhhlib_orf_parts.db')
            cursor_orf_parts = db_orf_parts.cursor()
            try:
                cursor_orf_parts.execute('CREATE TABLE parts(id INTEGER PRIMARY KEY, Full TEXT, Frame1 TEXT, CDR1 TEXT, Frame2 TEXT, CDR2 TEXT, Frame3 TEXT, CDR3 TEXT, Frame4 TEXT, CDS TEXT, DNA TEXT)')
                db_orf_parts.commit()
            except:
                pass
            status = ''
            while True:
                if status == 'complete':
                    break    
                list_DNA_seqs = []
                for n in range(0,100000):   
                    linesX4 = [f.readline().strip()] + [f.readline().strip()] + [f.readline().strip()] + [f.readline().strip()]
                    if len(linesX4[0]) > 4:
                        list_DNA_seqs += [linesX4[1]]
                    else:
                        status = 'complete' 
                        break
                pool = multiprocessing.Pool()
                result = pool.map(DNA_orf_split, list_DNA_seqs)
                for i in result:
                    cursor_orf_parts.execute('INSERT INTO parts(Full, Frame1, CDR1, Frame2, CDR2, Frame3, CDR3, Frame4, CDS, DNA) VALUES(?,?,?,?,?,?,?,?,?,?)', (i[0],i[1],i[2],i[3],i[4],i[5],i[6],i[7],i[8],i[9]))
                db_orf_parts.commit()
        row_num += 1
        col_num = 0
        save_feature_count = feature_count(s)
        for i in save_feature_count:
            worksheet_output_feature_count.write(row_num, col_num, i)
            col_num += 1
    date_stamp = datetime.date.today()
    workbook_output_feature_count.save(f'feature_count_{date_stamp}.xls')

    pool_aa_profile = multiprocessing.Pool()
    pool_aa_profile.map(aa_profile, list_sample)
