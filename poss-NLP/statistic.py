import operator
import itertools
import os

#parse a single row
def parseRow(row_str):
    #row_str = '[arabidopsis, 40, 0.1, 221, 0.0355, 25.9453125, True, False, 22, 0.3342, 117.64453125, True, False, 0, 0, 0, , , 0, 0, 0, , ]'
    row = row_str.strip().rstrip('\n')[1:-1]
    row_lst = row.split(', ')
    #print(row_lst)
    #print('%s-%s'%(row_lst[0],row_lst[1])) #key
    k = '%s-%s'%(row_lst[0],row_lst[1])
    #print([row_lst[7],float(row_lst[3]),float(row_lst[4]),float(row_lst[5])])#value1
    v1 = (row_lst[7],float(row_lst[3]),float(row_lst[4]),float(row_lst[5]))
    #print([row_lst[12],float(row_lst[8]),float(row_lst[9]),float(row_lst[10])])#value2
    v2 = (row_lst[12],float(row_lst[8]),float(row_lst[9]),float(row_lst[10]))
    return k,v1,v2,row_lst[2]

#parse a file
def parseFile(nlp_file):
    global record_dict
    fp = open(nlp_file, encoding='utf-8')
    line = fp.readline()
    #print(line)
    while (line != ''):
        line.strip()  # remove the left and right blanks
        if line.startswith('#'):  # skip the notation
            line = fp.readline()
            continue
        if line == '\n':  # skip the empty line
            line = fp.readline()
            continue
        #print(line)
        rowResult = parseRow(line)
        if rowResult[0] in record_dict:
            if rowResult[1][0] == 'False':
                record_dict[rowResult[0]]['S-in'].append(rowResult[1][1:4])
            else:
                record_dict[rowResult[0]]['S-over'].append(rowResult[1][1:4])
            if rowResult[2][0] == 'False':
                record_dict[rowResult[0]]['O-in'].append(rowResult[2][1:4])
            else:
                record_dict[rowResult[0]]['O-over'].append(rowResult[2][1:4])
        else: #initialize
            tmp_dict = {'statistic': []} #instore the statics
            if rowResult[1][0] == 'False':
                tmp_dict['S-in'] = [rowResult[1][1:4], ]
                tmp_dict['S-over'] = []
            else:
                tmp_dict['S-in'] = []
                tmp_dict['S-over'] = [rowResult[1][1:4], ]
            if rowResult[2][0] == 'False':
                tmp_dict['O-in'] = [rowResult[2][1:4],]
                tmp_dict['O-over'] = []
            else:
                tmp_dict['O-in'] = []
                tmp_dict['O-over'] = [rowResult[2][1:4],]
            record_dict[rowResult[0]] = tmp_dict
        line = fp.readline()
    fp.close()

#compute the average of values in the time limit
def calResult():
    for k,v in record_dict.items():
        record_dict[k]['statistic'] = dict()
        record_dict[k]['statistic']['CNT_S_in'] = len(v['S-in']) #the count of all groups
        record_dict[k]['statistic']['CNT_O_in'] = len(v['O-in'])

        record_dict[k]['statistic']['AVG_S_NLR'] = round( sum(list(map(lambda x: x[0], record_dict[k]['S-in']))) / record_dict[k]['statistic']['CNT_S_in'] ,2)
        record_dict[k]['statistic']['AVG_S_Time'] = round( sum(list(map(lambda x: x[1], record_dict[k]['S-in']))) / record_dict[k]['statistic']['CNT_S_in'] ,4)
        record_dict[k]['statistic']['AVG_S_Memory'] = round( sum(list(map(lambda x: x[2], record_dict[k]['S-in']))) / record_dict[k]['statistic']['CNT_S_in'] ,2)
        if record_dict[k]['statistic']['CNT_O_in'] == 0:
            record_dict[k]['statistic']['AVG_O_NLR'] = -1
            record_dict[k]['statistic']['AVG_O_Time'] = -1
            record_dict[k]['statistic']['AVG_O_Memory'] = -1
        else:
            record_dict[k]['statistic']['AVG_O_NLR'] = round( sum(list(map(lambda x: x[0], record_dict[k]['O-in']))) / record_dict[k]['statistic']['CNT_O_in'] ,2)
            record_dict[k]['statistic']['AVG_O_Time'] = round( sum(list(map(lambda x: x[1], record_dict[k]['O-in']))) / record_dict[k]['statistic']['CNT_O_in'] ,4)
            record_dict[k]['statistic']['AVG_O_Memory'] = round( sum(list(map(lambda x: x[2], record_dict[k]['O-in']))) / record_dict[k]['statistic']['CNT_O_in'] ,2)

        if record_dict[k]['statistic']['CNT_S_in'] != 30: #normally, all specific is inside the timelimit
            print('The count of %s is wrong.'%k)
        # if record_dict[k]['statistic']['CNT_O_in'] * 2 < record_dict[k]['statistic']['CNT_S_in']:
        #     print(k)
        #     print(record_dict[k])

#get the nodes used in the graph
def getNodes(): #(40,163)  (80,276)  (120,422)  (160,557)  (200,670)  (240,835)  (280,976)  (320,1098)
    global nodes_dict
    Nodes_S_NLR = {'mammalian': [], 'fission': [], 'budding': [], 'arabidopsis': [], 'thelper': [], 'tcrNet': []}
    Nodes_S_Time = {'mammalian': [], 'fission': [], 'budding': [], 'arabidopsis': [], 'thelper': [], 'tcrNet': []}
    Nodes_S_Memory = {'mammalian': [], 'fission': [], 'budding': [], 'arabidopsis': [], 'thelper': [], 'tcrNet': []}
    Nodes_O_NLR = {'mammalian': [], 'fission': [], 'budding': [], 'arabidopsis': [], 'thelper': [], 'tcrNet': []}
    Nodes_O_Time = {'mammalian': [], 'fission': [], 'budding': [], 'arabidopsis': [], 'thelper': [], 'tcrNet': []}
    Nodes_O_Memory = {'mammalian': [], 'fission': [], 'budding': [], 'arabidopsis': [], 'thelper': [], 'tcrNet': []}

    for k in record_dict.keys():
        #print(k)
        GRN_name = k.split('-')[0]
        E_cnt = int(k.split('-')[1])
        Nodes_S_NLR[GRN_name].append( (E_cnt,record_dict[k]['statistic']['AVG_S_NLR'] ))
        Nodes_S_Time[GRN_name].append( (E_cnt,record_dict[k]['statistic']['AVG_S_Time']))
        Nodes_S_Memory[GRN_name].append( (E_cnt,record_dict[k]['statistic']['AVG_S_Memory']))
        if record_dict[k]['statistic']['CNT_O_in'] * 2 >= record_dict[k]['statistic']['CNT_S_in']:
            Nodes_O_NLR[GRN_name].append((E_cnt, record_dict[k]['statistic']['AVG_O_NLR']))
            Nodes_O_Time[GRN_name].append((E_cnt, record_dict[k]['statistic']['AVG_O_Time']))
            Nodes_O_Memory[GRN_name].append((E_cnt, record_dict[k]['statistic']['AVG_O_Memory']))
    nodes_dict['S_NLR'] = Nodes_S_NLR
    nodes_dict['S_Time'] = Nodes_S_Time
    nodes_dict['S_Memory'] = Nodes_S_Memory
    nodes_dict['O_NLR'] = Nodes_O_NLR
    nodes_dict['O_Time'] = Nodes_O_Time
    nodes_dict['O_Memory'] = Nodes_O_Memory
    # print(nodes_dict['S_NLR'])
    # print(nodes_dict['S_Time'])
    # print(nodes_dict['S_Memory'])
    # print(nodes_dict['O_NLR'])
    # print(nodes_dict['O_Time'])
    # print(nodes_dict['O_Memory'])

#get the timeout records
def getTimeoutParts():
    dict_default = {200:0, 240:0, 280:0, 320:0, 360:0, 400:0, 440:0, 480:0}
    dict_timeout = {'mammalian':dict_default.copy(), 'fission':dict_default.copy(), 'budding':dict_default.copy(), 'arabidopsis':dict_default.copy(), 'thelper':dict_default.copy(), 'tcrNet':dict_default.copy() }
    for k in record_dict.keys():
        #if record_dict[k]['statistic']['CNT_O_in'] * 2 < record_dict[k]['statistic']['CNT_S_in'] :
        if record_dict[k]['statistic']['CNT_O_in'] < record_dict[k]['statistic']['CNT_S_in']:
            print(k,record_dict[k]['statistic']['CNT_O_in'],record_dict[k]['statistic']['CNT_S_in'] - record_dict[k]['statistic']['CNT_O_in'],record_dict[k]['statistic']['CNT_S_in'])
            #print(k.split('-')[0],k.split('-')[1],record_dict[k]['statistic']['CNT_S_in'] - record_dict[k]['statistic']['CNT_O_in'])
            dict_timeout[k.split('-')[0]][int(k.split('-')[1])] = record_dict[k]['statistic']['CNT_S_in'] - record_dict[k]['statistic']['CNT_O_in']
    print(dict_timeout)
    for v in dict_timeout.values():
        print(' & '.join(list(map(lambda x: str(x), v.values()))) )



def printMammalian():
    global record_dict
    fp = open(log_dir + file_name, encoding='utf-8')
    line = fp.readline()
    #print(line)
    while (line != ''):
        line.strip()  # remove the left and right blanks
        if line.startswith('#'):  # skip the notation
            line = fp.readline()
            continue
        if line == '\n':  # skip the empty line
            line = fp.readline()
            continue
        #print(line)
        rowResult = parseRow(line)
        #print(rowResult)
        #('mammalian-40', ('False', 140.0, 0.0379, 25.87109375), ('False', 19.0, 0.4192, 111.2578125), '0.0')
        E_cnt = int(rowResult[0].split('-')[1])
        if rowResult[3] == '0.0':
            #print('\multirow{4}{*}{%s} & 0 & '%E_cnt)
            s1 = '\multirow{4}{*}{%s} & 0 & '%E_cnt
        else:
            #print(' & %s & '%rowResult[3])
            s1 = ' & %s & '%rowResult[3]
        #print(' %d & %s & %s  &' % (int(rowResult[1][1]), rowResult[1][2], str(round(float(rowResult[1][3]), 2))))
        s2 = ' %d & %s & %s  &' % (int(rowResult[1][1]), rowResult[1][2], str(round(float(rowResult[1][3]), 2)))
        if rowResult[2][0] == 'False':
            #print(' %d & %s & %s \\\\ \n \cline{2-8} ' % (int(rowResult[2][1]), rowResult[2][2], str(round(float(rowResult[2][3]), 2))))
            s3 = ' %d & %s & %s \\\\ \n ' % (int(rowResult[2][1]), rowResult[2][2], str(round(float(rowResult[2][3]), 2)))
        else:
            #print(' -1 & -1 & -1 \\\\ \n \cline{2-8} ')
            s3 = ' -1 & -1 & -1 \\\\ \n  '
        if rowResult[3] == '0.8':
            s4 = '\hline '
        else:
            s4 = '\cline{2-8}  '
        print(s1,s2,s3,s4)
        #break
        line = fp.readline()
    fp.close()


if __name__ == "__main__":

    file_name = 'specific30.txt' #10&320&0.1&test.txt   10&320&0.1.txt  30&320&0.1.txt    30&480&0.1.txt Specific.txt
    #file_name = 'multiB.log'
    root_dir = os.getcwd()
    log_dir = root_dir + '/log/'
    record_dict = dict()
    nodes_dict = dict()

    parseFile(log_dir + file_name)
    #print(record_dict)
    #printMammalian()

    calResult()
    getNodes()
    print(nodes_dict)
    getTimeoutParts()
