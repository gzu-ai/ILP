# -*- coding: UTF-8 -*-

# To generate the random examples, ie., (I,O)s, from a normal logic program P

# To run as follows:
# <this-file> <program>

import numpy as np
import os  # for file reading and writting
import sys, getopt  # for argv argument
import random
import copy
# import pdb
# import math  # for ceil function
# import itertools
# import re
# from clingo.control import Control
# from guppy import hpy
# from multiprocessing import Pool, Manager
import math
import time
from functools import wraps
import memory_profiler
import json
from subprocess import Popen, PIPE
import re


#Define this decorator to measure the consumption of time and memory. This measure can cross different processes.
def fn_memory_time(function):
    @wraps(function)
    def function_inner(*args, **kwargs):
        m0 = memory_profiler.memory_usage()[0]
        t0 = time.time()
        m1,result = memory_profiler.memory_usage((function,args,kwargs),retval=True,include_children=True,max_iterations=1,max_usage=True)# its unit is in MiB
        t1 = time.time()
        timeCost = str(round(t1 - t0,4))
        memoryCost = str(m1 - m0)
        if silent_level < 6:
            print("Total running cost of %s is %s seconds and %s MB." %(function.__name__, timeCost, memoryCost))
        return timeCost, memoryCost, result
        #return result
    return function_inner

class Trans:
    input = dict()
    output = dict()

    def __init__(self, dicI, dicO):
        self.input = dicI
        self.output = dicO

class Examples:
    atom_ids = []
    atom_names = []
    Exps = []
    Is = []
    NaiveSpecificSol = None
    LastSpecificSol = None


    def printE(self):
        for e in self.Exps:
            str_input = ';'.join(list('(%s, %s)' % (key, str(value)) for key, value in e.input.items()))
            str_output = ';'.join(list('(%s, %s)' % (key, str(value)) for key, value in e.output.items()))
            str_trans = '{%s} -> {%s} ' % (str_input, str_output)
            print(str_trans)


    # suppose the I is the same.
    # the structure of O is the same. Namely, all atoms ['1', '2', '3'] occur in the list
    # flag_comp in {-1,0,1} representing {less,equal,greater}
    # Compare two set of examples. return True if meeting the desired relationship; otherwise, return False.
    def examples_comp(self, E2, atoms, flag_comp):  # examples2 is also an instance of Example
        examples1 = self.Exps
        examples2 = E2.Exps
        extra_evidence = False
        judgement = True
        for i in range(len(examples1)):
            for a in atoms:
                v1 = examples1[i].output[a]
                v2 = examples2[i].output[a]
                if flag_comp == 0: #judgement is enough to judge their equivalence
                    if v1 != v2:
                        judgement = False
                        break
                else:
                    if (v1 > v2) - (v1 < v2) == -flag_comp:
                        judgement = False
                        break
                    elif (v1 > v2) - (v1 < v2) == flag_comp:
                        extra_evidence = True
        if flag_comp != 0:#extra_evidence means desired greater exists. judgement means others are all equal
            judgement = judgement and extra_evidence
        if silent_level < 5:
            if judgement:
                print('matched!')
            else:
                print('The two examples is not desired!')
                print('************examples1*************')
                self.printE()
                print('************examples2**********')
                E2.printE()
        return judgement


class cRule:
    head = '0'
    pos = set() #set of strings
    neg = set() #set of strings
    necessity = 0

    '''
    #it's unnecessary to impletly define this construction function unless you want to change something
    # initialize an empty instance
    def __init__(self):
        self.head = '0'
        self.pos = set()
        self.neg = set()
        self.necessity = 0
    '''

    # to build the rule from a given string in rule format
    def refreshRuleFromAttributes(self, head_str, pos_set, neg_set, necessity_float):
        self.head = head_str
        self.pos = pos_set
        self.neg = neg_set
        self.necessity = necessity_float

    # to build the rule from a given string in rule format
    def refreshRuleFromLine(self, rule_str):  # (1 :- 2,not 1.,0.4)
        self.pos = set()  # initialize the attribution since many instances created consequently
        self.neg = set()
        poss_rule = rule_str.strip().rstrip('\n')[1:-1]
        r_and_n = poss_rule.split('.,')
        self.necessity = float(r_and_n[1].strip())
        h_and_b = r_and_n[0].split(':-')
        self.head = h_and_b[0].strip()
        if h_and_b[1].count(',') > 0:  # many bodies
            list_body = h_and_b[1].split(',')
        elif len(h_and_b[1].strip()) > 0:  # contain only one element #one body
            list_body = [h_and_b[1]]
        else:  # no body
            return
        for str_atom in list_body:
            str_atom = str_atom.strip()
            if str_atom.startswith('not'):  # negative body
                self.neg.add(str_atom[3:].strip())
            else:
                self.pos.add(str_atom.strip())
        #return self

    def printR(self):
        pos_body = ','.join(self.pos)
        neg_body = ','.join(list(map(lambda x: 'not %s' % x, self.neg)))
        if pos_body != '' and neg_body != '':
            pos_body = pos_body + ','
        str_tmp = '(%s :- %s%s.,%s)' % (self.head, pos_body, neg_body, str(self.necessity))
        print(str_tmp)

    # return a random subset of atom_set
    def genSubset(self, atoms_set):
        atoms_list = list(atoms_set)
        new_set_len = random.randint(0, len(atoms_list))
        if new_set_len == 0:
            return set()
        else:
            new_list = atoms_list[0:new_set_len]
            return set(new_list)

    def specificRule(self, universe_set, magnification):
        additional_atoms = universe_set.difference(self.pos).difference(self.neg)
        pos_addition = self.genSubset(additional_atoms)
        left_additional_atoms = additional_atoms.difference(pos_addition)
        neg_addition = self.genSubset(left_additional_atoms)
        self.pos = self.pos.union(pos_addition)
        self.neg = self.neg.union(neg_addition)
        self.necessity = round(random.randint(1, round(self.necessity * magnification)) / magnification, 2)

    #compute beta applicable  {'1': 0.5,'2': 0.8} is an interpretation
    def applicable(self, interp): # return 0 or a fraction
        a_input = frozenset(interp.keys())
        if self.pos.issubset(a_input) and self.neg.isdisjoint(a_input):
            beta_lst = []
            beta_lst.append(self.necessity)
            for e in self.pos:
                beta_lst.append(interp[e])
            return min(beta_lst)
        else:
            return 0

    # checking if the body is satisfied by an interpretation interp
    # represented in set
    def satisfy_body(self, interp):#removable
        for m in self.pos:
            if not (m in interp):
                return False
        for m in self.neg:
            if m in interp:
                return False
        return True



class cNLP:
    atom_names = list()
    atom_ids = list() # a list of strings which starts from 1 and ends with num_of_atoms
    str_atom_map = '' #store the map from names to ids
    num_of_atoms = 0
    Rules = []
    RomdomBackground = [] #a list of Rules
    RandomIs = []
    RandomExps = None #an instance of Examples
    magnification = 100 #set it 100. Accordingly, the necessity of formula contains only two decimal places
    sTransCnt = 0

    #Assume we only use self.Ruses when cNLP() is called without practical construction function refreshProgramFromFile followed. Then this __init__ reset self.Rules.
    def __init__(self):
        self.Rules = []
        self.magnification = 100

    def refreshProgramFromAttributes(self,atom_ids,str_atom_map):
        self.atom_ids = atom_ids
        self.str_atom_map = str_atom_map

    # read the normal logic program from the nlp_file
    def refreshProgramFromFile(self, nlp_file):
        self.Rules = [] #instances can share attributions in python. So a new object should initialize its attributions explicitly or implicitly.
        fp = open(nlp_file)
        line = fp.readline()
        while (line != ''):
            line.strip()  # remove the left and right blanks
            if line == '':#skip the empty line
                line = fp.readline()
                continue
            if line.startswith('#'):  # read the real node name and its id
                line = line.lstrip('#').lstrip().rstrip('\n')
                line1 = line.replace(',', '')
                if line1.isdigit():#atom_ids
                    #atoms_str = line.split(',')
                    self.atom_ids = line.split(',')
                    #self.atom_ids = list(map(lambda x: int(x), atoms_str))
                    self.num_of_atoms = len(self.atom_ids)
                else:#atom_names
                    self.atom_names = line.split(',')
                line = fp.readline()
                continue;
            rule = cRule()
            rule.refreshRuleFromLine(line)
            self.Rules.append(rule)
            line = fp.readline()
        fp.close()
        self.str_atom_map = '# ' + ','.join(self.atom_names) + '\n# ' + ','.join(self.atom_ids) + '\n'
        return self


    def program_extend(self,list_Rules): #extend a list of rules to return a new cNLP. Itself does not change
        newNLP = copy.deepcopy(self)
        newNLP.Rules.extend(list_Rules)
        return newNLP

    def printP(self):
        for r in self.Rules:
            r.printR()


    def print_B(self):
        for r in self.RomdomBackground:
            r.printR()

    def writeProgramTxt(self,P_list,P_path):
        str_all_rules = self.str_atom_map
        for rule in P_list:
            pos_body = ','.join(rule.pos)
            neg_body = ','.join(list(map(lambda x: 'not %s' % x, rule.neg)))
            if pos_body != '' and neg_body != '':
                pos_body = pos_body + ','
            str_tmp = '(%s :- %s%s.,%s)' % (rule.head, pos_body, neg_body, str(rule.necessity))
            str_all_rules = str_all_rules + str_tmp + '\n'
        with open(P_path, 'w') as file_e:
            file_e.write(str_all_rules)

    def writeBackgroundTxt(self, list_B):  # [('cdh1', frozenset({'cdc20'}), frozenset(), 0.9),('ubch10', frozenset({'cdh1', 'ubch10'}), frozenset(), 0.7)]
        backgrounds_txt_path = backgrounds_txt_dir + '/%s_%d_%d.lp' %(GRN_name,num_of_Is,num_of_Bs)
        self.writeProgramTxt(list_B,backgrounds_txt_path)

    def writeItselfTxt(self,txt_path):
        self.writeProgramTxt(self.Rules, txt_path)

    def writeBackgroundASP(self, list_B, t_str):
        background_asp_path = background_asp_dir + '/%s_%d_%d+%d.lp' % (GRN_name, num_of_Is, num_of_Bs, 1)
        str_all_rules = ''
        i = 0
        while i < len(list_B):
            rule = list_B[i]
            i = i + 1
            if t_str == 'optimalBound' or t_str == 'optimalLeast':
                pos_body = ''.join(list(map(lambda x: 'bg_body_pos(%d,%s).\n' % (i, x), rule.pos)))
                neg_body = ''.join(list(map(lambda x: 'bg_body_neg(%d,%s).\n' % (i, x), rule.neg)))
                head = 'bg_head(%d,%s).\n' % (i, rule.head)
                facts_temp = head + pos_body + neg_body + 'bg_necessity(%d,%d).\n' % (i, round(rule.necessity * self.magnification))
            else:
                pos_body = ''.join(list(map(lambda x: 'bgRulePos(%s,%s).\n' % ('b' + str(i), x), rule.pos)))
                neg_body = ''.join(list(map(lambda x: 'bgRuleNeg(%s,%s).\n' % ('b' + str(i), x), rule.neg)))
                head = 'bgRuleHead(%s,%s,%d).\n' % ('b' + str(i), rule.head, round(rule.necessity * self.magnification))
                facts_temp = head + pos_body + neg_body
            str_all_rules = str_all_rules + facts_temp
        with open(background_asp_path, 'w') as file_e:
            file_e.write(str_all_rules)
        return str_all_rules



    def writeSeparatedBackgroundASP(self, list_B, t_str):
        # put these asp strings into a dictionary
        dic_atom_asp = dict()
        for k in range(1, self.num_of_atoms + 1):  # initialize dic_atom_asp
            dic_atom_asp[k] = ''
        i = 0
        while i < len(list_B):
            rule = list_B[i]
            i = i + 1
            head = int(rule.head)
            if t_str == 'optimalBoundSep' or t_str == 'optimalLeastSep':
                pos_body = ''.join(list(map(lambda x: 'bg_body_pos(%d,%s).\n' % (i, x), rule.pos)))
                neg_body = ''.join(list(map(lambda x: 'bg_body_neg(%d,%s).\n' % (i, x), rule.neg)))
                necessity = 'bg_necessity(%d,%d).\n' % (i, round(rule.necessity * self.magnification))
            else:
                pos_body = ''.join(list(map(lambda x: 'bgRulePos(%s,%s).\n' % ('b' + str(i), x), rule.pos)))
                neg_body = ''.join(list(map(lambda x: 'bgRuleNeg(%s,%s).\n' % ('b' + str(i), x), rule.neg)))
                necessity = 'bgRuleHead(%s,%d).\n' % ('b' + str(i), round(rule.necessity * self.magnification))
            str_rule = pos_body + neg_body + necessity
            dic_atom_asp[head] = dic_atom_asp[head] + str_rule
        # write this dictionary into separated files
        for k in range(1, self.num_of_atoms + 1):
            background_asp_path = background_asp_dir + '/%s_%d_%d+%d.lp' % (GRN_name, num_of_Is, num_of_Bs, k)
            with open(background_asp_path, 'w') as file_b:
                file_b.write(dic_atom_asp[k])
        return dic_atom_asp

    '''
        It only generate about 10% specific rules as background rule. Namely, num_of_Bs =  len(self.Rules) / 10 < len(self.Rules)
    '''
    def genBackground(self, ratio):
        global num_of_Bs
        num_of_Bs = math.ceil(len(self.Rules) * ratio)
        self.RomdomBackground = []
        list_idxes = random.sample(range(0, len(self.Rules)), num_of_Bs)
        for i in list_idxes:
            new_rule = copy.deepcopy(self.Rules[i])
            new_rule.specificRule(set(self.atom_ids),self.magnification)
            self.RomdomBackground.append(new_rule)


    # compute T_P^d
    # Here we assume the logic program is normal
    # I: an interpretation in dict format, e.g., {'1': 0.3,'3': 0.4}
    def consequence(self, I_dic):  # return another interpretation including possibilistic atoms whose necessity is 0
        O_dic = dict()
        for k in self.atom_ids:
            O_dic[k] = 0
        for r in self.Rules:
            if O_dic.get(r.head) < r.applicable(I_dic):
                O_dic[r.head] = r.applicable(I_dic)
        return O_dic

    def test_consequence(self,interp):
        print(self.consequence(interp))

    def calTransitions(self,Is):
        exampers = list()
        for possInput in Is:
            possOutput = self.consequence(possInput)
            trans = Trans(possInput,possOutput)
            exampers.append(trans)
        return exampers

    def genExamples(self):
        self.RandomExps = Examples()
        self.RandomExps.Exps = self.calTransitions(self.RandomIs)
        self.RandomExps.Is = self.RandomIs
        self.RandomExps.atom_ids = self.atom_ids
        self.RandomExps.atom_names = self.atom_names

    def singleTransCnt(self):
        total = 0
        for trans in self.RandomExps.Exps:
            for v in trans.output.values():
                if v > 0:
                    total = total + 1
        self.sTransCnt = total




    def writeExamplesTxt(self,list_E): #{(a, 0.3); (s, 0.4)} -> {(s, 0.1); (b, 0.4)}
        examples_txt_path = examples_txt_dir + '/%s_%d.lp' %(GRN_name,num_of_Is)
        str_all_exampers = self.str_atom_map
        for e in list_E:
            str_input = ';'.join(list('(%s, %s)'%(key,str(value)) for key, value in e.input.items()))
            str_output = ';'.join(list('(%s, %s)' % (key, str(value)) for key, value in e.output.items()))
            str_trans = '{%s} -> {%s} '%(str_input,str_output)
            #print(str_trans)
            str_all_exampers = str_all_exampers + str_trans + '\n'
        with open(examples_txt_path, 'w') as file_e:
            file_e.write(str_all_exampers)



    def writeExamplesASP(self,list_E, t_str): #{(a, 0.3); (s, 0.4)} -> {(s, 0.1); (b, 0.4)}
        examples_asp_path = examples_asp_dir + '/%s_%d+%d.lp' % (GRN_name, num_of_Is, 1)
        str_all_exampers = ''
        idx = 1
        for i in range(0, len(list_E)):
            facts_temp = ''
            trans = list_E[i]
            for k_i, v_i in trans.input.items():
                if t_str == 'optimalBound' or t_str == 'optimalLeast':
                    facts_temp = facts_temp + 'input(%d .. %d,%s,%d).\n' % (idx, idx + self.num_of_atoms - 1, k_i, round(v_i * self.magnification))
                else:
                    facts_temp = facts_temp + 'input(%s,%s,%d).\n' % ('e' + str(i + 1), k_i, round(v_i * self.magnification))
            for k_i, v_i in trans.output.items():
                if t_str == 'optimalBound' or t_str == 'optimalLeast':
                    facts_temp = facts_temp + 'output(%d,%s,%d).\n' % (idx, k_i, round(v_i * self.magnification))
                    idx = idx + 1
                else:
                    if v_i > 0:
                        facts_temp = facts_temp + 'output(%s,%s,%d).\n' % ('e' + str(i + 1), k_i, round(v_i * self.magnification))
            str_all_exampers = str_all_exampers + facts_temp
        with open(examples_asp_path, 'w') as file_e:
            file_e.write(str_all_exampers)
        return str_all_exampers



    def writeSeparatedExamplesASP(self,list_E, t_str): #[{(a, 0.3); (s, 0.4)} -> {(s, 0.1); (b, 0.4)}  ]，list_E的输出部分包含了0
        # put these asp strings into a dictionary
        dic_atom_asp = dict()
        for i in range(1, self.num_of_atoms + 1): #initialize dic_atom_asp
            dic_atom_asp[i] = ''
        facts_input = '' #the public input facts
        for i in range(1, len(list_E) + 1):
            trans = list_E[i - 1]
            for k_i, v_i in trans.input.items():
                if t_str == 'optimalBoundSep' or t_str == 'optimalLeastSep':
                    facts_input = facts_input + 'input(%d,%s,%d).\n' % (i, k_i, round(v_i * self.magnification))
                else:
                    facts_input = facts_input + 'input(%s,%s,%d).\n' % ('e' + str(i), k_i, round(v_i * self.magnification))
            for j in range(1, self.num_of_atoms + 1):
                outN = trans.output[str(j)]
                if t_str == 'optimalBoundSep' or t_str == 'optimalLeastSep':
                    facts_output = 'output(%d,%d).\n' % (i, round(outN * self.magnification))
                else:
                    facts_output = 'output(%s,%d).\n' % ('e' + str(i), round(outN * self.magnification))
                dic_atom_asp[j] = dic_atom_asp[j] + facts_output
        for i in range(1, self.num_of_atoms + 1): #append facts_input into each elements of dic_atom_asp
            dic_atom_asp[i] = dic_atom_asp[i] + facts_input
        #write this dictionary into separated files
        for i in range(1, self.num_of_atoms + 1):
            examples_asp_path = examples_asp_dir + '/%s_%d+%d.lp' % (GRN_name, num_of_Is, i)
            with open(examples_asp_path, 'w') as file_e:
                file_e.write(dic_atom_asp[i])
        return dic_atom_asp

    def writeAuxiliaryASP(self, t_str):
        pauxiliary_path = auxiliary_dir + '/%s.lp' % GRN_name
        facts = 'universe(1 .. %s). \n' %self.num_of_atoms
        if t_str[0:12] == 'optimalBound':
            facts = facts + 'top(%d). \n' % self.magnification
        with open(pauxiliary_path, 'w') as file_e:
            file_e.write(facts)
        return facts



    # get the possibilistic interpretation from an binary  integer  representation
    def getPossInterp(self, integ):
        I_dic = {}
        for i in range(1, self.num_of_atoms + 1):
            if integ & 0b1 == 1:
                I_dic[str(i)] = round(random.random(),2)
            integ = integ >> 1
        return I_dic


    # generate num_of_Is random I for the <I,J> examples
    def genRandomIs(self, num_of_Is):
        self.RandomIs = list()
        list_integer = list(np.random.randint(0, 2**self.num_of_atoms + 1, num_of_Is, "int64"))
        for i in list_integer:
            possI = self.getPossInterp(i)
            self.RandomIs.append(possI)

    # deciding if the set J1 is a subset of J2 that are represented in integers
    def subset(self, J1, J2):
        if J1 & J2 == J1:
            return True
        else:
            return False

    # compute the number bodies that are satisfied by an interpretation I
    # which is represented in set
    def satisfied_bodys(self, Inp):
        s = 0
        sat_heads = set()
        for rule in self.Rules:
            if rule.head in sat_heads: continue
            if rule.satisfy_body(Inp):
                s = s + 1
                sat_heads.add(rule.head)
        return s

class cGRN:
    str_atom_map = ''
    NLP = None

    def genNecessityList(self,sampleNo):
        mu = 0.7
        sigma = 0.15
        np.random.seed(0)
        samples = np.random.normal(mu, sigma, sampleNo)
        #print(samples)
        i = 0
        while i < len(samples):
            if samples[i] > 1:
                samples[i] = 1
            elif samples[i] <= 0:
                samples[i] = 0.01
            else:
                samples[i] = round(samples[i], 2)
            i = i + 1
        # print(samples)
        return samples

    def refreshNecessities(self):
        necessities = self.genNecessityList(len(self.NLP.Rules))
        i = 0
        while i < len(necessities):
            self.NLP.Rules[i].necessity = necessities[i]
            i = i + 1


    def getRuleFromGRN(self,rule_str):
        rule = cRule()
        rule.pos = set()
        rule.neg = set()
        rule.necessity = 0
        [h, b] = rule_str.split(":-")
        rule.head = h
        if b.count(',') > 0:
            lb = b.split(',')
            lbint = list(map(int, lb))
        else:
            lbint = [int(b)]
        for atom in lbint:
            if atom > 0:
                rule.pos.add(str(atom))
            if atom < 0:
                rule.neg.add(str(abs(atom)))
        return rule



class CMDsClingo:
    cmd_list = list()
    atom_ids = list()
    magnification = 1
    resultsDict = dict()
    modelsDict = dict()
    solsList = list() # a list of cRules
    timeoutFlag = False
    optimumFlag = True
    ruleCnt = 0


    def __init__(self, lst, atom_ids, magnification):
        self.cmd_list = lst
        self.atom_ids = atom_ids
        self.magnification = magnification
        timeoutFlag = False
        optimumFlag = True

    @fn_memory_time
    def exeCMDsParallel(self,t_str):
        ResultList = dict()
        totalParts = 1
        if t_str[-3:] == 'Sep':
            totalParts = len(self.atom_ids)
        for i in range(0,totalParts):
            cmd = self.cmd_list[i]
            child = Popen(cmd, stdin=PIPE, stdout=PIPE, stderr=PIPE, shell=True)
            ResultList[str(i+1)] = child
        for k,v in ResultList.items():
            self.resultsDict[k] = v.communicate()[0].decode("gbk")


    def preParseModels(self, res):
        pattern = '^Models\s*:\s*(\d+)(\+?)\s'
        searchObj = re.search(pattern, res, re.M)
        num = searchObj.group(1)
        add = searchObj.group(2)
        return num, add

    def preParseAnswer(self, res, t_str):
        if t_str[0:8] == 'specific':
            pattern = 'Answer:\s*\d*\s*([A-Za-z_,()0-9 ]*)\s*SATISFIABLE'
        else:
            pattern = 'Answer:\s*\d*\s*([A-Za-z_,()0-9 ]*)\s*Optimization:'
        searchObj = re.search(pattern, res, re.M)
        AS = searchObj.group(1)
        return AS

    def preParseOptimum(self, res):
        pattern = '\s+Optimum\s*:\s*(\w+)\s'
        searchObj = re.search(pattern, res, re.M)
        opt = searchObj.group(1)
        return opt

    def preParseCount(self, res):
        pattern = '\s+Optimum\s*:\s*(\d+)\s'
        searchObj = re.search(pattern, res, re.M)
        opt = searchObj.group(1)
        return opt



    def parseResult(self, resultK,t_str):
        resultDict = dict()
        resultDict['error'] = False
        resultDict['AS'] = ''
        resultDict['optmum'] = True
        resultDict['timeout'] = False
        #print('********resultK******', resultK)
        if resultK == '':
            #print('******ERROR********')
            resultDict['error'] = True
        else:
            #print('others')
            modelInf = self.preParseModels(resultK)
            if modelInf[0] != '0':
                #print('******StableModel********')
                resultDict['AS'] = self.preParseAnswer(resultK,t_str)
                if t_str[0:8] == 'specific':
                    resultDict['optmum'] = False
                    resultDict['countR'] = resultDict['AS'].count('lastRuleHead(')
                else:
                    ptInf = self.preParseOptimum(resultK)
                    if ptInf != 'yes':
                        # print('******Not OPTIMUM********')
                        resultDict['optmum'] = False
            if modelInf[1] == '+':
                #print('******Timeout********')
                resultDict['timeout'] = True
        return resultDict

    def getModels(self,t_str):
        for k,v in self.resultsDict.items():
            self.modelsDict[k] = self.parseResult(v,t_str)
        self.timeoutFlag = list(map(lambda x: x['timeout'], self.modelsDict.values())).count(True) > 0
        self.optimumFlag = not list(map(lambda x: x['optmum'], self.modelsDict.values())).count(False) > 0
        def sumR(lst):
            total = 0
            for i in lst:
                total = total + i
            return total
        self.ruleCnt = sumR(list(map(lambda x: x['countR'], self.modelsDict.values())))



    def parseModel(self, stableModel, head_atom, magnification,t_str):
        def addRule(rule_head, rule_pos, rule_neg, rule_N):
            nonlocal LastOptimalSol
            crule = cRule()
            crule.refreshRuleFromAttributes(rule_head, rule_pos, rule_neg, rule_N)
            LastOptimalSol.append(crule)
        LastOptimalSol = list()  # a list of Rules whose head is head_atom
        model = stableModel.split(' ')
        if t_str[0:12] == 'optimalBound':
            for rule in list(filter(lambda x: x.startswith('last_rule('), model)):
                parameters = rule[10:-1].split(',')
                index_start = 10 + len(parameters[0])
                if t_str[-3:] == 'Sep':
                    rule_head = head_atom
                    rule_N = round(int(parameters[1]) / magnification, 2)
                else:
                    rule_head = str(parameters[1])
                    rule_N = round(int(parameters[2]) / magnification, 2)
                rule_pos = frozenset(list(map(lambda x: x[index_start:-1], list(filter(lambda x: x.startswith('rule_pos(%s,' % parameters[0]), model)))))
                rule_neg = frozenset(list(map(lambda x: x[index_start:-1], list(filter(lambda x: x.startswith('rule_neg(%s,' % parameters[0]), model)))))
                addRule(rule_head, rule_pos, rule_neg, rule_N)
        else:
            for rule in list(filter(lambda x: x.startswith('lastRuleHead('), model)):
                id = rule[13:-1]
                parameters = id.split(',')
                index_start = len('lastRulePos(') + len(id) + 1
                if t_str[-3:] == 'Sep':
                    rule_head = head_atom
                    rule_N = round(int(parameters[1]) / magnification, 2)
                else:
                    rule_head = str(parameters[1])
                    rule_N = round(int(parameters[2]) / magnification, 2)
                rule_pos = frozenset(list(map(lambda x: x[index_start:-1], list(filter(lambda x: x.startswith('lastRulePos(%s,' % id), model)))))
                rule_neg = frozenset(list(map(lambda x: x[index_start:-1], list(filter(lambda x: x.startswith('lastRuleNeg(%s,' % id), model)))))
                addRule(rule_head, rule_pos, rule_neg, rule_N)
        return LastOptimalSol


    def getRules(self,t_str):  # use a list of cRules to store the solutions
        self.solsList = list()
        if t_str[-3:] == 'Sep':
            for id in self.atom_ids:
                self.solsList.extend(self.parseModel(self.modelsDict[id]['AS'], id, self.magnification, t_str))
        else:
            self.solsList.extend(self.parseModel(self.modelsDict['1']['AS'],'', self.magnification, t_str))

def record_singleton(type_str,time,ram,len,match,timeout,stCnt,lenCnt):
    global record_list
    record_list.extend([type_str, len, time, ram, match, timeout, stCnt, lenCnt])


def gen_E_B_H_check_new_compare():
    #公用部分
    global search_type_list
    primary_program_path = primary_program_dir + '/%s.lp' % GRN_name
    NLP = cNLP()
    NLP.refreshProgramFromFile(primary_program_path)
    NLP.genRandomIs(num_of_Is)
    NLP.genExamples()
    NLP.genBackground(Bs_ratio)

    NLP.singleTransCnt()
    sT_Cnt = NLP.sTransCnt

    LastSol = None
    ExamplesASP = None
    BackgroundASP = None
    ExamplesASPdic = None
    BackgroundASPdic = None
    AuxiliaryASP = None
    match_flag = None
    timeout_flag = None
    rule_Cnt = 0


    def facts(t_str):
        nonlocal ExamplesASP, BackgroundASP, AuxiliaryASP, ExamplesASPdic, BackgroundASPdic, AuxiliaryASP,LastSol
        if t_str[-3:] == 'Sep':
            ExamplesASPdic = NLP.writeSeparatedExamplesASP(NLP.RandomExps.Exps, t_str)
            BackgroundASPdic = NLP.writeSeparatedBackgroundASP(NLP.RomdomBackground, t_str)
        else:
            ExamplesASP = NLP.writeExamplesASP(NLP.RandomExps.Exps, t_str)
            BackgroundASP = NLP.writeBackgroundASP(NLP.RomdomBackground, t_str)
        AuxiliaryASP = NLP.writeAuxiliaryASP(t_str)
        LastSol = cNLP()
        LastSol.refreshProgramFromAttributes(NLP.atom_ids, NLP.str_atom_map)

    def getOptimalParallel(t_str):
        nonlocal timeout_flag
        nonlocal rule_Cnt
        global decoratorList
        def getCMDs(t_str):
            global type_file_dict
            cmd_list = list()
            exe_path = 'clingo' # 'c:/Users/86188/.vscode/extensions/ffrankreiter.answer-set-programming-language-support-0.7.0/clingo_win.exe'
            encoding_path = root_dir + '/%s.lp' % type_file_dict[t_str]
            parameter_path = '--outf=0 --verbose=1 --parallel-mode=2,compete --time-limit=3600 --models 0 --opt-mode=opt --quiet=1' #--time-limit=1800
            totalParts = 1
            if t_str[-3:] == 'Sep':
                totalParts = len(NLP.atom_ids)
            for k in range(1, totalParts + 1): #如果是sep则有多个文件，否则只有一个。但是格式一样
                background_asp_path = background_asp_dir + '/%s_%d_%d+%d.lp' % (GRN_name, num_of_Is, num_of_Bs, k)
                examples_asp_path = examples_asp_dir + '/%s_%d+%d.lp' % (GRN_name, num_of_Is, k)
                pauxiliary_path = auxiliary_dir + '/%s.lp' % GRN_name
                cmd_str = '%s %s %s %s %s %s ' % (exe_path, encoding_path, background_asp_path, examples_asp_path, pauxiliary_path, parameter_path)
                cmd_list.append(cmd_str)
            return cmd_list

        commandList = getCMDs(t_str)
        if silent_level < 3:
            print(commandList)
        myCMDS = CMDsClingo(commandList, NLP.atom_ids, NLP.magnification)
        decoratorList = myCMDS.exeCMDsParallel(t_str)  # the memory consumption will contain it's children's part because of  include_children=True along memory_profiler.memory_usage()
        myCMDS.getModels(t_str)

        # myCMDS.getRules(t_str)
        # if silent_level < 4:
        #     print(myCMDS.modelsDict)
        #     #print(myCMDS.solsList)
        # LastSol.Rules = myCMDS.solsList

        timeout_flag = str(myCMDS.timeoutFlag)
        rule_Cnt = myCMDS.ruleCnt

    def transMatch():
        nonlocal LastSol, match_flag
        H = LastSol.program_extend(NLP.RomdomBackground)
        E2 = Examples()
        E2.Exps = H.calTransitions(NLP.RandomIs)
        match_flag = NLP.RandomExps.examples_comp(E2, NLP.atom_ids, 0)

    def recordResult(t_str):
        global decoratorList
        record_singleton(t_str, decoratorList[0], decoratorList[1], len(LastSol.Rules), str(match_flag), timeout_flag,sT_Cnt, rule_Cnt)
        if silent_level < 4:
            if t_str[0:7] == 'optimal':
                LastSol.writeItselfTxt(hypotheses_dir + '/optimal' + '/%s_%d_%d.lp' % (GRN_name, num_of_Is, num_of_Bs))
            else:
                LastSol.writeItselfTxt(hypotheses_dir + '/specific' + '/%s_%d_%d.lp' % (GRN_name, num_of_Is, num_of_Bs))



    def computeMatchRecord(t_str):
        facts(t_str)
        getOptimalParallel(t_str)

        # transMatch()
        recordResult(t_str)


    # 私有部分

    types = search_type_list.split('+')
    for type in types:
        computeMatchRecord(type)

    #公用部分
    def intermediatePrint():
        print('***************newProgram***************')
        NLP.printP()
        print('***************RandomIs***************')
        print(NLP.RandomIs)
        print('***************RandomExps***************')
        NLP.RandomExps.printE()
        print('***************RomdomBackground***************')
        NLP.print_B()
        print('******LastSol******')
        LastSol.printP()


    def intermediateWrite():
        NLP.writeExamplesTxt(NLP.RandomExps.Exps)
        NLP.writeBackgroundTxt(NLP.RomdomBackground)

    if silent_level < 1:
        intermediatePrint()
    if silent_level < 2:
        intermediateWrite()



def printExecutionSpecification(pName):
    global type_file_dict
    print("Instruction: %s <GRN_name> <num_of_Is> <Bs_ratio> <search_type_list>  <silent_level> <file_name>" %pName)
    print('The value range of these parameters are as follows:')
    print('GRN_name in {mammalian, fission, budding, arabidopsis, thelper, tcrNet}')
    print('num_of_Is in {40, 80, 120, 160, 200, 240, 280, 320}')
    print('Bs_ratio in {0, 0.1, 0.4, 0.8}')
    print('search_type_list is a list concatenated by plus, such as "optimalBound+optimalConstraintR", of searching type.')
    print('The collection of candidate types is %s.' % str(list(type_file_dict.keys())))
    print('silent_level in {0, 1, 2, 3, 4, 5, 6}')
    print('parameter file_name assigns the place where log the experimental record')
    '''
    6 + write necessary results
    5 + print time;memory
    4 + print check result
    3 + print modelsDict for debug.
    2 + write ASP; print commandList for debug otherwise.
    1 + intermediateWrite
    0 + intermediatePrint
    '''

def externalCall():
    global GRN_name,num_of_Is,Bs_ratio,search_type_list,silent_level,file_name,type_file_dict
    wrongParameter = False

    def checkSearchTypeList(t_list):
        types = t_list.split('+')
        if set(types).issubset(set(type_file_dict.keys())):
            return True
        else:
            return False


    if len(sys.argv) != 7:
        printExecutionSpecification(sys.argv[0])
        sys.exit(0)

    if sys.argv[1] not in {'mammalian', 'fission', 'budding', 'arabidopsis', 'thelper', 'tcrNet'}:
        print('Instruction:GRN_name in {mammalian, fission, budding, arabidopsis, thelper, tcrNet}')
        wrongParameter = True
    if sys.argv[2] not in {'40', '80', '120', '160', '200', '240', '280', '320'}:
        print('Instruction:num_of_Is in {40, 80, 120, 160, 200, 240, 280, 320}')
        wrongParameter = True
    if sys.argv[3] not in {'0', '0.1', '0.4', '0.8'}:
        print('Instruction:Bs_ratio in {0, 0.1, 0.4, 0.8}')
        wrongParameter = True
    if not checkSearchTypeList(sys.argv[4]):
        print('search_type_list is a list concatenated by plus, such as "optimalBound+optimalConstraintR", of searching type.')
        print('The collection of candidate types is %s.' % str(list(type_file_dict.keys())))
        wrongParameter = True
    # if sys.argv[4] not in {'specific', 'optimal_batch', 'optimal_sep', 'optimal_parallel', 'overall', 'compare'}:
    #     print('Instruction:search_type_list in {specific, optimal_batch, optimal_sep, optimal_parallel, overall, compare}')
    #     wrongParameter = True
    if sys.argv[5] not in {'0', '1', '2', '3', '4', '5', '6'}:
        print('Instruction:silent_level in {0, 1, 2, 3, 4, 5, 6}')
        wrongParameter = True

    if wrongParameter:
        print("Instruction: %s <GRN_name> <num_of_Is> <Bs_ratio> <search_type_list>  <silent_level>" % sys.argv[0])
        sys.exit(0)

    GRN_name = sys.argv[1]
    num_of_Is = int(sys.argv[2])
    Bs_ratio = round(float(sys.argv[3]), 1)
    search_type_list = sys.argv[4]
    silent_level = int(sys.argv[5])
    file_name = sys.argv[6]

if __name__ == "__main__":
    GRN_name = 'thelper'
    num_of_Is = 3000 #{40, 80, 120, 160, 200, 240, 280, 320, 360, 400, 440, 480}
    Bs_ratio = 0.1 #the actual num_of_Bs is determined by Bs_ratio * (NLP.Rules)
    #search_type_list = 'optimalBound+optimalConstraintR+optimalConstraintSP+specificTrimed+specificNaive+optimalApplicableR+' \
    #                  'optimalBoundSep+optimalConstraintRSep+optimalConstraintSPSep+specificTrimedSep+specificNaiveSep' #refer to type_file_dict
    #search_type_list = 'specificTrimedSep+optimalConstraintRSep'
    #search_type_list = 'optimalBound+optimalConstraintR'
    #search_type_list = 'optimalBound+optimalApplicableR+optimalConstraintR'
    #search_type_list = 'optimalBound+optimalLeast+optimalLeastSep+optimalApplicableR+optimalApplicableSP'
    #search_type_list = 'optimalApplicableR+optimalApplicableA+optimalApplicableS+optimalApplicableSP' #optimalApplicableSP
    #search_type_list = 'optimalBoundSep+optimalConstraintRSep'
    #search_type_list = 'specificNaiveSep+specificTrimedSep'
    #search_type_list = 'specificTrimed+specificTrimedSep'
    #search_type_list = 'specificTrimedSep+optimalLeastSep'
    #search_type_list = 'optimalBoundSep'
    search_type_list = 'specificTrimedSep'
    silent_level = 5
    file_name = 'log.log'

    type_file_dict = {'optimalBound':'MHSC','optimalLeast':'MHSC-cut-uniform','optimalConstraintR':'R-ILTP','optimalConstraintSP':'SP-ILTP','specificTrimed':'ILTP',
                      'specificNaive':'ILTP-naive','optimalApplicableR':'R-ILTP-app','optimalApplicableA':'A-ILTP-app','optimalApplicableS':'S-ILTP-app',
                      'optimalApplicableSP':'SP-ILTP-app', 'optimalBoundSep':'MHSC-head', 'optimalLeastSep':'MHSC-cut-uniform-head',
                      'optimalConstraintRSep':'R-ILTP-head','optimalConstraintSPSep':'SP-ILTP-head','specificTrimedSep':'ILTP-head',
                      'specificNaiveSep':'ILTP-naive-head'}
    record_list = None

    #comment this function if doing internal debug or invocation
    #externalCall()

    # GRN_name  is an element from GRN_name_list
    GRN_name_list = ['thelper', 'tcrNet'] #'mammalian', 'fission', 'budding', 'arabidopsis', 'thelper', 'tcrNet'
    #Is_amount_list = [40, 80, 120, 160, 200, 240, 280, 320, 360, 400, 440, 480] #[10, 20, 40, 80, 160, 640]  # the num_of_Is is an element from Is_amount_list
    Is_amount_list = [500, 1000, 1500, 2000, 2500, 3000, 3500, 4000, 4500] #  [500, 1500, 2000, 2500] [3000, 3500, 4000, 4500, 5000] [720, 800]
    Bs_ratio_list = [0.1, 0.4, 0.8]
    num_of_Bs = 0 
    times = 10  # run 10 times and get the average. 2 fro margin.


    root_dir = os.getcwd()
    #print(root_dir)
    source_GRN_dir = root_dir + '/GRN'
    primary_program_dir = root_dir + '/primary_program'
    examples_txt_dir = root_dir + '/examples_txt'
    examples_asp_dir = root_dir + '/examples_asp'
    backgrounds_txt_dir = root_dir + '/background_txt'
    background_asp_dir = root_dir + '/background_asp'
    auxiliary_dir = root_dir + '/auxiliary'
    hypotheses_dir = root_dir + '/hypotheses'    

    def exe():
        gen_E_B_H_check_new_compare()
        with open(root_dir + '/log/' + file_name, 'a+', encoding='utf-8') as file:
            file.write(json.dumps(record_list).replace('\"', '') + '\n')
        file.close()


    # record_list = [GRN_name, num_of_Is, Bs_ratio]
    # exe()

    #repeat
    # for i in range(times):
    #     record_list = [GRN_name, num_of_Is, Bs_ratio]
    #     exe()

    #tcrNet
    for i in Is_amount_list:
        num_of_Is = i
        for i in range(times):
            record_list = [GRN_name, num_of_Is, Bs_ratio]
            exe()

    #optimalApplicableR+optimalApplicableSP
    # for g in GRN_name_list:
    #     GRN_name = g
    #     for i in range(times):
    #         record_list = [GRN_name, num_of_Is, Bs_ratio]
    #         exe()

    #specificTrimedSep + specificNaiveSep
    # for i in Is_amount_list:
    #     num_of_Is = i
    #     for g in GRN_name_list:
    #         GRN_name = g
    #         for i in range(times):
    #             record_list = [GRN_name, num_of_Is, Bs_ratio]
    #             exe()



    #optimalBoundSep+optimalConstraintRSep
    # for g in GRN_name_list:
    #     GRN_name = g
    #     for b in Bs_ratio_list:
    #         Bs_ratio = b
    #         #for i in range(times):
    #         record_list = [GRN_name, num_of_Is, Bs_ratio]
    #         exe()




