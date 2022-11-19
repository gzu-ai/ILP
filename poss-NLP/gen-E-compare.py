# -*- coding: UTF-8 -*-

# To generate the random examples, ie., (I,O)s, from a normal logic program P

# To run as follows:
# <this-file> <program>

import pdb
import numpy as np
import math  # for ceil function
import os  # for file reading and writting
import sys, getopt  # for argv argument
import random
import copy
import itertools
import re
from clingo.control import Control
import math
#from guppy import hpy
import time
from functools import wraps
import memory_profiler
from multiprocessing import Pool, Manager
import json
from subprocess import Popen, PIPE
import re


#define this decorator to measure the execution time of a function
def fn_timer(function):
    @wraps(function)
    def function_timer(*args, **kwargs):
        t0 = time.time()
        result = function(*args, **kwargs)
        t1 = time.time()
        if silent_level < 4:
            print("Total time running %s: %s seconds" %(function.__name__, str(t1 - t0)))
        return result
    return function_timer

#define this decorator to measure the memory consumption of a function
def fn_memory(function):
    @wraps(function)
    def function_memory(*args, **kwargs):
        if silent_level < 4:
            print ('Memory (Before): {} MB'.format(memory_profiler.memory_usage()))
        result = function(*args, **kwargs)
        if silent_level < 4:
            print('Memory (After): {} MB'.format(memory_profiler.memory_usage()))
        return result
    return function_memory

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

    def __init__(self):
        atom_ids = []
        atom_names = []
        Exps = []
        Is = []

    # read the possibilistic normal logic program from the nlp_file
    def refreshExampleFromFile(self, nlp_file):
        fp = open(nlp_file)
        line = fp.readline()
        while (line != ''):
            line.strip()  # remove the left and right blanks
            IO = line.rstrip('\n')
            if IO == '':#skip the empty line
                line = fp.readline()
                continue
            if IO.startswith('#'):#skip the str map
                str_map = IO.lstrip('#').lstrip()
                line1 = str_map.replace(',', '')
                if line1.isdigit():  # atom_ids
                    self.atom_ids = str_map.split(',')
                else:  # atom_names
                    self.atom_names = str_map.split(',')
                line = fp.readline()
                continue
            I_and_O = IO.split('->')
            I_dic = dict()
            I_atoms = I_and_O[0].strip().lstrip('{').rstrip('}').split(';')
            # logging.debug(I_atoms)
            for atom in I_atoms:
                if atom == '':
                    break
                matchObj = re.match("^\((\w+), (.+)\)$", atom.strip(), re.I)
                I_dic[matchObj.group(1)] = float(matchObj.group(2))
            O_dic = dict()
            O_atoms = I_and_O[1].strip().lstrip('{').rstrip('}').split(';')
            for atom in O_atoms:
                if atom == '':
                    break
                matchObj = re.match("^\((\w+), (.+)\)$", atom.strip(), re.I)
                O_dic[matchObj.group(1)] = float(matchObj.group(2))
            trans = Trans(I_dic, O_dic)
            self.Exps.append(trans)
            self.Is.append(I_dic)
            line = fp.readline()
        fp.close()

    def printE(self):
        for e in self.Exps:
            str_input = ';'.join(list('(%s, %s)' % (key, str(value)) for key, value in e.input.items()))
            str_output = ';'.join(list('(%s, %s)' % (key, str(value)) for key, value in e.output.items()))
            str_trans = '{%s} -> {%s} ' % (str_input, str_output)
            print(str_trans)

    def SpecificRule(self,I_dic,id,id_N): #self.SpecificRule(trans.input,id,O_dic[id])
        body_positive = frozenset( key for key, value in I_dic.items() if float(value) >= id_N)
        body_negative = frozenset(set(self.atom_ids).difference(frozenset(I_dic.keys()))) #list list
        rule = cRule()
        rule.refreshRuleFromAttributes(id,body_positive,body_negative,id_N) #refreshRuleFromAttributes(self, head_str, pos_set, neg_set, necessity_float):
        return rule


    def computeNaiveSpecificSol(self):
        self.NaiveSpecificSol = cNLP()
        self.NaiveSpecificSol.str_atom_map = '# ' + ','.join(self.atom_names) + '\n# ' + ','.join(self.atom_ids) + '\n'
        for trans in self.Exps:
            for O_key,O_value in trans.output.items():
                if O_value > 0:
                    rule = self.SpecificRule(trans.input,O_key,O_value)
                    self.NaiveSpecificSol.Rules.append(rule)

    #the efficientcy is significantly improved than computeLastSpecificSol, since the small loop and tiny filter.
    def computeLastSpecificSol(self, B_rules):  # self.NaiveSpecificSol is a NLP. B_rules is a list of rules.
        self.LastSpecificSol = cNLP()
        self.LastSpecificSol.atom_ids = self.atom_ids
        self.LastSpecificSol.str_atom_map = '# ' + ','.join(self.atom_names) + '\n# ' + ','.join(self.atom_ids) + '\n'
        group_result = itertools.groupby(self.NaiveSpecificSol.Rules, lambda x: (x.head, x.pos, x.neg))
        for key, group in group_result:
            max_necessity = max(list(x.necessity for x in group))
            if len(list(filter(lambda x: x.head == key[0] and x.pos == key[1] and x.neg == key[2] and x.necessity >= max_necessity, B_rules))) == 0:
                new_rule = cRule()
                new_rule.refreshRuleFromAttributes(key[0], key[1], key[2], max_necessity)
                self.LastSpecificSol.Rules.append(new_rule)

    #@fn_memory
    #@fn_timer
    @fn_memory_time
    def computeLastSpecificSol_bk(self,B_rules): #self.NaiveSpecificSol is a NLP. B_rules is a list of rules.
        self.LastSpecificSol  = cNLP()
        self.LastSpecificSol.atom_ids = self.atom_ids
        self.LastSpecificSol.str_atom_map = '# ' + ','.join(self.atom_names) + '\n# ' + ','.join(self.atom_ids) + '\n'
        for rule in self.NaiveSpecificSol.Rules:
            great_inter_P = list(filter(lambda x: x.head == rule.head and x.pos == rule.pos and x.neg == rule.neg and x.necessity > rule.necessity, self.NaiveSpecificSol.Rules))
            great_inter_B = list(filter(lambda x: x.head == rule.head and x.pos == rule.pos and x.neg == rule.neg and x.necessity >= rule.necessity, B_rules))
            # rule is the greatest in P; rule is also greater than any b in B; the same rule has not been added into P_new
            if len(great_inter_P) == 0 and len(great_inter_B) == 0 and rule not in self.LastSpecificSol.Rules:
                self.LastSpecificSol.Rules.append(rule)

    def writeLastSpecificSol(self):
        self.LastSpecificSol.writeItselfTxt(hypotheses_dir + '/specific' + '/%s_%d_%d.lp' % (GRN_name, num_of_Is, num_of_Bs))


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
                print('DISmatched!')
                def printDISmatched():
                    print('The two examples is not desired!')
                    print('************examples1*************')
                    self.printE()
                    print('************examples2**********')
                    E2.printE()
                if silent_level < 1:
                    printDISmatched()
        return judgement

    def examples_comp_bk(self, examples2, atoms, flag_comp): #examples2 is a list of trans
        examples1 = self.Exps
        agree = False
        for i in range(len(examples1)):
            for a in atoms:
                v1 = examples1[i].output[a]
                v2 = examples2[i].output[a]
                if flag_comp == 0:
                    if v1 != v2:
                        return False
                else:
                    if (v1 > v2) - (v1 < v2) == -flag_comp:
                        return False
                    elif (v1 > v2) - (v1 < v2) == flag_comp:
                        agree = True
        if flag_comp != 0:
            return agree
        return True


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

    def writeBackgroundASP(self, list_B):
        background_asp_path = background_asp_dir + '/%s_%d_%d.lp' %(GRN_name,num_of_Is,num_of_Bs)
        str_all_rules = ''
        i = 0
        while i < len(list_B):
            rule = list_B[i]
            i = i + 1
            pos_body = ''.join(list(map(lambda x: 'bg_body_pos(%d,%s).\n' % (i, x), rule.pos)))
            neg_body = ''.join(list(map(lambda x: 'bg_body_neg(%d,%s).\n' % (i, x), rule.neg)))
            head = 'bg_head(%d,%s).\n' % (i, rule.head)
            facts_temp = head + pos_body + neg_body + 'bg_necessity(%d,%d).\n' % (i, round(rule.necessity * self.magnification))
            str_all_rules = str_all_rules + facts_temp
        if silent_level < 3:
            with open(background_asp_path, 'w') as file_e:
                file_e.write(str_all_rules)
        return str_all_rules

    def writeSeparatedBackgroundASP(self, list_B):
        # put these asp strings into a dictionary
        dic_atom_asp = dict()
        for k in range(1, self.num_of_atoms + 1):  # initialize dic_atom_asp
            dic_atom_asp[k] = ''
        i = 0
        while i < len(list_B):
            rule = list_B[i]
            i = i + 1
            head = int(rule.head)
            pos_body = ''.join(list(map(lambda x: 'bg_body_pos(%d,%s).\n' % (i, x), rule.pos)))
            neg_body = ''.join(list(map(lambda x: 'bg_body_neg(%d,%s).\n' % (i, x), rule.neg)))
            necessity = 'bg_necessity(%d,%d).\n' % (i, round(rule.necessity * self.magnification))
            str_rule = pos_body + neg_body + necessity
            dic_atom_asp[head] = dic_atom_asp[head] + str_rule
        if silent_level < 3 or search_type == 'compare':
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

    def writeExamplesASP(self,list_E): #{(a, 0.3); (s, 0.4)} -> {(s, 0.1); (b, 0.4)}
        examples_asp_path = examples_asp_dir + '/%s_%d.lp' %(GRN_name,num_of_Is)
        str_all_exampers = ''
        idx = 1
        for i in range(0, len(list_E)):
            facts_temp = ''
            trans = list_E[i]
            for k_i, v_i in trans.input.items():
                facts_temp = facts_temp + 'input(%d .. %d,%s,%d).\n' % (idx, idx + self.num_of_atoms - 1, k_i, round(v_i * self.magnification))
            for j in range(1, self.num_of_atoms + 1):
                outN = trans.output[str(j)]
                facts_temp = facts_temp + 'output(%d,%d,%d).\n' % (idx, j, round(outN * self.magnification))
                idx = idx + 1
            str_all_exampers = str_all_exampers + facts_temp
        if silent_level < 3:
            with open(examples_asp_path, 'w') as file_e:
                file_e.write(str_all_exampers)
        return str_all_exampers

    def writeSeparatedExamplesASP(self,list_E): #[{(a, 0.3); (s, 0.4)} -> {(s, 0.1); (b, 0.4)},  ]
        # put these asp strings into a dictionary
        dic_atom_asp = dict()
        for i in range(1, self.num_of_atoms + 1): #initialize dic_atom_asp
            dic_atom_asp[i] = ''
        facts_input = '' #the public input facts
        for i in range(1, len(list_E) + 1):
            trans = list_E[i - 1]
            for k_i, v_i in trans.input.items():
                facts_input = facts_input + 'input(%d,%s,%d).\n' % (i, k_i, round(v_i * self.magnification))
            for j in range(1, self.num_of_atoms + 1):
                outN = trans.output[str(j)]
                facts_output = 'output(%d,%d).\n' % (i, round(outN * self.magnification))
                dic_atom_asp[j] = dic_atom_asp[j] + facts_output
        for i in range(1, self.num_of_atoms + 1): #append facts_input into each elements of dic_atom_asp
            dic_atom_asp[i] = dic_atom_asp[i] + facts_input
        #write this dictionary into separated files
        for i in range(1, self.num_of_atoms + 1):
            examples_asp_path = examples_asp_dir + '/%s_%d+%d.lp' % (GRN_name, num_of_Is, i)
            if silent_level < 3 or search_type == 'compare':
                with open(examples_asp_path, 'w') as file_e:
                    file_e.write(dic_atom_asp[i])
        return dic_atom_asp

    def writeAuxiliaryASP(self):
        pauxiliary_path = auxiliary_dir + '/%s.lp' % GRN_name
        facts = 'top(%d). \n' %self.magnification
        facts = facts + 'universe(1 .. %s). \n' %self.num_of_atoms
        if silent_level < 3 or search_type == 'compare':
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


class ASP:
    Examples = ''
    Background = ''
    Auxiliary = ''
    Encoding = ''
    StableModel = ''
    LastOptimalSol = list() #a list of rules

    def __init__(self, str_Examples, str_Background, str_Auxiliary):
        self.Examples = str_Examples
        self.Background = str_Background
        self.Auxiliary = str_Auxiliary

    def getEncoding(self,file_path):
        with open(file_path, encoding="utf-8") as f:
            self.Encoding = f.read()
        f.closed

    def getStableModel(self):
        facts = self.Examples +  self.Background + self.Auxiliary
        ctl = Control(["0", "--opt-mode=opt", "--parallel-mode=2,compete"])
        ctl.add("instance", [], facts)
        ctl.ground([("instance", [])])
        ctl.add("base", [], self.Encoding)
        ctl.ground([("base", [])])
        with ctl.solve(yield_=True) as handler:
            for model in handler:
                self.StableModel = str(model.symbols(shown=True))[1:-1]

    def parseTotalModel(self,magnification):
        self.LastOptimalSol = list()
        model = self.StableModel.split(', ')
        for rule in list(filter(lambda x: x.startswith('last_rule('), model)):
            parameters = rule[10:-1].split(',')
            index_start = 10 + len(parameters[0])
            rule_head = str(parameters[1])
            rule_pos = frozenset(list(map(lambda x: x[index_start:-1], list(filter(lambda x: x.startswith('rule_pos(%s,' % parameters[0]), model)))))
            rule_neg = frozenset(list(map(lambda x: x[index_start:-1], list(filter(lambda x: x.startswith('rule_neg(%s,' % parameters[0]), model)))))
            rule_N = round(int(parameters[2]) /  magnification, 2)
            crule = cRule()
            crule.refreshRuleFromAttributes(rule_head, rule_pos, rule_neg, rule_N)
            self.LastOptimalSol.append(crule)

    def parseSeparatedModel(self,head_atom,magnification):
        self.LastOptimalSol = list() #a list of Rules whose head is head_atom
        model = self.StableModel.split(', ')
        for rule in list(filter(lambda x: x.startswith('last_rule('), model)):
            parameters = rule[10:-1].split(',')
            index_start = 10 + len(parameters[0])
            rule_head = head_atom
            rule_pos = frozenset(list(map(lambda x: x[index_start:-1], list(filter(lambda x: x.startswith('rule_pos(%s,' % parameters[0]), model)))))
            rule_neg = frozenset(list(map(lambda x: x[index_start:-1], list(filter(lambda x: x.startswith('rule_neg(%s,' % parameters[0]), model)))))
            rule_N = round(int(parameters[1]) /  magnification, 2)
            crule = cRule()
            crule.refreshRuleFromAttributes(rule_head, rule_pos, rule_neg, rule_N)
            self.LastOptimalSol.append(crule)


    def computeLastOptimal(self,head_atom,magnification,assigned_dir):
        if head_atom == '':
            file_path = assigned_dir + '/MHSC.lp'
        else:
            file_path = assigned_dir + '/MHSC_head.lp'
        self.getEncoding(file_path)
        self.getStableModel()
        #print(self.StableModel)
        if head_atom == '':
            self.parseTotalModel(magnification)
        else:
            self.parseSeparatedModel(head_atom,magnification)




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

    #read GRN, generate NLP, write NLP
    def refreshGRNFromFile(self, GRN_name):
        self.NLP = cNLP()
        file_name = '/%s.lp'%GRN_name
        nlp_file = source_GRN_dir + file_name
        fp = open(nlp_file)
        line = fp.readline()
        while (line != ''):
            line.strip()  # remove the left and right blanks
            if line == '':  # skip the empty line
                line = fp.readline()
                continue
            if line.startswith('#'):  # read the real node name and its id
                self.str_atom_map = self.str_atom_map + line
                line = fp.readline()
                continue
            rule = self.getRuleFromGRN(line)
            self.NLP.Rules.append(rule)
            line = fp.readline()
        fp.close()
        self.NLP.str_atom_map = self.str_atom_map
        self.refreshNecessities()
        self.NLP.writeItselfTxt(primary_program_dir + file_name)



class CMDsClingo:
    cmd_list = list()
    atom_ids = list()
    magnification = 1
    resultsDict = dict()
    modelsDict = dict()
    solsList = list() # a list of cRules
    timeoutFlag = False
    optimumFlag = True


    def __init__(self, lst, atom_ids, magnification):
        self.cmd_list = lst
        self.atom_ids = atom_ids
        self.magnification = magnification
        timeoutFlag = False
        optimumFlag = True

    @fn_memory_time
    def exeCMDsParallel(self):
        ResultList = dict()
        for i in range(0,len(self.atom_ids)):
            cmd = self.cmd_list[i]
            child = Popen(cmd, stdin=PIPE, stdout=PIPE, stderr=PIPE, shell=True)
            ResultList[str(i+1)] = child
        for k,v in ResultList.items():
            self.resultsDict[k] = v.communicate()[0].decode("gbk")

    def exeCMDsParallel_bk(self):
        ResultList = list()
        for cmd in self.cmd_list:
            child = Popen(cmd, stdin=PIPE, stdout=PIPE, stderr=PIPE, shell=True)
            ResultList.append(child)
        self.resultsList = list(map(lambda x: x.communicate()[0].decode("gbk"), ResultList))

    def preParseModels(self, res):
        pattern = '^Models\s*:\s*(\d+)(\+?)\s'
        searchObj = re.search(pattern, res, re.M)
        num = searchObj.group(1)
        add = searchObj.group(2)
        return num, add

    def preParseAnswer(self, res):
        pattern = 'Answer:\s*\d*\s*([A-Za-z_,()0-9 ]*)\s*Optimization:'
        searchObj = re.search(pattern, res, re.M)
        AS = searchObj.group(1)
        return AS

    def preParseOptimum(self, res):
        pattern = '\s+Optimum\s*:\s*(\w+)\s'
        searchObj = re.search(pattern, res, re.M)
        opt = searchObj.group(1)
        return opt

    def parseResult(self, resultK):
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
                resultDict['AS'] = self.preParseAnswer(resultK)
                ptInf = self.preParseOptimum(resultK)
                if ptInf != 'yes':
                    #print('******Not OPTIMUM********')
                    resultDict['optmum'] = False
            if modelInf[1] == '+':
                #print('******Timeout********')
                resultDict['timeout'] = True
        return resultDict

    def getModels(self):
        for k,v in self.resultsDict.items():
            self.modelsDict[k] = self.parseResult(v)
        self.timeoutFlag = list(map(lambda x: x['timeout'], self.modelsDict.values())).count(True) > 0
        self.optimumFlag = not list(map(lambda x: x['optmum'], self.modelsDict.values())).count(False) > 0

    def getModels_bk(self):
        for r in self.resultsList:
            self.modelsList.append(self.parseResult(r))
        self.timeoutFlag = list(map(lambda x: x['timeout'], self.modelsList)).count(True) > 0
        self.optimumFlag = not list(map(lambda x: x['optmum'], self.modelsList)).count(False) > 0

    def parseSeparatedModel(self, stableModel, head_atom, magnification):
        LastOptimalSol = list()  # a list of Rules whose head is head_atom
        model = stableModel.split(' ')
        for rule in list(filter(lambda x: x.startswith('last_rule('), model)):
            parameters = rule[10:-1].split(',')
            index_start = 10 + len(parameters[0])
            rule_head = head_atom
            rule_pos = frozenset(list(map(lambda x: x[index_start:-1], list(filter(lambda x: x.startswith('rule_pos(%s,' % parameters[0]), model)))))
            rule_neg = frozenset(list(map(lambda x: x[index_start:-1], list(filter(lambda x: x.startswith('rule_neg(%s,' % parameters[0]), model)))))
            rule_N = round(int(parameters[1]) / magnification, 2)
            crule = cRule()
            crule.refreshRuleFromAttributes(rule_head, rule_pos, rule_neg, rule_N)
            LastOptimalSol.append(crule)
        return LastOptimalSol

    def getRules(self):  # use a list of cRules to store the solutions
        for id in self.atom_ids:
            self.solsList.extend(self.parseSeparatedModel(self.modelsDict[id]['AS'], id, self.magnification))

    def getRules_bk(self): #use a list of cRules to store the solutions
        for i in range(0,len(self.atom_ids)):
            self.solsList.extend(self.parseSeparatedModel(self.modelsDict[i]['AS'],self.atom_ids[i],self.magnification))




def readWriteNLPs():
    global GRN_name_list
    for name in GRN_name_list:
        GRN = cGRN()
        GRN.refreshGRNFromFile(name)


def gen_E_B():
    primary_program_path = primary_program_dir + '/%s.lp'%GRN_name
    NLP = cNLP()
    NLP.refreshProgramFromFile(primary_program_path)
    NLP.printP()
    NLP.genRandomIs(num_of_Is)
    print(NLP.RandomIs)
    NLP.genExamples()
    print(NLP.RandomExps)
    NLP.writeExamplesTxt(NLP.RandomExps.Exps)
    NLP.writeExamplesASP(NLP.RandomExps.Exps)
    NLP.genBackground(Bs_ratio)
    NLP.print_B()
    NLP.writeBackgroundTxt(NLP.RomdomBackground)
    NLP.writeBackgroundASP(NLP.RomdomBackground)
    NLP.writeAuxiliaryASP()

def classify():
    hypotheses_path = hypotheses_dir + '/optimal/%s_%d_%d.lp' % (GRN_name, num_of_Is, num_of_Bs)
    backgrounds_txt_path = backgrounds_txt_dir + '/%s_%d_%d.lp' % (GRN_name, num_of_Is, num_of_Bs) #if putting clingo solver in this program too, NLP.RomdomBackground can be directly used.
    examples_txt_path = examples_txt_dir + '/%s_%d.lp' %(GRN_name,num_of_Is)
    NLP = cNLP()
    NLP.refreshProgramFromFile(hypotheses_path)
    #NLP.printP()
    NLP2 = cNLP()
    NLP2.refreshProgramFromFile(backgrounds_txt_path)
    #NLP2.printP()
    H = NLP.program_extend(NLP2.Rules)
    H.printP()
    E1 = Examples()
    E1.refreshExampleFromFile(examples_txt_path)
    #E2 = H.calTransitions(E1.Is)
    E2 = Examples()
    E2.Exps = H.calTransitions(E1.Is)
    E1.examples_comp(E2, H.atom_ids, 0)



def subprocessTask(k_str,E_str,B_str,A_str,M_int, root_dir):
    aspProblem = ASP(E_str,B_str,A_str)
    aspProblem.computeLastOptimal(k_str, M_int, root_dir)
    #time.sleep(20)
    return aspProblem.LastOptimalSol


def record_singleton(type_str,time,ram,len,match):
    global record_dict
    record_dict['Time_%s' % type_str] = time
    record_dict['RAM_%s' % type_str] = ram
    record_dict['P_%s' % type_str] = len
    record_dict['match_%s' % type_str] = match
    if float(time) > 10800 :  # timeLimit_value 30 * 60 = 1800
        timeout = True
    else:
        timeout = False
    record_dict['timeout_%s' % type_str] = str(timeout)

#integrate the above functions together
def gen_E_B_H_check_singleton(): #search_type in {specific, optimal_batch, optimal_sep, optimal_parallel}
    primary_program_path = primary_program_dir + '/%s.lp'%GRN_name
    NLP = cNLP()
    NLP.refreshProgramFromFile(primary_program_path)
    NLP.genRandomIs(num_of_Is)
    NLP.genExamples()
    NLP.genBackground(Bs_ratio)

    if search_type == 'specific':
        @fn_memory_time
        def computeSpecificSol():
            NLP.RandomExps.computeNaiveSpecificSol()
            NLP.RandomExps.computeLastSpecificSol(NLP.RomdomBackground)
        decoratorList = computeSpecificSol()
        LastSol = NLP.RandomExps.LastSpecificSol
        if silent_level < 4:
            NLP.RandomExps.writeLastSpecificSol()
    else:
        if search_type == 'optimal_batch':
            ExamplesASP = NLP.writeExamplesASP(NLP.RandomExps.Exps)
            BackgroundASP = NLP.writeBackgroundASP(NLP.RomdomBackground)
        else:
            ExamplesASPdic = NLP.writeSeparatedExamplesASP(NLP.RandomExps.Exps)
            BackgroundASP = NLP.writeSeparatedBackgroundASP(NLP.RomdomBackground)
        AuxiliaryASP = NLP.writeAuxiliaryASP()
        LastSol = cNLP()
        LastSol.refreshProgramFromAttributes(NLP.atom_ids,NLP.str_atom_map)

        if search_type == 'optimal_batch':
            @fn_memory_time
            def computeOptimalBatched():
                aspProblem = ASP(ExamplesASP, BackgroundASP, AuxiliaryASP)
                aspProblem.computeLastOptimal('', NLP.magnification, root_dir)
                return aspProblem.LastOptimalSol
            decoratorList = computeOptimalBatched()
            LastSol.Rules = decoratorList[2]
        elif search_type == 'optimal_sep':
            @fn_memory_time
            def computeOptimalIncrementally():
                tmp_list = list()
                for k in range(1, NLP.num_of_atoms + 1):
                    aspProblem = ASP(ExamplesASPdic[k], BackgroundASP[k], AuxiliaryASP)
                    aspProblem.computeLastOptimal(str(k),NLP.magnification, root_dir)
                    tmp_list.extend(aspProblem.LastOptimalSol)
                return tmp_list
            decoratorList = computeOptimalIncrementally()
            LastSol.Rules = decoratorList[2]
        elif search_type == 'optimal_parallel':
            @fn_memory_time
            def computeOptimalParallelIncrementally():
                print('Here is a test!')
                applyResultList = list()
                pool = Pool(NLP.num_of_atoms)
                for k in range(1, NLP.num_of_atoms + 1):
                    p_result = pool.apply_async(func=subprocessTask, args=(str(k), ExamplesASPdic[k], BackgroundASP[k], AuxiliaryASP, NLP.magnification, root_dir))
                    applyResultList.append(p_result)
                print('%d processes is running!'%len(pool._cache))
                pool.close()
                pool.join()
                programsList = list(map(lambda x: x.get(), applyResultList))
                lastProgram = list(itertools.chain(*programsList))
                return lastProgram


            decoratorList = computeOptimalParallelIncrementally()
            LastSol.Rules = decoratorList[2]
        if silent_level < 4:
            LastSol.writeItselfTxt(hypotheses_dir + '/optimal' + '/%s_%d_%d.lp' % (GRN_name, num_of_Is, num_of_Bs))


    H = LastSol.program_extend(NLP.RomdomBackground)
    E2 = Examples()
    E2.Exps = H.calTransitions(NLP.RandomIs)
    match_flag = NLP.RandomExps.examples_comp(E2, NLP.atom_ids, 0)

    if search_type == 'specific':
        record_singleton('specific', decoratorList[0], decoratorList[1], len(LastSol.Rules), str(match_flag))
    else:
        record_singleton(search_type[search_type.index('_') + 1:], decoratorList[0], decoratorList[1], len(LastSol.Rules), str(match_flag))


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


#run the 4 functions sequently toward to a same data source
def gen_E_B_H_check_overall(): #search_type == overall
    primary_program_path = primary_program_dir + '/%s.lp'%GRN_name
    NLP = cNLP()
    NLP.refreshProgramFromFile(primary_program_path)
    NLP.genRandomIs(num_of_Is)
    NLP.genExamples()
    NLP.genBackground(Bs_ratio)

    LastSol = None
    ExamplesASP = None
    BackgroundASP = None
    ExamplesASPdic = None
    BackgroundASPdic = None
    AuxiliaryASP = None
    global record_dict

    def getSpecific():
        @fn_memory_time
        def computeSpecificSol():
            NLP.RandomExps.computeNaiveSpecificSol()
            NLP.RandomExps.computeLastSpecificSol(NLP.RomdomBackground)

        decoratorList = computeSpecificSol()
        nonlocal LastSol
        LastSol = NLP.RandomExps.LastSpecificSol
        record_singleton('specific', decoratorList[0], decoratorList[1], len(LastSol.Rules), '')
        if silent_level < 4:
            NLP.RandomExps.writeLastSpecificSol()

    def prepareOptimal():
        nonlocal ExamplesASP,BackgroundASP,ExamplesASPdic,BackgroundASPdic,AuxiliaryASP

        if search_type == 'optimal_batch':
            ExamplesASP = NLP.writeExamplesASP(NLP.RandomExps.Exps)
            BackgroundASP = NLP.writeBackgroundASP(NLP.RomdomBackground)
        else:
            ExamplesASPdic = NLP.writeSeparatedExamplesASP(NLP.RandomExps.Exps)
            BackgroundASPdic = NLP.writeSeparatedBackgroundASP(NLP.RomdomBackground)
        AuxiliaryASP = NLP.writeAuxiliaryASP()
        LastSol = cNLP()
        LastSol.refreshProgramFromAttributes(NLP.atom_ids, NLP.str_atom_map)

    def getOptimalBatch():
        #nonlocal ExamplesASP,BackgroundASP,AuxiliaryASP
        @fn_memory_time
        def computeOptimalBatched():
            aspProblem = ASP(ExamplesASP, BackgroundASP, AuxiliaryASP)
            aspProblem.computeLastOptimal('', NLP.magnification, root_dir)
            return aspProblem.LastOptimalSol

        decoratorList = computeOptimalBatched()
        LastSol.Rules = decoratorList[2]
        record_singleton('batch', decoratorList[0], decoratorList[1], len(LastSol.Rules), '')

    def getOptimalSep():
        @fn_memory_time
        def computeOptimalIncrementally():
            tmp_list = list()
            for k in range(1, NLP.num_of_atoms + 1):
                aspProblem = ASP(ExamplesASPdic[k], BackgroundASPdic[k], AuxiliaryASP)
                aspProblem.computeLastOptimal(str(k), NLP.magnification, root_dir)
                tmp_list.extend(aspProblem.LastOptimalSol)
            return tmp_list

        decoratorList = computeOptimalIncrementally()
        LastSol.Rules = decoratorList[2]
        record_singleton('sep', decoratorList[0], decoratorList[1], len(LastSol.Rules), '')

    def getOptimalParallel():

        @fn_memory_time
        def computeOptimalParallelIncrementally():
            applyResultList = list()
            pool = Pool(NLP.num_of_atoms)
            for k in range(1, NLP.num_of_atoms + 1):
                p_result = pool.apply_async(func=subprocessTask, args=(str(k), ExamplesASPdic[k], BackgroundASPdic[k], AuxiliaryASP, NLP.magnification, root_dir))
                applyResultList.append(p_result)
            pool.close()
            pool.join()
            programsList = list(map(lambda x: x.get(), applyResultList))
            lastProgram = list(itertools.chain(*programsList))
            return lastProgram

        decoratorList = computeOptimalParallelIncrementally()
        LastSol.Rules = decoratorList[2]
        record_singleton('parallel', decoratorList[0], decoratorList[1], len(LastSol.Rules), '')


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
        if search_type != 'specific' and silent_level < 4:
            LastSol.writeItselfTxt(hypotheses_dir + '/optimal' + '/%s_%d_%d.lp' % (GRN_name, num_of_Is, num_of_Bs))

    def resultCheck():
        H = LastSol.program_extend(NLP.RomdomBackground)
        E2 = Examples()
        E2.Exps = H.calTransitions(NLP.RandomIs)
        match_flag = NLP.RandomExps.examples_comp(E2, NLP.atom_ids, 0)
        if search_type == 'specific':
            record_dict['match_specific'] = str(match_flag)
        else:
            record_dict['match_%s'%search_type[search_type.index('_') + 1:]] = str(match_flag)


    def printWriteCheck():
        if silent_level < 1:
            intermediatePrint()
        if silent_level < 2:
            intermediateWrite()
        resultCheck()


    search_type = 'specific'
    getSpecific();    printWriteCheck()
    search_type = 'optimal_parallel'
    prepareOptimal();    getOptimalParallel();    printWriteCheck()
    search_type = 'optimal_sep'
    prepareOptimal();    getOptimalSep();    printWriteCheck()
    search_type = 'optimal_batch'
    prepareOptimal();    getOptimalBatch();    printWriteCheck()

    '''
    #the original respective function logic
    if search_type == 'specific':
        getSpecific()
    else:
        prepareOptimal()
        if search_type == 'optimal_batch':
            getOptimalBatch()
        elif search_type == 'optimal_sep':
            getOptimalSep()
        elif search_type == 'optimal_parallel':
            getOptimalParallel()
    printWriteCheck()
    '''


#run only specific and optimal_parallel toward to a same data source via using cmd implementation
def gen_E_B_H_check_compare(): #search_type == compare
    primary_program_path = primary_program_dir + '/%s.lp'%GRN_name
    NLP = cNLP()
    NLP.refreshProgramFromFile(primary_program_path)
    NLP.genRandomIs(num_of_Is)
    NLP.genExamples()
    NLP.genBackground(Bs_ratio)

    LastSol = None
    ExamplesASPdic = None
    BackgroundASPdic = None
    AuxiliaryASP = None
    global record_dict

    def getSpecific():
        @fn_memory_time
        def computeSpecificSol():
            NLP.RandomExps.computeNaiveSpecificSol()
            NLP.RandomExps.computeLastSpecificSol(NLP.RomdomBackground)

        decoratorList = computeSpecificSol()
        nonlocal LastSol
        LastSol = NLP.RandomExps.LastSpecificSol
        record_singleton('specific', decoratorList[0], decoratorList[1], len(LastSol.Rules), '')
        NLP.RandomExps.writeLastSpecificSol()

    def prepareOptimal():
        nonlocal  ExamplesASPdic,BackgroundASPdic,AuxiliaryASP

        ExamplesASPdic = NLP.writeSeparatedExamplesASP(NLP.RandomExps.Exps)
        BackgroundASPdic = NLP.writeSeparatedBackgroundASP(NLP.RomdomBackground)
        AuxiliaryASP = NLP.writeAuxiliaryASP()
        LastSol = cNLP()
        LastSol.refreshProgramFromAttributes(NLP.atom_ids, NLP.str_atom_map)

    def getOptimalParallel():
        def getCMDs():
            cmd_list = list()
            exe_path = 'clingo'#'c:/Users/86188/.vscode/extensions/ffrankreiter.answer-set-programming-language-support-0.4.2/clingo_win.exe'
            encoding_path = root_dir + '/MHSC_head.lp'
            parameter_path = '--outf=0 --verbose=1 --parallel-mode=2,compete --time-limit=10800 --models 0 --opt-mode=opt --quiet=1' # timeLimit_value
            for k in range(1,len(NLP.atom_ids) + 1):
                background_asp_path = background_asp_dir + '/%s_%d_%d+%d.lp' % (GRN_name, num_of_Is, num_of_Bs, k)
                examples_asp_path = examples_asp_dir + '/%s_%d+%d.lp' % (GRN_name, num_of_Is, k)
                pauxiliary_path = auxiliary_dir + '/%s.lp' % GRN_name
                cmd_str = '%s %s %s %s %s %s ' % (exe_path, encoding_path, background_asp_path, examples_asp_path, pauxiliary_path, parameter_path)
                cmd_list.append(cmd_str)
            return cmd_list

        commandList = getCMDs()
        if silent_level < 3:
            print(commandList)
        myCMDS = CMDsClingo(commandList, NLP.atom_ids, NLP.magnification)
        decoratorList = myCMDS.exeCMDsParallel()  # the memory consumption will contain it's children's part because of  include_children=True along memory_profiler.memory_usage()
        myCMDS.getModels()
        myCMDS.getRules()
        if silent_level < 4:
            print(myCMDS.modelsDict)
            #print(myCMDS.solsList)
        LastSol.Rules = myCMDS.solsList
        record_singleton('parallel', decoratorList[0], decoratorList[1], len(LastSol.Rules), '')
        record_dict['timeout_parallel'] = str(myCMDS.timeoutFlag)


    def getOptimalParallel_bk():
        @fn_memory_time
        def computeOptimalParallelIncrementally():
            applyResultList = list()
            pool = Pool(NLP.num_of_atoms)
            for k in range(1, NLP.num_of_atoms + 1):
                p_result = pool.apply_async(func=subprocessTask, args=(str(k), ExamplesASPdic[k], BackgroundASPdic[k], AuxiliaryASP, NLP.magnification, root_dir))
                applyResultList.append(p_result)
            pool.close()
            pool.join()
            programsList = list(map(lambda x: x.get(), applyResultList))
            lastProgram = list(itertools.chain(*programsList))
            return lastProgram
        decoratorList = computeOptimalParallelIncrementally()
        LastSol.Rules = decoratorList[2]
        record_singleton('parallel', decoratorList[0], decoratorList[1], len(LastSol.Rules), '')



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
        LastSol.writeItselfTxt(hypotheses_dir + '/optimal' + '/%s_%d_%d.lp' % (GRN_name, num_of_Is, num_of_Bs))

    tmp_search_type = ''
    def resultCheck():
        H = LastSol.program_extend(NLP.RomdomBackground)
        E2 = Examples()
        E2.Exps = H.calTransitions(NLP.RandomIs)
        match_flag = NLP.RandomExps.examples_comp(E2, NLP.atom_ids, 0)
        if tmp_search_type == 'specific':
            record_dict['match_specific'] = str(match_flag)
        else:
            record_dict['match_parallel'] = str(match_flag)

    def printWriteCheck():
        if silent_level < 1:
            intermediatePrint()
        intermediateWrite()
        resultCheck()

    tmp_search_type = 'specific'
    getSpecific();    printWriteCheck()
    tmp_search_type = 'optimal_parallel'
    prepareOptimal();    getOptimalParallel();    printWriteCheck()


def printExecutionSpecification(pName):
    print("Instruction: %s <GRN_name> <num_of_Is> <Bs_ratio> <search_type>  <silent_level> <file_name>" %pName)
    print('The value range of these parameters are as follows:')
    print('GRN_name in {mammalian, fission, budding, arabidopsis, thelper, tcrNet}')
    print('num_of_Is in {40, 80, 120, 160, 200, 240, 280, 320, 100, 10000}')
    print('Bs_ratio in {0, 0.1, 0.4, 0.8}')
    print('search_type in {specific, optimal_batch, optimal_sep, optimal_parallel, overall, compare}')
    print('silent_level in {0, 1, 2, 3, 4, 5, 6}')
    print('parameter file_name assigns the place where log the experimental record')
    '''
    6 + write necessary results
    5 + print time;memory
    4 + print check result
    3 + write LastSol when search_type is not compare; print modelsDict for debug otherwise.
    2 + write ASP when search_type is not compare; print commandList for debug otherwise.
    1 + intermediateWrite when search_type is not compare
    0 + intermediatePrint + printDISmatched
    '''

def externalCall():
    global GRN_name,num_of_Is,Bs_ratio,search_type,silent_level,file_name
    wrongParameter = False
    if len(sys.argv) != 7:
        printExecutionSpecification(sys.argv[0])
        sys.exit(0)

    if sys.argv[1] not in {'mammalian', 'fission', 'budding', 'arabidopsis', 'thelper', 'tcrNet'}:
        print('Instruction:GRN_name in {mammalian, fission, budding, arabidopsis, thelper, tcrNet}')
        wrongParameter = True
    if sys.argv[2] not in {'2','40', '80', '120', '160', '200', '240', '280', '320', '360', '400', '440', '480', '520', '100', '10000'}:
        print('Instruction:num_of_Is in {2, 40, 80, 120, 160, 200, 240, 280, 320, 360, 400, 440, 480, 520, 100, 10000}')
        wrongParameter = True
    if sys.argv[3] not in {'0', '0.1', '0.4', '0.8'}:
        print('Instruction:Bs_ratio in {0, 0.1, 0.4, 0.8}')
        wrongParameter = True
    if sys.argv[4] not in {'specific', 'optimal_batch', 'optimal_sep', 'optimal_parallel', 'overall', 'compare'}:
        print('Instruction:search_type in {specific, optimal_batch, optimal_sep, optimal_parallel, overall, compare}')
        wrongParameter = True
    if sys.argv[5] not in {'0', '1', '2', '3', '4', '5', '6'}:
        print('Instruction:silent_level in {0, 1, 2, 3, 4, 5, 6}')
        wrongParameter = True

    if wrongParameter:
        print("Instruction: %s <GRN_name> <num_of_Is> <Bs_ratio> <search_type>  <silent_level>" % sys.argv[0])
        sys.exit(0)

    GRN_name = sys.argv[1]
    num_of_Is = int(sys.argv[2])
    Bs_ratio = round(float(sys.argv[3]), 1)
    search_type = sys.argv[4]
    silent_level = int(sys.argv[5])
    file_name = sys.argv[6]

if __name__ == "__main__":

    GRN_name = 'mammalian'
    num_of_Is = 80
    Bs_ratio = 0 #the actual num_of_Bs is determined by Bs_ratio * (NLP.Rules)
    search_type = 'optimal_parallel' #search_type in {specific, optimal_batch, optimal_sep, optimal_parallel, overall, compare}
    silent_level = 2
    file_name = 'log.log'

    #comment this function if doing internal debug or invocation
    externalCall()

    GRN_name_list = ['mammalian', 'fission', 'budding', 'arabidopsis', 'thelper', 'tcrNet']  # GRN_name  is an element from GRN_name_list
    Is_amount_list = [10, 20, 40, 80, 160, 640]  # the num_of_Is is an element from Is_amount_list
    num_of_Bs = 0
    times = 10  # run 10 times and get the average



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

    record_dict = dict()
    record_dict['GRN'] = GRN_name
    record_dict['E'] = num_of_Is
    record_dict['Bs_ratio'] = Bs_ratio
    record_dict['P_specific'] = 0
    record_dict['Time_specific'] = 0
    record_dict['RAM_specific'] = 0
    record_dict['match_specific'] = ''
    record_dict['timeout_specific'] = ''
    record_dict['P_parallel'] = 0
    record_dict['Time_parallel'] = 0
    record_dict['RAM_parallel'] = 0
    record_dict['match_parallel'] = ''
    record_dict['timeout_parallel'] = ''
    record_dict['P_sep'] = 0
    record_dict['Time_sep'] = 0
    record_dict['RAM_sep'] = 0
    record_dict['match_sep'] = ''
    record_dict['timeout_sep'] = ''
    record_dict['P_batch'] = 0
    record_dict['Time_batch'] = 0
    record_dict['RAM_batch'] = 0
    record_dict['match_batch'] = ''
    record_dict['timeout_batch'] = ''


    if search_type == 'overall':
        gen_E_B_H_check_overall()
    elif search_type == 'compare':
        gen_E_B_H_check_compare()
    else:
        gen_E_B_H_check_singleton()

    #print(record_dict)
    with open(root_dir + '/log/' + file_name, 'a+', encoding='utf-8') as file:
        alldata = list(record_dict.values())
        file.write(json.dumps(alldata).replace('\"', '') + '\n')
    file.close()

    '''
    if search_type == 'specific':
        gen_E_B_SpecificH_check()
    elif search_type == 'optimal_batch':
        gen_batch_E_B_OptimalH_check()
    elif search_type == 'optimal_sep':
        gen_sep_E_B_OptimalH_check()
    elif search_type == 'optimal_parallel':
        gen_parallel_E_B_OptimalH_check()
    '''
    #gen_E_B()
    #test2()
    #readWriteNLPs()
    #classify()



