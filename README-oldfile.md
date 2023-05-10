Inductive learning programing of AI lab. Generally, we use gen-E-compare.py implimented in Python 3.6 to test the major features. But the key lies in MHSC_head.lp solved by clingo for a preferable solution.

introduction of directories:
Directory GRN	stores the original genetical regulation networks.
Directory primary_program	stores the modified GRNs. Namely, the poss-NLPs after adding necessity to the rules in GRN.
Directory auxiliary stores the auxiliary files.
Directory  examples_txt stores the examples in the form of tuple.
Directory  examples_asp stores the examples in the form of predicates.
Directory background_txt stores the background in the form of rules.
Directory background_asp stores the background in the form of predicates.
Directory hypotheses stores the result. Namely, the learned poss-NLP.
Directory log is the suggested directory to place log.
Directory gen_bash is used to store gen.sh.


introduction of files:
MHSC_V0.lp is a ASP solving the task with three preferences.
MHSC.lp is a ASP solving the task with one preference only.
MHSC_head.lp is a ASP solving the task with one preference regardless of heads.
gen-E-compareV0.py invokes MHSC_V0.lp.
gen-E-compare.py invokes MHSC.lp or MHSC_head.lp.
gen.sh constructs a number of commands and executes them sequentially.
statistic.py does statistics for experiments.

One can perform this experiment via Linux command 
python3.6 -f -g -e -b -t -s -l
such as the suggested "python3.6 gen-E-compare.py budding $E 0.1 compare 5 fact.log"
Options:
-f the executing file. Namely, gen-E-compare.py|gen-E-compareV0.py
-g the name of a GRN. Namely,  mammalian|fission|budding|arabidopsis|thelper|tcrNet
-e the amount of the examples. Namely,  2|40|80|120|160|200|240|280|320|360|400|440|480|520|100|10000
-b the ratio of the background. Namely,  0|0.1|0.4|0.8
-t the type of the processing. Namely, specific|optimal_batch|optimal_sep|optimal_parallel|overall|compare
-s the silent level. The higher of this value, the less contents printed in the result. This parameter is mainly used for debugging.  Namely, 0|1|2|3|4|5|6
-l the name of the log.

The result is in the form 
[name, |E|, |ratioB|, NLR1, time1, memory1, matchFlag1, timeoutFlag1, NLR2, time2, memory2, matchFlag2, timeoutFlag2, NLR3, time3, memory3, matchFlag3, timeoutFlag3, NLR4, time4, memory4, matchFlag4, timeoutFlag4]
such as 
[fission, 40, 0.1, 118, 0.2994, 34.28125, True, False, 18, 2.0539, 442.99609375, True, False, 18, 1.3712, 74.91796875, True, False, 18, 1.9915, 160.4296875, True, False].
The first three element name, |E|, |ratioB| respectively stores -g, -e, -b in the parameters.
The left element concludes the results of four types of the processing corresponding to the parameter  specific|optimal_batch|optimal_sep|optimal_parallel.
The result matchFlag is True if the learned program explains the examples.
The result timeoutFlag is True if timeout occurs.
