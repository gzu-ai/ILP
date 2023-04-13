
%试试优化求最值的地方

%%%% get the candidate solution space%%%%%%%%%
{lastRuleHead(I, H, N)} :- output(I, H, N).
{rulePos(I, H, N, A)} :- output(I, H, N), input(I, A, NA), NA >= N.
{ruleNeg(I, H, N, A)} :- universe(A), output(I, H, N), not input(I, A, _). 


%%%% get the psuedo solution space%%%%%%%%%
%get the concomitant outputs
coOutput(I, A, 0) :- input(I, _, _), universe(A), not output(I, A, _).
%get the elements of psuedo solution spaces
limitation(I, H, N) :- coOutput(I, H, N).
limitation(I, H, N) :- output(I, H, N).
limitPos(I, H, N, A) :- limitation(I, H, N), input(I, A, NA), NA > N.
limitNeg(I, H, N, A) :- universe(A), limitation(I, H, N), not input(I, A, _). 

%%%% no rule occures in a  psuedo solution space%%%%%%%%%%%%%
limitStrictlySubsume(I, J, H, N1, N2) :- lastRuleHead(I, H, N1), limitation(J, H, N2), I != J, N1 > N2.
notLimitSubsetPos(I, J, H, N1, N2) :- limitStrictlySubsume(I, J, H, N1, N2), rulePos(I, H, N1, A), not limitPos(J, H, N2, A).
notLimitSubsetNeg(I, J, H, N1, N2) :- limitStrictlySubsume(I, J, H, N1, N2), ruleNeg(I, H, N1, A), not limitNeg(J, H, N2, A).
:- limitStrictlySubsume(I, J, H, N1, N2), not notLimitSubsetPos(I, J, H, N1, N2), not notLimitSubsetNeg(I, J, H, N1, N2).


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %given background covers a single transition
% %the negative body part of applicable requirements is unsatisfied. Namely, a intersection exists.
% notBodyNegSat(R,I) :- bgRuleNeg(R,A),input(I,A,_).
% %the positive body part of applicable requirements is unsatisfied. Namely, not a subset.
% notBodyPosSat(R,I) :- bgRulePos(R,A),output(I,_,_),not input(I,A,_).
% %R is applicable in I, since all the three applicable requirements satisfy.
% applicable(R, I, H, NR,NI) :- bgRuleHead(R,H, NR),output(I,H,NI),not notBodyNegSat(R,I),not notBodyPosSat(R,I).
% %R cover <I,{(H,NI)}>
% coverB(R, I, H, NI) :- applicable(R, I, H, NR, NI),#min{N2 : input(I,A,N2),bgRulePos(R,A);NR} = BETA , BETA == NI.

%existed rule from example covers a single transition which corresponds to no rule at last
needCoverTran(I, H, N) :- output(I, H, N), not lastRuleHead(I, H, N).
notBodyNegSatE(I, H, N, J) :- lastRuleNeg(I, H, N, A),needCoverTran(J,_,_),input(J,A,_).
notBodyPosSatE(I, H, N, J) :- lastRulePos(I, H, N, A),needCoverTran(J,_,_),not input(J,A,_).
applicableE(I, J, H, N, NI) :- lastRuleHead(I, H, N),needCoverTran(J,H,NI), not notBodyNegSatE(I, H, N, J),not notBodyPosSatE(I, H, N,J), N >= NI.
%coverE(I, J, H, N, NI) :- applicableE(I, J, H, N, NI),#min{N2 : input(J,A,N2),lastRulePos(I, H, N, A);N} = BETA , BETA == NI.


% coverE(I, J, H, N, NI) :- applicableE(I, J, H, N, NI),#min{N2 : input(J,A,N2),lastRulePos(I, H, N, A)} = BETA.

%coverE(I, J, H, N, NI) :-applicableE(I, J, H, N, NI), N == NI.


coverE(I, J, H, N, NI) :- applicableE(I, J, H, N, NI), #min{N2 : input(J,A,N2),lastRulePos(I, H, N, A); N3:lastRuleHead(I, H, N3)}  == NI.

%:- applicableE(I, J, H, N, NI), lastRulePos(I, H, N, A), input(J,A,N2), N2 < NI.

% coverE(I, J, H, N, NI) :- applicableE(I, J, H, N, NI), lastRulePos(I, H, N, A), input(J,A,N2), N2 == NI.
% coverE(I, J, H, N, NI) :- applicableE(I, J, H, N, NI), N == NI.

%a single transition must has corresponding rule at last if it is not covered by any existed rule.
% :- output(I, H, N), not lastRuleHead(I, H, N), not coverB(_, I, H, N), not coverE(_, I, H, _, N).
:- needCoverTran(I, H, N), not coverE(_, I, H, _, N).

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%leave the unsubsumed rules  as the last rules only%%%%%%%%%%
%the head must exist, whereas the body does not.
lastRulePos(I, H, N, A) :- lastRuleHead(I, H, N), rulePos(I, H, N, A).
lastRuleNeg(I, H, N, A) :- lastRuleHead(I, H, N), ruleNeg(I, H, N, A).

#show lastRuleHead/3.
#show lastRulePos/4.
#show lastRuleNeg/4.
%内容与ILTP一致
%#show subsumed/3.
%#show strictSubsumes/5.
%#show mutualSubsumes/5.
%#show winMutualSubsumes/3.
% #show limitStrictlySubsume/5.
% #show rule_pos/4.
% #show rule_neg/4.
%#show last_rule/3.
% #show input/3.
% #show output/3.
% #show coverE/5.
% #show applicableE/5.

% get the three measures %
%the last program has C rules.
last_cardinality(C) :- C = #count{(I, H, N): lastRuleHead(I, H, N)}.
%the amount of all atoms appearing in the program.
%atoms_amount(Pos_cnt + Neg_cnt + Head_cnt) :- Pos_cnt = #count{(I, H, N, A): lastRulePos(I, H, N, A) },Neg_cnt = #count{(I, H, N, A): lastRuleNeg(I, H, N, A) },last_cardinality(Head_cnt).
%the count of atoms in the program
%atoms_cnt(Atoms_cnt) :- Atoms_cnt = #count{A: lastRulePos(I, H, N, A); A: lastRuleNeg(I, H, N, A); H:lastRuleHead(J, H, N2)}.

% point out the optimal directions %
%Change the two into weak constraints
%less rules as possible. This measure is more significant than atoms_amount whose priority is lower.
#minimize{C@3: last_cardinality(C)}.

%holistic less atoms including repetition (shorter rules) as possible
%#minimize{C@2: atoms_amount(C)}.
%holistic less atoms excluding repetition as possible
%#minimize{C@1: atoms_cnt(C)}.

%:- lastRuleHead(e2,6,56).

% input(e1,3,58).
% input(e1,5,59).
% input(e1,7,83).
% input(e1,8,87).
% input(e1,10,97).
% output(e1,6,75).

% input(e2,1,48).
% input(e2,2,43).
% input(e2,5,14).
% input(e2,7,42).
% input(e2,8,56).
% input(e2,10,12).
% output(e2,6,56).

% lastRuleHead(e1,6,75).
% lastRuleHead(e2,6,56).
% lastRulePos(e1,6,75,8).



% input(e1,3,27).
% input(e1,4,86).
% input(e1,6,43).
% input(e1,9,84).
% output(e1,3,27).
% output(e1,4,27).
% output(e1,9,43).
% input(e2,2,77).
% input(e2,4,62).
% input(e2,5,28).
% input(e2,7,38).
% input(e2,9,90).
% input(e2,10,77).
% output(e2,4,28).
% output(e2,7,62).
% output(e2,8,57).
% input(e3,1,32).
% input(e3,4,97).
% input(e3,5,5).
% input(e3,7,46).
% input(e3,8,65).
% output(e3,2,32).
% output(e3,6,65).
% output(e3,7,67).


% universe(1 .. 10). 
%top(100). 
