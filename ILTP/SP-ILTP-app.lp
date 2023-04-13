%%%% get the candidate solution space%%%%%%%%%
1 {lastRuleHead(J, H, N); trimedRuleHead(J, H, N)} 1 :- output(J, H, N).
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
limitStrictlySubsume(I, J, H, N1, N2) :- lastRuleHead(I, H, N1), limitation(J, H, N2), N1 > N2.
notLimitSubsetPos(I, J, H, N1, N2) :- limitStrictlySubsume(I, J, H, N1, N2), lastRulePos(I, H, N1, A), not limitPos(J, H, N2, A).
notLimitSubsetNeg(I, J, H, N1, N2) :- limitStrictlySubsume(I, J, H, N1, N2), lastRuleNeg(I, H, N1, A), not limitNeg(J, H, N2, A).
:- limitStrictlySubsume(I, J, H, N1, N2), not notLimitSubsetPos(I, J, H, N1, N2), not notLimitSubsetNeg(I, J, H, N1, N2).


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%given background covers a single transition
%the negative body part of applicable requirements is unsatisfied. Namely, a intersection exists.
notBodyNegSat(R,I) :- bgRuleNeg(R,A),input(I,A,_).
%the positive body part of applicable requirements is unsatisfied. Namely, not a subset.
notBodyPosSat(R,I) :- bgRulePos(R,A),output(I,_,_),not input(I,A,_).
%R is applicable in I, since all the three applicable requirements satisfy.
applicable(R, I, H, NR,NI) :- bgRuleHead(R,H, NR),output(I,H,NI),not notBodyNegSat(R,I),not notBodyPosSat(R,I).
%R cover <I,{(H,NI)}>
coverB(R, I, H, NI) :- applicable(R, I, H, NR, NI),#min{N2 : input(I,A,N2),bgRulePos(R,A);NR} = BETA , BETA == NI.

%existed rule from example covers a single transition
notBodyNegSatE(I, H, N, J) :- lastRuleNeg(I, H, N, A),input(J,A,_), I != J.
notBodyPosSatE(I, H, N, J) :- lastRulePos(I, H, N, A),output(J,_,_),not input(J,A,_), I != J.
applicableE(I, J, H, N, NI) :- lastRuleHead(I, H, N),output(J,H,NI),not notBodyNegSatE(I, H, N, J),not notBodyPosSatE(I, H, N,J), I != J.
coverE(I, J, H, N, NI) :- applicableE(I, J, H, N, NI),#min{N2 : input(J,A,N2),lastRulePos(I, H, N, A);N} = BETA , BETA == NI.


covered(I, H, N) :- coverB(_, I, H, N).
covered(I, H, N) :- coverE(_, I, H, _, N).

%a single transition must has corresponding rule at last if it is not covered by any existed rule.
:- trimedRuleHead(I, H, N), not covered(I, H, N).

%%%%leave the unsubsumed rules  as the last rules only%%%%%%%%%%
%the head must exist, whereas the body does not.
%{lastRuleHead(J, H, N)} :- output(J, H, N).
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
lastCardinality(C) :- C = #count{(I, H, N): lastRuleHead(I, H, N)}.
%the amount of all atoms appearing in the program.
atomsAmount(PosCnt + NegCnt + HeadCnt) :- PosCnt = #count{(I, H, N, A): lastRulePos(I, H, N, A) },NegCnt = #count{(I, H, N, A): lastRuleNeg(I, H, N, A) }, HeadCnt = #count{(I, H, N): lastRuleHead(I, H, N)}.
%the count of atoms in the program
atomsCnt(AtomsCnt) :- AtomsCnt = #count{A: lastRulePos(I, H, N, A); A: lastRuleNeg(I, H, N, A); H:lastRuleHead(J, H, N2)}.

% point out the optimal directions %
%Change the two into weak constraints
%less rules as possible. This measure is more significant than atomsAmount whose priority is lower.
#minimize{C@3: lastCardinality(C)}.

%holistic less atoms including repetition (shorter rules) as possible
#minimize{C@2: atomsAmount(C)}.
%holistic less atoms excluding repetition as possible
#minimize{C@1: atomsCnt(C)}.
