%changed two places
% get the support %
hypotheses(I, N) :- output(I, N), N > 0.
{rule_pos(I,A)} :- hypotheses(I,_), output(I,  N), input(I, A, NA), NA >= N.
{rule_neg(I,A)} :- hypotheses(I,_), universe(A), not input(I, A, _). 

% get the limitation %
limitation(I, N) :- output(I, N).
limitPos(I, A) :- output(I,  N), input(I, A, NA), NA > N.
limitNeg(I, A) :- limitation(I,_), universe(A), not input(I, A, _).  

%analyze the four conditions of limitation which causes incompatibility, in order to set necessity constraints %
%hypothesis I and limitation J share the same head H
hypoLimitHead(I,J) :- hypotheses(I,_),limitation(J,_).
%positive body of I is not a subset of its corresponding limit_pos J. Namely, condition 1 of incompatibility unsatisfied
limitPosUnsat(I,J) :- rule_pos(I,A),hypoLimitHead(I,J),not limitPos(J,A).
limitNegUnsat(I,J) :- rule_neg(I,A),hypoLimitHead(I,J),not limitNeg(J,A).
%if the body part, viz. condition 1 and 2, of incompatibility is satisfied, the N of I must be constrainted by J. otherwise, free N.
necessityLimited(I,J) :- hypoLimitHead(I,J),not limitPosUnsat(I,J),not limitNegUnsat(I,J).
%it is not allowed that three conditions of incompatibility are inevitably satisfied simultaneously. 
%This constraint makes the solving faster, although the premised coherency makes this check meaningless.
:- necessityLimited(I,J),limitation(J,N1),hypotheses(I,N2),N1 < N2.


% cognize combinable rules  %
%hypothesis I1 and hypothesis J2 share the same head H. Rule out reflexivity via an order I1 != I2.
hypoHead(I1,I2) :- hypotheses(I1,_),hypotheses(I2,_),I1 != I2.
%I1 has some positive atom beyond the counterpart of I2.
notinPos(I1,I2) :- hypoHead(I1,I2),rule_pos(I1,A),not rule_pos(I2,A).
notinNeg(I1,I2) :- hypoHead(I1,I2),rule_neg(I1,A),not rule_neg(I2,A).
%the positive body of I1 is equivalent to the counterpart of I2. Rule out symmetry via an order I1 < I2.
atomPosCoincide(I1,I2) :- hypoHead(I1,I2),not notinPos(I1,I2), not notinPos(I2,I1), I1 < I2.
atomNegCoincide(I1,I2) :- hypoHead(I1,I2),not notinNeg(I1,I2), not notinNeg(I2,I1), I1 < I2.
%body I1 and body I2 are combinable.
atomCombinable(I1,I2) :- atomPosCoincide(I1,I2), atomNegCoincide(I1,I2).
% figure out the group identification and its boundaries %
%the body of I1 acts as an identification of group in which I2 is a member of it.
groupIdMember(I1,I2) :- atomCombinable(I1,I2), not atomCombinable(_,I1).
%I is a group id
groupId(I) :- groupIdMember(I,_).
%I2 has been covered by some group.
coveredRule(I2):- groupIdMember(I1,I2).
%In group I1, get the holistic lower bound L which meets all boundaries
ruleCombinedHolistic(I1,L) :- groupId(I1), hypotheses(I1,N1),  L = #max{N1; N2: groupIdMember(I1, I2), hypotheses(I2,N2)}.

%for those single hypothesis I, get the L and U
ruleCombinedHolistic(I,N) :- hypotheses(I,N), not groupId(I),not coveredRule(I).
%any holistic rule should satisfy U >= L
%:- ruleCombinedHolistic(_,L,U),L > U.


% take background into consideration %
%hypothesis I and background rule J share the same head H
hypoBgHead(I,J) :- hypotheses(I,_),bg_head(J).
%whether the hypothesis I in P  is removable because of the rule J in background B
combinablePB(I,J) :- hypoBgHead(I,J),ruleCombinedHolistic(I,L),bg_necessity(J,N), L <= N.
%at least one atom in the positive body of I belonging to P is outside the counterpart of J belonging to B
posPNotinB(I,J) :- combinablePB(I,J),rule_pos(I,A),not bg_body_pos(J,A).
%at least one atom in the positive body of J belonging to B is outside the counterpart of I belonging to P
posBNotinP(I,J) :- combinablePB(I,J),not rule_pos(I,A),bg_body_pos(J,A).
negPNotinB(I,J) :- combinablePB(I,J),rule_neg(I,A),not bg_body_neg(J,A).
negBNotinP(I,J) :- combinablePB(I,J),not rule_neg(I,A),bg_body_neg(J,A).
%the hypothesis I in P  is removed because of the rule J in background
combinedPB(I,J) :- combinablePB(I,J),not posPNotinB(I,J),not posBNotinP(I,J),not negPNotinB(I,J),not negBNotinP(I,J).


% settle the last rules and its cardinality, as well as atoms count %
%the last rule I whose head is H has the necessity L
last_rule(I,L) :- ruleCombinedHolistic(I,L),hypotheses(I,_),not combinedPB(I,_).


% get the three measures %
%the last program has C rules.count once for the entire group, while count once for every rule outside a group.
last_cardinality(C) :- C = #count{I: last_rule(I,_)}.
%the amount of all atoms appearing in the program.
%atoms_amount(Pos_cnt + Neg_cnt + Head_cnt) :- Pos_cnt = #count{(I,A): rule_pos(I,A),last_rule(I,_,_) },Neg_cnt = #count{(I,A): rule_neg(I,A),last_rule(I,_,_)},last_cardinality(Head_cnt).
%the count of atoms in the program
%atoms_cnt(Atoms_cnt) :- Atoms_cnt = #count{A: last_rule(I,_,_),rule_pos(I,A); A: last_rule(I,_,_),rule_neg(I,A); H:last_rule(I,_,_),hypotheses(I,H)}.


% point out the optimal directions %
%Change the two into weak constraints
%less rules as possible. This measure is more significant than atoms_amount whose priority is lower.
#minimize{C@3: last_cardinality(C)}.
%holistic less atoms including repetition (shorter rules) as possible
%#minimize{C@2: atoms_amount(C)}.
%holistic less atoms excluding repetition as possible
%#minimize{C@1: atoms_cnt(C)}.


#show rule_pos/2.
#show rule_neg/2.
#show last_rule/2.





% universe(1 .. 10). 


% output(e1,0).
% output(e2,17).
% output(e3,71).
% output(e4,0).
% output(e5,0).
% output(e6,0).
% output(e7,0).
% output(e8,42).
% output(e9,8).
% output(e10,0).
% output(e11,0).
% output(e12,0).
% output(e13,0).
% output(e14,0).
% output(e15,43).
% output(e16,29).
% output(e17,64).
% output(e18,0).
% output(e19,13).
% output(e20,22).
% output(e21,0).
% output(e22,15).
% output(e23,37).
% output(e24,28).
% output(e25,0).
% output(e26,0).
% output(e27,25).
% output(e28,0).
% output(e29,0).
% output(e30,0).
% output(e31,67).
% output(e32,0).
% output(e33,59).
% output(e34,31).
% output(e35,0).
% output(e36,67).
% output(e37,0).
% output(e38,0).
% output(e39,0).
% output(e40,0).
% input(e1,1,94).
% input(e1,2,4).
% input(e1,5,99).
% input(e1,8,7).
% input(e1,9,9).
% input(e1,10,85).
% input(e2,2,3).
% input(e2,3,42).
% input(e2,4,91).
% input(e2,5,47).
% input(e2,6,17).
% input(e2,7,68).
% input(e2,9,14).
% input(e2,10,40).
% input(e3,1,39).
% input(e3,2,65).
% input(e3,4,73).
% input(e3,5,62).
% input(e3,6,83).
% input(e3,7,71).
% input(e3,8,47).
% input(e3,9,49).
% input(e4,1,91).
% input(e4,6,16).
% input(e4,7,24).
% input(e4,8,40).
% input(e4,10,99).
% input(e5,2,4).
% input(e5,4,39).
% input(e5,6,45).
% input(e5,8,61).
% input(e5,10,86).
% input(e6,3,26).
% input(e6,4,58).
% input(e6,5,86).
% input(e6,6,74).
% input(e6,9,48).
% input(e6,10,89).
% input(e7,3,60).
% input(e7,6,23).
% input(e7,7,32).
% input(e8,2,3).
% input(e8,3,94).
% input(e8,4,42).
% input(e8,7,62).
% input(e8,9,49).
% input(e8,10,38).
% input(e9,2,12).
% input(e9,3,41).
% input(e9,4,20).
% input(e9,6,35).
% input(e9,7,8).
% input(e9,8,1).
% input(e9,9,60).
% input(e10,2,71).
% input(e11,5,2).
% input(e11,6,14).
% input(e11,7,47).
% input(e11,8,35).
% input(e11,9,39).
% input(e11,10,78).
% input(e12,1,21).
% input(e12,2,68).
% input(e12,6,97).
% input(e12,7,62).
% input(e12,9,34).
% input(e12,10,59).
% input(e13,1,36).
% input(e13,2,100).
% input(e13,3,0).
% input(e13,6,49).
% input(e13,7,89).
% input(e13,8,53).
% input(e13,10,96).
% input(e14,2,77).
% input(e14,3,70).
% input(e14,5,72).
% input(e14,6,72).
% input(e14,10,5).
% input(e15,1,2).
% input(e15,2,19).
% input(e15,4,43).
% input(e15,5,80).
% input(e15,7,94).
% input(e15,8,100).
% input(e15,9,86).
% input(e16,2,57).
% input(e16,3,80).
% input(e16,4,29).
% input(e16,5,54).
% input(e16,7,57).
% input(e16,8,70).
% input(e16,10,7).
% input(e17,2,67).
% input(e17,3,89).
% input(e17,4,84).
% input(e17,5,9).
% input(e17,6,64).
% input(e17,7,77).
% input(e17,8,64).
% input(e17,9,62).
% input(e17,10,51).
% input(e18,2,26).
% input(e18,3,80).
% input(e18,5,11).
% input(e18,6,94).
% input(e18,7,88).
% input(e18,8,14).
% input(e19,3,91).
% input(e19,4,99).
% input(e19,5,61).
% input(e19,6,13).
% input(e19,7,43).
% input(e19,8,69).
% input(e19,9,42).
% input(e20,2,82).
% input(e20,3,32).
% input(e20,4,22).
% input(e20,7,23).
% input(e20,8,80).
% input(e20,9,66).
% input(e20,10,87).
% input(e21,4,60).
% input(e21,5,38).
% input(e21,6,45).
% input(e21,8,58).
% input(e22,1,97).
% input(e22,2,45).
% input(e22,7,15).
% input(e22,8,9).
% input(e23,2,1).
% input(e23,4,43).
% input(e23,5,98).
% input(e23,6,96).
% input(e23,7,37).
% input(e23,8,46).
% input(e23,10,67).
% input(e24,4,28).
% input(e24,9,99).
% input(e24,10,50).
% input(e25,1,62).
% input(e25,2,6).
% input(e25,4,25).
% input(e25,6,61).
% input(e25,10,80).
% input(e26,1,11).
% input(e26,3,46).
% input(e26,5,35).
% input(e26,10,9).
% input(e27,1,84).
% input(e27,2,60).
% input(e27,3,4).
% input(e27,4,25).
% input(e27,5,94).
% input(e27,7,53).
% input(e27,8,11).
% input(e27,9,69).
% input(e28,1,79).
% input(e28,2,12).
% input(e28,3,56).
% input(e28,8,35).
% input(e28,10,23).
% input(e29,1,93).
% input(e29,2,94).
% input(e29,3,33).
% input(e29,4,69).
% input(e29,6,13).
% input(e29,8,15).
% input(e29,9,49).
% input(e30,3,23).
% input(e30,6,78).
% input(e30,7,90).
% input(e30,8,16).
% input(e30,9,45).
% input(e30,10,98).
% input(e31,1,81).
% input(e31,4,84).
% input(e31,8,81).
% input(e31,10,15).
% input(e32,1,51).
% input(e32,2,36).
% input(e32,6,43).
% input(e32,8,81).
% input(e32,9,42).
% input(e33,1,3).
% input(e33,2,8).
% input(e33,4,59).
% input(e33,5,3).
% input(e33,6,75).
% input(e33,7,88).
% input(e33,8,17).
% input(e33,9,96).
% input(e34,2,67).
% input(e34,4,31).
% input(e34,8,37).
% input(e34,9,44).
% input(e34,10,89).
% input(e35,1,98).
% input(e35,2,48).
% input(e35,3,19).
% input(e35,9,6).
% input(e35,10,16).
% input(e36,1,27).
% input(e36,3,7).
% input(e36,4,92).
% input(e36,5,75).
% input(e36,7,1).
% input(e37,3,4).
% input(e37,4,22).
% input(e37,5,4).
% input(e37,6,83).
% input(e37,10,66).
% input(e38,1,26).
% input(e38,2,89).
% input(e38,5,17).
% input(e38,6,48).
% input(e38,7,92).
% input(e38,8,43).
% input(e38,9,13).
% input(e39,1,100).
% input(e39,2,20).
% input(e39,6,100).
% input(e39,10,26).
% input(e40,3,16).
% input(e40,5,95).
% input(e40,6,90).
% input(e40,8,84).
