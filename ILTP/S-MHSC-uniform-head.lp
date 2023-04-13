#defined bgRuleHead/2.
#defined bgRulePos/2.
#defined bgRuleNeg/2.

% get the support %
hypotheses(I, N) :- output(I,  N), N > 0.
supportPos(I, A) :- hypotheses(I,  N), input(I, A, NA), NA >= N.
supportNeg(I,  A) :- hypotheses(I,  N), universe(A), not input(I, A, _). 
supportMargin(I,  A) :- hypotheses(I,  N), input(I, A, NA), NA == N.

%get the concomitant outputs
coOutput(I,  0) :- input(I, _, _), not output(I,  _).

% get the limitation %
limitation(I,  N) :- coOutput(I,  N).
limitation(I,  N) :- output(I,  N).
limitPos(I,  A) :- limitation(I, N), input(I, A, NA), NA > N.
limitNeg(I,  A) :- universe(A), limitation(I,  N), not input(I, A, _). 

% set the body atoms of the rule candidate. generate only the subset of supporting set  %
{rulePos(I,  A)} :- supportPos(I,  A).
{ruleNeg(I,  A)} :- supportNeg(I,  A).

%analyze the four conditions of limitation which causes incompatibility, in order to set necessity constraints %
limitStrictlySubsume(I, J) :- hypotheses(I, _),limitation(J,  _).
notLimitSubsetPos(I, J) :- limitStrictlySubsume(I, J), rulePos(I,  A), not limitPos(J,  A).
notLimitSubsetNeg(I, J) :- limitStrictlySubsume(I, J), ruleNeg(I,  A), not limitNeg(J,  A).
mayLimited(I, J) :- limitStrictlySubsume(I, J), not notLimitSubsetPos(I, J), not notLimitSubsetNeg(I, J).
:- mayLimited(I, J), hypotheses(I,  N1),limitation(J,  N2), N1 > N2.


% set the two boundaries of necessity of the rule candidate %
%if lower_bound exists
lower(I,  N) :- hypotheses(I,  N).
% set the type of necessity of the rule candidate %
%if positive body separates from supporting margin, the necessity is fixed.
unfixed(I) :- hypotheses(I,  _),rulePos(I,  A),supportMargin(I,  A). %a range
fixed(I) :- hypotheses(I,  _),not unfixed(I).%a value

%if the necessity is fixed, the upper bound is equal to the support necessity. So it is N
upper(I,  N) :- fixed(I), hypotheses(I,  N).

%if the necessity is unfixed, and limited, then the only possible value of its upper bound is limit_necessity
%Here, a necessity could be constrainted by many limitation. So we choose the least.
upper(I, N_min) :- unfixed(I),N_min = #min{N: mayLimited(I, J),limitation(J, N)}.

% cognize combinable rules  %
%hypothesis I1 and hypothesis J2 share the same head H. Rule out reflexivity via an order I1 != I2.
hypoHead(I1,I2) :- hypotheses(I1, N1),hypotheses(I2,N2),I1 != I2.
%I1 has some positive atom beyond the counterpart of I2.
notinPos(I1,I2) :- hypoHead(I1,I2),rulePos(I1,A),not rulePos(I2,A).
notinNeg(I1,I2) :- hypoHead(I1,I2),ruleNeg(I1,A),not ruleNeg(I2,A).
%the positive body of I1 is equivalent to the counterpart of I2. Rule out symmetry via an order I1 < I2.
atomPosCoincide(I1,I2) :- hypoHead(I1,I2),not notinPos(I1,I2), not notinPos(I1,I2), I1 < I2.
atomNegCoincide(I1,I2) :- hypoHead(I1,I2),not notinNeg(I1,I2), not notinNeg(I1,I2), I1 < I2.
%body I1 and body I2 are combinable.
atomCombinable(I1,I2) :- atomPosCoincide(I1,I2), atomNegCoincide(I1,I2).
% figure out the group identification and its boundaries %
%the body of I1 acts as an identification of group in which I2 is a member of it.
groupIdMember(I1,I2) :- atomCombinable(I1,I2), not atomCombinable(_,I1).
%I is a group id
groupId(I) :- groupIdMember(I,_).
%I2 has been covered by some group.
coveredRule(I):- groupIdMember(_,I).
%In group I1, get the holistic lower bound L and upper bound U
ruleCombinedHolistic(I1,L,U) :- groupId(I1), 
L = #max{L2: groupIdMember(I1,I2),lower(I2,L2); L1: lower(I1,L1)},
U = #min{U2: groupIdMember(I1,I2),upper(I2,U2); U1: upper(I1,U1)}.
%for those single hypothesis I, get the L and U
ruleCombinedHolistic(I,L,U) :- hypotheses(I,_),not groupId(I),not coveredRule(I), lower(I,L), upper(I,U).
%any holistic rule should satisfy U >= L
:- ruleCombinedHolistic(_,L,U),L > U.

% take background into consideration %
%hypothesis I and background rule J share the same head H
hypoBgHead(I,R) :- hypotheses(I,_),bgRuleHead(R,_).
%whether the hypothesis I in P  is removable because of the rule J in background B
combinablePB(I,R) :- hypoBgHead(I,R), ruleCombinedHolistic(I,L,U), bgRuleHead(R, NR), NR >= L, NR <= U.

%at least one atom in the positive body of I belonging to P is outside the counterpart of J belonging to B
posPNotinB(I,R) :- combinablePB(I,R),rulePos(I,A),not bgRulePos(R,A).
%at least one atom in the positive body of J belonging to B is outside the counterpart of I belonging to P
posBNotinP(I,R) :- combinablePB(I,R),not rulePos(I,A),bgRulePos(R,A).
negPNotinB(I,R) :- combinablePB(I,R),ruleNeg(I,A),not bgRuleNeg(R,A).
negBNotinP(I,R) :- combinablePB(I,R),not ruleNeg(I,A),bgRuleNeg(R,A).
%the hypothesis I in P  is removed because of the rule J in background
combinedPB(I,R) :- combinablePB(I,R),not posPNotinB(I,R),not posBNotinP(I,R),not negPNotinB(I,R),not negBNotinP(I,R).

% settle the last rules and its cardinality, as well as atoms count %
%the last rule I whose head is H has the necessity L
lastRuleHead(I,L) :- ruleCombinedHolistic(I,L,_),hypotheses(I,_),not combinedPB(I,_).
lastRulePos(I, N, A) :- lastRuleHead(I, N), rulePos(I,A).
lastRuleNeg(I, N, A) :- lastRuleHead(I, N), ruleNeg(I,A).



#show lastRuleHead/2.
#show lastRulePos/3.
#show lastRuleNeg/3.
%#show atomCombinable/5.


% get the three measures %
%the last program has C rules.
%lastCardinality(C) :- C = #count{(I, H): lastRuleHead(I, H)}.
%the amount of all atoms appearing in the program.
atomsAmount(PosCnt + NegCnt + HeadCnt) :- PosCnt = #count{(I, N, A): lastRulePos(I, N, A) },NegCnt = #count{(I, N, A): lastRuleNeg(I, N, A) }, HeadCnt = #count{(I, H): lastRuleHead(I, H)}.
%the count of atoms in the program
%atomsCnt(AtomsCnt) :- AtomsCnt = #count{A: lastRulePos(I, N, A); B: lastRuleNeg(I, N, B); H:lastRuleHead(J, H)}.

% point out the optimal directions %
%Change the two into weak constraints
%less rules as possible. This measure is more significant than atomsAmount whose priority is lower.
%#minimize{C@3: lastCardinality(C)}.

%holistic less atoms including repetition (shorter rules) as possible
#minimize{C@2: atomsAmount(C)}.
%holistic less atoms excluding repetition as possible
%#minimize{C@1: atomsCnt(C)}.


% universe(1 .. 10). 

% output(e1,0).
% output(e2,0).
% output(e3,0).
% output(e4,59).
% input(e1,2,61).
% input(e1,3,1).
% input(e1,4,51).
% input(e1,6,4).
% input(e1,9,91).
% input(e1,10,90).
% input(e2,2,61).
% input(e2,6,10).
% input(e2,8,95).
% input(e3,5,37).
% input(e3,7,37).
% input(e3,10,42).
% input(e4,1,59).
% input(e4,6,19).

% bgRulePos(b1,7).
% bgRuleNeg(b1,9).
% bgRuleNeg(b1,2).
% bgRuleNeg(b1,5).
% bgRuleNeg(b1,10).
% bgRuleNeg(b1,1).
% bgRuleNeg(b1,8).
% bgRuleHead(b1,1).

