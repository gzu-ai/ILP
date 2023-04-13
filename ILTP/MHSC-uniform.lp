#defined bgRuleHead/3.
#defined bgRulePos/2.
#defined bgRuleNeg/2.

% get the support %
hypotheses(I, H, N) :- output(I, H, N), N > 0.
supportPos(I, H, A) :- hypotheses(I, H, N), input(I, A, NA), NA >= N.
supportNeg(I, H, A) :- hypotheses(I, H, N), universe(A), not input(I, A, _). 
supportMargin(I, H, A) :- hypotheses(I, H, N), input(I, A, NA), NA == N.

%get the concomitant outputs
%coOutput(I, A, 0) :- input(I, _, _), universe(A), not output(I, A, _).

% get the limitation %
%limitation(I, H, N) :- coOutput(I, H, N).
limitation(I, H, N) :- output(I, H, N).
limitPos(I, H, A) :- limitation(I, H, N), input(I, A, NA), NA > N.
limitNeg(I, H, A) :- universe(A), limitation(I, H, N), not input(I, A, _). 

% set the body atoms of the rule candidate. generate only the subset of supporting set  %
{rulePos(I, H, A)} :- supportPos(I, H, A).
{ruleNeg(I, H, A)} :- supportNeg(I, H, A).

%analyze the four conditions of limitation which causes incompatibility, in order to set necessity constraints %
limitStrictlySubsume(I, J, H) :- hypotheses(I, H, _),limitation(J, H, _).
notLimitSubsetPos(I, J, H) :- limitStrictlySubsume(I, J, H), rulePos(I, H, A), not limitPos(J, H, A).
notLimitSubsetNeg(I, J, H) :- limitStrictlySubsume(I, J, H), ruleNeg(I, H, A), not limitNeg(J, H, A).
mayLimited(I, J, H) :- limitStrictlySubsume(I, J, H), not notLimitSubsetPos(I, J, H), not notLimitSubsetNeg(I, J, H).
:- mayLimited(I, J, H), hypotheses(I, H, N1),limitation(J, H, N2), N1 > N2.


% set the two boundaries of necessity of the rule candidate %
%if lower_bound exists
lower(I, H, N) :- hypotheses(I, H, N).
% set the type of necessity of the rule candidate %
%if positive body separates from supporting margin, the necessity is fixed.
unfixed(I, H) :- hypotheses(I, H, _),rulePos(I, H, A),supportMargin(I, H, A). %a range
fixed(I, H) :- hypotheses(I, H, _),not unfixed(I, H).%a value

%if the necessity is fixed, the upper bound is equal to the support necessity. So it is N
upper(I, H, N) :- fixed(I, H), hypotheses(I, H, N).

%if the necessity is unfixed, and limited, then the only possible value of its upper bound is limit_necessity
%Here, a necessity could be constrainted by many limitation. So we choose the least.
upper(I, H, N_min) :- unfixed(I, H),N_min = #min{N: mayLimited(I, J, H),limitation(J, H, N)}.
%if the necessity is unfixed, and not limited, then set the default upper bound top.
%upper(I, H, Y) :- unfixed(I, H), top(Y).

% cognize combinable rules  %
%hypothesis I1 and hypothesis J2 share the same head H. Rule out reflexivity via an order I1 != I2.
hypoHead(I1,I2,H) :- hypotheses(I1, H, N1),hypotheses(I2,H,N2),I1 != I2.
%I1 has some positive atom beyond the counterpart of I2.
notinPos(I1,I2,H) :- hypoHead(I1,I2,H),rulePos(I1,H,A),not rulePos(I2,H,A).
notinNeg(I1,I2,H) :- hypoHead(I1,I2,H),ruleNeg(I1,H,A),not ruleNeg(I2,H,A).
%the positive body of I1 is equivalent to the counterpart of I2. Rule out symmetry via an order I1 < I2.
atomPosCoincide(I1,I2,H) :- hypoHead(I1,I2,H),not notinPos(I1,I2,H), not notinPos(I1,I2,H), I1 < I2.
atomNegCoincide(I1,I2,H) :- hypoHead(I1,I2,H),not notinNeg(I1,I2,H), not notinNeg(I1,I2,H), I1 < I2.
%body I1 and body I2 are combinable.
atomCombinable(I1,I2,H) :- atomPosCoincide(I1,I2,H), atomNegCoincide(I1,I2,H).
% figure out the group identification and its boundaries %
%the body of I1 acts as an identification of group in which I2 is a member of it.
groupIdMember(I1,I2,H) :- atomCombinable(I1,I2,H), not atomCombinable(_,I1,H).
%I is a group id
groupId(I,H) :- groupIdMember(I,_,H).
%I2 has been covered by some group.
coveredRule(I,H):- groupIdMember(_,I,H).
%In group I1, get the holistic lower bound L and upper bound U
ruleCombinedHolistic(I1,H,L,U) :- groupId(I1,H), 
L = #max{L2: groupIdMember(I1,I2,H),lower(I2,H,L2); L1: lower(I1,H,L1)},
U = #min{U2: groupIdMember(I1,I2,H),upper(I2,H,U2); U1: upper(I1,H,U1)}.
%for those single hypothesis I, get the L and U
ruleCombinedHolistic(I,H,L,U) :- hypotheses(I,H,_),not groupId(I,H),not coveredRule(I,H), lower(I,H,L), upper(I,H,U).
%any holistic rule should satisfy U >= L
:- ruleCombinedHolistic(_,_,L,U),L > U.

% take background into consideration %
%hypothesis I and background rule J share the same head H
hypoBgHead(I,R,H) :- hypotheses(I,H,_),bgRuleHead(R,H,_).
%whether the hypothesis I in P  is removable because of the rule J in background B
combinablePB(I,R,H) :- hypoBgHead(I,R,H), ruleCombinedHolistic(I,H,L,U), bgRuleHead(R,H, NR), NR >= L, NR <= U.

%at least one atom in the positive body of I belonging to P is outside the counterpart of J belonging to B
posPNotinB(I,R,H) :- combinablePB(I,R,H),rulePos(I,H,A),not bgRulePos(R,A).
%at least one atom in the positive body of J belonging to B is outside the counterpart of I belonging to P
posBNotinP(I,R,H) :- combinablePB(I,R,H),not rulePos(I,H,A),bgRulePos(R,A).
negPNotinB(I,R,H) :- combinablePB(I,R,H),ruleNeg(I,H,A),not bgRuleNeg(R,A).
negBNotinP(I,R,H) :- combinablePB(I,R,H),not ruleNeg(I,H,A),bgRuleNeg(R,A).
%the hypothesis I in P  is removed because of the rule J in background
combinedPB(I,R,H) :- combinablePB(I,R,H),not posPNotinB(I,R,H),not posBNotinP(I,R,H),not negPNotinB(I,R,H),not negBNotinP(I,R,H).

% settle the last rules and its cardinality, as well as atoms count %
%the last rule I whose head is H has the necessity L
lastRuleHead(I,H,L) :- ruleCombinedHolistic(I,H,L,_),hypotheses(I,H,_),not combinedPB(I,_,H).
lastRulePos(I, H, N, A) :- lastRuleHead(I, H, N), rulePos(I,H, A).
lastRuleNeg(I, H, N, A) :- lastRuleHead(I, H, N), ruleNeg(I,H, A).



#show lastRuleHead/3.
#show lastRulePos/4.
#show lastRuleNeg/4.
%#show atomCombinable/5.


% get the three measures %
%the last program has C rules.count once for the entire group, while count once for every rule outside a group.
lastCardinality(C) :- C = #count{(I,H): lastRuleHead(I,H,_)}.


% point out the optimal directions %
%Change the two into weak constraints
%less rules as possible. This measure is more significant than atoms_amount whose priority is lower.
#minimize{C@3: lastCardinality(C)}.

% universe(1 .. 10). 

% input(1 .. 10,2,63).
% input(1 .. 10,3,61).
% input(1 .. 10,5,44).
% input(1 .. 10,6,5).
% input(1 .. 10,8,88).
% output(1,1,0).
% output(2,2,0).
% output(3,3,5).
% output(4,4,0).
% output(5,5,0).
% output(6,6,75).
% output(7,7,0).
% output(8,8,0).
% output(9,9,5).
% output(10,10,0).
% input(11 .. 20,2,56).
% input(11 .. 20,4,28).
% input(11 .. 20,6,20).
% input(11 .. 20,8,60).
% input(11 .. 20,10,86).
% output(11,1,0).
% output(12,2,0).
% output(13,3,0).
% output(14,4,28).
% output(15,5,0).
% output(16,6,60).
% output(17,7,0).
% output(18,8,57).
% output(19,9,0).
% output(20,10,0).
% input(21 .. 30,3,74).
% input(21 .. 30,8,93).
% input(21 .. 30,9,86).
% output(21,1,0).
% output(22,2,0).
% output(23,3,74).
% output(24,4,0).
% output(25,5,0).
% output(26,6,75).
% output(27,7,0).
% output(28,8,0).
% output(29,9,32).
% output(30,10,0).
% input(31 .. 40,2,74).
% input(31 .. 40,3,73).
% input(31 .. 40,5,53).
% input(31 .. 40,6,32).
% input(31 .. 40,8,96).
% input(31 .. 40,10,67).
% output(31,1,0).
% output(32,2,0).
% output(33,3,0).
% output(34,4,0).
% output(35,5,32).
% output(36,6,75).
% output(37,7,0).
% output(38,8,57).
% output(39,9,32).
% output(40,10,0).
% input(41 .. 50,2,90).
% input(41 .. 50,3,89).
% input(41 .. 50,6,27).
% input(41 .. 50,7,66).
% input(41 .. 50,9,62).
% input(41 .. 50,10,29).
% output(41,1,0).
% output(42,2,0).
% output(43,3,0).
% output(44,4,0).
% output(45,5,0).
% output(46,6,0).
% output(47,7,0).
% output(48,8,29).
% output(49,9,27).
% output(50,10,0).
% input(51 .. 60,2,78).
% input(51 .. 60,4,28).
% input(51 .. 60,5,19).
% input(51 .. 60,6,43).
% input(51 .. 60,7,72).
% input(51 .. 60,8,63).
% output(51,1,0).
% output(52,2,0).
% output(53,3,0).
% output(54,4,0).
% output(55,5,0).
% output(56,6,63).
% output(57,7,28).
% output(58,8,0).
% output(59,9,0).
% output(60,10,0).
% input(61 .. 70,1,78).
% input(61 .. 70,2,75).
% input(61 .. 70,4,65).
% input(61 .. 70,5,44).
% input(61 .. 70,8,10).
% input(61 .. 70,9,63).
% input(61 .. 70,10,87).
% output(61,1,0).
% output(62,2,78).
% output(63,3,0).
% output(64,4,0).
% output(65,5,0).
% output(66,6,10).
% output(67,7,65).
% output(68,8,57).
% output(69,9,0).
% output(70,10,0).
% input(71 .. 80,5,77).
% input(71 .. 80,9,98).
% output(71,1,0).
% output(72,2,0).
% output(73,3,0).
% output(74,4,0).
% output(75,5,72).
% output(76,6,0).
% output(77,7,0).
% output(78,8,0).
% output(79,9,32).
% output(80,10,0).