% universe(1 .. 3). 


% input(e3,2,50). 
% output(e3,3,30).

% bgRuleHead(b1, 3, 40).
% bgRulePos(b1, 2).
% bgRuleNeg(b1, 1).
% bgRuleNeg(b1, 3).

universe(1 .. 2). 


input(e3,1,60). 
output(e3,2,60).

bgRuleHead(b1, 2, 70).
% bgRulePos(b1, 2).
% bgRuleNeg(b1, 1).
% bgRuleNeg(b1, 3).

%get the concomitant outputs
coOutput(I, A, 0) :- input(I, _, _), universe(A), not output(I, A, _).
anyOutput(I, A, N) :- coOutput(I, A, N).
anyOutput(I, A, N) :- output(I, A, N).

%the negative body part of applicable requirements is unsatisfied. Namely, a intersection exists.
notBodyNegSat(R,I) :- bgRuleNeg(R,A),input(I,A,_).
%the positive body part of applicable requirements is unsatisfied. Namely, not a subset.
notBodyPosSat(R,I) :- bgRulePos(R,A),anyOutput(I,_,_),not input(I,A,_).
%R is applicable in I, since all the three applicable requirements satisfy.
applicable(R, I, H, NR,NI) :- bgRuleHead(R,H, NR),anyOutput(I,H,NI),not notBodyNegSat(R,I),not notBodyPosSat(R,I).
%I is incompatible w.r.t. R after taking BETA into consideration.
incompatible(R, I, H) :- applicable(R, I, H, NR, NI),#min{N2 : input(I,A,N2),bgRulePos(R,A);NR} = BETA , BETA > NI.
%if existing incompatible situation, E is inconsistent w.r.t. B.
incosistent :- incompatible(_,_,_).
%if not inconsistent, than consistent.
consistent :- not incosistent.

#show incosistent/0.
#show consistent/0.

% #show incompatible/3.
% #show applicable/5.
% #show notBodyPosSat/2.