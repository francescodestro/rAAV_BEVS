%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% states
% x(1) = target cell (#/mL)
%
% virions
% x(2) = BV rep (#/mL) 
% x(3) = BV cap (#/mL) 
% x(4) = BV GOI (#/mL)
%
% cells infected by one virus
% x(5) = rep
% x(6) = cap
% x(7) = goi
%
% cells infected by two viruses
% x(8) = rep_cap
% x(9) = rep_GOI
% x(10) = cap_goi
% 
% cells infected by three viruses
% x(11) = triple-infected
%
% content of infected cells (for i=1:7)
% x(12+(i-1)*22) = bound repBV
% x(13+(i-1)*22) = bound capBV
% x(14+(i-1)*22) = bound GOIBV
% x(15+(i-1)*22) = nuclear rep DNA
% x(16+(i-1)*22) = nuclear cap DNA
% x(17+(i-1)*22) = nuclear GOI DNA
% x(18+(i-1)*22) = rep78 mRNA
% x(19+(i-1)*22) = rep52 mRNA
% x(20+(i-1)*22) = VP mRNA
% x(21+(i-1)*22) = GOI mRNA
% x(22+(i-1)*22) = GOI number of copies
% x(23+(i-1)*22) = rep78 protein conc.
% x(24+(i-1)*22) = rep52 protein conc.
% x(25+(i-1)*22) = empty capsids conc.
% x(26+(i-1)*20)) = GOI protein
% x(27+(i-1)*22) = filled capsids conc.
% 
% content of dead cells
% x(28+(i-1)*22)) = number of dead triple-infected cells
% x(29+(i-1)*22)) = number of GOI copies in dead cells
% x(30+(i-1)*22)) = rep78 protein conc.
% x(31+(i-1)*22)) = rep52 protein conc.
% x(32+(i-1)*22)) = empty capsids conc.
% x(33+(i-1)*22)) = filled capsids conc.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%