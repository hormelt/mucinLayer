% lipidcomp.m
% calculates the volumes of stock solutions to use for lipid mixtures of desired composition
%
% Raghuveer Parthasarathy Oct. 28, 2002
% last modified July 2, 2008 -- added ZB-2200-8; streamlined calc.
%

clear all



% -------------------------------------------------------------------------


fs = 'lipidcomp.m --- calculations for lipid mixtures.  RP 10/02';
disp(fs);

% table of molecular weights --- g/mol

j=1;  % counter
lipid(j) = {'DMPC'};   mw(j) = 677.95; j = j+1;
lipid(j) = {'Egg-PC'}; mw(j) = 760.08; j = j+1;
lipid(j) = {'DPPC'};   mw(j) = 734.05; j = j+1;
lipid(j) = {'DOPC'};   mw(j) = 785.6;  j = j+1;
lipid(j) = {'DODAP'};  mw(j) = 648.06; j = j+1;
lipid(j) = {'DOTAP'};  mw(j) = 698.55; j = j+1;
lipid(j) = {'DOPS'};   mw(j) = 810.04; j = j+1;
lipid(j) = {'DOPA'};   mw(j) = 722.96; j = j+1;
lipid(j) = {'DOPE'};   mw(j) = 744.04; j = j+1;
lipid(j) = {'DMPS'};   mw(j) = 701.85; j = j+1;

lipid(j) = {'cholesterol'}; mw(j) = 386.66; j = j+1;
lipid(j) = {'sphingomyelin (egg)'}; mw(j) = 703.03; j = j+1;

lipid(j) = {'Plant-PI'};  mw(j) = 857.04;  j = j+1;
lipid(j) = {'PI(4)P'};    mw(j) = 1001.18; j = j+1;
lipid(j) = {'PI(4,5)P2'}; mw(j) = 1098.19; j = j+1;
lipid(j) = {'GM1'};       mw(j) = 1545;    j = j+1;

lipid(j) = {'TR-DHPE'};  mw(j) = 1381.85; j = j+1;
lipid(j) = {'MB-DHPE'};  mw(j) = 944.14;  j = j+1;
lipid(j) = {'16:0 12:0 NBD-PC'}; mw(j) = 857.05; j = j+1;
lipid(j) = {'16:0  6:0 NBD PC'}; mw(j) = 771.4;  j = j+1;
lipid(j) = {'NBD-PS (head)'}; mw(j) = 985.2; j = j+1;
lipid(j) = {'NBD-PS (tail 18:1 12:0)'}; mw(j) = 901.1; j = j+1;

lipid(j) = {'16:0 biotinyl-cap-PE'}; mw(j) = 1053.4; j = j+1;
lipid(j) = {'18:1 biotinyl-cap-PE'}; mw(j) = 1105.48; j = j+1;

lipid(j) = {'Trehalose Dimycolate (TDM)'}; mw(j) = 2636.5; j = j+1;
lipid(j) = {'ZB-2200-15'}; mw(j) = 790.54; j = j+1;
lipid(j) = {'ZB-2200-8'}; mw(j) = 594.33; j = j+1;

%---  deleted
%lipid(j) = {'NBD-PG'};  mw(j) = 867.96; j = j+1;
%lipid(j) = {'gramicidin A'}; mw(j) = 1882.3; j = j+1;
%lipid(j) = {'DMEPC'}; mw(j) = 742.46; j = j+1;
%lipid(j) = {'DOPG'}; mw(j) = 797.04; j = j+1;
%lipid(j) = {'DNP-cap-PE'}; mw(j) = 993.20; j = j+1;
%lipid(j) = {'20:1 PC (cis)'}; mw(j) = 841.65; j = j+1;
%lipid(j) = {'13:0 PC'}; mw(j) = 649.47; j = j+1;
%lipid(j) = {'BODIPY-TR PI(4,5)P_2'}; mw(j) = 1730.81; j = j+1;
%lipid(j) = {'DOEPC'};  mw(j) = 850.64; j = j+1;
%lipid(j) = {'18:1 cap PE'}; mw(j) = 879.18; j = j+1;
%lipid(j) = {'16:0 cap PE'}; mw(j) = 805.13; j = j+1;

fs = ' '; disp(fs); disp(fs);

nlist = size(mw,2);
fs = sprintf('\t #   %s   \t %s', 'lipid', 'MW (g/mol)');
disp(fs);
for i=1:nlist,
   fs = sprintf('\t %d  %s   \t %.2f', i, char(lipid(i)), mw(i));
   disp(fs);
end
fs = sprintf('More lipids (not shown) -- see lipidcomp.m');
disp(fs);


% ---------------------------------------------------------------------
% Desired composition

fs = ' '; disp(fs); disp(fs); disp(fs);

nlipids = input('   Enter the number of components:  ');
nlipids = round(nlipids);  % make sure an integer!

W_T = input('   Enter the total lipid mass desired (mg):  ');
n = zeros(1, nlipids);
p = zeros(1, nlipids);
stock = zeros(1, nlipids);
for j=1:nlipids
   nt = input('   Enter the lipid number (from list) [99 for new lipid]:  ');
      if (nt==99)
         % User wants a lipid not on the list --- add
         nlist = nlist + 1;
         nt = nlist;
         temp = input('      * Enter name of new lipid:   ', 's');
         lipid(nlist) = {temp};
         mw(nlist) = input('      * Enter molecular weight of new lipid (g/mol):   ');
      end
      n(j) = nt;
      p(j) = input('     Enter the mole percent desired for this lipid: ');
          % note p(i) is percent, not fraction
      stock(j) = ... 
         input('     Enter the solution concentration of this lipid (mg/mL): ');
end
fs = ' '; disp(fs); disp('   Composition complete');
if (abs(sum(p)-100.0) > 0.01)
   fs = '***************************************************';
   disp(fs);
   fs = '*** ERROR:  sum of percentages does not equal 100!';
   disp(fs);
   fs = input('   Press Enter to continue --- recommend Ctrl-C to abort.');
   disp(fs);
end

% -------------------
% Calculation

fs = ' ';
for i=1:3;
   disp(fs);
end
mol_T = W_T ./ (10.0 * sum(p.*mw(n)));  % total number of moles (note W_T in mg)
mol = p.*mol_T/100; % moles, not needed, but will calculate and display
w = mol_T .* (p/100.0) .* mw(n) * 1000.0; % mass of constituent i, in mg
vol = 1000*w./stock; % volume of stock to use, in uL

disp('    ');
fs = sprintf('pct \t lipid \t Mol. wt. (g/mol) \t mol \t weight (mg) \t stock (mg/mL) \t vol (uL)');
disp(fs);
for i=1:nlipids,
   fs = sprintf('%.1f \t %s \t %.2f \t %.3e \t %.2e \t %.2f \t %.1f', ...
      p(i), char(lipid(n(i))), mw(n(i)), mol(i), w(i), stock(i), vol(i));
   disp(fs);
end
