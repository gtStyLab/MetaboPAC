function BST_penalty = pen_BSTfit_hynne(nT,absolute_concMatrix,Vcalc)
% Calculate penalty for fitting BST equation to metabolite and flux data.

% Set up linear equations
Av2 = [ones(nT,1), log(absolute_concMatrix(1:nT,2)), log(absolute_concMatrix(1:nT,8)), log(absolute_concMatrix(1:nT,17))];
Av3 = [ones(nT,1), log(absolute_concMatrix(1:nT,2)), log(absolute_concMatrix(1:nT,21))];
Av4 = [ones(nT,1), log(absolute_concMatrix(1:nT,8)), log(absolute_concMatrix(1:nT,12))];
Av5 = [ones(nT,1), log(absolute_concMatrix(1:nT,12)), log(absolute_concMatrix(1:nT,20)), log(absolute_concMatrix(1:nT,21))];
Av6 = [ones(nT,1), log(absolute_concMatrix(1:nT,6)), log(absolute_concMatrix(1:nT,7)), log(absolute_concMatrix(1:nT,19))];
Av7 = [ones(nT,1), log(absolute_concMatrix(1:nT,6)), log(absolute_concMatrix(1:nT,7))];
Av8 = [ones(nT,1), log(absolute_concMatrix(1:nT,6)), log(absolute_concMatrix(1:nT,9)), log(absolute_concMatrix(1:nT,18)), log(absolute_concMatrix(1:nT,22))];
Av9 = [ones(nT,1), log(absolute_concMatrix(1:nT,5)), log(absolute_concMatrix(1:nT,11)), log(absolute_concMatrix(1:nT,18)), log(absolute_concMatrix(1:nT,21))];
Av10 = [ones(nT,1), log(absolute_concMatrix(1:nT,5)), log(absolute_concMatrix(1:nT,11))];
Av11 = [ones(nT,1), log(absolute_concMatrix(1:nT,14))];
Av12 = [ones(nT,1), log(absolute_concMatrix(1:nT,1)), log(absolute_concMatrix(1:nT,22))];
Av13 = [ones(nT,1), log(absolute_concMatrix(1:nT,4)), log(absolute_concMatrix(1:nT,15))];
Av14 = [ones(nT,1), log(absolute_concMatrix(1:nT,15))];
Av15 = [ones(nT,1), log(absolute_concMatrix(1:nT,7)), log(absolute_concMatrix(1:nT,9)), log(absolute_concMatrix(1:nT,22))];
Av16 = [ones(nT,1), log(absolute_concMatrix(1:nT,3)), log(absolute_concMatrix(1:nT,16))];
Av17 = [ones(nT,1), log(absolute_concMatrix(1:nT,16))];
Av18 = [ones(nT,1), log(absolute_concMatrix(1:nT,1)), log(absolute_concMatrix(1:nT,10))];
Av19 = [ones(nT,1), log(absolute_concMatrix(1:nT,10))];
Av20 = [ones(nT,1), log(absolute_concMatrix(1:nT,10)), log(absolute_concMatrix(1:nT,13))];
Av21 = [ones(nT,1), log(absolute_concMatrix(1:nT,13))];
Av22 = [ones(nT,1), log(absolute_concMatrix(1:nT,8)), log(absolute_concMatrix(1:nT,21))];
Av23= [ones(nT,1), log(absolute_concMatrix(1:nT,21))];
Av24 = [ones(nT,1), log(absolute_concMatrix(1:nT,5)), log(absolute_concMatrix(1:nT,20)), log(absolute_concMatrix(1:nT,21))];
A = blkdiag(Av2, Av3, Av4, Av5, Av6, Av7, Av8, Av9, Av10, Av11, Av12, Av13, Av14, Av15, Av16, Av17, Av18, Av19, Av20, Av21, Av22, Av23, Av24);
b = reshape(log(Vcalc(:,2:24)),[numel(Vcalc(:,2:24)),1]);

% Solve for linearized best fit
x = A\b;

% Switch a terms out of log-space
x([1 5 8 11 15 19 22 27 32 35 37 40 43 45 49 52 54 57 59 62 64 67 69]) = exp(x([1 5 8 11 15 19 22 27 32 35 37 40 43 45 49 52 54 57 59 62 64 67 69]));

% Calculate predicted BST flux
BST_flux(:,1) = x(1).*absolute_concMatrix(1:nT,2).^x(2).*absolute_concMatrix(1:nT,8).^x(3).*absolute_concMatrix(1:nT,17).^x(4);
BST_flux(:,2) = x(5).*absolute_concMatrix(1:nT,2).^x(6).*absolute_concMatrix(1:nT,21).^x(7);
BST_flux(:,3) = x(8).*absolute_concMatrix(1:nT,8).^x(9).*absolute_concMatrix(1:nT,12).^x(10);
BST_flux(:,4) = x(11).*absolute_concMatrix(1:nT,12).^x(12).*absolute_concMatrix(1:nT,20).^x(13).*absolute_concMatrix(1:nT,21).^x(14);
BST_flux(:,5) = x(15).*absolute_concMatrix(1:nT,6).^x(16).*absolute_concMatrix(1:nT,7).^x(17).*absolute_concMatrix(1:nT,19).^x(18);
BST_flux(:,6) = x(19).*absolute_concMatrix(1:nT,6).^x(20).*absolute_concMatrix(1:nT,7).^x(21);
BST_flux(:,7) = x(22).*absolute_concMatrix(1:nT,6).^x(23).*absolute_concMatrix(1:nT,9).^x(24).*absolute_concMatrix(1:nT,18).^x(25).*absolute_concMatrix(1:nT,22).^x(26);
BST_flux(:,8) = x(27).*absolute_concMatrix(1:nT,5).^x(28).*absolute_concMatrix(1:nT,11).^x(29).*absolute_concMatrix(1:nT,18).^x(30).*absolute_concMatrix(1:nT,21).^x(31);
BST_flux(:,9) = x(32).*absolute_concMatrix(1:nT,5).^x(33).*absolute_concMatrix(1:nT,11).^x(34);
BST_flux(:,10) = x(35).*absolute_concMatrix(1:nT,14).^x(36);
BST_flux(:,11) = x(37).*absolute_concMatrix(1:nT,1).^x(38).*absolute_concMatrix(1:nT,22).^x(39);
BST_flux(:,12) = x(40).*absolute_concMatrix(1:nT,4).^x(41).*absolute_concMatrix(1:nT,15).^x(42);
BST_flux(:,13) = x(43).*absolute_concMatrix(1:nT,15).^x(44);
BST_flux(:,14) = x(45).*absolute_concMatrix(1:nT,7).^x(46).*absolute_concMatrix(1:nT,9).^x(47).*absolute_concMatrix(1:nT,22).^x(48);
BST_flux(:,15) = x(49).*absolute_concMatrix(1:nT,3).^x(50).*absolute_concMatrix(1:nT,16).^x(51);
BST_flux(:,16) = x(52).*absolute_concMatrix(1:nT,16).^x(53);
BST_flux(:,17) = x(54).*absolute_concMatrix(1:nT,1).^x(55).*absolute_concMatrix(1:nT,10).^x(56);
BST_flux(:,18) = x(57).*absolute_concMatrix(1:nT,10).^x(58);
BST_flux(:,19) = x(59).*absolute_concMatrix(1:nT,10).^x(60).*absolute_concMatrix(1:nT,13).^x(61);
BST_flux(:,20) = x(62).*absolute_concMatrix(1:nT,13).^x(63);
BST_flux(:,21) = x(64).*absolute_concMatrix(1:nT,8).^x(65).*absolute_concMatrix(1:nT,21).^x(66);
BST_flux(:,22) = x(67).*absolute_concMatrix(1:nT,21).^x(68);
BST_flux(:,23) = x(69).*absolute_concMatrix(1:nT,5).^x(70).*absolute_concMatrix(1:nT,20).^x(71).*absolute_concMatrix(1:nT,21).^x(72);

% Calculate error between predicted flux using BST and inferred flux
BST_error(1) = real(sqrt(sum(sum((Vcalc(1:nT,2) - BST_flux(:,1)).^2)))) + abs(imag(sqrt(sum(sum((Vcalc(1:nT,2) - BST_flux(:,1)).^2)))));
BST_error(2) = real(sqrt(sum(sum((Vcalc(1:nT,3) - BST_flux(:,2)).^2)))) + abs(imag(sqrt(sum(sum((Vcalc(1:nT,3) - BST_flux(:,2)).^2)))));
BST_error(3) = real(sqrt(sum(sum((Vcalc(1:nT,4) - BST_flux(:,3)).^2)))) + abs(imag(sqrt(sum(sum((Vcalc(1:nT,4) - BST_flux(:,3)).^2)))));
BST_error(4) = real(sqrt(sum(sum((Vcalc(1:nT,5) - BST_flux(:,4)).^2)))) + abs(imag(sqrt(sum(sum((Vcalc(1:nT,5) - BST_flux(:,4)).^2)))));
BST_error(5) = real(sqrt(sum(sum((Vcalc(1:nT,6) - BST_flux(:,5)).^2)))) + abs(imag(sqrt(sum(sum((Vcalc(1:nT,6) - BST_flux(:,5)).^2)))));
BST_error(6) = real(sqrt(sum(sum((Vcalc(1:nT,7) - BST_flux(:,6)).^2)))) + abs(imag(sqrt(sum(sum((Vcalc(1:nT,7) - BST_flux(:,6)).^2)))));
BST_error(7) = real(sqrt(sum(sum((Vcalc(1:nT,8) - BST_flux(:,7)).^2)))) + abs(imag(sqrt(sum(sum((Vcalc(1:nT,8) - BST_flux(:,7)).^2)))));
BST_error(8) = real(sqrt(sum(sum((Vcalc(1:nT,9) - BST_flux(:,8)).^2)))) + abs(imag(sqrt(sum(sum((Vcalc(1:nT,9) - BST_flux(:,8)).^2)))));
BST_error(9) = real(sqrt(sum(sum((Vcalc(1:nT,10) - BST_flux(:,9)).^2)))) + abs(imag(sqrt(sum(sum((Vcalc(1:nT,10) - BST_flux(:,9)).^2)))));
BST_error(10) = real(sqrt(sum(sum((Vcalc(1:nT,11) - BST_flux(:,10)).^2)))) + abs(imag(sqrt(sum(sum((Vcalc(1:nT,11) - BST_flux(:,10)).^2)))));
BST_error(11) = real(sqrt(sum(sum((Vcalc(1:nT,12) - BST_flux(:,11)).^2)))) + abs(imag(sqrt(sum(sum((Vcalc(1:nT,12) - BST_flux(:,11)).^2)))));
BST_error(12) = real(sqrt(sum(sum((Vcalc(1:nT,13) - BST_flux(:,12)).^2)))) + abs(imag(sqrt(sum(sum((Vcalc(1:nT,13) - BST_flux(:,12)).^2)))));
BST_error(13) = real(sqrt(sum(sum((Vcalc(1:nT,14) - BST_flux(:,13)).^2)))) + abs(imag(sqrt(sum(sum((Vcalc(1:nT,14) - BST_flux(:,13)).^2)))));
BST_error(14) = real(sqrt(sum(sum((Vcalc(1:nT,15) - BST_flux(:,14)).^2)))) + abs(imag(sqrt(sum(sum((Vcalc(1:nT,15) - BST_flux(:,14)).^2)))));
BST_error(15) = real(sqrt(sum(sum((Vcalc(1:nT,16) - BST_flux(:,15)).^2)))) + abs(imag(sqrt(sum(sum((Vcalc(1:nT,16) - BST_flux(:,15)).^2)))));
BST_error(16) = real(sqrt(sum(sum((Vcalc(1:nT,17) - BST_flux(:,16)).^2)))) + abs(imag(sqrt(sum(sum((Vcalc(1:nT,17) - BST_flux(:,16)).^2)))));
BST_error(17) = real(sqrt(sum(sum((Vcalc(1:nT,18) - BST_flux(:,17)).^2)))) + abs(imag(sqrt(sum(sum((Vcalc(1:nT,18) - BST_flux(:,17)).^2)))));
BST_error(18) = real(sqrt(sum(sum((Vcalc(1:nT,19) - BST_flux(:,18)).^2)))) + abs(imag(sqrt(sum(sum((Vcalc(1:nT,19) - BST_flux(:,18)).^2)))));
BST_error(19) = real(sqrt(sum(sum((Vcalc(1:nT,20) - BST_flux(:,19)).^2)))) + abs(imag(sqrt(sum(sum((Vcalc(1:nT,20) - BST_flux(:,19)).^2)))));
BST_error(20) = real(sqrt(sum(sum((Vcalc(1:nT,21) - BST_flux(:,20)).^2)))) + abs(imag(sqrt(sum(sum((Vcalc(1:nT,21) - BST_flux(:,20)).^2)))));
BST_error(21) = real(sqrt(sum(sum((Vcalc(1:nT,22) - BST_flux(:,21)).^2)))) + abs(imag(sqrt(sum(sum((Vcalc(1:nT,22) - BST_flux(:,21)).^2)))));
BST_error(22) = real(sqrt(sum(sum((Vcalc(1:nT,23) - BST_flux(:,22)).^2)))) + abs(imag(sqrt(sum(sum((Vcalc(1:nT,23) - BST_flux(:,22)).^2)))));
BST_error(23) = real(sqrt(sum(sum((Vcalc(1:nT,24) - BST_flux(:,23)).^2)))) + abs(imag(sqrt(sum(sum((Vcalc(1:nT,24) - BST_flux(:,23)).^2)))));

BST_penalty = nansum(BST_error);