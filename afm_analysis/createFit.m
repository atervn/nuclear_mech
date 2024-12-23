function [fitresult, gof] = createFit(e7, e8)
%CREATEFIT(E7,E8)
%  Create a fit.
%
%  Data for 'untitled fit 1' fit:
%      X Input : e7
%      Y Output: e8
%  Output:
%      fitresult : a fit object representing the fit.
%      gof : structure with goodness-of fit info.
%
%  See also FIT, CFIT, SFIT.

%  Auto-generated by MATLAB on 21-Oct-2023 21:00:25


%% Fit: 'untitled fit 1'.
[xData, yData] = prepareCurveData( e7, e8 );

% Set up fittype and options.
ft = fittype( 'power1' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Lower = [-Inf 1.5];
opts.StartPoint = [0.0666674079653832 1.21987193756246];
opts.Upper = [Inf 1.5];

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

% Plot fit with data.
% figure( 'Name', 'untitled fit 1' );
% h = plot( fitresult, xData, yData );
% legend( h, 'e8 vs. e7', 'untitled fit 1', 'Location', 'NorthEast', 'Interpreter', 'none' );
% % Label axes
% xlabel( 'e7', 'Interpreter', 'none' );
% ylabel( 'e8', 'Interpreter', 'none' );
% grid on


