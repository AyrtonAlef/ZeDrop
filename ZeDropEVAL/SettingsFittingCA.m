function [UseCircFit,n_circfit,PlotCircleImages,UsePolyFit,n_polyfit,poly_ord,PlotPolynomialImages,UsePolyFitDropen,n_polyfitDropen,poly_ordDropen,PlotPolynomialDropenImages,UseEllipseFit,PlotEllipticImages,UseMaskMethod,fillBlankHorLines,marginMaskBaselineDet,maskSize,PlotMaskImages] = SettingsFittingCA(UseCircFit,n_circfit,PlotCircleImages,UsePolyFit,n_polyfit,poly_ord,PlotPolynomialImages,UsePolyFitDropen,n_polyfitDropen,poly_ordDropen,PlotPolynomialDropenImages,UseEllipseFit,PlotEllipticImages,UseMaskMethod,fillBlankHorLines,marginMaskBaselineDet,maskSize,PlotMaskImages)
%SETTINGSFITTINGCA Enables the user to edit fitting settings
%   INPUT/OUTPUT:
% UseCircFit - Option to enable circle fitting (0 = false, 1 = true)
% n_circfit - number of fitting points for circle fitting (Suggested: 100)
% PlotCircleImages - Option to plot circle fitting results (0 = false, 1 = true)
% UsePolyFit - Option to enable polynomial fitting (0 = false, 1 = true)
% n_polyfit - Number of fitting points for polynomial fitting (Suggested: 10)
% poly_ord - Polynomial order for polynomial fitting (Suggested: 4)
% PlotPolynomialImages - Option to plot polynomial fitting images (0 = false, 1 = true)
% UsePolyFitDropen - Option to enable polynomial fitting from Dropen (0 = false, 1 = true)
% n_polyfitDropen - Number of fitting points for polynomial fitting from Dropen (Suggested: 100)
% poly_ordDropen - polynomial order for polynomial fitting from Dropen (Suggested: (i) for CA < 60° -> poly_ord = 2 and (ii) for CA > 60° -> poly_ord = 3)
% PlotPolynomialDropenImages - Option to plot polynomial fitting (Dropen) images (0 = false, 1 = true)
% UseEllipseFit - Option to enable double sided elliptical fitting (0 = false, 1 = true)
% PlotEllipticImages - Option to plot double sided elliptical fitting images (0 = false, 1 = true)
% UseMaskMethod - Option to enable mask method (0 = false, 1 = true)
% fillBlankHorLines = 1; %Option for filling blank horizontal lines (1 = enable and 0 = disable). Option for Akbari and Antonini 2021 code (BaselineDetOpt == 3).
% marginMaskBaselineDet = 10; %Minimum pixel distance (in y coordinate) between points along the detected edge to be considered part of the drop profile (Suggestion = jum_dist). Option for Akbari and Antonini 2021 code (BaselineDetOpt == 3). 
% maskSize = 20; %Mask size (neighbour index)
% PlotMaskImages - Option to plot mask method images (0 = false, 1 = true)

%Fitting method settings
while 1
    %Current settings
    fprintf('-------------------------------------------------------------- \n');
    fprintf('- Current fitting settings: \n');
    fprintf(strcat("UseCircFit = ",num2str(UseCircFit)," \n"));
    fprintf(strcat("n_circfit = ",num2str(n_circfit)," \n"));
    fprintf(strcat("PlotCircleImages = ",num2str(PlotCircleImages)," \n"));
    fprintf(strcat("UsePolyFit = ",num2str(UsePolyFit)," \n"));
    fprintf(strcat("n_polyfit = ",num2str(n_polyfit)," \n"));
    fprintf(strcat("poly_ord = ",num2str(poly_ord)," \n"));
    fprintf(strcat("PlotPolynomialImages = ",num2str(PlotPolynomialImages)," \n"));
    fprintf(strcat("UsePolyFitDropen = ",num2str(UsePolyFitDropen)," \n"));
    fprintf(strcat("n_polyfitDropen = ",num2str(n_polyfitDropen)," \n"));
    fprintf(strcat("poly_ordDropen = ",num2str(poly_ordDropen)," \n"));
    fprintf(strcat("PlotPolynomialDropenImages = ",num2str(PlotPolynomialDropenImages)," \n"));
    fprintf(strcat("UseEllipseFit = ",num2str(UseEllipseFit)," \n"));
    fprintf(strcat("PlotEllipticImages = ",num2str(PlotEllipticImages)," \n"));
    fprintf(strcat("UseMaskMethod = ",num2str(UseMaskMethod)," \n"));
    fprintf(strcat("fillBlankHorLines = ",num2str(fillBlankHorLines)," \n"));
    fprintf(strcat("marginMaskBaselineDet = ",num2str(marginMaskBaselineDet)," \n"));
    fprintf(strcat("maskSize = ",num2str(maskSize)," \n"));
    fprintf(strcat("PlotMaskImages = ",num2str(PlotMaskImages)," \n"));
    fprintf('-------------------------------------------------------------- \n');

    %Edit Settings
    fprintf('- Edit settings... \n')
    fprintf('Fitting method settings: \n')
    fprintf('1. Circle fitting; \n');
    fprintf('2. Polynomial fitting; \n');
    fprintf('3. Polynomial fitting (Dropen); \n');
    fprintf('4. Ellipse fitting; \n');
    fprintf('5. Mask method; \n');
    ansFittingMethod = input('Select a fitting method option [0 to exit]: ');
    if isempty(ansFittingMethod)
        ansFittingMethod = 0;
    end
    if ansFittingMethod == 0
        break;
    elseif ansFittingMethod == 1 %Circle fitting settings
        fprintf('Circle fitting settings... \n');
        ansCircFit = input('Do you like to enable circle fitting (Y/N) [Y]? ','s');
        if isempty(ansCircFit)
            ansCircFit = 'Y';
        end
        if ansCircFit == 'Y'
            UseCircFit = 1; %Enable circle fitting
            n_circfit = input('Number of fitting points [100]: ');
            if isempty(n_circfit)
                n_circfit = 100;
            end
            ansPlotCircFit = input('Do you like to plot circle fitting results (Y/N) [Y]? ','s');
            if isempty(ansPlotCircFit)
                ansPlotCircFit = 'Y';
            end
            if ansPlotCircFit == 'Y'
                PlotCircleImages = 1;
            else
                PlotCircleImages = 0;
            end
        else
            UseCircFit = 0; %Disable circle fitting
            PlotCircleImages = 0; %Disabel circle fitting images
        end
        
    elseif ansFittingMethod == 2 %Polynomial fitting settings
        fprintf('Polynomial fitting settings... \n');
        ansPolyFit = input('Do you like to enable polynomial fitting (Y/N) [Y]? ','s');
        if isempty(ansPolyFit)
            ansPolyFit = 'Y';
        end
        if ansPolyFit == 'Y'
            UsePolyFit = 1; %Enable polynomial fitting
            n_polyfit = input('Number of fitting points [10]: ');
            if isempty(n_polyfit)
                n_polyfit = 10;
            end
            poly_ord = input('Polynomial order [4]: ');
            if isempty(poly_ord)
                poly_ord = 4;
            end
            ansPlotPolyFit = input('Do you like to plot polynomial fitting results (Y/N) [Y]? ','s');
            if isempty(ansPlotPolyFit)
                ansPlotPolyFit = 'Y';
            end
            if ansPlotPolyFit == 'Y'
                PlotPolynomialImages = 1;
            else
                PlotPolynomialImages = 0;
            end
        else
            UsePolyFit = 0; %Disable polynomial fitting
            PlotPolynomialImages = 0; %Disable polynomial fitting images
        end
        
    elseif ansFittingMethod == 3 %Polynomial fitting (Dropen) settings
        fprintf('Polynomial fitting (Dropen) settings... \n');
        ansPolyFitDropen = input('Do you like to enable polynomial fitting (from Dropen) (Y/N) [Y]? ','s');
        if isempty(ansPolyFitDropen)
            ansPolyFitDropen = 'Y';
        end
        if ansPolyFitDropen == 'Y'
            UsePolyFitDropen = 1; %Enable polynomial fitting (Dropen)
            n_polyfitDropen = input('Number of fitting points [100]: ');
            if isempty(n_polyfitDropen)
                n_polyfitDropen = 100;
            end
            poly_ordDropen = input('Polynomial order [3]: ');
            if isempty(poly_ordDropen)
                poly_ordDropen = 3;
            end
            ansPlotPolyFitDropen = input('Do you like to plot polynomial fitting (Dropen) results (Y/N) [Y]? ','s');
            if isempty(ansPlotPolyFitDropen)
                ansPlotPolyFitDropen = 'Y';
            end
            if ansPlotPolyFitDropen == 'Y'
                PlotPolynomialDropenImages = 1;
            else
                PlotPolynomialDropenImages = 0;
            end
        else
            UsePolyFitDropen = 0; %Disable polynomial fitting (Dropen)
            PlotPolynomialDropenImages = 0; %Disable polynomial fitting images
        end
        
    elseif ansFittingMethod == 4 %Ellipse fitting settings
        fprintf('Ellipse fitting settings... \n');
        ansEllipseFit = input('Do you like to enable double sided elliptical fitting (Y/N) [Y]? ','s');
        if isempty(ansEllipseFit)
            ansEllipseFit = 'Y';
        end
        if ansEllipseFit == 'Y'
            UseEllipseFit = 1; %Enable ellipse fitting
            ansPlotEllipseFit = input('Do you like to plot double sided elliptical fitting results (Y/N) [Y]? ','s');
            if isempty(ansPlotEllipseFit)
                ansPlotEllipseFit = 'Y';
            end
            if ansPlotEllipseFit == 'Y'
                PlotEllipticImages = 1;
            else
                PlotEllipticImages = 0;
            end
        else
            UseEllipseFit = 0; %Disable ellipse fitting
            PlotEllipticImages = 0; %Disable ellipse fitting images
        end
        
    elseif ansFittingMethod == 5 %Mask method settings
        fprintf('Mask method settings... \n');
        ansMaskMethod = input('Do you like to enable mask method (Y/N) [Y]? ','s');
        if isempty(ansMaskMethod)
            ansMaskMethod = 'Y';
        end
        if ansMaskMethod == 'Y'
            UseMaskMethod = 1; %Enable mask method
            ansfillBlankHorLines = input('Fill blank horizontal lines (Y/N) [Y]? ','s');
            if isempty(ansfillBlankHorLines)
                ansfillBlankHorLines = 'Y';
            end
            if ansfillBlankHorLines == 'N'
                fillBlankHorLines = 0;
            else
                fillBlankHorLines = 1;
            end
            marginMaskBaselineDet = input('Minimum pixel distance (y coordinate) along the detected edge to be considered part of the drop profile [10]: ');
            if isempty(marginMaskBaselineDet )
                marginMaskBaselineDet  = 10;
            end
            maskSize = input('Mask size (neighbour index) [20]: ');
            if isempty(maskSize)
                maskSize  = 20;
            end
            ansPlotMaskMethod = input('Do you like to plot mask method results (Y/N) [Y]? ','s');
            if isempty(ansPlotMaskMethod)
                ansPlotMaskMethod = 'Y';
            end
            if ansPlotMaskMethod == 'Y'
                PlotMaskImages = 1;
            else
                PlotMaskImages = 0;
            end
        else
            UseMaskMethod = 0; %Disable mask method
            PlotMaskImages = 0; %Disable mask method images
        end
        
    end
end
end

