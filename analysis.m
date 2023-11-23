
nSims = 108;
nCases = 12;

values = {};

means = [];
stds = [];

ind = 1;
figure(2);
hold on
tempValues = [];
for i = 1:nSims
        
    a = dlmread(append(".\afm_fit_results_w\afm_sim_finalw_",num2str(i),".csv"));
    
    if length(a(:,1)) > 10
    
        plot(a(:,1),a(:,2),'-k')
        hold on

        [fitresult, gof] = createFit(a(:,1), a(:,2));

        tempValues(end+1) = fitresult.a*3/4/sqrt(3.31e-6)*(1-0.5^2);

        ind = ind + 1;

        if ind == nCases+1

            values = [values {tempValues}];

            means(end+1) = mean(tempValues);
            stds(end+1) = std(tempValues);

            ind = 1;
            tempValues = [];
        end
    end

    
end

hold off