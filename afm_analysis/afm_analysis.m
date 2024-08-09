
nSims = 288;
nPerCases = 288;
values = {};

means = [];
stds = [];

ind = 1;
figure(4);
hold on
tempValues = [];
for i = 1:nSims
        
    a = dlmread(append(".\AFM_INF_multi_sensitivity_",sprintf('%0*d',3,i),".csv"));
    
    if length(a(:,1)) > 1
    
        plot(a(:,1),a(:,2),'-r')
        % pause(1)
        hold on

        [fitresult, gof] = createFit(a(:,1), a(:,2));

        tempValues(end+1) = fitresult.a*3/4/sqrt(3.31e-6)*(1-0.5^2);

        ind = ind + 1;

        if ind == nPerCases+1

            values = [values {tempValues}];

            means(end+1) = mean(tempValues);
            stds(end+1) = std(tempValues);

            ind = 1;
            tempValues = [];
        end
    end

    
end

hold off