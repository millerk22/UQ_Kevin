% Plot high/low confidence MNIST digits
dpHighConfidence = './visualization/MNIST/higher_confidence/';
dpLowConfdience  = './visualization/MNIST/lower_confdience/';
% use subtightplot
subplot = @(m,n,p) subtightplot (m, n, p, [0.01 0.05], [0.1 0.01], [0.1 0.01]);
digitList        = [1, 4, 7, 9];
highCandidate    = {[135,   151,    183,    232],
                    [603,   613,    686,    774],
                    [1125,  1166,   1200,   1394],
                    [1622,  1699,   1729,   1975]
                    };
lowCandidate    =  {[224,   439,    345,    305],
                    [680,   914,    835,    917],
                    [1163,  1205,   1319,   1321],
                    [1598,  1758,   1797,   1834]
                    };
fig = figure;
for  idigit = 1: numel(digitList)
    for i =1:4
        fn = sprintf('digit_%d_image_%d.fig', digitList(idigit),...
            highCandidate{idigit}(i));
        subplot_ = subplot(4, 4, (idigit - 1) * 4 + i);
        fig_     = openfig([dpHighConfidence, fn]);
        ax       = gca;
        fig__     = get(ax, 'children');
        fclose(fig_)
        copyobj(fig__, subplot_);
    end
end