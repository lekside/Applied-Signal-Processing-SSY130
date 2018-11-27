function f_handle = plot_FT_fs(X, name, fs, filename, figpos)

    N = length(X);
    f_handle = figure('Color','white','Position',figpos);
    sb1 = subplot(2,1,1);
    plot( linspace(0,1-1/N,N)*fs ,abs(X), 'Marker','o','MarkerFaceColor','red', 'MarkerSize', 3); grid on;
    ylabel(['magnitude[ ',name,'(k) ]'])
    xlabel('f [Hz]')
    sb2 = subplot(2,1,2);
    plot( linspace(0,1-1/N,N)*fs ,angle(X)*180/pi, 'Marker','o','MarkerFaceColor','red', 'MarkerSize', 3); grid on;
    ylabel(['angle[ ',name,'(k) ] [º]'])
    xlabel('f [Hz]')
    
    linkaxes([sb1, sb2], 'x');
    
    if ~strcmp(name,'')
        if ~ exist(fullfile(pwd,'images'),'dir'), mkdir images; end
        set(gca,'LooseInset',get(gca,'TightInset'))
        saveas(gcf, fullfile(pwd,'images/',filename),'epsc')
    end
end

