function f_handle = plot_FT(X, name, figpos)

    N = length(X);
    f_handle = figure('Color','white','Position',figpos);
    sb1 = subplot(2,1,1);
    plot( linspace(0,1-1/N,N) ,abs(X), 'Marker','o','MarkerFaceColor','red', 'MarkerSize', 3); grid on;
    ylabel(['magnitude[ ',name,'(k) ]'])
    xlabel('w / w_s')
    sb2 = subplot(2,1,2);
    plot( linspace(0,1-1/N,N) ,angle(X)*180/pi, 'Marker','o','MarkerFaceColor','red', 'MarkerSize', 3); grid on;
    ylabel(['angle[ ',name,'(k) ] [º]'])
    xlabel('w / w_s')
    
    linkaxes([sb1, sb2], 'x');
end

