function f_handle = plot_data_time(z, name, figpos)

    N = length(z);
    f_handle = figure('Color','white','Position',figpos);
    plot(real(z))
    sb1 = subplot(2,1,1);
    plot( linspace(0,N,N) ,real(z), 'Marker','o','MarkerFaceColor','red', 'MarkerSize', 3); grid on;
    ylabel(['real[ ',name,'(n) ]'])
    xlabel('N / f_s [s]')
    sb2 = subplot(2,1,2);
    plot( linspace(0,N,N) ,imag(z), 'Marker','o','MarkerFaceColor','red', 'MarkerSize', 3); grid on;
    ylabel(['im[ ',name,'(k) ]'])
    xlabel('N / f_s [s]')
    
    linkaxes([sb1, sb2], 'x');
end

