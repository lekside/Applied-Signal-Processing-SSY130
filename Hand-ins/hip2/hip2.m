%NO_PFILE
% HIP2

% Perform the following steps:
%   1) In student_sols.m, update the student_id variable as described
%   there.
%
%   2) In student_sols.m, complete all the partially-complete functions.
%   (Indicated by a '%TODO: ...' comment). Note that all the functions you
%   need to complete are located in student_sols.m (and only in
%   student_sols.m). You can test these functions by running this file,
%   which will apply a self-test to your functions. When all functions pass
%   the self-test, a unique password will be printed in the terminal. Be
%   sure to include this password in your submission.
%
%   3) Now that the functions in student_sols.m are completed, continue
%   working with this file. Notably, your finished functions will be used
%   to evaluate the behavior of the assignment.
%
% -------------------------------------------------------------------------
%                    Note on function handles
% -------------------------------------------------------------------------
% In this file, we will make use of function handles. A function handle is
% a variable that refers to a function. For example:
%
% x = @plot
%
% assigns a handle to the plot function to the variable x. This allows to
% for example do something like
%
% x(sin(linspace(0,2*pi)))
%
% to call the plot function. Usefully for you, there exist function handles
% to all the functions you've written in student_sols.m. See below for
% exactly how to call them in this assignment.
%
%
% Final note: files with a .p extension are intentionally obfusticated
% (they cannot easily be read). These files contain the solutions to the
% tasks you are to solve (and are used in order to self-test your code).
% Though it is theoretically possible to break into them and extract the
% solutions, doing this will take you *much* longer than just solving the
% posed tasks =)

% Do some cleanup
clc
close all
clear variables
format short eng

% Perform all self-tests of functions in student_sol.m
% % % % % % apply_tests();

% Load student-written functions
funs = student_sols();

% Call your function to get the generated filter coefficients
h = funs.gen_filter();

% Here are some sample plots to illustrate the behavior of your filter.
% Feel free to modify, re-use, or completely remove the following lines.



%% Plot the filter coefficiencts and magnitude/phase response

% figure('Color','white');
% stem(h);
% title('Filter coefficients');
% 
% figure('Color','white');
% N_fft = 1e3;    %Zero-pad FFT for increased frequency resolution
% plot(abs(fft(h, N_fft)));
% title('Filter magnitude response');
% xlabel('A frequency unit (which?)');
% ylabel('|H|');
% 
% figure('Color','white');
% plot(unwrap(angle(fft(h, N_fft))));
% title('Filter phase response');
% xlabel('A frequency unit (which?)');
% ylabel('arg(H)');