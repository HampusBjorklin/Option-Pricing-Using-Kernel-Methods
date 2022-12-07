function q = rate_of_convergence(err, stepsize)
q = NaN;
if length(err) ~= length(stepsize)
    disp('Error and step size vector not equally long.')

else
    q = zeros(1, length(err)-1);
    for i = 2:length(err)
        q(i) = log(err(i)/err(i-1))/log(stepsize(i)/stepsize(i-1));
    end
    q = mean(q);
    loglog(stepsize(2:end), err(2:end));
    hold on
    loglog(stepsize([2, length(err)]), err([2, length(err)]));
end
