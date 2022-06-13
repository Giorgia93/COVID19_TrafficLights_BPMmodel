function new_cov = get_cov(old_cov, elig_pop, p1, verbose)
% Function that recalculates vaccination coverage per age group when given
% a new national coverage. The same proportion of eligible unvaccinated
% people in each group is assumed to get vaccinated to reach the new
% coverage.
% INPUTS:
% - old_cov = row vector of old coverage proportions per age group
% - elig_pop = column vector of eligible population counts per age group
% - p1 = new desired average coverage (as a proportion)
% - verbose = boolean parameter. If =1 prints out new coverage per age gr.
% OUTPUT:
% - new_cov = row vector of new coverage proportions per age group


elig_vax = old_cov' .* elig_pop;
elig_unvax = (1 - old_cov') .* elig_pop;
q = (p1 * sum(elig_pop) - sum(elig_vax)) / sum(elig_unvax);
newvax = q * elig_unvax;
% p1_check = sum(elig_vax + newvax) ./ sum(elig_pop);

new_cov = ((elig_vax + newvax) ./ elig_pop)';
new_cov(elig_pop == 0) = 0;

if verbose == 1
    fprintf("Old cov: %i%%\t%i%%\t%i%%\t%i%%\t%i%%\t%i%%\t%i%%\t%i%%\t%i%%\t%i%%\t%i%%\t%i%%\t%i%%\t%i%%\t%i%%\t%i%%\n", floor(100*old_cov'))
    fprintf("New cov: %i%%\t%i%%\t%i%%\t%i%%\t%i%%\t%i%%\t%i%%\t%i%%\t%i%%\t%i%%\t%i%%\t%i%%\t%i%%\t%i%%\t%i%%\t%i%%\n", floor(100*new_cov'))
    fprintf("New national coverage = %i%% of %i+ vaccinated\n\n", round(100*sum(elig_pop .* new_cov') / sum(elig_pop)));
end


end